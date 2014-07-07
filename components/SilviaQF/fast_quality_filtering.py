#!/usr/bin/env python

import os
import sys
import math
import re
from datetime import datetime
import gzip

def usage():
    """This function displays the usage of the script."""

    print "\nUSAGE: ./thisprogram.py input_file[.fastq] output_folder_path\n"


def check_input_path(path):
    """This function creates a folder if it does not already exist."""

    if os.path.exists(path) == False:
        os.system("mkdir " + path)

class ParseFastQ(object):
    """This class returns a read-by-read fastq parser. It can be used as:
    parser.next()
    - OR -
    for rec in parser:
        ... do something with rec ...
    rec is a tuple: (seqHeader,seqStr,qualHeader,qualStr).
    """

    def __init__(self, file_path):
        self._file = gzip.open(file_path) #, 'rU')
        self._current_line_number = 0

    def __iter__(self):
        return self

    def next(self):
        """This function reads the next element, parses it and does some
        controls. It returns a tuple: (seqHeader,seqStr,qualHeader,qualStr).
        """

        # Get next four lines.
        elem_list = []
        counter = 0
        while counter < 4:
            line = self._file.readline()
            self._current_line_number += 1
            if line:
                elem_list.append(line.strip('\n'))
            else:
                elem_list.append(None)
            counter += 1

        # Check lines for expected form.
        trues = [bool(x) for x in elem_list].count(True)
        nones = elem_list.count(None)

        # Check for acceptable end of file.
        if nones == 4:
            raise StopIteration

        # Make sure we got 4 full lines of data.
        assert trues == 4, \
            "ERROR: It looks like I encountered a premature EOF or empty line.\
            Please check the file near line %s" % (self._current_line_number)

        # Make sure we are in the correct "register".
        assert elem_list[0].startswith("@"), \
            "ERROR: The 1st line in fastq element does not start with '@'.\
            Please check FastQ file and try again."
        assert elem_list[2].startswith("+"), \
            "ERROR: The 3rd line in fastq element does not start with '+'.\n\
            Please check FastQ file and try again."

        # Return fastQ data as tuple.
        return tuple(elem_list)

def collect_min_max_ascii():
    """This function gets the numerical value corresponding to min and max
    ASCII character contained in the quality lines of the input FASTQ file.
    """

    parser = ParseFastQ(inp_file) 
    min_ascii = "NULL"; max_ascii = "NULL"
    read_n = 0
    for record in parser:
        read_n += 1
        if read_n >= 1000000:
           break
        min_c = min(record[3]); min_n = ord(min_c)
        max_c = max(record[3]); max_n = ord(max_c)
        if min_n <= min_ascii or min_ascii == "NULL":
            min_ascii = min_n
        if max_n >= max_ascii or max_ascii == "NULL":
            max_ascii = max_n
    return min_ascii, max_ascii

def detect_codification():
    """This function detects the quality codification used in the FASTQ file:

    SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS..................................
    ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX...
    ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII...
    .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ...
    LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL.................................
    !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijk
    |                         |    |        |                              |   
    33                        59   64       73                            104  

    S - Sanger        Phred+33,  raw reads typically (0, 40)
    X - Solexa        Solexa+64, raw reads typically (-5, 40)
    I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
    J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
       with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
    L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)

    (taken from http://en.wikipedia.org/wiki/FASTQ_format#Encoding).
    For simplicity, S and L are joined as 'Sanger', I and J as 'Illumina' and
    X as 'Solexa'. In the most ambiguous case (64 < min, max < 73), the reads
    qualities will be assigned to 'Sanger'.
    """

    if max_ascii > 74:
        if min_ascii < 64:
            encoding = "Solexa"
            offset = 64
        else:
            encoding = "Illumina"
            offset = 64
    else:
        encoding = "Sanger"
        offset = 33
    return encoding, offset

def phred_quality_filter():
    """This function filters the reads contained in the input file assuming 
    that their qualities are encoded in Phred scores. It returns a fasta file.
    """

    out_file = open(out_folder + "/" + infRoot + ".fa", "w")
    parser = ParseFastQ(inp_file) 
    tot_reads = 0
    good_reads = 0

    for record in parser:
        tot_reads += 1
        q_sum = 0
        curr_l = 0
        for i in range(len(record[3])):
            q_sum += ord(record[3][i]) - offset
            q_mean = q_sum / float(i + 1)
            if q_mean >= min_quality:
                curr_l = i + 1
        if record[1][:curr_l + 1].count("N") <= max_n:
            if curr_l >= min_length:
                out_file.write('>' + record[0][1:] + '\n' + 
                record[1][:curr_l + 1] + '\n')
                good_reads += 1
    out_file.close()

    print 'Codification type: ', encoding
    print 'Original number of reads: ', tot_reads
    print 'Reads that passed the quality filter: ', good_reads, "(% 0.2f" % \
    ((good_reads / float(tot_reads)) * 100), "%)"

    return tot_reads, good_reads

def solexa_quality_filter():
    """This function filters the reads contained in the input file assuming 
    that their qualities are encoded in Solexa scores. The formula to convert
    Solexa into Phred comes from 'Cock et al., Nucleic Acids Research 2010'.
    It returns a fasta file.
    """

    out_file = open(out_folder + "/" + infRoot + ".fa", "w") #.fa
    parser = ParseFastQ(inp_file)
    tot_reads = 0
    good_reads = 0

    for record in parser:
        tot_reads += 1
        q_sum = 0
        curr_l = 0
        for i in range(len(record[3])):
            q_sum += 10 * math.log10(10 ** ((ord(record[3][i]) - offset) / 
            10.0) + 1)    # The quality is converted from Solexa to Phred. 
            q_mean = q_sum / float(i + 1)
            if q_mean >= min_quality:
                curr_l = i + 1
        if record[1][:curr_l + 1].count("N") <= max_n:
            if curr_l >= min_length:
                out_file.write('>' + record[0][1:] + '\n' +
                record[1][:curr_l + 1] + '\n')
                good_reads += 1
    out_file.close()

    print 'Codification type:', encoding
    print 'Original number of reads:', tot_reads
    print 'Reads that passed the quality filter:', good_reads, "(% 0.2f" % \
    ((good_reads / float(tot_reads)) * 100), "%)"

    return tot_reads, good_reads

###############################################################################
# MAIN #
###############################################################################

if __name__ == '__main__':

    T1 = datetime.now()

    # Read input parameters and set fixed ones
    if len(sys.argv) != 3:
        usage()           # Display usage and exit
        sys.exit(2)
    
    inp_file = sys.argv[1]
    [path, filename] = os.path.split(inp_file)
    infRoot = re.sub('\.fastq', '', filename)
    infRoot = re.sub('\.gz$', '', infRoot)

    out_folder = sys.argv[2]
   
    max_n = 2             # Max amount of Ns in the filtered read sequence
    min_length = 25       # Min length for the filtered read
    min_quality = 20      # Min average quality along the filtered read

    # Check if out_folder exists, otherwise create it
    check_input_path(out_folder)

    # Detect quality scores codification, filter reads and create a fasta file
    min_ascii, max_ascii = collect_min_max_ascii()
    encoding, offset = detect_codification()

    T2 = datetime.now()
    print T2 - T1, "time for encoding detection"

    if encoding == "Sanger" or encoding == "Illumina":
        tot_reads, good_reads = phred_quality_filter()
    elif encoding == "Solexa":
        tot_reads, good_reads = solexa_quality_filter()

    T3 = datetime.now()
    print T3 - T2, "time for quality filtering"
    print T3 - T1, "total time"
    
