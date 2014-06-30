#!/usr/bin/env python

import os
import sys
import math
from datetime import datetime
import re
import gzip

def usage():
    """This function displays the usage of the script."""

    print "\nUSAGE: ./thisprogram.py input_file[.fastq] output_folder_path plot_before_filepath(.pdf) plot_after_filepath(.pdf) log_file\n"


def foobar(file_path):

    gzipped = True
    try:
        fin = gzip.open(file_path)
        fin.readline()
    except IOError:
        gzipped = False
    fin.close()

    if gzipped:
        fin = gzip.open(file_path)
    else:
        fin = open(file_path)

    elem_list = []
    for line in fin:
        elem_list.append(line.rstrip())
        if len(elem_list) == 4:
            yield tuple(elem_list)
            elem_list = []

    if len(elem_list) != 0:
        print "check file, perhaps it is corrupted!"
        print elem_list
    fin.close()


class ParseFastQ(object):
    """This class returns a read-by-read fastq parser. It can be used as:
    parser.next()
    - OR -
    for rec in parser:
        ... do something with rec ...
    rec is a tuple: (seqHeader,seqStr,qualHeader,qualStr).
    """

    def __init__(self, file_path):
        self._file = open(file_path, 'rU')
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


def collect_min_max_ascii(inp_file):
    """This function gets the numerical value corresponding to min and max
    ASCII character contained in the quality lines of the input FASTQ file.
    """

    #parser = ParseFastQ(inp_file)
    min_ascii = "NULL"; max_ascii = "NULL"
    read_n = 0
    for record in foobar(inp_file):
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

def detect_codification(min_ascii, max_ascii):
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

def update_read_length_frequencies(curr_l, read_length_frequencies):
    """This function takes as input the length of the current read and updates
    the corresponding length frequency in the input dictionary.
    """

    try:
        read_length_frequencies[curr_l] += 1
    except KeyError:
        read_length_frequencies[curr_l] = 1

def phred_quality_filter(inp_file, out_file, offset, encoding):
    """This function filters the reads contained in the input file assuming 
    that their qualities are encoded in Phred scores. It returns a fasta file.
    """
    max_n = 2             # Max amount of Ns in the filtered read sequence
    min_length = 25       # Min length for the filtered read
    min_quality = 20      # Min average quality along the filtered read

    #parser = ParseFastQ(inp_file)
    tot_reads = 0
    good_reads = 0

    tot_read_length_frequencies = {}
    good_read_length_frequencies = {}

    o = open(out_file, 'w')

    for record in foobar(inp_file):
        tot_reads += 1
        q_sum = 0
        curr_l = 0
        update_read_length_frequencies(len(record[3]),
                                       tot_read_length_frequencies)
        for i in xrange(len(record[3])):
            q_sum += ord(record[3][i]) - offset
            q_mean = q_sum / float(i + 1)
            if q_mean >= min_quality:
                curr_l = i + 1
        if record[1][:curr_l].count("N") <= max_n:
            if curr_l >= min_length:
                o.write(">%s\n%s\n" % (record[0][1:], record[1][:curr_l]))
                good_reads += 1
                update_read_length_frequencies(curr_l,
                                               good_read_length_frequencies)
    o.close()

    print 'Codification type: ', encoding
    print 'Original number of reads: ', tot_reads
    print 'Reads that passed the quality filter: ', good_reads, "(% 0.2f" % \
    ((good_reads / float(tot_reads)) * 100), "%)"

    return tot_reads, good_reads, tot_read_length_frequencies, \
           good_read_length_frequencies

def solexa_quality_filter(inp_file, out_file, offset, encoding):
    """This function filters the reads contained in the input file assuming 
    that their qualities are encoded in Solexa scores. The formula to convert
    Solexa into Phred comes from 'Cock et al., Nucleic Acids Research 2010'.
    It returns a fasta file.
    """
    max_n = 2             # Max amount of Ns in the filtered read sequence
    min_length = 25       # Min length for the filtered read
    min_quality = 20      # Min average quality along the filtered read

    #parser = ParseFastQ(inp_file)
    tot_reads = 0
    good_reads = 0

    tot_read_length_frequencies = dict()
    good_read_length_frequencies = dict()

    o = open(out_file, 'w')

    for record in foobar(inp_file):
        tot_reads += 1
        q_sum = 0
        curr_l = 0
        update_read_length_frequencies(len(record[3]),
                                       tot_read_length_frequencies)
        for i in xrange(len(record[3])):
            q_sum += 10 * math.log10(10 ** ((ord(record[3][i]) - offset) / 
            10.0) + 1)    # The quality is converted from Solexa to Phred. 
            q_mean = q_sum / float(i + 1)
            if q_mean >= min_quality:
                curr_l = i + 1
        if record[1][:curr_l].count("N") <= max_n:
            if curr_l >= min_length:
                o.write('>%s\n%s\n' % (record[0][1:], record[1][:curr_l]))
                good_reads += 1
                update_read_length_frequencies(curr_l,
                                               good_read_length_frequencies)
    o.close()

    print 'Codification type:', encoding
    print 'Original number of reads:', tot_reads
    print 'Reads that passed the quality filter:', good_reads, "(% 0.2f" % \
    ((good_reads / float(tot_reads)) * 100), "%)"

    return tot_reads, good_reads, tot_read_length_frequencies, \
           good_read_length_frequencies


def plot_read_length_frequencies(plot_before, plot_after, tot_read_length_frequencies,
                                 good_read_length_frequencies):
    """This function takes as input the dictionary of the read length
    frequencies before and after the quality filter and plots their distribution
    through an R script.
    """

    # create the x and y vectors of tot reads for R
    x_tot = list(); y_tot = list()
    for k,v in sorted(tot_read_length_frequencies.iteritems()):
        x_tot.append(k)
        y_tot.append(v)

    # create the x and y vectors of good reads for R
    x_good = list(); y_good = list()
    for k,v in sorted(good_read_length_frequencies.iteritems()):
        x_good.append(k)
        y_good.append(v)

    r_script_name = os.path.join(os.path.split(out_file)[0], "r_script.R")
    r_script = open(r_script_name, "w")
    r_script.write("x_tot=c(" + ",".join(["%s" % e for e in x_tot]) + ")\n")
    r_script.write("y_tot=c(" + ",".join(["%s" % e for e in y_tot]) + ")\n")
    r_script.write("x_good=c(" + ",".join(["%s" % e for e in x_good]) + ")\n")
    r_script.write("y_good=c(" + ",".join(["%s" % e for e in y_good]) + ")\n")
    r_script.write("min_x=min(x_tot, x_good)\n")
    r_script.write("max_x=max(x_tot, x_good)\n")
    r_script.write("min_y=min(y_tot, y_good)\n")
    r_script.write("max_y=max(y_tot, y_good)\n")
    r_script.write("y_tot_revcumsum=rev(cumsum(rev(y_tot)))\n")
    r_script.write("y_good_revcumsum=rev(cumsum(rev(y_good)))\n")

    # commands to plot the read length distribution before quality filtering
    r_script.write("pdf('" + plot_before + "')\n")
    r_script.write("par(mar=c(5,6,4,2)+0.1)\n")
    r_script.write("plot(x_tot, y_tot_revcumsum, main='Read length distribution before \
QC', xlab='read length', ylab='number of reads with at least read length\\n', log='y', type='o', pch=19, \
col='blue', xlim=c(min_x,max_x), ylim=c(min_y,max_y), las=1)\n")

    # commands to plot the read length distribution after quality filtering
    r_script.write("pdf('" + plot_after + "')\n")
    r_script.write("par(mar=c(5,6,4,2)+0.1)\n")
    r_script.write("plot(x_good, y_good_revcumsum, main='Read length distribution after \
QC', xlab='read length', ylab='number of reads with at least read length\\n', log='y', type='o', pch=19, \
col='blue', xlim=c(min_x,max_x), ylim=c(min_y,max_y), las=1)\n")

    r_script.write("dev.off()\n\n")
    r_script.close()

    # execute the R script
    os.system("R --vanilla < " + r_script_name)


###############################################################################
# MAIN #
###############################################################################

if __name__ == '__main__':

    T1 = datetime.now()

    # Read input parameters and set fixed ones
    #if len(sys.argv) != 3:
    #    usage()           # Display usage and exit
    #    sys.exit(2)
    #inp_file = sys.argv[1]

    inp_file = sys.argv[1]

    out_file = sys.argv[2]
    plot_before = sys.argv[3]
    plot_after = sys.argv[4]
    log_file = sys.argv[5]

    logf = open(log_file, 'w')

    # Check if out_folder exists, otherwise create it

    # Detect quality scores codification, filter reads and create a fasta file
    min_ascii, max_ascii = collect_min_max_ascii(inp_file)
    encoding, offset = detect_codification(min_ascii, max_ascii)

    T2 = datetime.now()
    print T2 - T1, "time for encoding detection"
    logf.write("%s time for encoding detection\n" %(T2-T1))

    if encoding == "Sanger" or encoding == "Illumina":
        tot_reads, good_reads, tot_read_length_frequencies, \
        good_read_length_frequencies = phred_quality_filter(inp_file,
                                       out_file, offset, encoding)
    elif encoding == "Solexa":
        tot_reads, good_reads, tot_read_length_frequencies, \
        good_read_length_frequencies = solexa_quality_filter(inp_file,
                                       out_file, offset, encoding)
    #plot_read_length_frequencies(out_folder, tot_read_length_frequencies, good_read_length_frequencies)
    plot_read_length_frequencies(plot_before, plot_after, tot_read_length_frequencies,
                                 good_read_length_frequencies)

    T3 = datetime.now()
    print T3 - T2, "time for quality filtering"
    logf.write("%s time for quality filtering\n" %(T3-T2))
    print T3 - T1, "total time"
    logf.write("%s total time" %(T3-T1))

    logf.close()
