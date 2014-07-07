#!/usr/bin/env python

###############################################################################
# SOME COMMENT...
###############################################################################

import os
import getopt
import sys
from random import sample

def usage():
    """This function displays the usage of the script."""

    print "\nUSAGE: ./thisprogram.py -i alignment_file[.aln] -o output_file\n"

def check_input_file(name_file):
    """This function ensures that the input file path is properly defined."""

    if os.path.exists(name_file) == False:
        print "\nFile named '" + name_file + "' does not exist. Exiting."
        sys.exit(2)

class Column:
    """This class creates an object containing all the informations about a
    column coming from a genomic alignment.
    """

    def __init__(self, orgs, seq):
        self._seq = seq
        self._orgs = orgs
        self._sorted_orgs = self.sort_orgs(orgs)
        self._sorted_seq = self.sort_seq(self._sorted_orgs, orgs, seq)
        self._sorted_pattern = self.build_pattern(self._sorted_seq)

    def build_pattern(self, sorted_seq):
        if sorted_seq[0] == "-":
            pattern = "-"
            for lett in sorted_seq[1:]:
                if lett == sorted_seq[0]:
                    pattern += "-"
                else:
                    pattern += "0"
        else:
            pattern = "1"
            for lett in sorted_seq[1:]:
                if lett == sorted_seq[0]:
                    pattern += "1"
                elif lett == "-":
                    pattern += "-"
                else:
                    pattern += "0"
        return pattern

    def sort_orgs(self, orgs):
        sorted_orgs = sorted(orgs)
        return sorted_orgs

    def sort_seq(self, sorted_orgs, orgs, seq):
        new_seq = []
        pairs = dict(zip(orgs, seq))

        for s in sorted_orgs:
           new_seq.append(pairs[s])
        return new_seq

def read_input_alignments():
    """This function takes a single alignment at the time in the form of a
    dictionary ( ex: {'>rheMac2': ['-', 'A', ...], '>>mm9': ['A', 'T', ...]} )
    and stores its columns as single objects into the dictionary 'all_columns'.
    Ex: all_columns = {'>>mm9_>hg18':{'10':[<col1>, <col2>,...], '1-':[...]}, 
                       '>>mm9_>rheMac2_hg18':{...}}.
    """

    inp = open(aln_file)
    all_columns = {}
    aln = {}
    while True:
        line = inp.readline()
        if not line or line == "\n":
            break
        # beginning of a new alignment
        if line.startswith(">>") == True and aln != {}:
            store_columns(aln, all_columns)
            aln = {}
        aln[line[:line.find("_")]] = list((inp.readline()).strip("\n"))
    # process the last alignment of the file
    store_columns(aln, all_columns)
    inp.close()
    return all_columns

def store_columns(aln, all_columns, ref = '>>hg19'):
    """This function takes as input an alignment like the following
       Ex: {'>rheMac2': ['-', 'A', ...], '>>mm9': ['A', 'T', ...]}
    and creates a list of objects, one per each column of the alignment.
    Each object is added to the whole set of columns, 'all_columns', that is
    passed to the upper function 'read_input_alignments()'.
    """

    orgs = aln.keys()
    for i in range(len(aln[ref])):
        seq = []
        for o in orgs:
            seq.append(aln[o][i])
        # create an object of the current alignment column
        col = Column(orgs, seq)
        # join the sorted orgs names of the column in a string
        key = "_".join(col._sorted_orgs)
        # update all_columns adding the new col
        try:
            all_columns[key][col._sorted_pattern].append(col)
        except KeyError:
            try:
                all_columns[key][col._sorted_pattern] = [col]
            except KeyError:
                all_columns[key] = {col._sorted_pattern: [col]}
    return all_columns

def shuffle_input_alignments():
    """This function takes a single alignment at the time and passes it to the
    function 'exchange_columns()', which substitutes each column with another
    having the same alignment pattern. The output file contains the shuffled
    version of the input alignments.
    """

    inp = open(aln_file)
    out = open(out_file, "w")
    aln = {}
    while True:
        line = inp.readline()
        if not line or line == "\n":
            break
        # beginning of a new alignment
        if line.startswith(">>") == True and aln != {}:
            new_aln = exchange_columns(aln, all_columns)
            # write the shuffled alignment in 'out_file'
            out.write(">>hg19\n" + "".join(new_aln[">>hg19"]) + "\n")
            for key in sorted(new_aln.keys()):
                if key != ">>hg19":
                    out.write(key + "\n")
                    out.write("".join(new_aln[key]) + "\n")
            aln = {}
        aln[line[:line.find("_")]] = list((inp.readline()).strip("\n"))
    # process the last alignment of the file
    new_aln = exchange_columns(aln, all_columns)
    # write the shuffled alignment in 'out_file'
    out.write(">>hg19\n" + "".join(new_aln[">>hg19"]) + "\n")
    for key in sorted(new_aln.keys()):
        if key != ">>hg19":
            out.write(key + "\n")
            out.write("".join(new_aln[key]) + "\n")
    inp.close()
    out.close()

def exchange_columns(aln, all_columns, ref = '>>hg19'):
    """This function takes as input an alignment like the following
       Ex: {'>rheMac2': ['-', 'A', ...], '>>mm9': ['A', 'T', ...]}
    and creates a list of objects, one per each column of the alignment.
    Then, for each column, it looks in the dictionary 'all_columns' for another
    column having the same organisms and alignment pattern. Once found, the new
    column is added to the returned list and removed from 'all_columns'. If a
    substitute column cannot be found, it is simply kept in the same position.
    """

    orgs = aln.keys()       # not sorted list of alignment organisms
    new_cols = []
    for i in range(len(aln[ref])):
        seq = []
        for o in orgs:
            seq.append(aln[o][i])
        # create an object of the current alignment column
        col = Column(orgs, seq)
        # join the sorted orgs names of the column in a string
        key = "_".join(col._sorted_orgs)
        # get the column's sorted pattern (ex: '10-1')
        pat = col._sorted_pattern
        # skip the substitution if there are no other columns with same pattern
        if len(all_columns[key][pat]) == 1 and \
        all_columns[key][pat][0]._sorted_seq == col._sorted_seq:
            new_cols.append(col)
            all_columns[key][pat].pop()
            continue
        # search a random column with same pattern but different nuclotides
        else:
            pat_col = len(all_columns[key][pat])
            for rnd in sample(xrange(pat_col), pat_col):
                if all_columns[key][pat][rnd]._sorted_seq != col._sorted_seq:
                    break
            new_cols.append(all_columns[key][pat][rnd])
        # remove the substituted column from 'all_column'
        all_columns[key][pat].pop(rnd)
    # convert the list of objects into an alignment dictionary
    new_aln = {} ; i = 0
    for org in sorted(orgs):
        for obj in range(len(new_cols)):
            try:
                new_aln[org].append(new_cols[obj]._sorted_seq[i])
            except KeyError:
                new_aln[org] = [new_cols[obj]._sorted_seq[i]]
        i += 1
    return new_aln


if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:")
        if len(sys.argv) == 1:
            usage()
            sys.exit(2)
        if len(sys.argv) > 2 and len(sys.argv) < 5:
            print "\nOne or more input arguments is missing. Exiting."
            usage()
            sys.exit(2)    
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)
    for o, a in opts:
        if o == "-i":
            aln_file = a
            check_input_file(aln_file)
        elif o == "-o":
            out_file = a
        else:
            assert False, "unhandled option"

    all_columns = read_input_alignments()
    shuffle_input_alignments()
