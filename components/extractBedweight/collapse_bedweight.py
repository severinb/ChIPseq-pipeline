#!/usr/bin/env python

import fileinput

def main():
    """
    This codes reads from stdin.
    Input data need to be sorted by 5' ends (pioSortBed9 -s5).
    For reads starting at the same position (same 5' end) weights in 5th column get summed.
    For + strand reads I print chrom, read_start, read_start +1, ., strand
    For - strand reads I print chrom, read_end-1, read_end, ., strand
    Output is again sorted by 5' end.
    """

    init = True

    for line in fileinput.input():
        l = line.split()
        if init:
            curr_chrom, curr_start, curr_end = l[:3]
            curr_score = float(l[4])
            curr_strand = l[5]
            init = False
        else:
            if l[0] == curr_chrom:
                if l[5] == curr_strand:
                    if curr_strand == '+':
                        if l[1] == curr_start:
                            curr_score += float(l[4])
                        else:
                            print '\t'.join([curr_chrom, curr_start, str(int(curr_start) + 1), '.', str(curr_score), curr_strand])
                            curr_start, curr_end = l[1:3]
                            curr_score = float(l[4])
                            curr_strand = l[5]
                    else:
                        if l[2] == curr_end:
                            curr_score += float(l[4])
                        else:
                            print '\t'.join([curr_chrom, str(int(curr_end) - 1), curr_end, '.', str(curr_score), curr_strand])
                            curr_start, curr_end = l[1:3]
                            curr_score = float(l[4])
                            curr_strand = l[5]
                else:
                    if curr_strand == '+':
                        print '\t'.join([curr_chrom, curr_start, str(int(curr_start) + 1), '.', str(curr_score), curr_strand])
                        curr_start, curr_end = l[1:3]
                        curr_score = float(l[4])
                        curr_strand = l[5]
                    else:
                        print '\t'.join([curr_chrom, str(int(curr_end) - 1), curr_end, '.', str(curr_score), curr_strand])
                        curr_start, curr_end = l[1:3]
                        curr_score = float(l[4])
                        curr_strand = l[5]
            else:                
                if curr_strand == '+':
                    print '\t'.join([curr_chrom, curr_start, str(int(curr_start) + 1), '.', str(curr_score), curr_strand])
                    curr_chrom, curr_start, curr_end = l[:3]
                    curr_score = float(l[4])
                    curr_strand = l[5]
                else:
                    print '\t'.join([curr_chrom, str(int(curr_end) - 1), curr_end, '.', str(curr_score), curr_strand])
                    curr_chrom, curr_start, curr_end = l[:3]
                    curr_score = float(l[4])
                    curr_strand = l[5]

    if curr_strand == '+':
        print '\t'.join([curr_chrom, curr_start, str(int(curr_start) + 1), '.', str(curr_score), curr_strand])
    else:
        print '\t'.join([curr_chrom, str(int(curr_end) - 1), curr_end, '.', str(curr_score), curr_strand])

if __name__ == '__main__':
    main()
