#!/usr/bin/env python

import component_skeleton.main
import os, sys
from pylab import *

def annotate_human_mouse(peaks, annotfile, outfile, r=500):

    # parse annotations:
    #chrom  strand  begin   end     length  cluster clusterType     motevoRegion    CpGclass(motevoRegion)  clusterMembers  genesOfNotNecessarilyAllClusterMembers[member|symbol|gid|name]
    # chr1    +       69090   69090   1       hg19_v1_chr1_+_69090_69090      transcript      hg19_chr1_68590_69590_+ at      NM_001005484    NM_001005484|OR4F5|79501|olfactory receptor, family 4, subfamily F, member 5
    # chr1    +       367639  367658  20      hg19_v1_chr1_+_367639_367658    transcripts     hg19_chr1_367139_368158_+       at      BC137547 BC137568 NM_001005221 NM_001005224 NM_001005277        BC137547|OR4F3|26683|olfactory receptor, family 4, subfamily F, member 3        BC137568|OR4F3|26683|olfactory receptor, family 4, subfamily F, member 3        NM_001005221|OR4F29|729759|olfactory receptor, family 4, subfamily F, member 29 NM_001005224|OR4F3|26683|olfactory receptor, family 4, subfamily F, member 3    NM_001005277|OR4F16|81399|olfactory receptor, family 4, subfamily F, member 16)

    # key: chrom_start_end, value: [transcripts, genes_of_transcripts ], i.e. promoter_coordinates: [transcripts, genes]
    # it's a dict of dicts. first key chromosome: second key: coords: values: [transcripts, genes_of_transcripts ], i.e. promoter_coordinates: [transcripts, genes]
    annot_dict = {}

    for line in open(annotfile):
        if line.startswith('#'):
            continue

        t = line.strip().split('\t')

        try:
            annot_dict[t[0]]['%s_%s' %(t[2], t[3])] = [t[9], t[10:]]
        except KeyError:
            annot_dict[t[0]] = {}
            annot_dict[t[0]]['%s_%s' %(t[2], t[3])] = [t[9], t[10:]]

        #annot_dict['%s_%s_%s' %(t[0], t[2], t[3])] = [t[9], t[10:]]


    # find for each peak whether it is in a promoter. Extend promoters by +-r bps.
    # chr1    231473620       231473733       reg1000000.p1   20.250  +
    # chr19   54605995        54606123        reg1000001.p1   19.279  +
    # chr7    112580129       112580232       reg1000008.p2   18.527  +

    annotated = 0
    nonannotated = 0
    mto = 0

    o = open(outfile, 'w')
    o.write('#chrom\tstart\tstop\tpeakID\tZ-score\tstrand\tpromoter_members\tgenes[member|symbol|gid|name]\n')

    for line in open(peaks):
        t = line.strip().split()
        chrom = t[0]
        start = int(t[1])
        stop = int(t[2])

        found = []
        try:
            for p in annot_dict[chrom]:
                pt = p.split('_')

                pstart = int(pt[0])
                pstop = int(pt[1])

                if start >= pstart-r and stop <= pstop+r:
                    found.append(annot_dict[chrom][p])
        except KeyError:
            print chrom
            nonannotated += 1
            continue

        if len(found) > 0:
            annotated += 1
            o.write('%s\t%s\t%s\t%s\t%s\t%s\t' %(chrom, start, stop, t[3], t[4], t[5]))
            for i in found:
                o.write('%s ' %i[0])
            for i in found:
                for j in i[1]:
                    o.write('\t%s' %j)
            o.write('\n')
            if len(found)>1:
                mto +=1
        else:
            nonannotated += 1


    o.close()


    print annotated
    print nonannotated
    print mto

    return annotated, nonannotated, mto



def annotate_drosophila(peaks, annotfile, outfile, r=500):


    # parse annotations:
    # compared to human and mouse I don't have promoters for drosophila. Thus maybe look for the closest gene. Or Gene inside some range of peak
    # #Chromosome Name        Gene Start (bp) Gene End (bp)   Associated Gene Name    Strand  Ensembl Gene ID
    # 2L      7529    9484    CG11023 1       FBgn0031208
    # 2L      9839    21372   l(2)gl  -1      FBgn0002121
    # 2L      21919   25151   Ir21a   -1      FBgn0031209
    # 2L      25402   59242   Cda5    -1      FBgn0051973

    # key: chrom_start_end, value: [gene name, Ensembl gene id ]
    # one dict for genes on plus strand and another one for the minus strand genes.
    annot_dict_plus = {}
    annot_dict_minus = {}

    for line in open(annotfile):
        if line.startswith('#'):
            continue

        t = line.strip().split('\t')
        if int(t[4]) > 0:
            annot_dict_plus['%s,%s,%s' %(t[0], t[1], t[2])] = [t[3], t[5]]
        else:
            annot_dict_minus['%s,%s,%s' %(t[0], t[1], t[2])] = [t[3], t[5]]


    # find for each peak whether it is in a promoter. Extend promoters by +-r bps.
    # chr1    231473620       231473733       reg1000000.p1   20.250  +
    # chr19   54605995        54606123        reg1000001.p1   19.279  +
    # chr7    112580129       112580232       reg1000008.p2   18.527  +


    annotated = 0

    nonannotated = 0
    mto = 0 #more than one annotation

    o = open(outfile, 'w')
    o.write('#chrom\tstart\tstop\tpeakID\tZ-score\tstrand\tgene_name\tgene_ensembl_id\n')

    for line in open(peaks):
        t = line.strip().split()
        chrom = t[0]
        start = int(t[1])
        stop = int(t[2])
        strand = t[5]

        found = []
        # First search the genes on the positive strand of the annotation file. Do this separately, because I have to extend the gene coordinates differently, depending on whether there on plus or minus strand.
        for p in annot_dict_plus:
            pt = p.split(',')
            pchrom = pt[0]
            pstart = int(pt[1])

            promoter_start = pstart-2*r
            promoter_stop = pstart + r #allow the peak to overlap the gene
            
            if chrom != 'chr' + pchrom:
                continue

            if start >= promoter_start and stop <= promoter_stop:
                found.append(annot_dict_plus[p])

        # Second go over the annotations of the minus strand
        for p in annot_dict_minus:
            pt = p.split(',')
            pchrom = pt[0]
            pstop = int(pt[2])

            promoter_start = pstop - r  #allow the peak to overlap the gene
            promoter_stop = pstop + 2*r

            if chrom != pchrom:
                continue

            if start >= promoter_start and stop <= promoter_stop:
                found.append(annot_dict_minus[p])

        # Now check how many annotations there were and write it to output file
        if len(found) > 0:
            annotated += 1
            o.write('%s\t%s\t%s\t%s\t%s\t%s\t' %(chrom, start, stop, t[3], t[4], t[5]))
            for i in found:
                o.write('%s ' %i[0])
            for i in found:
                for j in i[1]:
                    o.write('\t%s' %j)
            o.write('\n')

            if len(found)>1:
                mto +=1
        else:
            nonannotated += 1

    o.close()


    print annotated
    print nonannotated
    print mto

    return annotated, nonannotated, mto


def execute(cf):

    peaks = cf.get_input("peaks")
    outfile = cf.get_output("peakAnnotations")
    log_file = cf.get_output("log_file")

    genome = cf.get_parameter("genome", "string")
    annotfile = cf.get_parameter("annotationFile", "string")


    if genome == 'hg19' or genome == 'mm9':
        annotated, nonannotated, mto = annotate_human_mouse(peaks, annotfile, outfile)
    elif genome == 'dm3':
        annotated, nonannotated, mto = annotate_drosophila(peaks, annotfile, outfile)

    f = open(log_file, 'w')
    f.write('%i peaks were annotated out of %i peaks.\n%i peaks were annotated to more than one promoter.\n' %(annotated, annotated+nonannotated, mto))
    f.close()

    return 0


component_skeleton.main.main(execute)
