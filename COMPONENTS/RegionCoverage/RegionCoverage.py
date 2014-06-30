#!/usr/bin/env python
import sys
from string import *
import subprocess
import gzip
import os, re
import component_skeleton.main
from datetime import datetime

def regions2Length(peaksfile, length):
    """Makes all peaks to equal length. Adds bp symmetrically to each peak.
    """

    f = open(peaksfile, 'r')
    peakseqlen = os.path.join('intermediate', outfile + '_eqlenpeaks')
    o = open(peakseqlen ,'w')
    eof = False

    while not eof:
        line = f.readline()
        if line:
            peak = line.split()
            start = int(peak[1])
            end = int(peak[2])
            peaklength = end - start
            diff = length - peaklength
            if not diff == 0:
                newstart = start - int(diff/2)
                newend = end + int(diff/2)
                o.write('\t'.join([peak[0], str(newstart), str(newend), '+']) + '\n')
            else:
                o.write('\t'.join([peak[0], str(start), str(end), '+']) + '\n')
        else:
            eof = True
            break

    f.close()
    o.close()

    return peakseqlen


def toFragLength(readfile, ind, fraglen, intermediate):
    """This function extends reads (about 35 bp) to actual fragment length (calculated by Nick's FragmentLength component).
    """

    try:
        f = gzip.open(readfile)
        f.readline()
        f.close()
        f = gzip.open(readfile)
    except IOError:
        f.close()
        f = open(readfile)

    outfraglenfile = os.path.join(intermediate, 'reads2fraglength.shifted.bedweight' + str(ind))
    o = open(outfraglenfile, 'w')
    eof = False

    shift = int(fraglen/2)

    for line in f:
        read = line.strip().split()
        fraglen_read = []
        fraglen_read.append(read[0]) #chromosome name
        #center = int((int(read[1])+int(read[2]))/2)
        #start = center-shift
        #end = center+shift
        strand = read[5]
        if strand == '+':
            start = int(read[1]) - shift
            end  = int(read[1]) + shift
        elif strand == '-':
            start = int(read[2]) - shift
            end = int(read[2]) + shift
        if start < 0:                     #use for negative starts 0
            fraglen_read.append(str(0))
            if end <= 0:                    #use for negative or 0 end positions 1
                fraglen_read.append(str(1))
            else:
                fraglen_read.append(str(end))
        else:
            fraglen_read.append(str(start))
            fraglen_read.append(str(end))
        fraglen_read.append(read[3]) #description something
        fraglen_read.append(read[4]) #weight
        fraglen_read.append(read[5]) #strand
        o.write('\t'.join(fraglen_read) + '\n')

    f.close()
    o.close()
    
    print '\nDone: %s' %outfraglenfile
    return outfraglenfile


def intersectBed(reads, regions, intermediate, ind, BedToolsPath):
    """Runs intersectBed from bedtools. 
    """

    outfile = os.path.join(intermediate, 'RegionReadIntersection.' + str(ind))
    #'/import/bc2/home/nimwegen/GROUP/bedtools/intersectBed -a %s -b %s -wb -wa > %s' %(regions, reads, outfile),
    #proc = subprocess.Popen('%s intersect -a %s -b %s -wb -wa > %s' %(os.path.join(BedToolsPath, 'bedtools'), regions, reads, outfile),
    proc = subprocess.Popen('%s intersect -a %s -b %s -wb -wa > %s' %(BedToolsPath, regions, reads, outfile),
                             stdout=subprocess.PIPE,
                             stderr= subprocess.PIPE,
                             shell=True
                            )
    stdout_value, stderr_value = proc.communicate()
    print stdout_value
    print stderr_value

    if proc.poll() > 0:
        print '\tstderr:', repr(stderr_value.rstrip())
        return -1
    else:
        print '\nRead Region intersection done.'
        return (0, outfile)



def weightedCoverage(intersect_replicate, intermediate, ind):
    """Takes file with peaks and intersecting reads and computes weighted bp resolution coverage.
    """

    ###put all reads and regions into a dictionary. keys are regions, values are list of intersecting reads.
    #format of intersect file:
    #chr17   57863500        57864750        +       chr17   57863402        57863538        sq3524693       1       +
    #chr17   57863500        57864750        +       chr17   57863416        57863552        sq1764563       1       +
    #second file has to be exactly in this format. The first file can have further columns

    region_reads_dict = {}  #dictionary: keys(region): values([(read_start, read_end, read_weight)])
    f = open(intersect_replicate, 'r')

    for line in f:
        tmp = line.strip().split()
        if tmp[-5] < tmp[1]:   #Put overlapping reads to non overlapping ones. This means that reads that are lower or higher than the region get truncated ends to region start or end posititon respectively.
            tmp[-5] = tmp[1]
        if tmp[-4] > tmp[2]:
            tmp[-4] = tmp[2]
        region_position = ','.join(tmp[:4])
        read_position = tmp[-5:-3]  #don't take chromsome, because reads always match region (intersectBed did it).
        try:
            region_reads_dict[region_position].append((read_position[0], read_position[1], tmp[-2]))
        except KeyError:
            region_reads_dict[region_position] = [(read_position[0], read_position[1], tmp[-2])]            

    f.close()

    ###convert list of intersecting reads into bp resolved coverage: Prints them to a file of the form:
    ###region_position	bp	coverage
    cov_replicate = os.path.join(intermediate, 'regionCoverage.' + str(ind))
    o = open(cov_replicate, 'w')

    for region in region_reads_dict.keys():
        reg_tmp = region.split(',')
        reg_chrom = reg_tmp[0]
        reg_start = int(reg_tmp[1])
        reg_end = int(reg_tmp[2])
        region_dict = {}             #keys (1:reg_end - reg_start +1): values (summed weights of reads) or (bp: coverage)
        for i in xrange(1, reg_end - reg_start + 1):
            region_dict[i] = 0

        for read in region_reads_dict[region]:
            for bp in xrange(max(1,int(read[0])-reg_start+1), max(1,int(read[1])-reg_start+1)):   #26.7.12 13:47 added max functions
                try:  #23.10.12 added this exception for reads that overlap the region. IntersectBed reports all reads that overlap at least by one bp. I think this also takes care of this max(1, read[0] - reg_start) etc stuff. 
                    region_dict[bp] += float(read[2])
                except KeyError:  
                    pass

        for bp in region_dict:
            line = '\t'.join([reg_chrom, 
                             str(reg_start), 
                             str(reg_end), 
                             reg_tmp[-1], 
                             str(bp) ,
                             str(region_dict[bp]), '\n'])
            o.write(line)

    o.close()
    
    return cov_replicate



def coverageAverage(replicates, intermediate):
    """calculates the average coverage for each bp of the replicates.
    """

    averagefile = os.path.join(intermediate, 'coverageAverage')
    o = open(averagefile, 'w')

    #a dictionary for region: bp : [count1, count2 ..]
    #e.g. {chr1,12000,13000,+ : {1 : [0.3, 4], 2: [1, 2}, ...}
    CovDictDict = {}

    #a line: chr1    228760000       228761000       +       1       2.26259889921
    for rep in replicates:
        for line in open(rep):
            t = line.split()
            region = ','.join(t[:4])
            bp = int(t[4])
            weight = float(t[5])
            try:
                CovDictDict[region][bp].append(weight)
            except KeyError:
                try:
                    CovDictDict[region][bp] = [weight]
                except KeyError:
                    CovDictDict[region] = {}
                    CovDictDict[region][bp] = [weight]

    #now write average
    for region in CovDictDict:
        reglist = region.split(',')
        for bp in CovDictDict[region]:
            weightlist = CovDictDict[region][bp]
            o.write('\t'.join(reglist + [str(bp), str(sum(weightlist)/len(weightlist))]))
            o.write('\n')

    o.close()

    return averagefile


def split_peaks(averagefile, outdir):
    """Stores each peak and it's coverage (average of replicates) into a separate file. (In the end it's a directory with 500 files. One for each of the top 500 meged peaks)
       This 'peaks' directory is then the input for the flexible window peak refiner.
    """

    os.mkdir(outdir)

    #Alternative splitting, maybe faster:
    currchr = ''
    currstart = 0
    currend = 0
    currhandle = None
    for line in open(averagefile):
        reg = line.split()
        if reg[0] != currchr or int(reg[1]) != currstart or int(reg[2]) != currend:
            if currhandle:
                currhandle.close()
            currhandle = open(os.path.join(outdir, '_'.join(reg[:3])), 'w')
            currhandle.write(line)
            currchr = reg[0]
            currstart = int(reg[1])
            currend = int(reg[2])
        else:
            currhandle.write(line)



def execute(cf):
    """This code should now be adapted for several replicates.
       FOR ANDURIL INCORPORATION: put all reads to fragment length already in composite components and give 
       the output files as a string (like for the non fragment length reads). String can be given as parameter
       or input for THIS component. (component name like, PeakCoverage or regionCoverage).
       How it computes coverage:
       -take reads from all replicates and put them to fragment length
       -pull out reads that intersect regions (intersectBed bedtools)
       -set overlapping reads positions to start and end position of region
       -sum up read weights over each region position (bp resolution)
    """

    read_filesstr = cf.get_parameter("read_files", "string")
    regions_file = cf.get_input("regions")
    outdir = cf.get_output("out_dir")
    intermediate = cf.get_output("intermediate")
    regionlength = cf.get_parameter("regionlength", "int")
    BedToolsPath = cf.get_parameter("BedToolsPath", "string")
    toFraglen = cf.get_parameter("toFraglength", "boolean")
    logfile = cf.get_output("log_file")

    T1 = datetime.now()

    os.mkdir(intermediate)

    read_fileslist = read_filesstr.split()

    if regionlength > 0:
        #activate next line to make all peaks to length "length"
        regions_file = regions2Length(regions_file, regionlength) #does set length of all peaks to constant length

    T2 = datetime.now()

    replog = [] #list of [replicatenumber, fraglen, fraglenTime, intersectBedTime, WeightTime, replicatenumber2, ....]

    ##process all replicates
    i = 0
    replicates = []
    for rf in read_fileslist:

        TL1 = datetime.now()

        i = i+1

        if toFraglen:
            ##First get Fragment length from the log of FragmentLength component
            fraglenRoot = os.path.join(os.path.split(rf)[0], '..', 'intermediate')
            interdir = os.listdir(fraglenRoot)
            for f in interdir:
                m = re.match('\S*.res$', os.path.join(fraglenRoot, f))
                if m:
                    fraglenF = m.group()
            fraglen = int(float(open(fraglenF).read().strip().split()[-1]))

            print '\nProcessing: %s with fragment length: %s' %(rf,fraglen)
            replog.append('Sample %s has fragment length %i' %(os.path.split(rf)[1], fraglen))
      
            ##Set reads to fragment length
            fraglenfile = toFragLength(rf, i, fraglen, intermediate)
        else:
            replog.append('')
            fraglenfile = rf

        TL2 = datetime.now()



        ##get region intersecting reads
        (number, intersect_replicate) = intersectBed(fraglenfile, regions_file, intermediate, i, BedToolsPath)

        TL3 = datetime.now()

        if number != 0:
            print 'non zero exit of intersectBed'
            sys.exit(1)
        else:
            cov_replicate = weightedCoverage(intersect_replicate, intermediate, i)
            TL4 = datetime.now()

            replicates.append(cov_replicate)
            os.remove(fraglenfile) #remove large intermediate file with shifted fragment length reads.

            replog.append('\tTime for setting reads to fragment length: %s' %str(TL2-TL1))
            replog.append('\tTime for intersecting reads with regions: %s' %str(TL3-TL2))
            replog.append('\tTime for computing the weighted coverage: %s' %str(TL4-TL3))

    T3 = datetime.now()

    ##average coverage over replicates and split single peaks into separate files. Only if there are replicates.
    if len(replicates) != 1:
        averagefile = coverageAverage(replicates, intermediate)

        for covfile in replicates:
            os.system('rm %s' %covfile)
    else:
        averagefile = replicates[0]


    T4 = datetime.now()

    split_peaks(averagefile, outdir)     # splits all the regions from regions_file into separate files (region_file).


    T5 = datetime.now()
    time = '\n'.join(replog + ['','Running Time for:', 
                               '\tSetting regions to length %s: %s' %(regionlength, str(T2-T1)),
                               '\tOverall replicate Processing: %s' %str(T3-T2),
                               '\tAveraging Coverage: %s' %str(T4-T3),
                               '\tSplitting peaks to separate files: %s' %str(T5-T4),
                               '\tOverall: %s' %str(T5-T1)])
    lf = open(logfile, 'w')
    lf.write(time)
    lf.close
                        

    return 0

component_skeleton.main.main(execute)

