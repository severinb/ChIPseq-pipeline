#!/usr/bin/env python
import component_skeleton.main
import subprocess
import os, re
from datetime import datetime
import gzip
from pylab import savefig
import sys 

def remove_multi_mappers(in_file, intermediate):

        #check whether input file is gzipped:
        try:
                gop = gzip.open(in_file)
                gop.readline()
                gzipped = True
                gop.close()
        except IOError:
                gop.close()
                gzipped = False

        outfile = os.path.join(intermediate, 'input_without_multi_mappers.bedweight')
        fw = open(outfile, 'w')

     	#filters out multi-mappers: only write lines that have weight of 1 or higher
        if gzipped:
                for line in gzip.open(in_file):
                        items=line.strip().split()
                        if float(items[4])>=1.0:
                                fw.write("%s\n"%line.strip())
        else:
                for line in open(in_file):
                        items=line.strip().split()
                        if float(items[4])>=1.0:
                                fw.write("%s\n"%line.strip())

        fw.close()

        return outfile

def sample_infile(in_file, intermediate):
    """
    This function makes the input file smaller by taking a subsample of the input reads.
    Reads are selected randomly, so that we end up with  a certain number of reads.
    In case the input file has less reads than final_size, we don't change the initial input file.
    """

    final_size = 5000000

    proc1 = subprocess.Popen('less %s | wc -l' %in_file, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
    r = proc1.communicate()

    n_lines = float(r[0].strip())
    print 'Number of reads after removing multimappers: ', int(n_lines)

    ratio = final_size/n_lines
    if ratio >= 0.9: #make it 0.9 because it also takes time to make the file smaller, so there is no use in doing it unless it really makes the file smaller
        return in_file

    outfile = os.path.join(intermediate, 'smaller_input_without_multi_mappers.bedweight')

    proc2 = subprocess.Popen("less %s | awk 'BEGIN {srand()} !/^$/ { if (rand() <= %f) print $0}' - > %s" %(in_file, ratio, outfile), stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
    r = proc2.communicate()

    return outfile


def shift_file(in_file, given_fraglen, out_dir):

        try:
                gop = gzip.open(in_file)
                gop.readline()
                gzipped = True
                gop.close()
        except IOError:
                gop.close()
                gzipped = False


        tmp_file = os.path.join(out_dir, 'tmp.shifted.bedweight')
        res_shifted = open(tmp_file, 'w')

        def shifting(pos, val, strand):
                if strand == '+':
                        return str(int(pos)+(int(val)/2))
                else:
                        return str(int(pos)-(int(val)/2))

        if gzipped:
                for a_line in gzip.open(in_file):
                        tmp = a_line.rstrip().split()
                        if int(shifting(tmp[1], given_fraglen, tmp[5])) > 0:
                                res_shifted.write('\t'.join([tmp[0],
                                                             shifting(tmp[1], given_fraglen, tmp[5]),
                                                             shifting(tmp[2], given_fraglen, tmp[5]),
                                                             tmp[3],
                                                             tmp[4],
                                                             tmp[5],
                                                             '\n',
                                                             ]))
        else:
                for a_line in open(in_file):
                        tmp = a_line.rstrip().split()
                        if int(shifting(tmp[1], given_fraglen, tmp[5])) > 0:
                                res_shifted.write('\t'.join([tmp[0],
                                                             shifting(tmp[1], given_fraglen, tmp[5]),
                                                             shifting(tmp[2], given_fraglen, tmp[5]),
                                                             tmp[3],
                                                             tmp[4],
                                                             tmp[5],
                                                             '\n',
                                                             ]))

        res_shifted.close()


        return tmp_file


def run(command):

    proc = subprocess.Popen (command,
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

    return 0


def execute(cf):
    in_dir = cf.get_input("in_dir")    #Optional, but in pipeline this is used by default
    in_file = cf.get_input("in_file")  #Optional. Usually this isn't used. But if this one is defined it used (instead of in_dir). I added this in-port to be able to run pipeline from here with an already existing bedweight file 
    out_dir = cf.get_output("out_dir")
    intermediate_out = cf.get_output("intermediate")
    log_file = cf.get_output("FragmentLength_log")
    plot_file = cf.get_output("FragmentLength_plot")
    repeatPath = cf.get_parameter("repeatPath", "string")
    pythonPATH = cf.get_parameter("pythonPATH", "string")
    perlPATH = cf.get_parameter("perlPATH", "string")
    FMIpath = cf.get_parameter("FMIpath", "string")
    given_fraglen = cf.get_parameter("fraglen", "float")

    T1 = datetime.now()
  
    os.mkdir(out_dir)
    #intermediate_out = os.path.join(os.path.split(out_dir)[0], 'intermediate')
    os.mkdir(intermediate_out)


    #If in_file is not defined (as anduril in-port), then use in_dir to read in file from there. This is the default way.
    if not in_file:
        in_file = os.path.join(in_dir, os.listdir(in_dir)[0])

    tmp_in_file = remove_multi_mappers(in_file, intermediate_out)
    #tmp_in_file = sample_infile(tmp_in_file1, intermediate_out)

    if given_fraglen < 0:
        print "find shift\n"
        #run shift finder
        FraglenCommand = pythonPATH + ' findmaxshift.py %s %s %s %s %s %s %s' % (tmp_in_file, out_dir, intermediate_out, log_file, plot_file, perlPATH, repeatPath)
        r1 = run(FraglenCommand)
        if r1 < 0:
            return -1
        # extract fragment length
        with open(tmp_in_file + '.res') as resfile:
            fraglen = float(resfile.readline().strip().split()[1])

        print 'Found fraglen: ', fraglen

        print 'Shift input file'
        tmp_file = shift_file(in_file, fraglen, out_dir)


    else: #shift reads with given fragment length
        print 'Use given fragment length %f' %given_fraglen

        # produce dummy plot, so that nobody complains.
        savefig(plot_file)

        # write given fraglen to output file so that following components will find it
        resout = open(os.path.join(intermediate_out, 'fraglen.res'), 'w')
        resout.write("FRAG_LENGTH: %f\n" %given_fraglen)
        resout.close()

        print 'Shift input file'
        tmp_file = shift_file(in_file, given_fraglen, out_dir)



    out_file = os.path.join(out_dir, 'outfile.gz') 


    T2 = datetime.now()
    print "\nSort File\n"
    #sort temoprary outfile and remove it
    SortCommand = '%s/soft/pioSortBed9 -s5 --input-file %s | gzip -9 > %s' %(FMIpath, tmp_file, out_file) #sort by 5' ends
    r2 = run(SortCommand)
    if r2 < 0:
        return -1

    T3 = datetime.now()

    os.system('rm %s' %' '.join([tmp_file, tmp_in_file]))

    T4 = datetime.now()
    time = '\nRunning time for fragment length finder:\n\t-Finding Shift and Shifting: %s\n\t-Sorting and gzipping bedweight: %s\n\t-Overall: %s\n' %(T2-T1, T3-T2, T4-T1)
    lf = open(log_file, 'a')
    lf.write(time)
    lf.close

    return 0

component_skeleton.main.main(execute)
