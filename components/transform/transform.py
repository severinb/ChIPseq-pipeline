#!/usr/bin/env python
import component_skeleton.main
import subprocess
import sys, re, os
from datetime import datetime
import gzip

def execute(cf):
    in_file = cf.get_input("in_file")
    out_file = cf.get_output("out_file") #this is a unique fasta file
    log_file = cf.get_output("transform_log")
    FMIpath = cf.get_parameter("FMIpath", "string")
    adaptorfile = cf.get_input("adaptor3automatic")
    adaptorString = cf.get_parameter("adaptor3", "string")
    perlPATH = cf.get_parameter("perlPATH", "string")
    sortTMPDIR = cf.get_parameter("sortTMPDIR", "string")

    T1 = datetime.now()

    #to be sure to use the right perl
    os.environ['PATH'] =  perlPATH + ':' + os.environ['PATH']

    #get adaptor: Prioritize manually defined adaptor
    if adaptorString == 'NONE':
        af = open(adaptorfile, 'r')
        adaptor3 = af.read()
        af.close()
    else:
        adaptor3 = adaptorString
    

    # check first whether input file is gzipped. If so, then produce an gunzipped version for transform input
    gzipped = True
    try:
        fin = gzip.open(in_file)
        fin.readline()
    except IOError:
        gzipped = False
    fin.close()

    if gzipped:
        tmpfile = os.path.join(os.path.split(out_file)[0], 'tmpfile')
        o = open(tmpfile, 'w')
        for line in gzip.open(in_file):
            o.write(line)
        o.close()
        in_file = tmpfile
    

    proc1 = subprocess.Popen (perlPATH + ' ' + FMIpath+'/soft/filterAdaptors.pl -L ' + log_file + ' -s -3 ' + adaptor3 + ' -i ' + in_file + ' -F fasta > ' + out_file, 
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              shell=True
                              )
    stdout_value, stderr_value = proc1.communicate()
    print stdout_value
    print stderr_value

    if proc1.poll() > 0:
	    print '\tstderr:', repr(stderr_value.rstrip())
	    return -1

    lf = open(log_file, 'a')
    lf.write('\nTop 10 tags before adaptor removal:\n')
    lf.close()


    #Do the topTags pipe
    #"cat " + os.path.join(in_dir, infilename) + " | grep -v \">\" | sort -T " + sortTMPDIR  + " | uniq -c | sort -T " + sortTMPDIR + " -gr | perl -we '" + 'my $i=0; while(<>) { /\s*(\d+)\s+(\S+)/; print ">seq" . $i++ . "\t" . $1 . "\n" . $2 . "\n"; last if $i==200}' + "' > " + os.path.join(out_dir, infRoot) + ".topTags"
 
    p1 = subprocess.Popen ("cat " + in_file,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           shell=True
                           )

    p2 = subprocess.Popen ("grep -v \">\"",
                           stdin=p1.stdout,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           shell=True
                           )

    p3 = subprocess.Popen ("sort -T " + sortTMPDIR,
                           stdin=p2.stdout,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           shell=True
                           )

    p4 = subprocess.Popen ("uniq -c",
                           stdin=p3.stdout,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           shell=True
                           )

    p5 = subprocess.Popen ("sort -T " + sortTMPDIR + " -gr",
                           stdin=p4.stdout,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           shell=True
                           )

    p6 = subprocess.Popen ("perl -we \'my $i=0; while(<>) { /\\s*(\\d+)\\s+(\\S+)/; print \">seq\" . $i++ . \"\\t\" . $1 . \"\\n\" . $2 . \"\\n\"; last if $i==200}\' > " + out_file + ".topTags" ,
                           stdin=p5.stdout,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           shell=True
                           )

    sout, serr = p6.communicate()
    print sout
    print serr

    proc3 = subprocess.Popen ("head -20 " + out_file + ".topTags >> " + log_file, 
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              shell=True
                              )

    stdout_value, stderr_value = proc3.communicate()
    print stdout_value
    print stderr_value

    if proc3.poll() > 0:
	    print '\tstderr:', repr(stderr_value.rstrip())
	    return -1	

    if gzipped:
        os.system('rm %s' %tmpfile)

		
    T2 = datetime.now()
    time = '\nRunning time for adaptor removal: ' + str(T2-T1) + '\n'
    lf = open(log_file, 'a')
    #lf.write(time)
    lf.close()

    return 0


component_skeleton.main.main(execute)
