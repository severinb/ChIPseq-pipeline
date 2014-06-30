#!/usr/bin/env python
import component_skeleton.main
import subprocess
import re, os, sys
from string import *
from datetime import datetime

def execute(cf):
    #in_dir = cf.get_input("in_dir")
    desc_dir = cf.get_input("desc_dir")
    out_dir = cf.get_output("out_dir")
    logfile = cf.get_output("extractWig_log")
    FMIid = cf.get_parameter("FMIid", "string")
    annoType = cf.get_parameter("annoType","string")
    mismatches = cf.get_parameter("mismatches", "int")
    width = cf.get_parameter("width", "int")
    FMIpath = cf.get_parameter("FMIpath", "string")
    perlPATH = cf.get_parameter("perlPATH", "string")

    T1 = datetime.now()

    #to be sure to use the right perl
    os.environ['PATH'] =  perlPATH + ':' + os.environ['PATH']

    os.system("mkdir %s" % (out_dir))

    outfile = os.path.join(out_dir, ''.join([FMIid, '.', annoType, '.w', str(width), '.m', str(mismatches), '.wig']))

    command = ' '.join([perlPATH,
                        FMIpath+'/soft/extractData.pl',
                        FMIid,
                        annoType,
                        'genome -f | ',
                        FMIpath+'/soft/frag2wig.pl',
                        '- -m %d' %mismatches,
                        '-n %s' %FMIid,
                        '-w %d' %width,
                        '| gzip -9 > %s.gz' %outfile 
                        ])

    print 'Extracting %s\n' %FMIid
    proc = subprocess.Popen (command,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True
                         )
    stdout_value, stderr_value = proc.communicate()
    print stdout_value
    print stderr_value

    if proc.poll() > 0:
	    print '\tstderr:', repr(stderr_value.rstrip())
	    return -1		

    T2 = datetime.now()
    time = 'Running time for wig extraction: ' + str(T2-T1) + '\n'
    lf = open(logfile, 'a')
    lf.write(time)
    lf.close

    return 0


component_skeleton.main.main(execute)
