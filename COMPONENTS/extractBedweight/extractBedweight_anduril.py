#!/usr/bin/env python

import subprocess
import re, os, sys
from string import *
from datetime import datetime
import component_skeleton.main

def execute(cf):
    out_dir = cf.get_output("out_dir")
    logfile = cf.get_output("extractBedweight_log")
    FMIid = cf.get_parameter("FMIid", "string")
    annoType = cf.get_parameter("annoType","string")
    mismatches = cf.get_parameter("mismatches", "int")
    FMIpath = cf.get_parameter("FMIpath", "string")
    perlPATH = cf.get_parameter("perlPATH", "string")

    T1 = datetime.now()

    #to be sure to use the right perl
    os.environ['PATH'] =  perlPATH + ':' + os.environ['PATH']

    os.system("mkdir %s" % (out_dir))

    outfile = os.path.join(out_dir,''.join([FMIid, '.', annoType, '.m', str(mismatches), '.bedweight.gz']))
    tmpoutfile = os.path.join(out_dir, 'tmp')

    ##Extraction of Bedweight
    command = ' '.join([perlPATH, FMIpath+'/soft/extractData.pl',
                        FMIid,
                        annoType,
                        'genome -f |',
                        FMIpath+'/soft.bc2/frag2bedweight.pl',
                        '- -m %s ' %mismatches,
                        ' > %s' %tmpoutfile
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
    print proc.poll()
    if proc.poll() > 0:
        print '\tstderr:', repr(stderr_value.rstrip())
        return -1



    command2 = ' '.join(['%s/soft/pioSortBed6 --input-file %s' %(FMIpath, tmpoutfile),
                        '| gzip -9 > %s' %outfile
                        ])

    T2 = datetime.now()

    ##Sorting of Bedweight
    print '\nStart Sorting\n'
    proc = subprocess.Popen (command2,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True
                             )

    stdout_value, stderr_value = proc.communicate()
    print stdout_value
    print stderr_value
    print proc.poll()
    if proc.poll() > 0:
        print '\tstderr:', repr(stderr_value.rstrip())
        return -1

    #remove tmp file
    os.system('rm %s' %tmpoutfile)

   
    T3 = datetime.now()
    time = 'Running time for:\n\t-bedweight extraction: %s\n\t-bedweight sorting: %s\n\t-overall: %s\n' %(T2-T1, T3-T2, T3-T1)
    lf = open(logfile, 'a')
    lf.write(time)
    lf.close	


    return 0

component_skeleton.main.main(execute)
