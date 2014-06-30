#!/usr/bin/env python
import component_skeleton.main
import subprocess
import re, os, sys
from string import *
from datetime import datetime

def execute(cf):

    out_dir = cf.get_output("out_dir")
    logfile = cf.get_output("log_file")
    FMIid = cf.get_parameter("FMIid", "string")
    annoType = cf.get_parameter("annoType","string")
    maxhits = cf.get_parameter("maxhits", "int")
    FMIpath = cf.get_parameter("FMIpath", "string")
    sortTMPDIR = cf.get_parameter("sortTMPDIR", "string")
    perlPATH = cf.get_parameter("perlPATH", "string")


    T1 = datetime.now()

    #to be sure to use the right perl
    os.environ['PATH'] =  perlPATH + ':' + os.environ['PATH']

    os.system("mkdir %s" %out_dir)
 
    outfile = os.path.join(out_dir,''.join([FMIid, '.', annoType, '.m', str(maxhits), '.bed.gz']))
    tmpoutfile = os.path.join(out_dir, 'tmp')

    command = ' '.join([perlPATH, FMIpath+'/soft/extractData.pl',
                        FMIid,
                        annoType,
                        'genome -f |',
                        FMIpath+'/soft/frag2bed.pl',
                        '- -m %i -s' %maxhits,
                        '> %s' %tmpoutfile
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



    command2 = ' '.join(['%s/soft/pioSortBed6 --input-file %s' %(FMIpath, tmpoutfile),
                        '| gzip -9 > %s' %outfile
                        ])

    T2 = datetime.now()

    ##Sorting of Bed
    print '\nStart Sorting\n'
    proc = subprocess.Popen (command2,
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

    #remove tmp file
    os.system('rm %s' %tmpoutfile)

   
    T3 = datetime.now()
    time = 'Running time for:\n\t-bed extraction: %s\n\t-bed sorting: %s\n\t-overall: %s\n' %(T2-T1, T3-T2, T3-T1)
    lf = open(logfile, 'a')
    lf.write(time)
    lf.close	


    return 0

component_skeleton.main.main(execute)
