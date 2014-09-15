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
    FMI_output_dir = cf.get_parameter("FMI_output_dir", "string")

    T1 = datetime.now()

    #to be sure to use the right perl
    os.environ['PATH'] =  perlPATH + ':' + os.environ['PATH']

    os.system("mkdir %s" % (out_dir))

    outfile = os.path.join(out_dir,''.join([FMIid, '.', annoType, '.m', str(mismatches), '.bedweight.gz']))
    tmpoutfile = os.path.join(out_dir, 'tmp')


    command1 = ' '.join([perlPATH, FMIpath+'/soft/extractData.pl',
                         FMIid,
                         annoType,
                         'genome',
                         FMI_output_dir,
                         '-f'])

    command2 = ' ' .join([FMIpath+'/soft.bc2/frag2bedweight.pl',
                          '- -m %s ' %mismatches])
                          #' > %s' %tmpoutfile])

    command3 = ' ' .join([FMIpath+'/soft/pioSortBed9',
                          '-s5', #sort by 5' end
                          '--input-file -' #read from stdin
                          ])

    command4 = './collapse_bedweight.py'

    command5 = 'gzip -9 > %s' %outfile


    print 'Extracting %s\n' %FMIid
 
    p1 = subprocess.Popen(command1,
                          stdout=subprocess.PIPE,
                          shell=True)

    p2 = subprocess.Popen(command2,
                          stdin=p1.stdout,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          shell=True)

    p3 = subprocess.Popen(command3,
                          stdin=p2.stdout,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          shell=True)

    p4 = subprocess.Popen(command4,
                          stdin=p3.stdout,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          shell=True)

    p5 = subprocess.Popen(command5,
                          stdin=p4.stdout,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          shell=True)

    p1.stdout.close()
    p2.stdout.close()
    p3.stdout.close()
    p4.stdout.close()
    stdout_value, stderr_value = p5.communicate()

    print 'p5 returncode', p5.returncode, p5.poll()
    print 'p4 returncode', p4.returncode, p4.poll()
    print 'p3 returncode', p3.returncode, p3.poll()
    print 'p2 returncode', p2.returncode, p2.poll()
    print 'p1 returncode', p1.returncode, p1.poll()

    print stdout_value
    print stderr_value

    if p5.poll() > 0:
        print 'p5 failed'
        return -1
    if p4.poll() > 0:
        print 'p4 failed'
        return -1
    if p3.poll() > 0:
        print 'p3 failed'
        return -1
    if p2.poll() > 0:
        print 'p2 failed'
        return -1
    if p1.poll() > 0:
        print 'p1 failed'
        return -1


    # command2 = ' '.join(['%s/soft/pioSortBed6 --input-file %s' %(FMIpath, tmpoutfile),
    #                      '| gzip -9 > %s' %outfile
    #                      ])


    # ##Sorting of Bedweight
    # print '\nStart Sorting\n'
    # proc = subprocess.Popen (command2,
    #                          stdout=subprocess.PIPE,
    #                          stderr=subprocess.PIPE,
    #                          shell=True
    #                          )

    # stdout_value, stderr_value = proc.communicate()
    # print stdout_value
    # print stderr_value
    # if proc.poll() > 0:
    #     print '\tstderr:', repr(stderr_value.rstrip())
    #     return -1

    # #remove tmp file
    # os.system('rm %s' %tmpoutfile)

   
    T2 = datetime.now()
    time = 'Running time %s\n' %(T2-T1)
    lf = open(logfile, 'a')
    lf.write(time)
    lf.close	


    return 0

component_skeleton.main.main(execute)
