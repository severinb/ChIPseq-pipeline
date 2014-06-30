#! /import/bc2/home/nimwegen/GROUP/local/unladen/bin/python

import component_skeleton.main
from string import *
import yaml
import random
import sys, os
import subprocess
from datetime import datetime


def run(cmd):
    """
    given a command, it runs it
    """
  
    proc = subprocess.Popen (cmd,                      
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                shell=True
                         )
    stdout_value, stderr_value = proc.communicate()
    print stdout_value
    print stderr_value

    if proc.poll() > 0:
        sys.stderr.write ( "\nError\n" )
        print '\tstderr:', repr(stderr_value.rstrip())
        return False
    else:
        return True


def execute(cf):
    infiles_string = cf.get_parameter("infiles_string", "string")
    out_dir = cf.get_output("out_dir")
    intermediate = cf.get_output("intermediate")
    logfile = cf.get_output("ReadCount_log")
    winsize = cf.get_parameter("winsize", "int")
    stepsize = cf.get_parameter("stepsize", "int")
    chrominfo = cf.get_parameter("chrominfo", "string")
    verbose = cf.get_parameter("verbose", "string")
    weight = cf.get_parameter("weight", "string")
    pythonPATH = cf.get_parameter("pythonPATH", "string")

    T1 = datetime.now()

    os.mkdir(intermediate)    


    params = {'INPUT FILE': infiles_string.split(),
              'WINDOW': winsize,
              'STEP': stepsize,
              'RESULT DIR': out_dir,
              'CHROMOSOME INFO': chrominfo
        }


    try:
        f = open(intermediate + '/params.yaml', "w")
        yaml.dump(params, f)
        f.flush()
    except IOError,e:
        print e
        print 'An error occured while making parameter file!'
        sys.exit(-1)


    def option(opt, val):
        if val:
            return '-%s' % opt
        else:
            return ''

 # '/import/bc2/home/nimwegen/GROUP/hseq_pipeline/saeed/ChIP_read_counter/read_count.py'
    cmd = ' '.join([pythonPATH,
                    'read_count.py',
                    '-y %s -w' % (intermediate + '/params.yaml'),
                    option('w',weight),
                    option('v',verbose)
                    ])
    run(cmd)

    T2 = datetime.now()
    time = 'Running time for read counter: ' + str(T2-T1) + '\n'
    lf = open(logfile, 'w')
    lf.write(time)
    lf.close

    return 0

component_skeleton.main.main(execute)
    
   
