#!/usr/bin/env python
import component_skeleton.main
import subprocess
import re, os
from string import *

def execute(cf):
    FMIid = cf.get_parameter("FMIid", "string")
    out_log = cf.get_output("numberMappableReads_log")
    annoType = cf.get_parameter("annoType", "string")
    FMIpath = cf.get_parameter("FMIpath", "string")
    perlPATH = cf.get_parameter("perlPATH", "string")


    """
    proc = subprocess.Popen (perlPATH + ' numberMappableReads.pl %s %s %s %s %s' % (FMIid, out_log, annoType, perlPATH, FMIpath),
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

    """

    #to be sure to use the right perl
    os.environ['PATH'] =  perlPATH + ':' + os.environ['PATH']

    command = ' '.join([perlPATH,
                        FMIpath+'/soft/extractData.pl -f',
                        FMIid,
                        annoType,
                        'genome |',
                        perlPATH,
                        FMIpath+'/soft.bc2/frag2totalGenomic.pl - >',
                        out_log])
    
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


    return 0


component_skeleton.main.main(execute)
