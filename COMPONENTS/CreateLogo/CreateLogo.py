#!/usr/bin/env python
import component_skeleton.main
import subprocess
import os, re
from string import *
from pylab import *


def renameWM(wm, logo):
    """
    Put the output logo path as name inside the WM file
    """

    tmpwm = os.path.join(os.path.split(logo)[0], 'tmp_wm')

    wmname = re.sub('.pdf$', '', logo)

    wmlines = open(wm).readlines()

    #NA /import/bc2/home/nimwegen/GROUP/ENCODE.ChIPseq/GM12878_myers_hudsonAlpha/BG_8_third/OUTPUT/ETS1_FgBg-runmotevoPG4_2/Logo
    wmlines[1] = 'NA %s\n' %wmname

    o = open(tmpwm, 'w')
    for i in wmlines:
        o.write(i)

    return tmpwm


def createLogo(mylogo_path, WM):
    
    proc = subprocess.Popen('%s -n -a -c -p -Y -F PDF -f %s' %(mylogo_path, WM),
                             stdout=subprocess.PIPE,
                             stderr= subprocess.PIPE,
                             shell=True
                            )

    stdout_value, stderr_value = proc.communicate()
    print stdout_value
    print stderr_value
    

def execute(cf):

    ##Ports and parameters
    wm = cf.get_input("WM")
    llog = cf.get_input("llog")
    auclog = cf.get_input("auclog")

    logo = cf.get_output("Logo")
    log = cf.get_output("log_file")
    AUC_plot = cf.get_output("sens_spec")

    mylogo_path = cf.get_parameter("mylogo_path", "string")


    tmpwm = renameWM(wm, logo)
    createLogo(mylogo_path, tmpwm)


    ##extract sequence likelihood and AUC for this WM from log files of ComputeLikelihood and WMQuality components
    wmnameInit = wm

    ldict = {}
    for line in open(llog):
        t = line.strip().split()
        ldict[t[1]] = t[2:]

    lvals = ldict[wmnameInit]


    aucdict = {}
    for line in open(auclog):
        t = line.strip().split()
        aucdict[t[1]] = [t[0], t[2]]

    aucval = aucdict[wmnameInit][1]
    old_wmname = aucdict[wmnameInit][0].rstrip(':')

    os.system('cp %s %s' %(os.path.join( os.path.split(auclog)[0], 'intermediate', old_wmname + '.pdf'), AUC_plot))


    o = open(log, 'w')
    text = '\n'.join(['- Original path: %s' %wmnameInit,
                      '- Sequence log-likelihood of %s at an optimal prior of %s' %(lvals[0], lvals[1]),
                      '- Area under precision recall curve: %s' %(aucval)])
    o.write(text)
    o.close()


    return 0


component_skeleton.main.main(execute)
                                                                 
