#!/usr/bin/env python
import component_skeleton.main
import subprocess
import os, re
from string import *
from pylab import *

def reverse_and_rename(wmFile, logo_rev):
    """ Read swiss regulon style WM file """

    M = False
    ok = False
    lines = wmFile.read().splitlines()
    for line in lines:
        if not (line[0:2]=="//" or line==""):
            if ok:
                elm = np.reshape(np.array(map(float,line.split()[1:5])),(1,4))
                if M is False:
                    M = elm
                else:
                    M = np.vstack((M,elm))
            else:
                if line[0:2]=='NA':
                    na,name = line.split()
                elif line[0:2]=='PO' or line[0:2]=='P0': # WM start after PO/0 line
                    ok = True
                else:
                    pass


    # build the reverse complement matrix
    Mrev = np.vstack((M[::-1,3],M[::-1,2],M[::-1,1],M[::-1,0])).T

    tmpwm = os.path.join(os.path.split(logo_rev)[0], 'tmp_wm_rev')
    wmname = 'Logo_rev'

    o = open(tmpwm, 'w')
    o.write('//\nNA %s\nP0\tA\tC\tG\tT\n' %(wmname))
    for j, i in enumerate(Mrev):
        o.write('\t'.join([str(j+1).zfill(2)] + map(str,i)))
        o.write('\n')

    o.write('//\n')

    return tmpwm


def renameWM(wm, logo):
    """
    Put the output logo path as name inside the WM file
    """

    tmpwm = os.path.join(os.path.split(logo)[0], 'tmp_wm')

    #wmname = os.path.split(re.sub('.pdf$', '', logo))[1]
    wmname = 'Logo'

    wmlines = open(wm).readlines()

    #NA /import/bc2/home/nimwegen/GROUP/ENCODE.ChIPseq/GM12878_myers_hudsonAlpha/BG_8_third/OUTPUT/ETS1_FgBg-runmotevoPG4_2/Logo
    wmlines[1] = 'NA %s\n' %wmname

    o = open(tmpwm, 'w')
    for i in wmlines:
        o.write(i)

    return tmpwm


def createLogo(mylogo_path, WM):

    # so that Logo gets printed to the right place with a nicer name...
    pwd = os.getcwd()
    os.chdir(os.path.split(WM)[0])

    proc = subprocess.Popen('%s -n -a -c -p -Y -F PDF -f %s' %(mylogo_path, WM),
                             stdout=subprocess.PIPE,
                             stderr= subprocess.PIPE,
                             shell=True
                            )

    stdout_value, stderr_value = proc.communicate()
    print stdout_value
    print stderr_value

    os.chdir(pwd)


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
    tmpwm_rev = reverse_and_rename(open(wm), os.path.join(os.path.split(logo)[0], 'Logo_rev.pdf'))

    createLogo(mylogo_path, tmpwm)
    createLogo(mylogo_path, tmpwm_rev)

    ##extract sequence likelihood and AUC for this WM from log files of ComputeLikelihood and WMQuality components
    wmnameInit = wm

    ldict = {}
    for line in open(llog):
        t = line.strip().split()
        ldict[t[0]] = t[1]

    lvals = ldict[wmnameInit]


    aucdict = {}
    for line in open(auclog):
        t = line.strip().split()
        aucdict[t[1]] = [t[0], t[2]]

    aucval = aucdict[wmnameInit][1]
    old_wmname = aucdict[wmnameInit][0].rstrip(':')

    os.system('cp %s %s' %(os.path.join( os.path.split(auclog)[0], 'intermediate', old_wmname + '.pdf'), AUC_plot))
    os.system('cp %s %s' %(os.path.join( os.path.split(auclog)[0], 'intermediate', old_wmname + '.png'), AUC_plot.rstrip('.pdf')+'.png'))


    o = open(log, 'w')
    text = '\n'.join(['- Motif name: %s' %os.path.split(wmnameInit)[1],
                      '- Enrichment score: %s' %(lvals),
                      '- Area under precision recall curve: %s' %(aucval)])
    o.write(text)
    o.close()


    return 0


component_skeleton.main.main(execute)
                                                                 
