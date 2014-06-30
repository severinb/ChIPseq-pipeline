#!/usr/bin/env python
import component_skeleton.main
import subprocess
import os, re
from string import *
from datetime import datetime

def prepareInputAlns(alns, wmlen, tmpfile):
    """
    This function sorts out alignments where the reference species sequence is shorter than WM length when Ns were subtracted.
    """

    o = open(tmpfile, 'w')

    writeAln = False #whether to write the current alignment
    openAln = False #the alignment header that is open at the moment

    tooshort = 0 #count of how many sequences were too short
    goodalns = 0 #count of good alignments

    for line in open(alns):
        if line.startswith('>'):
            if line.startswith('>>'):
                writeAln = False
                openAln = line
                continue
            else:
                if writeAln:
                    o.write(line)
                    continue
                else:
                    continue
        else:
            if writeAln:
                o.write(line)
                continue
            elif openAln:
                if (len(line.strip()) - line.count('N')) <= wmlen:
                    writeAln = False
                    openAln = False
                    tooshort += 1
                    continue
                else:
                    writeAln = True
                    o.write(openAln)
                    o.write(line)
                    goodalns += 1
                    continue
            else:
                continue


    return tooshort, goodalns


def execute(cf):
    alns = cf.get_input("infile") #result file from AlignPeaks spltting from SplitAlignments (typically even half)
    out_file = cf.get_output("out_file")
    tracked_file = cf.get_output("tracked_file")
    report = cf.get_output("report")
    interm = cf.get_output("intermediate")

    weightMat1 = cf.get_output("WeightMatrix1")
    logo1 = cf.get_output("Logo1")
    weightMat2 = cf.get_output("WeightMatrix2")
    logo2 = cf.get_output("Logo2")
    logfile = cf.get_output("log_file") 

    PGpath = cf.get_parameter("PhyloGibbsPATH", "string")
    mylogo_path = cf.get_parameter("mylogo_path", "string")
    markovorder = cf.get_parameter("markovorder", "int")
    wmlen = cf.get_parameter("WindowLength", "int")
    numberTFBS = cf.get_parameter("numberWindows", "int")
    numberMotives = cf.get_parameter("numberColours", "int")
    alignOrder = cf.get_parameter("AlignmentOrder", "int")
    genome = cf.get_parameter("genome", "string")

    T1 = datetime.now()

    os.mkdir(interm)

    ##Sort out 'N's to '-'
    #text = open(alnsIN).read()
    #alns = os.path.join(interm,'alignments_noN')
    #o = open(alns, 'w')
    #o.write(re.sub('N','-',text))
    #o.close()

    tmpfile = os.path.join(interm, 'tmpinfile')
    tooshort, numalns = prepareInputAlns(alns, wmlen, tmpfile)

    #get number of windows:
    if numberTFBS < 0:
        numberTFBS = int(0.7*numalns)

    ##Run PhyloGibbs

    #hardcoded parameters (maybe change this later)
    #wmlen = str(20)
    if genome == "hg19":
        tree = "((((hg19:0.96756,rheMac2:0.94394):0.90646,mm9:0.702855):0.97955,(bosTau6:0.82968,(equCab2:0.89787,canFam2:0.86039):0.98962):0.96777):0.85554,monDom5:0.65318)" #this is the proximity tree for hg19
    if genome == "hg18":
        tree = "((((hg18:0.96756,rheMac2:0.94394):0.90646,mm9:0.702855):0.97955,(bosTau3:0.82968,(equCab1:0.89787,canFam2:0.86039):0.98962):0.96777):0.85554,monDom4:0.65318)" #proximity tree for hg18 (same as hg19 just species renamed!)
    if genome == "mm9":
        tree = "((((hg19:0.96756,rheMac2:0.94394):0.90646,mm9:0.702855):0.97955,(bosTau7:0.82968,(equCab2:0.89787,canFam2:0.86039):0.98962):0.96777):0.85554,monDom5:0.65318)" #this is the proximity tree for hg19
    if genome == "dm3":
        tree = "((((((dm3:0.9427,droSim1:0.9277):0.9598,(droYak2:0.9012,droEre2:0.8985):0.9474):0.8869,droAna3:0.6859):0.9305,dp4:0.8723):0.9408,droWil1:0.5851):0.9801,((droVir3:0.8220,droMoj3:0.7749):0.9296,droGri2:0.7475):0.7139)"

    if alignOrder == 0:  #no alignments are considered
        treeArgument = ""
    else:
        treeArgument = '-L \"' + tree + '\"'
        
    command = ' '.join([PGpath,
                        '-f', tmpfile,
                        '-D', str(alignOrder),
                        #'-L\"%s\"' %tree,
                        treeArgument,
                        '-m', str(wmlen),
                        '-N', str(markovorder),
                        '-o', out_file,
                        '-t', tracked_file,
                        '-y', str(numberTFBS),
                        '-z', str(numberMotives),
                        '>', report])
    
    proc = subprocess.Popen(command,
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


    ##get wm from phylogibbs out_file

    f = open(out_file)
    lines = f.readlines()

    indices = [i for (i,x) in enumerate(lines) if x == '-------- Weight matrix for this motif (absolute base counts)---------\n']
    logopaths = [logo1, logo2]
    wmpaths = [weightMat1, weightMat2]

    for i in [0,1]:
        #sometimes phylogibbs gives just one WM without complaining. For this case I just give a WM with 1 everywhere to not get problems
        try:
            WM = lines[indices[i]+1:indices[i]+int(wmlen)+5]
        except IndexError:
            WM = ['//\n', 'NA\n', 'PO\tA\tC\tG\tT\tcons\tinf\n', '01\t1\t1\t1\t1\tN\t0.001\n', '02\t1\t1\t1\t1\tN\t0.001\n', '03\t1\t1\t1\t1\tN\t0.001\n', '04\t1\t1\t1\t1\tN\t0.001\n', '//\n']

        WM[1] = 'NA ' + re.sub('.pdf','',logopaths[i]) + '\n'

        wm = open(wmpaths[i], 'w')

        for line in WM:
            wm.write(line)

        wm.close()


        ##Produce Logo for WM
    
        proc = subprocess.Popen('%s -n -a -c -p -Y -F PDF -f %s' %(mylogo_path, wmpaths[i]),
                                stdout=subprocess.PIPE,
                                stderr= subprocess.PIPE,
                                shell=True
                                )
    
        stdout_value, stderr_value = proc.communicate()
        print stdout_value
        print stderr_value


    T2 = datetime.now()
    time = 'Running time for PhyloGibbs: ' + str(T2-T1) + '\n'
    text = '\n'.join(['PhyloGibbs parameters:',
                      '\t-number of used input alignments/sequences: %s' %numalns,
                      '\t-number of too short input alignments/sequences: %s' %tooshort,
                      '\t-alignment order: %s' %alignOrder,
                      '\t-window length: %s' %wmlen,
                      '\t-number of windows: %s' %numberTFBS,
                      '\t-number of colours: %s' %numberMotives,
                      '\t-markov order: %s' %markovorder])

    lf = open(logfile, 'w')
    lf.write(time)
    lf.write(text)
    lf.close
                        

    return 0


component_skeleton.main.main(execute)
                                                                 
