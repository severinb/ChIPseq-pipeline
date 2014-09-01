#!/usr/bin/env python
import component_skeleton.main
import subprocess
import os, re
from string import *
import ROC2
import datetime, time
import signal
from pylab import *
from scipy import interpolate

def giveMotevoParamFile(mode, genome, wmlen, inter_dir, tag, minpostwm, wmdiff, aligned, ufemodel_path, ATfreq, GCfreq):
    """
    Returns a parameter file for motevo.
    """

    ##UFE_models from genome_dict are not used anymore
    #UFEmodel_hg19 is UFE model for mammal species
    genome_dict = {}
    genome_dict['hg19'] = ['((((hg19:0.032973,rheMac2:0.057695):0.09821,mm9:0.352605):0.020666,(bosTau6:0.186713,(equCab2:0.107726,canFam2:0.150374):0.010431):0.032764):0.156024,monDom5:0.425899);', '/import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/MotEvo_v1.0/UFEmodels/UFEmodel_hg19']
    genome_dict['hg18'] = ['((((hg18:0.032973,rheMac2:0.057695):0.09821,mm9:0.352605):0.020666,(bosTau3:0.186713,(equCab1:0.107726,canFam2:0.150374):0.010431):0.032764):0.156024,monDom4:0.425899);', '/import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/MotEvo_v1.0/UFEmodels/UFE_mammals']
    #genome_dict['dm3'] = ['((((((dm3:0.059,droSim1:0.075):0.041,(droYak2:0.104,droEre2:0.107):0.054):0.120,droAna3:0.377):0.072,dp4:0.397):0.061,droWil1:0.536):0.020,((droVir3:0.196,droMoj3:0.255):0.073,droGri2:0.291):0.337);', '/import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/MotEvo_v1.0/UFEmodels/UFEmodel_dm3']
    genome_dict['dm3'] = ['((((((dm3:0.059,droSim1:0.075):0.041,(droYak2:0.104,droEre2:0.107):0.054):0.120,droAna3:0.377):0.072,dp4:0.397):0.061,droWil1:0.536):0.020,((droVir3:0.196,droMoj3:0.255):0.073,droGri2:0.291):0.337);', '/import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/MotEvo_v1.0/UFEmodels/dm3UFEparallel/UFEmodel_dm3']
    genome_dict['mm9'] = ['((((hg19:0.032973,rheMac2:0.057695):0.09821,mm9:0.352605):0.020666,(bosTau7:0.186713,(equCab2:0.107726,canFam2:0.150374):0.010431):0.032764):0.156024,monDom5:0.425899);', '/import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/MotEvo_v1.0/UFEmodels/UFEmodel_mm9']

    minpostWM = minpostwm

    priorfilepath = os.path.join(inter_dir, 'priors_' + tag)

    print '\nCreate motevo parameter file %s' %tag
    print 'aligned ', aligned
    if aligned:
        motevo_params = '\n'.join(['refspecies %s' %genome,
                                   'TREE %s' %genome_dict[genome][0],
                                   'Mode %s' %mode,
                                   'minposteriorWM %s' %minpostWM,
                                   'wmdiff %f' %wmdiff,
                                   'EMprior %s' %1,
                                   'priordiff %s' %0.01,
                                   'UFEwmprior %s' %200,
                                   'UFEwmfile %s' %ufemodel_path, #%genome_dict[genome][1],
                                   'UFEwmlen %s' %wmlen,
                                   'markovorderBG %s' %1,
                                   'bgprior %s' %0.99,
                                   'bg A %s' %ATfreq,
                                   'bg T %s' %ATfreq,
                                   'bg G %s' %GCfreq,
                                   'bg C %s' %GCfreq,
                                   'restrictparses %s' %0,
                                   'priorfile %s' %priorfilepath,
                                   'minposterior %s' %0.2])
    else:
        motevo_params = '\n'.join(['refspecies %s' %genome,
                                   'TREE (%s: 1)' %genome,
                                   'Mode %s' %mode,
                                   'minposteriorWM %s' %minpostWM,
                                   'wmdiff %f' %wmdiff,
                                   'EMprior %s' %1,
                                   'priordiff %s' %0.01,
                                   'markovorderBG %s' %1,
                                   'bgprior %s' %0.99,
                                   # 'bg A %s' %0.25,
                                   # 'bg T %s' %0.25,
                                   # 'bg G %s' %0.25,
                                   # 'bg C %s' %0.25,
                                   'restrictparses %s' %0,
                                   'priorfile %s' %priorfilepath,
                                   'minposterior %s' %0.2])
        
    params_path = os.path.join(inter_dir, 'motevo_WMref_params_' + tag)
    pf = open(params_path, 'w')
    pf.write(motevo_params)
    return params_path

    
def runMotevoWM(motevo_path, alignments, params, WM, interm, tag):
    """
    runs Motevo (process is called with timeout function for WM refinement)
    """
    
    pwd = os.getcwd()
    os.chdir(interm)

    print '\nrun Motevo %s' %tag
    stime = datetime.datetime.now()
    print stime
    proc = subprocess.Popen(motevo_path + ' %s %s \"%s\"' %(alignments, params, WM),  #> report%s %(tag)
                            stdout=subprocess.PIPE,
                            stderr= subprocess.PIPE,
                            shell=True
                            )

    while proc.poll() == None:
        print proc.poll()
        time.sleep(10)
        now = datetime.datetime.now()
        if (now - stime).seconds > 600:
            os.kill(proc.pid, signal.SIGKILL)
            os.waitpid(-1, os.WNOHANG)
            print '\nMotevo weight matrix refinement did not converge.\n'
            return 0
        
    print proc.stderr.read()
    print proc.stdout.read()
    os.chdir(pwd)
    
    if proc.poll() > 0:
        print '\nMotevo weight matrix refinement not successful.\n'
        return -1
    else:
        print '\nMotevo weight matrix refinement converged.\n'
        return 1



def createLogo(mylogo_path, WM, logo_dir):

    pwd = os.getcwd()
    os.chdir(logo_dir)

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
    """
    This function/component refines a given weight matrix.
    It also makes a logo of the refined WM.
    """

    ##Ports and parameters
    odd_path = cf.get_input("odd_file") #result file from AlignPeaks spltting from SplitAlignments (typically even half)
    even_path = cf.get_input("even_file")
    ufemodel_path = cf.get_input("UFEmodel") #optional
    basefreqs = cf.get_input("BaseFrequencies") #optional
    old_wm = cf.get_input("WM")
    WM2 = cf.get_input("WM2")
    interm = cf.get_output("intermediate")
    logo = cf.get_output("Logo")
    refWM = cf.get_output("refWM")

    genome = cf.get_parameter("genome", "string")
    minpostwm = cf.get_parameter("minpostwm", "float")
    WMdiff = cf.get_parameter("wmdiff", "float")
    mylogo_path = cf.get_parameter("mylogo_path", "string")
    motevo_path = cf.get_parameter("motevo_path", "string")
    aligned = cf.get_parameter("aligned", "boolean")

    if not aligned:
        ufemodel_path = None

        ATfreq = None
        GCfreq = None

    else:
        #read base frequencies
        f = open(basefreqs)
        ATfreq = float(f.readline().strip().split()[1])
        GCfreq = float(f.readline().strip().split()[1])
        f.close()


    os.mkdir(interm)
    
    wmlen = len(open(old_wm).readlines())-4
    

    if not WM2: #only do the refinement if there is no WM2 given
        #WM refinement on even. Save wms.update in out_dir
        #start WMref. If it takes longer than some time(2 min), then kill Motevo and restart with higher wmdiff (because sometime it doesn't converge)
        tag = 'wmRef'
        #(WMREFparams_path, sf) = giveMotevoParamFile('WMREF', genome, wmlen, interm, tag, minpostwm, WMdiff)
        #runMotevo(motevo_path, even_path, WMREFparams_path, old_wm, interm, tag)
        converged = 0  #tag for whether motevo WM refinement converged
        while not converged:
            print '\nRunning motevo WM refinement with wmdiff = %s' %WMdiff
            WMREFparams_path = giveMotevoParamFile('WMREF', genome, wmlen, interm, tag, minpostwm, WMdiff, aligned, ufemodel_path, ATfreq, GCfreq)
            converged = runMotevoWM(motevo_path, even_path, WMREFparams_path, old_wm, interm, tag)
            WMdiff *= 10
            if converged == -1:
                return -1

        #Create output WM with logo name in it. Change to logo_dir when running mylogo:
        logo_dir, logo_name = os.path.split(logo)
        wmlines = open(os.path.join(interm,'wms.updated')).readlines()
        #wmlines[1] = 'NA ' + re.sub('.pdf','',logo) + '\n'
        wmlines[1] = 'NA ' +  re.sub('.pdf','',logo_name) + '\n'

        #test whether one WM row consists just of zeros. If this is the case take a dummy WM with just ones in it
        lines = [line.strip().split()[1:-2] for line in wmlines[3:-1]]
        matrix = []
        for j in lines:
            matrix += j

        if sum(map(float, matrix)) == 0.0:
            wrong = 1
        else:
            wrong = 0
    
        o = open(refWM, 'w')
        for line in wmlines:
            if wrong:
                o.write(re.sub('0.000','1',line))
            else:
                o.write(line)
        o.close()

        wmlen2 = wmlen

    else:
        wmlen2 = len(open(WM2).readlines())-4
        os.system('touch ' + logo) #create output files, so that anduril does not complain
        os.system('touch ' + refWM)
        refWM = WM2 #then put WM2 to refWM, so that further processing works.
        wrong = 0 #unimportant variable, can be set to true or false, makes no difference.


    print 'create updated Logo'
    #createLogo(old_wm)
    createLogo(mylogo_path, refWM, logo_dir)

    return 0


component_skeleton.main.main(execute)
                                                                 
