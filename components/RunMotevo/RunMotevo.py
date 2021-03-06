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

    sitefilepath = '/scratch/motevo_sites_' + tag + '_%i' %randint(10000)
    priorfilepath = os.path.join(inter_dir, 'priors_' + tag)

    if mode == 'WMREF':
        print '\nCreate motevo parameter file %s' %tag
        print 'aligned ', aligned
        if aligned:
            motevo_params = '\n'.join(['refspecies %s' %genome,
                                       'TREE %s' %genome_dict[genome][0],
                                       'Mode %s' %mode,
                                       'minposteriorWM %s' %minpostWM,
                                       'wmdiff %f' %wmdiff,
                                       'EMprior %s' % 1,
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
                                       'sitefile %s' %sitefilepath,
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
                                       'bg A %s' %0.25,
                                       'bg T %s' %0.25,
                                       'bg G %s' %0.25,
                                       'bg C %s' %0.25,
                                       'restrictparses %s' %0,
                                       'sitefile %s' %sitefilepath,
                                       'priorfile %s' %priorfilepath,
                                       'minposterior %s' %0.2])

        params_path = os.path.join(inter_dir, 'motevo_WMref_params_' + tag)
        pf = open(params_path, 'w')
        pf.write(motevo_params)
        return (params_path, sitefilepath)

    if mode == 'TFBS':
        print '\nCreate motevo parameter file %s' %tag
        print 'aligned', aligned
        if aligned:
            motevo_params = '\n'.join(['refspecies %s' %genome,
                                       'TREE %s' %genome_dict[genome][0],
                                       'Mode %s' %mode,
                                       'EMprior %s' %1,
                                       'priordiff %s' %0.05,
                                       'UFEwmprior %s' %200,
                                       'UFEwmfile %s' %ufemodel_path, #%genome_dict[genome][1],
                                       'UFEwmlen %s' %wmlen,
                                       'UFEprint %s' %0,
                                       'markovorderBG %s' %1,
                                       'bgprior %s' %0.99,
                                       'bg A %s' %ATfreq,
                                       'bg T %s' %ATfreq,
                                       'bg G %s' %GCfreq,
                                       'bg C %s' %GCfreq,
                                       'restrictparses %s' %0,
                                       'sitefile %s' %sitefilepath,
                                       'priorfile %s' %priorfilepath,
                                       'printsiteals %s' %0,
                                       'minposterior %f' %0.0])
        else:
            motevo_params = '\n'.join(['refspecies %s' %genome,
                                       'TREE (%s: 1)' %genome,
                                       'Mode %s' %mode,
                                       'EMprior %s' %1,
                                       'priordiff %s' %0.05,
                                       'markovorderBG %s' %1,
                                       'bgprior %s' %0.99,
                                       'bg A %s' %0.25,
                                       'bg T %s' %0.25,
                                       'bg G %s' %0.25,
                                       'bg C %s' %0.25,
                                       'restrictparses %s' % 0,
                                       'sitefile %s' %sitefilepath,
                                       'priorfile %s' %priorfilepath,
                                       'printsiteals %s' %0,
                                       'minposterior %f' %0.0])            

        params_path = os.path.join(inter_dir, 'motevo_TFBS_params_' + tag)
        pf = open(params_path, 'w')
        pf.write(motevo_params)
        return (params_path, sitefilepath)    

    
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


    
def runMotevo(motevo_path, alignments, params, WM, interm, tag):
    """
    runs Motevo
    """
    
    pwd = os.getcwd()
    os.chdir(interm)

    print '\nrun Motevo %s' %tag
    proc = subprocess.Popen(motevo_path + ' %s %s \"%s\"' %(alignments, params, WM),  #> report%s %(tag)
                            stdout=subprocess.PIPE,
                            stderr= subprocess.PIPE,
                            shell=True
                            )

    stdout_value, stderr_value = proc.communicate()
    print stdout_value
    print stderr_value
    os.chdir(pwd)
    
    if proc.poll() > 0:
        print '\tstderr:', repr(stderr_value.rstrip())
        return -1
    else:
        return 0


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


def ROCanalysis(non_sites_true, non_sites_bg, ref_sites_true, ref_sites_bg, ROCplot, odd_path, bgPeaks):
    """
    This function produces a ROC plot with functions from the ROC.py module.
    """

    nonref_trueSitesPost_dict = ROC2.get_peak_posteriors(non_sites_true, odd_path)
    nonref_bgSitesPost_dict = ROC2.get_peak_posteriors(non_sites_bg, bgPeaks)
    ref_trueSitesPost_dict = ROC2.get_peak_posteriors(ref_sites_true, odd_path)
    ref_bgSitesPost_dict = ROC2.get_peak_posteriors(ref_sites_bg, bgPeaks)


    nonref_senslist,nonref_ppvlist = ROC2.TP_P_FN_N(nonref_trueSitesPost_dict, nonref_bgSitesPost_dict)
    ref_senslist, ref_ppvlist = ROC2.TP_P_FN_N(ref_trueSitesPost_dict, ref_bgSitesPost_dict)

    ROC2.plot_ROC(nonref_senslist, nonref_ppvlist, ref_senslist, ref_ppvlist, ROCplot)

    #compute area under the ROC curves
    def giveAUC(x, y):
        ##first fit both lists with a spline
        base = linspace(0,1,len(x))
        
        #ytck = interpolate.splrep(base, y, s=0)
        #y = interpolate.splev(base, ytck, der=0)

        #xtck = interpolate.splrep(base, x, s=0)
        #x = interpolate.splev(base, xtck, der=0)

        # print min(x), min(y)
        # print max(x), max(y)
        # ROC2.plot_ROC(x,y,x,y,os.path.join(os.path.split(ROCplot)[0],'spline_%s.pdf' %rand()))

        ##trapezoid integration
        #get averages between two points in y (ppv) list (y starts with 0 to 1)
        yavg = [mean(y[i:i+1]) for i in arange(len(y)-1)]
        #get differences between two points in (x starts with 1 to 0)
        xd = [x[i] - x[i+1] for i in arange(len(x)-1)]

        AUC = dot(xd,yavg)

        return AUC
    
    
    #return mean(nonref_senslist) + mean(nonref_ppvlist), mean(ref_senslist) + mean(ref_ppvlist)
    return giveAUC(nonref_senslist, nonref_ppvlist), giveAUC(ref_senslist, ref_ppvlist)

def execute(cf):
    """
    This function/component refines a given weight matrix, then predicts sites for old and refined WM and true sites and on background peaks.
    With this it makes a positive predictive value - sensitivity plot.
    It also makes a logo of the refined WM.
    """

    ##Ports and parameters
    odd_path = cf.get_input("odd_file") #result file from AlignPeaks spltting from SplitAlignments (typically even half)
    even_path = cf.get_input("even_file")
    ufemodel_path = cf.get_input("UFEmodel") #optional
    basefreqs = cf.get_input("BaseFrequencies") #optional
    old_wm = cf.get_input("WM")
    WM2 = cf.get_input("WM2")
    bgPeaks = cf.get_input("shuffledPeaks")
    interm = cf.get_output("intermediate")
    ROCplot = cf.get_output("sens_ppv")
    logo = cf.get_output("Logo")
    refWM = cf.get_output("refWM")
    qual_file_nonref = cf.get_output("nonref_WM_Quality")
    qual_file_ref = cf.get_output("ref_WM_Quality")

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
            (WMREFparams_path, sf) = giveMotevoParamFile('WMREF', genome, wmlen, interm, tag, minpostwm, WMdiff, aligned, ufemodel_path, ATfreq, GCfreq)
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


    ##start with ROC analysis
    #get sitefile for nonrefined WM TFBS on odd and background
    tag = 'tfbsNonOdd'
    (TFBSnon_paramsOdd,sitesNonOdd) = giveMotevoParamFile('TFBS', genome, wmlen, interm, tag, None, None, aligned, ufemodel_path, ATfreq, GCfreq)
    runMotevo(motevo_path, odd_path, TFBSnon_paramsOdd, old_wm, interm, tag)
    tag = 'tfbsNonBg'
    (TFBSnon_paramsBg,sitesNonBg) = giveMotevoParamFile('TFBS', genome, wmlen, interm, tag, None, None, aligned, ufemodel_path, ATfreq, GCfreq)
    runMotevo(motevo_path, bgPeaks, TFBSnon_paramsBg, old_wm, interm, tag)

    #get sitefile for refined WM TFBS
    tag = 'tfbsRefOdd'
    (TFBSref_paramsOdd,sitesRefOdd) = giveMotevoParamFile('TFBS', genome, wmlen2, interm, tag, None, None, aligned, ufemodel_path, ATfreq, GCfreq)
    runMotevo(motevo_path, odd_path, TFBSref_paramsOdd, refWM, interm, tag)
    tag = 'tfbsRefBg'
    (TFBSref_paramsBg,sitesRefBg) = giveMotevoParamFile('TFBS', genome, wmlen2, interm, tag, None, None, aligned, ufemodel_path, ATfreq, GCfreq)
    runMotevo(motevo_path, bgPeaks, TFBSref_paramsBg, refWM, interm, tag)

    #ROC analysis
    print 'create ROC plot'
    non_ref_Quality, ref_Quality = ROCanalysis(sitesNonOdd, sitesNonBg, sitesRefOdd, sitesRefBg, ROCplot, odd_path, bgPeaks)

    print "Non refined motif quality: ", non_ref_Quality
    if not wrong:
        print "refined motif quality: ", ref_Quality
    else:
        ref_Quality = 0
        print "refined motif quality (unsuccessful refinement): ", ref_Quality


    #write qualities to a file
    qnonref = open(qual_file_nonref, 'w')
    qref = open(qual_file_ref, 'w')
    qnonref.write(str(non_ref_Quality))
    qref.write(str(ref_Quality))
    qnonref.close()
    qref.close()

    print 'create updated Logo'
    #createLogo(old_wm)
    createLogo(mylogo_path, refWM, logo_dir)

    #clean up:
    os.system('rm %s %s %s %s' %(sitesNonOdd, sitesNonBg, sitesRefOdd, sitesRefBg))

    return 0


component_skeleton.main.main(execute)
                                                                 
