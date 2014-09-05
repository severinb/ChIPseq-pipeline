import subprocess
import os
import datetime
import time
import re
import uuid
from concatenate_motifs import concatenate

def cleanup(infile):
    """
    To delete the files in the scratch directory
    """
    cmd = "rm -f '%s'" % infile
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, \
                                stderr=subprocess.PIPE)
    result = proc.communicate()
    if proc.returncode:
        print result[1]
        print 'Problem with removing these file at RunMotevo: ' % infile
    return 0 


def create_motevo_param_file(param_filename, site_filename, prior_filename, genome, priordiff=0.05, minposterior=0.0, prior=None):    
    if not prior:
        prior = '0.99'
        EMprior = 1
    else:
        EMprior = 0
    param_file = open(param_filename, 'w')    
    param_file.write('\n'.join([
        'refspecies %s' % genome,
        'TREE (%s: 1)' % genome,
        'Mode TFBS',
        'EMprior %d' % EMprior,
        'priordiff %f' % priordiff,
        'markovorderBG 0',
        'bgprior %s' % prior,  # as an initial value for fitting the prior
        'bg A 0.25',
        'bg T 0.25',
        'bg G 0.25',
        'bg C 0.25',
        'restrictparses 0',
        'sitefile %s' % site_filename,
        'priorfile %s' % prior_filename,
        'minposterior %f' % minposterior,
        'printsiteals 0',
        ]))
    param_file.close()
    return 0


def concatenate_motifs(WMs, outdir, priors=None):
    motifNames = map(os.path.basename, WMs)
    newMotifName = '_'.join(motifNames)
    fname = os.path.join(outdir, newMotifName)  # the name of the new 'super' motif
    with open(fname, 'w') as outf:
        concatenate(WMs, priors, outf)  # from concatenate_motifs.py
    return fname


def loadAllPriors(priorFile):
    priors = {}
    if priorFile:
        priors = dict([(line.split()[0], \
                        line.split()[1]) for line in open(priorFile)])
    return priors


def randomLabel():
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(64))


def run_motevo(WM, sequences, interm_dir, genome, priorFile=None, minposterior=.0):
    motevo_path = '/import/bc2/home/nimwegen/GROUP/software/motevo_ver1.03/bin/motevo'
    priors = loadAllPriors(priorFile)
    stime = datetime.datetime.now()
    isConcatenatedMotif = False 
    if len(WM) > 1:
        WM = concatenate_motifs(WM, interm_dir, priors=priors)
        isConcatenatedMotif = True 
    else:
        WM = WM[0]
    
    motifName = os.path.basename(WM)
    # print '\nrunnig Motevo for %s' % motifName
    uniqueID = str(uuid.uuid4())
    TFname = os.path.basename(os.path.dirname(interm_dir)).replace('_FgBg-enrichmentScores_all_motifs', '')
    siteFilename = os.path.join('/scratch/', '%s_%s_%s.sites' % (motifName, TFname, uniqueID))
    priorFilename = os.path.join('/scratch/', '%s_%s_%s.priors' % (motifName, TFname, uniqueID))
    paramFilename = os.path.join('/scratch/', '%s_%s_%s.params' % (motifName, TFname, uniqueID))
    create_motevo_param_file(paramFilename, siteFilename, priorFilename, genome, \
                             minposterior=minposterior, prior=priors.get('background'))
        
    cmd = ' '.join([
        motevo_path,
        "\'%s\'" % sequences,
        "\'%s\'" % paramFilename,
        "\'%s\'" % WM ])
    # print cmd
    proc = subprocess.Popen(cmd,
                            stdout=subprocess.PIPE,
                            stderr= subprocess.PIPE,
                            shell=True)
    result = proc.communicate()   
    if proc.returncode:
        print result[0]
        print result[1]
        return None, None, None, None    
    if isConcatenatedMotif: # to remove the temporary motifs that are created via concatenation
        cleanup(WM)
    return siteFilename, priorFilename, paramFilename, WM


def extract_priors(prior_file):
    prior = [float(l.split()[1]) for l in open(prior_file) if re.search('^background', l)].pop()
    return prior

