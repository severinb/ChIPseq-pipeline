import subprocess
import os
import datetime
import time
import re

def create_motevo_param_file(param_filename, site_filename, prior_filename, genome, priordiff=0.05, minposterior=0.0, prior=None):
    
    if not prior:
        prior=0.99
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
        'bgprior %f' % prior,  # as an initial value for fitting the prior
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


def run_motevo(WM, sequences, interm_dir, genome, prior=None, minposterior=.0):
    motevo_path = '/import/bc2/home/nimwegen/GROUP/software/motevo_ver1.03/bin/motevo'
    stime = datetime.datetime.now()
    motifName = os.path.basename(WM)
    print '\nrunnig Motevo for %s' % motifName
    siteFilename = os.path.join(interm_dir, '%s.sites' % motifName)
    priorFilename = os.path.join(interm_dir, '%s.priors' % motifName)
    paramFilename = os.path.join(interm_dir, '%s.params' % motifName)
    create_motevo_param_file(paramFilename, siteFilename, priorFilename, genome, minposterior=minposterior, prior=prior)
    cmd = ' '.join([
        motevo_path,
        sequences,
        paramFilename,
        WM ])
    proc = subprocess.Popen(cmd,
                            stdout=subprocess.PIPE,
                            stderr= subprocess.PIPE,
                            shell=True)    
    while proc.poll() == None:
        # print proc.poll()
        time.sleep(10)
        now = datetime.datetime.now()
        if (now - stime).seconds > 600:
            os.kill(proc.pid, signal.SIGKILL)
            os.waitpid(-1, os.WNOHANG)
            print '\nMotevo weight matrix refinement did not converge.\n'
            return None, None
    # print proc.stderr.read()
    # print proc.stdout.read()
    if proc.poll() > 0:
        print '\nMotevo weight matrix refinement not successful.\n'
        return None, None
    else:
        print '\nMotevo weight matrix refinement converged.\n'
        return siteFilename, priorFilename
    return (siteFilename, priorFilename)


def extract_priors(prior_file):
    prior = [float(l.split()[1]) for l in open(prior_file) if re.search('^background', l)].pop()
    return prior

