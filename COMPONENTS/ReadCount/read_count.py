#! /import/bc2/home/nimwegen/GROUP/local/unladen/bin/python

import os
import optparse
from string import *
import warnings
from sys import exit
import subprocess
import gzip
import yaml
import contextlib
import sys
import logging


def option_parser():
    """ Makes a parser of input arguments to the program."""
    parser = optparse.OptionParser(usage='Usage: %prog <options>')

    parser.add_option('-y' , '--yaml' , help='YAML parameter file', \
                      dest='PARAM_FILE', action='store', type='string', nargs=1)

    parser.add_option('-w' , '--weight' , help='Is there any weight for the reads (Default False, means NO)', \
                      dest='IS_THERE_WEIGHT', action='store_true')

    parser.add_option('-v' , '--verbose' , help='Verbose mode', \
                      dest='VERBOSE', action='store_true')

    (opts,args) = parser.parse_args()
    if opts.PARAM_FILE is None:
        print 'options --yaml is mandatory\n'
        parser.print_help()
        exit(-1)

    return (opts,args) 


def YAML_parser (opts):
    """ given an optparse instance that have variable PARAM_FILE which is a
    parameter file in YAML format, it returns a dictionary containing the parameters"""

    in_file = opts.PARAM_FILE
    try:
        with contextlib.closing (open(in_file)) as params_file:
            params = yaml.load (params_file)            
    except IOError:
        os.syserr.write ("Error in openning file %s" \
                         % in_file )
        exit(-1)

    if 'INPUT FILE' not in params:
        print "In the parameter file, the tag \'INPUT FILE\' is missing, which is mandatory.\n"
        exit(-1)

    if 'CHROMOSOME INFO' not in params:
        print "In the parameter file \'CHROMOSOME INFO\' is missing..\n"
        exit(-1)

    if 'RESULT DIR' not in params:
        if opts.VERBOSE:
            print "The results will be copied in ./RESULT directory."
            print "To change the directory name, use tag \'RESULT DIR\' in the parameter file."
            params['RESULT DIR'] = 'RESULT'

    if 'RESULT FILE' not in params:
        params['RESULT FILE'] = 'results'

    if 'WINDOW' not in params:
        params['WINDOW'] = 1000

    if 'STEP' not in params:
        params['STEP'] = 500

    if params ['WINDOW'] <= 50 or params ['STEP'] <= 25:
        print "The length of the window/step is very short.\n"
        exit(-1)

    params['VERBOSE'] = opts.VERBOSE
    params['IS THERE WEIGHT'] = opts.IS_THERE_WEIGHT

    return params


def make_dir (dirname):
    """Given a name, it makes a new directory in the current directory"""
    proc = subprocess.Popen ('mkdir %s' % dirname,                      
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                shell=True
                         )
    stdout_value, stderr_value = proc.communicate()
    if proc.poll() > 0:
        if stderr_value.rfind('File exists') > 0:        
            sys.stderr.write ( "\nWarning in creating directory %s\n" % dirname)
            print '\tstderr:', repr(stderr_value.rstrip())
            print "The program continues, but it may remove/change file(s) that are already in %s directory" \
                  % dirname

        else:
         sys.stderr.write ( "\nError in creating directory %s\n" % dirname )
         print '\tstderr:', repr(stderr_value.rstrip())
         return -1
         
    return 0
        
def join_results(infile,outfile,rep):
    try:
        f = open (infile)
        res = open (outfile, 'w')
    except IOError, e:
        print e.message()
        return -1
    
    SUM = lambda x,y: float(x) + float(y)
    tmp = f.readline().rstrip().split()
    window = tmp[0:4] + ['0.00' for i in xrange(rep+1)]
    index = int(tmp.pop())
    window[3 + index] = tmp.pop()
    for line in f:
        tmp = line.rstrip().split()
        if window[0:2] != tmp[0:2]:
            window[-1] = str(reduce(SUM, window[4:-1])) + '\n'
            res.write('\t'.join(window))
            window = tmp[0:4] + ['0.00' for i in xrange(rep+1)]
        index = int(tmp.pop())
        window[3 + index] = tmp.pop()
    res.flush()
    f.close()
    return 0
            

if __name__ == '__main__':

    from time import clock
    start_time = clock()

    (opts,args) = option_parser()
    params = YAML_parser (opts)
    
    if make_dir (params ['RESULT DIR']):
        print  "Error in creating directory %s. Program quits!\n" \
              % params['RESULT DIR']
        exit(-1)

    from Count import *    # here the main things happen!
    if len(params['INPUT FILE']) > 1:  # if there are more than one replicate, it goes on multi-threading
        result = open (os.path.join(params['RESULT DIR'],"tmp"), 'w')        
        threads = [FileCount(params,params['INPUT FILE'][i],result, str(i+1)) \
                   for i in xrange(len(params['INPUT FILE']))]
        for t in threads: t.start()
        for t in threads: t.join()
        result.flush()

        proc = subprocess.Popen('sort -k 1,1 -k 2,2n -k 5,5n %s > %s' % (os.path.join(params ['RESULT DIR'], "tmp"), os.path.join(params ['RESULT DIR'],".results")),                      
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    shell=True
                                )
        stdout_value, stderr_value = proc.communicate()
        if proc.poll() > 0:
            print '\tstderr:', repr(stderr_value.rstrip())
            exit(-1)

        proc = subprocess.Popen ('rm %s' % (os.path.join(params ['RESULT DIR'], "tmp")),                      
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    shell=True
                             )
        stdout_value, stderr_value = proc.communicate()
        if proc.poll() > 0:
            print '\tstderr:', repr(stderr_value.rstrip())
            exit(-1)

        if params['VERBOSE']:
            print 'Merging the results for the input files into a single file'
        # print params['INPUT FILE']
        if join_results(os.path.join(params ['RESULT DIR'],".results"), \
                         os.path.join(params ['RESULT DIR'], params['RESULT FILE']), \
                         len(params['INPUT FILE'])):
            print '\tError in joining the results'
            exit(-1)

        proc = subprocess.Popen ('rm %s' % (os.path.join(params ['RESULT DIR'], ".results")),                      
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    shell=True
                             )
        stdout_value, stderr_value = proc.communicate()
        if proc.poll() > 0:
            print '\tstderr:', repr(stderr_value.rstrip())
            exit(-1)
    else:
        result = open(os.path.join(params['RESULT DIR'], "tmp"), 'w')
#        count_instance = Count( params,  params['INPUT FILE'][0] , result, '')
#        count_instance.count_all_chrom()
        threads = [FileCount(params,params['INPUT FILE'][i],result, '') \
                   for i in xrange(len(params['INPUT FILE']))]
        for t in threads: t.start()
        for t in threads: t.join()
        result.flush()
        os.system('sort -T /import/bc2/scratch/tmp/ -k 1,1 -k 2,2n -k 5,5n -u %s | awk \'{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$5}\' > %s' % \
                  (os.path.join(params ['RESULT DIR'], "tmp"), os.path.join(params['RESULT DIR'], "results")))
        os.system('rm %s' % os.path.join(params ['RESULT DIR'], "tmp"))
        # os.system('rm %s' %  os.path.join(params['RESULT DIR'], "tmp"))
    elapsed = (clock() - start_time)
    if opts.VERBOSE:
        print "The time took to finsih: ", elapsed

    
    
        


    
        
    
    
    
    
    

    

        
    

    


    
    



    



        
        
    

