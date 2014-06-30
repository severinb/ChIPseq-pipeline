#! /import/bc2/home/nimwegen/GROUP/local/unladen/bin/python


def option_parser():
    """ Makes a parser of input arguments to the program."""

    import optparse
    
    parser = optparse.OptionParser(usage='Usage: %prog <options>')

    parser.add_option('-i' , '--input' , help='A file which contains the input files', \
                      dest='FILE', action='store', type='string', nargs=1)

    parser.add_option('-d' , '--directory' , help='result directory', \
                      dest='DIR', action='store', type='string', nargs=1)

    parser.add_option('-c' , '--chrom' , help='chromosome info file', \
                      dest='CHR', action='store', type='string', nargs=1)    

    parser.add_option('-t' , '--weight' , help='Is there any weight for the reads (Default False, means NO)', \
                      dest='IS_THERE_WEIGHT', action='store_true')

    parser.add_option('-v' , '--verbose' , help='Verbose mode', \
                      dest='VERBOSE', action='store_true')

    parser.add_option('-w' , '--window' , help='window size', \
                      dest='WINDOW', action='store', type='int')
    
    parser.add_option('-s' , '--step' , help='step size', \
                      dest='STEP', action='store', type='int')


    (opts,args) = parser.parse_args()

    if opts.FILE is None:
        print 'options --input is mandatory\n'
        parser.print_help()
        exit(-1)    

    return (opts,args)


def run(cmd):
    """
    given a command, it runs
    """
    import subprocess
    proc = subprocess.Popen (cmd,                      
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                shell=True
                         )
    stdout_value, stderr_value = proc.communicate()
    if proc.poll() > 0:
        import sys
        sys.stderr.write ( "\nError in creating directory %s\n" % dirname )
        print '\tstderr:', repr(stderr_value.rstrip())
        return False
    else:
        return True


if __name__ == '__main__':
    (opts, args) = option_parser()
    from string import *
    params = {'INPUT FILE': [fnames.rstrip() for fnames in open(opts.FILE)],
              'WINDOW': opts.WINDOW,
              'STEP': opts.STEP,
              'RESULT DIR': opts.DIR,
              'CHROMOSOME INFO': opts.CHR
        }


    import yaml
    import random
    try:
        f = open('%s.params.yaml' % opts.FILE, "w")
        yaml.dump(params, f)
        f.flush()
    except IOError,e:
        print e
        print 'An error occured while making parameter file!'
        from sys import exit
        exit(-1)


    def option(opt, val):
        if val:
            return '-%s' % opt
        else:
            return ''
        
 # '/import/bc2/home/nimwegen/GROUP/hseq_pipeline/saeed/ChIP_read_counter/read_count.py'
    cmd = ' '.join(['./read_count.py',
                    '-y %s' % ('%s.params.yaml' % opts.FILE),
                    option('w',opts.IS_THERE_WEIGHT),
                    option('v',opts.VERBOSE),
                    ])
    print cmd
    # run(cmd)

    
                    
    
    

    
    
