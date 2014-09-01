#!/usr/bin/env python
import component_skeleton.main
import subprocess
import os
from datetime import datetime

def execute(cf):
    infile = cf.get_input("in_file")
    out_dir = cf.get_output("out_dir")
    fg_winlen = cf.get_parameter("fg_winlength", "int")
    bg_winlen = cf.get_parameter("bg_winlength", "int")
    plotfile = cf.get_output("Z_hist")
    perlPATH = cf.get_parameter("perlPATH", "string")
    R_LIBS_USER = cf.get_parameter("R_LIBS_USER", "string")
    noBG = cf.get_parameter("noBG", "boolean")
    logfile = cf.get_output("log_file")

    T1 = datetime.now()

    os.mkdir(out_dir)

    #infile = os.path.join(indir, 'all_counts')

    print 'fit sigma mu rho'

    if not R_LIBS_USER == '': 
        os.environ['R_LIBS_USER'] = R_LIBS_USER

    if noBG:
        # if there is no background sample given, just add a column to the input file with uniform background counts
        # #chr    start   stop    middle  fg_0    fg_1    fg_2    bg_0
        # chr1    2999250 2999750 2999500 0       0       0       1.000000
        # chr1    2999500 3000000 2999750 0       0       0       1.000000
        # chr1    2999750 3000250 3000000 1.500000        0       0       1.000000

        tmp_infile = os.path.join(os.path.split(out_dir)[0], 'tmp_infile')
        o = open(tmp_infile, 'w')

        for line in open(infile):
            if line.startswith('#'):
                o.write(line.strip() + '\tbg_0\n')
            else:
                o.write(line.strip() + '\t10.00000\n')

        o.close()

        command = perlPATH +' fit_sigma_mu_rho_with_replicatesS.pl %s %s %s %s %s' % (tmp_infile, fg_winlen, bg_winlen, out_dir, plotfile)

    else:

        command = perlPATH +' fit_sigma_mu_rho_with_replicatesS.pl %s %s %s %s %s' % (infile, fg_winlen, bg_winlen, out_dir, plotfile)


    #/import/bc2/soft/bin/perl5/perl
    proc = subprocess.Popen(command,
                            stdout=subprocess.PIPE,
                            stderr= subprocess.PIPE,
                            shell=True
                            )
    stdout_value, stderr_value = proc.communicate()
    print stdout_value
    print stderr_value

    if proc.poll() > 0:
        print '\tstderr:', repr(stderr_value.rstrip())
        return -1


    #add stats file (fitted mu, sigma and rho) to log_file
    f = open(os.path.join(out_dir, 'stats_get_zvals'))
    text = f.read()
    f.close()

    T2 = datetime.now()
    time = 'Running time for peak caller: ' + str(T2-T1) + '\n'
    lf = open(logfile, 'w')
    lf.write(text)
    lf.write(time)
    lf.close

    return 0

component_skeleton.main.main(execute)
