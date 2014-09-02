import subprocess
import re
from sys import stdout

def arguments():
    import argparse
    parser = argparse.ArgumentParser(description='Fits beta and prior, and also run calculate the enrichment scores')
    parser.add_argument('-w', '--wm',
                    action="store", dest="WMs", type=str, nargs='+')    
    parser.add_argument('-p', '--prior',
                    action="store", dest="prior", type=str)
    parser.add_argument('-o', '--outfile',
                    action="store", dest="outfile", type=str)    
    results = parser.parse_args()
    return results


def concatenate(WMs, priors, outf):
    cmd = 'cat ' + ' '.join(WMs)
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    if not priors:
        for line in proc.stdout: outf.write(line)
    else:
        for line in proc.stdout:
            outf.write(line)
            if re.search('^NA ', line):
                motifName = line.split()[-1].strip()
                outf.write('PW\t%s\n' % (priors[motifName]))
    return 0


def loadPriors(priorFile):
    priors = {}
    if args.prior:
        priors = dict([(line.split()[0], \
                          line.split()[1]) for line in open(priorFile)])
    return priors


def main():
    args = arguments()
    priors = loadPriors(args.prior)
    if not args.outfile:  # if no output filename is passes, simply write to the terminal
        concatenate(args.WMs, priors, stdout)
    else:   # otherwise write to the file
        with open(args.outfile, 'w') as outf:  
            concatenate(args.WMs, priors, outf)


if __name__ == '__main__':
    main()
