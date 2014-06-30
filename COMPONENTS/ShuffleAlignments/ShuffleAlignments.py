#!/usr/bin/env python
import component_skeleton.main
import subprocess
import os, re
from string import *
import random    

def getSpecies(infile):
    """
    This function extracts all species in the infile.
    The one species that is in lines starting with >> is the reference species and is first in the output string.
    """

    refspec = ''
    specs = []

    for line in open(infile):
        if line.startswith('>'):
            t = line.strip().split('_')
            if line.startswith('>>'):
                if re.sub('>>','',t[0]) != refspec:
                    refspec = re.sub('>>','',t[0])
                else:
                    continue
            else:
                if not re.sub('>','',t[0]) in specs:
                    specs.append(re.sub('>','',t[0]))
                else:
                    continue

    return refspec, ' '.join(specs)


def execute(cf):
    """
    This component runs phils (by default) or silvias peak shuffler to produce background.
    """
    
    infile = cf.get_input("infile") #result file from AlignPeaks or spltting from SplitAlignments (typically odd half)
    outfile = cf.get_output("outfile")
    mode = cf.get_parameter("mode", "string") #phil of silvia
    iterations = cf.get_parameter("iterations", "int")
    perlPATH = cf.get_parameter("perlPATH", "string")
    pythonPATH = cf.get_parameter("pythonPATH", "string")

    genome, specs = getSpecies(infile)

    print 'Genome: %s' %genome
    print 'Species: %s' %specs

    species_dict = {}
    species_dict['hg18'] = ['mm9', 'rheMac2', 'canFam2', 'bosTau3', 'equCab1', 'monDom4']
    species_dict['hg19'] = ['rheMac2', 'mm9', 'canFam2', 'bosTau6', 'equCab2', 'monDom5']
    species_dict['dm3']  = ['droSim1','droYak2','droEre2','droAna3','dp4','droWil1','droVir3','droMoj3','droGri2']   
    species_dict['mm9'] = ['hg19', 'rheMac2', 'canFam2', 'bosTau6', 'equCab2', 'monDom5']
    
    if mode == 'phil':
        print 'Taking Phil\'s shuffler\n'
        out = outfile+'_tmp'
        for i in range(iterations):
            proc = subprocess.Popen('%s shuffle_better_all_new.pl %s %s %s >> %s' %(perlPATH, infile, genome, specs, out), 
                                    stdout=subprocess.PIPE,
                                    stderr= subprocess.PIPE,
                                    shell=True
                                    )
            stdout_value, stderr_value = proc.communicate()
            print 'STDOUT: ', stdout_value
            print 'STDERR: ', stderr_value
                
            if proc.poll() > 0:
                print '\tstderr:', repr(stderr_value.rstrip())
                return -1
            
        #Rename shuffled peaks to give them an identity.
        f = open(out)
        o = open(outfile, 'w')
        for line in f:
            if line.startswith('>>'):
                line1 = '>>' + genome + '_' + str(random.randint(0,500000000)) + '_' + str(random.randint(0,500000000)) + '\n'
                o.write(line1)
            else:
                o.write(line)

        f.close()
        os.system('rm %s' %out)


    if mode == 'silvia':
        print 'Taking Silvia\'s shuffler\n'
        out = open(outfile+'_tmp', 'w')
        for i in range(iterations):
            print i+1
            tmp = outfile+'_tmp1'
            proc = subprocess.Popen('%s shuffling_silvia.py -i %s -o %s' %(pythonPATH, infile, tmp), 
                                    stdout=subprocess.PIPE,
                                    stderr= subprocess.PIPE,
                                    shell=True
                                    )
            stdout_value, stderr_value = proc.communicate()
            print 'STDOUT: ', stdout_value
            print 'STDERR: ', stderr_value
                
            if proc.poll() > 0:
                print '\tstderr:', repr(stderr_value.rstrip())
                return -1

            out.write(open(tmp).read())
            os.system('rm %s' %tmp)

        out.close()

      	#Rename shuffled peaks to give them an identity.
        f = open(outfile+'_tmp')
        o = open(outfile, 'w')
        for line in f:
            if line.startswith('>>'):
                line1 = '>>' + genome + '_' + str(random.randint(0,500000000)) + '_' + str(random.randint(0,500000000)) + '\n'
                o.write(line1)
            else:
                o.write(line)

        f.close()
        os.system('rm %s' %outfile+'_tmp')
    

    return 0


component_skeleton.main.main(execute)
                                                                 
