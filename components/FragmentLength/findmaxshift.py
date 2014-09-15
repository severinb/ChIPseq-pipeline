#! /import/bc2/home/nimwegen/GROUP/local/unladen/bin/python

import sys,os
import numpy as np
import re
#from scipy import interpolate
import subprocess
import gzip
import matplotlib.pyplot as plt


def smooth(x, window_len=11,window='hanning'):

	if x.ndim != 1:
		raise ValueError, "smooth only accepts 1 dimension arrays."
	
	if x.size<window_len:
		raise ValueError, "Input vector must be larger than window size"

	if window_len<3:
		return x
	
	if not window in ['flat','hanning',',hamming','bartlett','blackman']:
		raise ValueError, "invalid window type:  'flat', 'hanning', 'hamming', 'bartlett', 'blackman' "

	s=np.r_[x]
	if window== 'flat': 	#moving ave
		w=np.ones(window_len,'d')
	else:
		w=eval('np.'+window+'(window_len)')
	y=np.convolve(w/w.sum(),s,mode='same')
	return y

def derivative(x):

	newx=[]
	for i in range(len(x)-1):
		newx.append(x[i+1]-x[i])
	newx.append(newx[len(newx)-1])		#so it has same length as original vector

	return newx


def smooth_file(f,smoothdist=50):
	x=[]
	yo=[]
	y2=[]
	for line in f:
		items=line.strip().split()
		x.append(items[0])
		yo.append(float(items[1]))	

	y=np.array(yo)
	
	y2=smooth(y,smoothdist)
	y2prime=derivative(y2)
	y2doubleprime=derivative(y2prime)

	return (y2,y2prime,y2doubleprime)

def write_smoothed_file(V,fwname):

	fw=open(fwname,'w')
	i=0
	for v in V:
		fw.write("%i	%f\n"%(i,v))
		i+=1		 
	fw.close()
	return 1

def calculate_shift(y,y1,y2,skip=80,smoothstep=50):
				
	maxy=len(y)-smoothstep
	maxy1=len(y1)-smoothstep
	maxy2=len(y2)-smoothstep

	max1=max(y[skip:maxy])
	Z1=(max1-np.average(y[skip:maxy]))/np.std(y[skip:maxy])
	max2=-1

	c=skip
	minx=-99999
	for x in y[skip:maxy]:
		if x>minx:
			ymax=c
			minx=x
		c+=1	


	minx=99999
	c=skip
	

	for x in y2[skip:maxy2]:
		if x<minx:
			y2max=c
			minx=x
		c+=1
	max2tmp=0

	step=25
	for n in range(skip,maxy1):
		if y1[n]>=-1 and y1[n]<=1:
			pre=max(y1[n-step:n])
			post=min(y1[n:n+step])
			preder=np.average(y1[(n-step):n])
			postder=np.average(y1[n:(n+step)])
			if preder>0 and postder<0:
				if y[n]>max2tmp:
					max2=n	
					max2tmp=y[n]

#	max3=min(y2[:80])

	return (ymax,max2,y2max,Z1)



if __name__ == '__main__':

	if not(len(sys.argv))==8 and not len(sys.argv) == 9:
		sys.exit('usage: findmaxshift.py input-file outputdir intermediate_outdir logfile plotfile perlPath multi (repeatPath)')  

	infile = sys.argv[1]
        out_dir=sys.argv[2]
        interm_out=sys.argv[3]
        log_file = sys.argv[4]
        plot_file = sys.argv[5]
        perlPATH = sys.argv[6]
        multi = int(sys.argv[7])

        try:
                repeatpath = sys.argv[8]
                if not(os.path.exists(repeatpath)):
                        sys.exit('GENOME DIRECTORY FOR REPEATS - %s - DOES NOT EXIST'%repeatpath)
        except IndexError:
                repeatpath = ""


	reg_cor_outfile='%s.reg_cor' %infile

        #start Erik's shift auto correlation script
        #command = '%s shifts_repeats_nwk.pl %s %s %s %s' %(perlPATH, infile, reg_cor_outfile, multi, repeatpath)
        #command = '%s shifts_repeats_evn_sep2014.pl %s %s %s %s' %(perlPATH, infile, reg_cor_outfile, multi, repeatpath)
        command = '%s shifts_repeats_sb_12_9_2014.pl %s %s %s %s' %(perlPATH, infile, reg_cor_outfile, multi, repeatpath)
        proc = subprocess.Popen (command,
                                 stdout=subprocess.PIPE,
                                 stderr= subprocess.PIPE,
                                 shell=True
                                 )
        stdout_value, stderr_value = proc.communicate()
        print stdout_value
        print stderr_value

        if proc.poll() > 0:
                print '\tstderr:', repr(stderr_value.rstrip())
                sys.exit(-1)


        (y1,y2,y3)=smooth_file(open(reg_cor_outfile))						#calculates 3 potential shifts

	write_smoothed_file(y1,'%s.smoothed'%reg_cor_outfile)	


	minfrag=80 #minimum number of reads at beginning to ignore in calculating the maximum, set to 80
	(shift1,shift2,shift3,Z)=calculate_shift(y1,y2,y3,minfrag)


	fwout=open(log_file,'w')
	fwresout=open('%s.res' %infile,'w')

	if abs(shift1-minfrag)>10:
		outvalue=shift1
	elif shift2>0:
		outvalue=shift2
	elif shift3>0:
		outvalue=shift3
	else:
		outvalue=-1

	fwout.write("shiftd0    %f,     shiftd1 %f,     shiftd2 %f,     Zd0     %f - CHOSE: %f\n"%(shift1,shift2,shift3,Z,outvalue))
	fwout.close()
	fwresout.write("FRAG_LENGTH: %f\n"%outvalue)


	##added by Severin
	##plot smoothed correlation curve with outvalue
	x = range(len(y1))
   	m = int(max(y1))
	plt.plot(x, y1, 'k') #, label="correlation")
	plt.plot([int(outvalue)]*(m+m/10),range(m+m/10), 'r', label="Fragment size %s" %int(outvalue))
	plt.xlabel("Distance")
	plt.ylabel("Cross-Correlation")
	plt.legend()
	plt.savefig(plot_file)
        plt.savefig(plot_file.rstrip('.pdf')) #doesn't work?

