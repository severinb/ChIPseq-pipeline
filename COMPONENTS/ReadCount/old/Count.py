import warnings
import contextlib
import os
from Window import *
from Read import *
import threading
import optparse
import gzip
import re

class Count:

    def __init__ (self, params, in_file, result, tag):
        self.in_file = in_file
        self.params = params
        self.tag = tag
        Count.window_array_length = params['WINDOW'] / params['STEP']
        try:
            if re.search('\.gz$', in_file):
                self.read_file = gzip.open (in_file)
            else:
                self.read_file = open(in_file)
        except IOError:
            print "An error happened in openning file %s" \
                  % in_file
            
            exit(-1)

        Count.chromInfo = self.load_chromInfo(self.params['CHROMOSOME INFO'])
        print Count.chromInfo

        try:
            self.result = result
            # self.result = open ((os.path.join(self.params['RESULT DIR'],"results", 'r+')
            # self.result = open (os.path.join(self.params['RESULT DIR'],in_file), 'w')
        except IOError:
            print "Error in creating file result in %s" \
                  % os.path.join(self.params['RESULT DIR'],in_file)
        
        

    def load_chromInfo(self, in_file):
        """ This function loads the chromosome information into a dictionary"""

        try:
            with contextlib.closing (open(in_file)) as chrom_file:

                Count.chromInfo = dict([(line.split()[0], int(line.split()[1])) for line in  open (in_file) ]) #\
                           #if line.split()[0].rfind('random') < 0 ])

        except IOError:
            os.syserr.write ("Error in openning file %s" \
                             % in_file )
            exit(-1)

        return Count.chromInfo


    def count_all_chrom(self):
        """This method counts all the reads in all chromosomes in a sorted fashion """

        chroms = sorted(Count.chromInfo.keys())
        initial_read = ""
        for CHR in chroms:
            if self.params['VERBOSE']:
                print 'Working on chromosome %(chr)s in file %(file)s' \
                      % {'chr':CHR, 'file':self.in_file}
            
            if Count.chromInfo[CHR] < self.params['WINDOW']:
                WarningMessage = "The length of chromosome %s is lower than the window size %d" \
                                     % (CHR, params['WINDOW'])
                warnings.warn (WarningMessage, Warning)
                continue
            
            initial_read = self.read_count_in_chrom (CHR, initial_read)

            if initial_read == 'EOF':
                break
            
        return 1
           
    
    def read_count_in_chrom (self, Chr, initial_read):
        """This methods counts the number of reads in chromosome Chr.
        It, does not set the position in the read_file to zero"""

        window_array = [Window(Chr, -2, end=-1) \
                        for i in xrange(Count.window_array_length)]

        if initial_read != "":
            current_read = initial_read
        else:
            current_read = Read (self.read_file.readline(), is_weight_given=self.params['IS THERE WEIGHT'])
        
        while (current_read.give_chrom() == Chr):             
             index = 0             
             for window in window_array:
                 if current_read.give_end()  >  window.give_end():
                     if  window.give_count () > 0:
                         self.result.write ( '%s\t%s\n' % (window.info().rstrip(), self.tag))
                     window_start = (current_read.give_start() / self.params['WINDOW']) * self.params['WINDOW']
                     offset = index * self.params['STEP']
                     window.change_window_info ([Chr, window_start + offset, \
                                            min(window_start+self.params['WINDOW'] + offset, Count.chromInfo[Chr]), 0. ])

                 if current_read.give_start () >=  window.give_start () and \
                        current_read.give_start ()   <=  window.give_end ():
                     window.increase_read_count(current_read.give_weight())

                 index += 1

             read_line = self.read_file.readline()

             if read_line=="":                 
                 for window in window_array:
                     if window.give_count () > 0:
                         self.result.write ( '%s\t%s\n' % (window.info().rstrip(), self.tag))
                         
                 return 'EOF'

             current_read.change_info ( read_line, is_weight_given=self.params['IS THERE WEIGHT'] )

             while (current_read.give_end() > Count.chromInfo[current_read.give_chrom()] and \
                    current_read.give_chrom() == Chr) :

                 WarningMessage = "The read that starts at %d and ends %d is outside of the range of chromosome %s which has a length of %d!" \
                                 % (current_read.give_start(), current_read.give_end(), current_read.give_chrom(), Count.chromInfo[Chr])            
                 warnings.warn (WarningMessage, Warning)
                 read_line = self.read_file.readline()
                 if read_line != "":
                     current_read.change_info (read_line, is_weight_given=self.params['IS THERE WEIGHT'])
                 else:       
                     return 'EOF'

        for window in window_array:
             if  window.give_count () > 0:
                 self.result.write ( '%s\t%s\n' % (window.info ().rstrip(), self.tag))

        return current_read


class FileCount (threading.Thread):
    """This class allows us to make instances of Counter class and run them in a multi-threading fashion."""
    def __init__(self, params, in_file, result, tag):        
        self.in_file = in_file
        self.params = params
        condition = threading.Condition()
        threading.Thread.__init__(self, name=in_file, args=(condition,))        
        self.result = result
        self.tag = tag

    def run (self):
        try:
           count_instance = Count( self.params, self.in_file , self.result, self.tag)
           count_instance.count_all_chrom()
        except Exception,e:
            print "An error occured while counting the reads from file %s" \
                  % self.in_file, e

    @staticmethod
    def consumer(cond):
        """wait for the condition and use the resource"""
        t = threading.currentThread()
        with cond:
            cond.wait()

    def producer(cond):
        """set up the resource to be used by the consumer"""
        with cond:
            cond.notifyAll()
            
                

        
              
                 
             
             

        

        
        

        

        
        
    
