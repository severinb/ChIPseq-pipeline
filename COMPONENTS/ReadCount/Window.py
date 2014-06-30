
class Window:

    def __init__ (self, chrom, start, size):
        self.chrom = chrom
        self.start = start
        self.end = start + size
        self.read_count = 0.
        self.middle = (self.start+self.end)/2

    def __init__ (self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.read_count = 0.
        self.middle = (self.start+self.end)/2        

    def increase_read_count(self):   
        self.read_count += 1.

    def increase_read_count(self,w):      
        self.read_count += w

    def give_start(self):
        return self.start

    def give_end(self):
        return self.end

    def give_size(self):
        return (self.end - self-start)

    def give_middle(self):
        return self.middle

    def give_chrom(self):       
        return self.chrom

    def give_count(self):
        return self.read_count

    def change_count(self, f):
        try:
            self.read_count = f
        except Exception, e:
            print e
        return -1

    def change_chrom(self, Chr):
        self.chrom = Chr
        return self.chrom

    def change_start(self, start):
        self.start = start
        return self.start

    def change_end(self,end):
        self.end = end
        return self.end


    def change_window_info (self, window):     
        """
        This method take as an input a list of the following information:
        CHROM START END COUNT ...
        """
        try:
            self.chrom = window[0]
            self.start = int(window[1])
            self.end = int(window[2])
            self.read_count =  float(window[3])
            self.middle = (self.start+self.end)/2
        except IndexError:
            raise Exception ('NOTHING_PASSED')
        return 1

    def change_window_info_str (self, Str):     
        """
        This method take as an input string of the following information:
        CHROM START END COUNT ...
        """
        try:
            window = Str.rstrip().split()
            self.chrom = window[0]
            self.start = int(window[1])
            self.end = int(window[2])
            self.read_count =  float(window[4])
            self.middle = (self.start+self.end)/2
        except IndexError:
            raise Exception("EOF")
        
        return 0

    def change_window_info_window(self, window):
        self.chrom = window.give_chrom()
        self.start = window.give_start()
        self.end = window.give_end()
        self.middle = window.give_middle()
        self.read_count = window.give_count()
        
        return 0
        

    def info (self):
        window_info = \
                    '%s\t%d\t%d\t%d\t%.2f\n' \
                    % (self.chrom, self.start, self.end, \
                       self.middle, self.read_count )
        
        return window_info

    def info_no_count (self):
        window_info = \
                    '%s\t%d\t%d\t%d' \
                    % (self.chrom, self.start, self.end, \
                       self.middle )
        
        return window_info

    
    def __eq__(self, other):
        if hasattr(other, 'chrom') and hasattr(other, 'middle'):
            return other.chrom == self.chrom and other.middle == self.middle
        else:
            return False
                
    def __lt__(self, other):
        if hasattr(other, 'chrom') and hasattr(other, 'middle'):
            if other.chrom == self.chrom:
                return other.middle > self.middle
            else:
                return other.chrom > self.chrom            

    def __gt__(self, other):
        return not (self.__lt__(other) or self.__eq__(other))

    def __le__(self, other):
        return self.__eq__(other) or self.__lt__(other)

    def __ge__(self, other):
        return self.__eq__(other) or self.__gt__(other)        
    
    def __ne__(self, other):
        return not self.__eq__(other)

        
    
