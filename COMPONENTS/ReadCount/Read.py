

class Read:

    def __init__ (self, BED_line, is_weight_given):
        tmp = BED_line.split()
        self.chrom = tmp[0]
        self.start = int(tmp[1])
        self.end = int(tmp[2])
        self.strand = tmp[5].rstrip()
        if is_weight_given:
            self.weight = float(tmp[4])
        else:
            self.weight = 1.0

    def give_chrom(self):
        return self.chrom

    def give_strand():
        return self.strand    

    def give_start(self):
        return self.start

    def give_end(self):
        return self.end

    def give_weight(self):
        return self.weight

    def give_length(self):
        return (self.end - self.start)

    def change_info(self, BED_line, is_weight_given):
        tmp = BED_line.split()
        self.chrom = tmp[0]
        self.start = int(tmp[1])
        self.end = int(tmp[2])
        if is_weight_given:
            self.weight = float(tmp[4])
        else:
            self.weight = 1.0

        return 1

    def info(self):
        read_info = \
                    '%s\t%d\t%d\t%d\n' \
                    % (self.chrom, self.start, self.end, \
                       (self.weight))
        
        return read_info
    
