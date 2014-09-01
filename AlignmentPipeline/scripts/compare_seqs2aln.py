
import sys, getopt, os, shutil, re

###################

def main():
    # read arguments
    seq_dir = 'SEQS'
    aln_dir = 'ALIGNMENT'
    missing_dir = 'MISSING_SEQS'
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "s:a:")
    except getopt.GetoptError, err:
        print str(err)
        exit()
        
    for o, a in opts:
        if o == '-s':
            seq_dir = a
        elif o == '-a':
            aln_dir = a

    # check if missing_dir exists and create it otherwise

    if not os.path.isdir(missing_dir):
        os.mkdir(missing_dir)
                
    # get seq file names
    seq_files = os.listdir(seq_dir)
    aln_files = os.listdir(aln_dir)
    
    fout = open("%s/filelist" % missing_dir, 'w')
    for seq_file in seq_files:
        if not re.search('\.fna', seq_file):
            continue
        
        aln_file = seq_file.rstrip('fna') + 'aln'
        if aln_file not in aln_files:
            shutil.copy2(seq_dir+'/'+seq_file, missing_dir+'/'+seq_file)
            fout.write("%s\n" % seq_file)
    fout.close()
                
###################
                
if __name__ == '__main__':
    main()
