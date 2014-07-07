""" Remove dirs older than 7 days from scratch directory. """

import argparse
import os
import re
import sys
import time
import yaml

###################

def main():
    """ Main """
    parser = argparse.ArgumentParser(       
        description='Remove dirs older than 7 days from scratch directory.')

    parser.add_argument('--conf',
                        help='config file',
                        dest='config_file',
                        default='/import/wnz/home/crunch/PipeLine/web/config/crunch.conf')
    
    parser.add_argument('--log',
                        help='log file',
                        dest='log_file',
                        default='/import/wnz/home/crunch/log/sentinel.log')

    args = parser.parse_args()

    config = yaml.load(open(args.config_file))
    scratch_dir = config['scratch_dir']
    
    # get list of dirs from scratch
    dirs_in_scratch = [os.path.join(scratch_dir, x) for x in os.listdir(scratch_dir)
                       if os.path.isdir(os.path.join(scratch_dir, x))]
    # get list of data_xxxx dirs
    data_dirs = [x for x in dirs_in_scratch if re.search("data_\S{5}", x)]

    # get list of dirs older than 8 days
    for data_dir in data_dirs:
        if ((time.time() - os.path.getmtime(data_dir)) / (3600 * 24)) > 16:
            os.system("rm -rf %s" % data_dir)
            log = open(args.log_file, 'a')
            log.write("Deleted dir %s\n" % data_dir)
            log.close()



###################



###################

if __name__ == '__main__':
    main()
