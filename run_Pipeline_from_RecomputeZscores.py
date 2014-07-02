#! /usr/bin/env python

def arguments():
    import argparse
    parser = argparse.ArgumentParser(description="""Takes as input a YAML formated file,
    and initilizes the necessary files for running the CRUNCH pipeline. It then prints out
    the commands that are required for running the pipeline.
    Note that the version of the pipeline that is going to be run is not the 'whole' CRUNCH pipeline;
    namely steps that are involved with Bowtie and so on are skipped. Therefore, the current pipeline
    starts immediately after SelectPeaksMM component. """)

    parser.add_argument('-i', dest='param_file', action='store', required=True, \
                        help='The parameter file in YAML format.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    """    
    By providing a YAML file that contains the
    information required for running the pipeline, this code produces the necassary files,
    such as hosts.config and network.and files. 
    It then prints out commands for running the pipeline.    
    """

    import yaml, os
    arg = arguments()
    params = yaml.load( open(arg.param_file) )  # loading the parameters in the YAML file

    # following lines are copied, with a few changes, from the run_Pipeline.py file, written by Severin Berger
    repository_path = os.path.join(params['FMI_PATH'], 'samples')
    template_dir = params['ANDURIL_TEMPLATES_PATH']
    scripts_dir = os.path.join(os.path.split(template_dir)[0], 'scripts')

    ## Create queue dictionary
    # queue dict conatins queue types if cluster is used. Otherwise it only contains wrappers
    # (that are added in create_wrapper function) that do not submit to the queue (localrun.sh).
    queue_dict = {}          # queue_type: [queue_name, wrapper_path]
    if params['USE_CLUSTER']:
        queue_names = params['QUEUE_NAMES'].split()
        for q in queue_names:
            tmp = q.split('=')
            queue_dict[tmp[0]] = []
            queue_dict[tmp[0]].append(tmp[1])        
    
        ## Create list with all used output directories (Main/Project directory and FMI repository)
    outdirlist = []  #list of all output directories
    outdirlist.append(params['OUT_DIR'])

    ## Create Directory that contains all created files ans folders
    pipelineDir = os.path.abspath('PipelineFiles_from_RecomputeZscores')
    if not os.path.exists(pipelineDir):
        os.system('mkdir %s' % pipelineDir) 

    ## Call functions to create all files
    datafiles_csv = os.path.join(pipelineDir, 'datafiles.csv')
    csvf = open(datafiles_csv, 'a')
    csvf.write('filepath\tmode\tFMIid\tdesc\tWM\tformat\n')
    csvf.close()
    
        
