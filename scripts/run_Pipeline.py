#!/import/bc2/home/nimwegen/GROUP/local/bin/python
#/usr/bin/env python

import yaml
import sys
import os
import re
import subprocess
import shutil
from string import *

def create_csvFiles(IPkey, BGkey, flag):
    """create csv file with fg and bg files and return filepath. For samples that were already processed up to shifting, file is given to another csv file, so that processing isn't done again.
    """

    datafiles_csv = os.path.join(pipelineDir, 'datafiles.csv')
    csvf = open(datafiles_csv, 'a')
    #csvf.write('filepath\tmode\tFMIid\tdesc\tWM\n')

    existing_csv = os.path.join(pipelineDir, 'existingSamples.csv')
    existingcsv = open(existing_csv, 'a')
    #existingcsv.write('filepath\tmode\tWM\n')

    ##Process foreground files
    TFdict = params[IPkey]

    if TFdict:
        ##Loop over all TFs in IP_FILES. Write everything into datafiles.csv (path, modem, FMIid, description, WMpath) or existingSamples.csv
        for TF in TFdict:
            inlist = TFdict[TF].split()
            try:
                globalTFdict[TF] += inlist
            except Exception:
                globalTFdict[TF] = inlist

            #Do everything with descriptions: check whether they are there.
            if flag == 'fastq':
                try:
                    indescslist = params['IP_FASTQ_DESCRIPTION'][TF]
                except (TypeError, KeyError):
                    indescslist = []
            elif flag == 'fasta':
                try:
                    indescslist = params['IP_FASTA_DESCRIPTION'][TF]
                except (TypeError, KeyError):
                    indescslist = []
            else:
                indescslist = []


            for i, infile in enumerate(inlist):
                try:
                    desc = indescslist[i]
                    if not desc:
                        desc = "No sample description given"
                except (TypeError, IndexError):
                    desc = "No sample description given"
                   
                #get rid of all file extensions
                fileroot = os.path.split(infile)[1]
                fileFMIid = re.sub('\.gz$','',fileroot)
                fileFMIid = re.sub('\.fastq$','',fileFMIid)
                fileFMIid = re.sub('\.fa$','',fileFMIid)
                fileFMIid = re.sub('\.bedweight$','',fileFMIid)
                fileFMIid = re.sub('\.bed$','',fileFMIid)
                try:
                    rtag = params['RENAME_TAG']
                    fileFMIid = rtag + '_' + fileFMIid
                except (TypeError, KeyError):
                    pass

                outdir = os.path.join(repository_path, fileFMIid)
                outdirlist.append(outdir)

                if os.path.isdir(outdir):
                    print fileFMIid
                    print 'Warning: Some FG samples are already in the repository. Check all created pipeline files, especially datafiles.csv and existingSamples.csv and check repository/samples. Check whether fileRoot.bedweight.shifted.gz is present in existingSamples.csv. Otherwise edit and adapt exisitngSamples.csv and datafiles.csv.\n'
                    matched = False
                    for (p,d,f) in os.walk(outdir):
                        for fi in f:
                            #matchShiftfile = re.search(os.path.join(outdir,'\S+\d+-fraglen/out_dir',fileFMIid + '\S*.bedweight.shifted.gz$'), os.path.join(p,fi))
                            matchShiftfile = re.search(os.path.join(outdir,'\S+\d+-fraglen/out_dir/outfile.gz'), os.path.join(p,fi))
                            if matchShiftfile:
                                matched = True
                                existingcsv.write(matchShiftfile.group() + '\t' + TF + '\t' + wm + '\n')
                    if not matched: #Directory is present in FMI repository but shifted bedweight file not! Treat it as if sample wasn't present in FMI repository at all.
                        csvf.write(infile + '\t' + TF + '\t' + fileFMIid + '\t' + desc + '\t' + flag + '\n')
                else:
                    csvf.write(infile + '\t' + TF + '\t' + fileFMIid + '\t' + desc + '\t' + flag + '\n')


    ##Process background files
    bgfiles = params[BGkey]

    if bgfiles:
        bglist = bgfiles.split()
        global globalBGlist
        globalBGlist += bglist

        if flag == 'fastq':
            try:
                bgdescslist = params['BG_FASTQ_DESCRIPTION']
            except (TypeError, KeyError):
                bgdescslist = []
        elif flag == 'fasta':
            try:
                bgdescslist = params['BG_FASTA_DESCRIPTION']
            except (TypeError, KeyError):
                bgdescslist = []
        else:
            bgdescslist = []

        for i, bgfile in enumerate(bglist):
            try:
                desc = bgdescslist[i]
                if not desc:
                    desc = "No sample description given"
            except (TypeError, IndexError):
                desc = "No sample description given"

            #get rid of all file extensions
            fileroot = os.path.split(bgfile)[1]
            fileFMIid = re.sub('\.gz$','',fileroot)
            fileFMIid = re.sub('\.fastq$','',fileFMIid)
            fileFMIid = re.sub('\.fa$','',fileFMIid)
            fileFMIid = re.sub('\.bedweight$','',fileFMIid)
            fileFMIid = re.sub('\.bed$','',fileFMIid)
            try:
                rtag = params['RENAME_TAG']
                fileFMIid = rtag + '_' + fileFMIid
            except (TypeError, KeyError):
                pass

            outdir = os.path.join(repository_path, fileFMIid)
            outdirlist.append(outdir)
            
            if os.path.isdir(outdir):
                print fileFMIid
                print 'Warning: Some BG samples are already in the repository. Check all created pipeline files, especially datafiles.csv and existingSamples.csv and check repository/samples. Check whether fileRoot.bedweight.shifted.gz is present in existingSamples.csv. Otherwise edit and adapt exisitngSamples.csv and datafiles.csv.\n'
                matched = False
                for (p,d,f) in os.walk(outdir):
                    for fi in f:
                        #matchShiftfile = re.search(os.path.join(outdir,'\S+\d+-\S*fraglen/out_dir',fileFMIid + '\S*.bedweight.shifted.gz$'), os.path.join(p,fi))
                        matchShiftfile = re.search(os.path.join(outdir,'\S+\d+-\S*fraglen/out_dir/outfile.gz'), os.path.join(p,fi))
                        if matchShiftfile:
                            matched = True
                            existingcsv.write(matchShiftfile.group() + '\tBG\n')
                if not matched:
                    csvf.write(bgfile + '\tBG\t' + fileFMIid + '\t' + desc + '\t' + flag + '\n')
            else:
                csvf.write(bgfile + '\tBG\t' + fileFMIid + '\t'+ desc + '\t' + flag + '\n')

    csvf.close()
    existingcsv.close()
 
    return (os.path.abspath(datafiles_csv), os.path.abspath(existing_csv))


def create_network():
    """replace CSV_DATAFILE_PATH, ADAPTOR_3, ANNOTYPE, GENOME, IP_DESCFOLDER, BG_DESCFOLDER, FG_WINSIZE, BG_WINSIZE, STEPSIZE, CHROMINFOFILE_PATH, USER etc...
       in network.and template
       Return path of network.and file
    """

    if params['USE_CLUSTER']:
        ntemp = open(os.path.join(template_dir, 'network_template_cluster.and'),'r')
        #ntemp = open('network_template_cluster.v4.and','r')
        ntext = ntemp.read()  #network text
        ntemp.close()
        
    else:
        ntemp = open(os.path.join(template_dir, 'network_template.and'),'r')
        ntext = ntemp.read()  #network text
        ntemp.close()            

    if not params['DO_MOTIF_FINDING']:
        ntext = re.sub('PIPEFLAG', 'JUSTPEAKS', ntext)


    ntext = re.sub('CSV_DATAFILE_PATH', csvfile_path, ntext)

    ntext = re.sub('CSV_EXISTING_SAMPLES', existingcsv_path, ntext)

    ntext = re.sub('FMI_PATH', params['FMI_PATH'], ntext)

    ntext = re.sub('PERL_PATH', params['PERL_PATH'], ntext)

    ntext = re.sub('PYTHON_PATH', params['PYTHON_PATH'], ntext)

    ntext = re.sub('ALIGN_PIPE_PATH', params['ALIGN_PIPE_PATH'], ntext)

    ntext = re.sub('ANNOTYPE', params['ANNOTATION_TYPE'], ntext)

    ntext = re.sub('GENOME', params['GENOME'], ntext)

    ntext = re.sub('GENDIR_PATH', params['GENDIR_PATH'], ntext)

    ntext = re.sub('FG_WINSIZE', str(params['WINDOW']), ntext)

    ntext = re.sub('BG_WINSIZE', str(params['BACKGROUND_WINDOW']), ntext)

    ntext = re.sub('STEPSIZE', str(params['STEP']), ntext)

    ntext = re.sub('CHROMINFOFILE_PATH', params['CHROMOSOME_INFO'], ntext)

    ntext = re.sub('REPEAT_PATH', params['REPEAT_PATH'], ntext)

    ntext = re.sub('BEDTOOLS_PATH', params['BEDTOOLS_PATH'], ntext)

    ntext = re.sub('PHYLOGIBBS_PATH', params['PHYLOGIBBS_PATH'], ntext)

    ntext = re.sub('MYLOGO_PATH', params['MYLOGO_PATH'], ntext)

    ntext = re.sub('WEBLOGO_PATH', params['WEBLOGO_PATH'], ntext)

    ntext = re.sub('MOTEVO_PATH', params['MOTEVO_PATH'], ntext)

    ntext = re.sub('MOTEVOUFE_PATH', params['MOTEVOUFE_PATH'], ntext)

    ntext = re.sub('R_LIBS_USER_PATH', params['R_LIBS_USER'], ntext)

    ntext = re.sub('FDRVALUE', str(params['FDR']), ntext)

    ntext = re.sub('THREAD_NUMBER', str(params['THREAD_NUMBER']), ntext)

    try:
        project_l = params['PROJECTLEADER']
        ntext = re.sub('PROJECTLEADER', project_l, ntext)
    except (TypeError, KeyError):
        ntext = re.sub('PROJECTLEADER', '', ntext)

    ntext = re.sub('ANNOTATION_FILE', params['ANNOTATION_FILE'], ntext)

    ntext = re.sub('SORT_TMPDIR', params['SORT_TMPDIR'], ntext)

    ntext = re.sub('WMLIBRARY', params['WMLIBRARY'], ntext)

    ntext = re.sub('FMI_OUTPUT_DIR', FMI_output_dir, ntext)

    try:   
        ntext = re.sub('ADAPTERS_3', params['ADAPTOR'], ntext)
    except (TypeError, KeyError):
        ntext = re.sub('ADAPTERS_3', '', ntext)

    try:
        user = params['USER_NAME']
        ntext = re.sub('USER_NAME', user, ntext)
    except (TypeError, KeyError):
        ntext = re.sub('USER_NAME', '', ntext)

    net = open('network.and', 'w')
    net.write(ntext)
    net.close()
 
    return os.path.abspath('network.and')


def create_hosts_conf():
    """
    This funtion creates the hosts configuration file. Here the different queues and queue types are defined, so that anduril knows how to use the queue.
    This hosts.conf file also allows to use component specific execution directories.
    So for each replicate - by looping over the queue_dict - a set of hosts is defined with a particular execution directory (in FMI repository).
    """

    host_template = '\n'.join(["HostID = HOSTID",
                               "HostName = 127.0.0.1",
                               "RemoteExecutionDirectory = OUT_DIR",
                               "RemoteExecute = REMOTE_EXECUTE ${COMMAND}",
                               "CopyLocalRemote = rsync -a --exclude='**/.*' ${LOCAL_PATH} ${REMOTE_DIR}/",
                               "CopyRemoteLocal = rsync -a ${REMOTE_PATH} ${LOCAL_DIR}/",
                               "IsSharedFileSystem = SHARED",
                               "Wrapper = WRAPPER_PATH",
                               "PathMapping = /=/",
                               "\n"])

    hosts_conf = os.path.join(pipelineDir, 'hosts.conf')
    h = open(hosts_conf, 'w')

    #creating hosts config for non FMI components
    for queue in queue_dict:
        text = host_template

        text = re.sub('HOSTID', queue, text)
        text = re.sub('OUT_DIR', os.path.abspath(params['OUT_DIR']), text)
        text = re.sub('REMOTE_EXECUTE', os.path.join(scripts_dir, 'localrun.sh'), text)
        text = re.sub('SHARED', 'true', text)
        text = re.sub('WRAPPER_PATH', queue_dict[queue][1], text)

        h.write(text)


    h.close()


    return os.path.abspath(hosts_conf)


def create_wrappers():
    """
    This function produces wrappers for queue submission. For each queue type (long, verylong etc...) one wrapper is created.
    Additionally a wrapper is created that allows execution by anduril remote execution but without using the queue. This allows to specify the execution directory individually.
    If no queue should be used, just this noQueue wrapper is created. The noQueue wrapper is the localrun.sh script.
    """

    if params['USE_CLUSTER']:
        #all other queues
        for queue in queue_dict:
            temp = open(os.path.join(template_dir, 'wrapper_template.py'))
            text = temp.read()
            temp.close()

            text = re.sub('PYTHON_PATH', params['PYTHON_PATH'], text)
            try:
                project_l = params['PROJECTLEADER']
                text = re.sub('PROJECTLEADER', project_l, text)
            except (TypeError, KeyError):
                text = re.sub('PROJECTLEADER', '', text)

            text = re.sub('QUEUE_TYPE', queue_dict[queue][0], text)

            if queue.startswith('sjpn'):
                text = re.sub('ADDITIONAL_PARAMETERS', '-l sjpn=1', text) #call qsub with single job per node option
            elif queue.startswith('multi'):
                text = re.sub('ADDITIONAL_PARAMETERS', '-pe MT %i -l mem_total=100G' %params['THREAD_NUMBER'], text) #qsub for multicores with a lot of memory (needed for the mapping stuff).
            else:
                text = re.sub('ADDITIONAL_PARAMETERS', '', text)

            wrapper_filename = os.path.join(pipelineDir, 'wrapper_' + queue + '.py')
            o = open(wrapper_filename, 'w')
            o.write(text)
            o.close()
            os.chmod(wrapper_filename, 0755)

            queue_dict[queue].append(os.path.abspath(wrapper_filename))


    #noQueue
    wrapper_filename = os.path.join(pipelineDir, 'wrapper_noQueue')
    shutil.copyfile(os.path.join(scripts_dir, 'localrun.sh'), wrapper_filename)
    os.chmod(wrapper_filename, 0755)
    
    queue_dict['noQueue'] = ['noQueue', os.path.abspath(wrapper_filename)]


    return queue_dict


def start_with_bedweight():
    """
    Add existing bedweight files that are not in the FMI repository to existingSamples.csv.
    """
    existingcsv = open(existingcsv_path, 'a')

    ##first foreground samples:
    TFdict = params['IP_SHIFTEDBED_FILES']

    if TFdict:
        for TF in TFdict:

                
            if TFdict[TF]:
                inlist = TFdict[TF].split()
            else:
                break

            try:
                globalTFdict[TF] += inlist
            except KeyError:
                globalTFdict[TF] = inlist
                
            for bedweightFile in inlist:
                existingcsv.write(bedweightFile + '\t' + TF + '\n')

    ##Now background samples
    bgfiles = params['BG_SHIFTEDBED_FILES']

    if bgfiles:
        bglist = bgfiles.split()

        global globalBGlist
        globalBGlist += bglist

        for bedweightFile in bglist:
            existingcsv.write(bedweightFile + '\tBG\n')


if __name__ == '__main__':
    """
    Produces hosts.conf, network.and, datafiles.csv, wrapper files and qsub template and gives pipeline run command.

    Example config yaml file: params.v4.yaml
    """

    if not len(sys.argv) == 2:
        print '\nUsage: ./run_Pipeline.py config_file.yaml\n'
        sys.exit(0)


    # # load directory with all output html templates into jinja2:
    # pipeline_dir = os.path.split(os.path.split(os.path.realpath(__file__))[0])[0]
    # temp_env = jinja2.Environment(loader=jinja2.FileSystemLoader('%s' %(os.path.join(pipeline_dir, 'templates'))),
    #                               trim_blocks=True)


    ##read yaml config file
    configfile = sys.argv[1]
    cf = open(configfile)
    params = yaml.load(cf)
    cf.close()

    repository_path = os.path.join(params['FMI_PATH'], 'samples')
    template_dir = params['ANDURIL_TEMPLATES_PATH']
    scripts_dir = os.path.join(os.path.split(template_dir)[0], 'scripts')

    ##global variables that count foreground and background files
    globalTFdict = {}
    globalBGlist = []

    ##Create queue dictionary
    #queue dict conatins queue types if cluster is used. Otherwise it only contains wrappers (that are added in create_wrapper function) that do not submit to the queue (localrun.sh).
    queue_dict = {} #queue_type: [queue_name, wrapper_path]
    if params['USE_CLUSTER']:
        queue_names = params['QUEUE_NAMES'].split()
        for q in queue_names:
            tmp = q.split('=')
            queue_dict[tmp[0]] = []
            queue_dict[tmp[0]].append(tmp[1])


    ##Create list with all used output directories (Main/Project directory and FMI repository)
    outdirlist = []  #list of all output directories
    outdirlist.append(params['OUT_DIR'])

    ##Create Directory that contains all created files ans folders
    pipelineDir = 'PipelineFiles'
    if not os.path.exists(pipelineDir):
        os.system('mkdir %s' %pipelineDir)

    ##Create Main/Project output directory
    if not os.path.exists(params['OUT_DIR']):
        os.system('mkdir %s' %params['OUT_DIR'])

    ##Create FMI output directeory
    FMI_output_dir = os.path.join(os.getcwd(), 'OUTPUT_FMI')
    if not os.path.exists(FMI_output_dir):
        os.system('mkdir %s' %FMI_output_dir)
        os.system('mkdir %s/tmp' %FMI_output_dir)
    repository_path = FMI_output_dir

    ##Call functions to create all files
    datafiles_csv = os.path.join(pipelineDir, 'datafiles.csv')
    csvf = open(datafiles_csv, 'a')
    csvf.write('filepath\tmode\tFMIid\tdesc\tformat\n')
    csvf.close()

    existing_csv = os.path.join(pipelineDir, 'existingSamples.csv')
    existingcsv = open(existing_csv, 'a')
    existingcsv.write('filepath\tmode\n')
    existingcsv.close()

    (csvfile_path, existingcsv_path) = create_csvFiles('IP_FASTQ_FILES', 'BG_FASTQ_FILES', 'fastq')
    (csvfile_path, existingcsv_path) = create_csvFiles('IP_FASTA_FILES', 'BG_FASTA_FILES', 'fasta')
    (csvfile_path, existingcsv_path) = create_csvFiles('IP_BED_FILES', 'BG_BED_FILES', 'bed')
    
    start_with_bedweight() #shifted bedweight

    network_path = create_network()

    queue_dict = create_wrappers()

    hostsconf_path = create_hosts_conf()

    ##Print anduril command for running pipeline and also write it to a file
    command = ' '.join(['anduril run',
                        network_path,
                        '-c %s' %params['ANDURIL_COMPONENTS_PATH'],
                        '-d %s' %params['OUT_DIR'],
                        '--hosts %s' %os.path.join(pipelineDir, 'hosts.conf'),
                        '--perl-exec %s' %params['PERL_PATH'],
                        '--python-exec %s' %params['PYTHON_PATH'],
                        '--threads 8'])

    help_text = '\n'.join(['\n\nFurther information:',
                           '\t- Mappings for individual samples are here:',
                           '\t\t%s'%repository_path,
                           '\t- Log files for individual component instances are here: ./log/',
                           '\t- Log for complete run is here: ./log/_global',
                           '\t- It is advisable to run the second command above with nohup or screen',
                           '\t- To create report pages for looking at in the web browser run after Crunch has finished:',
                           '\t\t%s/make_output.py %s' %(scripts_dir, configfile),
                           '\t- Control the number of concurrent components using:',
                           '\t\t--threads NUMBER\n'])

    source_command = 'source %s' %os.path.join(os.path.split(template_dir)[0], 'config/crunch.env')

    oc = open('andurilCOMMAND','w')
    oc.write('1. Run:\n%s\n\n2. Run:\n' %source_command)
    oc.write(command)
    oc.write(help_text)
    oc.close()

    os.system('cat andurilCOMMAND')
