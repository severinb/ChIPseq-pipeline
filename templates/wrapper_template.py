#!PYTHON_PATH

import drmaa
import os, sys
import datetime


if __name__ == '__main__':

    """
    This is a new wrapper for anduril remote execution that replaces the bash wrapper from the anduril guys.
    The purpose of this is to be independent of SGE. (For this: how specific are the parameters given in JOB_PARAM (s.createJobTemplate()) for SGE?)
    It is a wrapper thus following things have to be filled in:
       -python to use for execution in header line
       -project_leader
       -queue_type (fs_long fs_long_hm)

    As far as I know with drmaa one is not able to define the program (python) with which it should run a script. --NO IT IS POSSIBLE, NO PROBLEM THERE
    Thus I could make all COMPONENT_WRAPPER_SCRIPT scripts executable with '/usr/bin/env python' and just add the good program (python) to the environment in this wrapper...
    """

    ##The command that gets executed by anduril-remote (using the eval command) looks like this:
    ##./wrapper.py /import/bc2/home/nimwegen/GROUP/local/bin/python /import/bc2/home/nimwegen/...../PipeLine140113/COMPONENTS/transform/transform_anduril.py /import/bc2/home/nimwegen/..../IggrabTNFAaBG/OUTPUT/NFKB_4-trans/_command
    ##wrapper.py PYTHON COMPONENT_WRAPPER_SCRIPT COMMAND_FILE

    #this could theoretically also be a perl path if I would use perl wrappers to run components...
    pythonpath = sys.argv[1]
    component_wrapper = sys.argv[2] #this one doesn't have to be executable
    commandfile = sys.argv[3]

    project_leader = 'PROJECTLEADER'
    queue_type = 'QUEUE_TYPE' #'fs_long'


    QueueFilesDir =  os.path.split(commandfile)[0]
    stderrpath = os.path.join(QueueFilesDir, 'job.stderr')
    stdoutpath = os.path.join(QueueFilesDir, 'job.stdout')


    #extract name of job
    f = open(commandfile)
    l = f.readline()
    f.close()
    jobname = l.strip().split('=')[1]


    # -w n is something that I added because there were these comlib errors (comlib error: got read error)
    # in JOB_params I added -shell yes so that the execution of the component together with the profiling scripts works.
    JOB_PARAM = '-q %s -P %s -e %s -o %s -j n -w n -N %s -cwd -V -b y -shell yes ADDITIONAL_PARAMETERS ' %(queue_type, project_leader, stderrpath, stdoutpath, jobname)
    timestats = '-f "# Real time                       : %E\\n# User time                       : %U\\n# System time                     : %S\\n# Percent of CPU this job got     : %P"'
    
    s = drmaa.Session()
    s.initialize()

    jt = s.createJobTemplate()
    jt.nativeSpecification = JOB_PARAM
    jt.remoteCommand = '/usr/bin/time'
    # monmem has a high number (1000TB) so that it never gets killed. Every 10 min it prints memory usage 
    jt.args = [timestats, pythonpath, component_wrapper, commandfile, '& \/import/bc2/soft/bin/monmem $! 1000000000 600 | grep -v "Memory threshold\|Start " 1>&2']


    T1 = datetime.datetime.now()
    print T1
    sys.stdout.flush()

    jobid = s.runJob(jt)
    print 'submitted', jobid
    sys.stdout.flush()

    retval = s.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)
    print 'Job: ' + str(retval.jobId) + ' finished with status ' + str(retval.hasExited) + ', ' + str(retval.exitStatus)
    sys.stdout.flush()

    print 'stdout:'
    sys.stdout.flush()
    os.system('cat %s' %stdoutpath)

    print 'stderr:'
    sys.stdout.flush()
    os.system('cat %s' %stderrpath)

    T2 = datetime.datetime.now()
    print T2
    print 'Running Time for %s: %s' %(jobname, str(T2-T1))
    sys.stdout.flush()
    

    s.deleteJobTemplate(jt)
    s.exit()

