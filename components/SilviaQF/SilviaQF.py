#!/usr/bin/env python
import component_skeleton.main
import subprocess
import re, os

def execute(cf):
    in_file = cf.get_input("in_file")
    out_file = cf.get_output("out_file")
    plot_before = cf.get_output("plot_before")
    plot_after = cf.get_output("plot_after")
    log_file = cf.get_output("log_file")
    pythonPATH = cf.get_parameter("pythonPATH", "string")

    proc = subprocess.Popen (pythonPATH + ' fast_quality_filtering_with_plot_2.py %s %s %s %s %s' % (in_file, out_file, plot_before, plot_after, log_file),
                             stdout=subprocess.PIPE,
                             stderr= subprocess.PIPE,
                             shell=True
                             )

    stdout_value, stderr_value = proc.communicate()

    print stdout_value
    print stderr_value

    if proc.poll() > 0:
        print '\tstderr:', repr(stderr_value.rstrip())
        return -1

    return 0


component_skeleton.main.main(execute)
