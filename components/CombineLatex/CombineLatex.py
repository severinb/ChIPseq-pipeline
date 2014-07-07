#!/usr/bin/env python
import component_skeleton.main
from string import *
import os
import re

def addFragment(frags, fragment):

    if fragment:
        l = open(fragment, 'r')
        frag = l.read()
        l.close()
        frags = frags + frag

    return frags

def execute(cf):
    l1 = cf.get_input('latex1')
    l2 = cf.get_input('latex2')
    l3 = cf.get_input('latex3')
    l4 = cf.get_input('latex4')
    l5 = cf.get_input('latex5')
    l6 = cf.get_input('latex6')
    l7 = cf.get_input('latex7')
    l8 = cf.get_input('latex8')
    l9 = cf.get_input('latex9')
    l10 = cf.get_input('latex10')
    l11 = cf.get_input('latex11')
    l12 = cf.get_input('latex12')
    l13 = cf.get_input('latex13')
    l14 = cf.get_input('latex14')
    l15 = cf.get_input('latex15')
    logstr = cf.get_parameter('LogStr', 'string')
    document = cf.get_output('document')
    pdfreport = cf.get_output('pdfreport')

    header = '\n'.join(['\\documentclass{article}', 
                        '\\usepackage{graphicx}',
                        '\\usepackage{listings}',
                        '\\usepackage[top=3cm, bottom=3cm, left=2cm, right=2cm]{geometry}',
                        '\\usepackage[section]{placeins}'
                        '\\usepackage[font=small]{caption}'
                        '\\lstset{language=, breaklines=true}',
                        '\\begin{document}\n'
                        ])

    footer = '\n\\end{document}'

    frags = ''

    if logstr:
        logfragments_paths = logstr.split()
        for fragment in logframents_paths:
            frags = addFragment(frags, fragment)

    lfrags = [l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15]
    for i in range(15):
        frags = addFragment(frags, lfrags[i])



    doc = open(document, 'w')
    doc.write(header + frags + footer)
    doc.close()

    os.system('pdflatex -jobname ' + re.sub('.pdf$', '', os.path.split(pdfreport)[1]) + ' -output-directory ' + os.path.split(pdfreport)[0] + ' ' + document)

    return 0


component_skeleton.main.main(execute)
