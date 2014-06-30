#!/usr/bin/env python
import component_skeleton.main
import os

def addLog(latex, logfile, caption):
    """
    Adds log file to latex script.
    """

    if logfile:
        l = open(logfile,'r')
        log = l.read()
        l.close()
        if caption:
            frag = '\n\\subsection*{' + caption + '} \n' + '\\begin{lstlisting}\n' + log + '\\end{lstlisting}\n' 
            latex = latex + frag
        else:
            frag = '\n\\subsection{}\n' + '\\begin{lstlisting}\n' + log + '\\end{lstlisting}\n'
            latex = latex + frag
        return latex
    else:
        return latex
    
def addPlot(latex, plot, caption):
    """
    Add plot to latex script.
    """
    
    if plot:
        if caption:
            frag = '\n\\begin{figure}[!h]\n\\centering\n\\graphicspath{{' + os.path.split(plot)[0] + '/}}\n\\includegraphics[width=0.4\\textwidth]{' + os.path.split(plot)[1] + '}\n\\caption{' + caption + '}\n\\end{figure}\n'
            latex = latex + frag 
        else:
            frag = '\n\\begin{figure}[!h]\n\\centering\n\\graphicspath{{' + os.path.split(plot)[0] + '/}}\n\\includegraphics[width=0.4\\textwidth]{' + os.path.split(plot)[1] + '}\n\\end{figure}\n'
            latex = latex + frag
        return latex
    else:
        return latex

def execute(cf):
    l1 = cf.get_input('log1')
    l2 = cf.get_input('log2')
    l3 = cf.get_input('log3')
    l4 = cf.get_input('log4')
    l5 = cf.get_input('log5')
    l6 = cf.get_input('log6')
    l7 = cf.get_input('log7')
    l8 = cf.get_input('log8')
    l9 = cf.get_input('log9')
    l10 = cf.get_input('log10')
    p1 = cf.get_input('plot1')
    p2 = cf.get_input('plot2')
    p3 = cf.get_input('plot3')
    p4 = cf.get_input('plot4')
    p5 = cf.get_input('plot5')
    p6 = cf.get_input('plot6')
    p7 = cf.get_input('plot7')
    p8 = cf.get_input('plot8')
    p9 = cf.get_input('plot9')
    p10 = cf.get_input('plot10')
    p11 = cf.get_input('plot11')
    p12 = cf.get_input('plot12')
    p13 = cf.get_input('plot13')
    p14 = cf.get_input('plot14')
    p15 = cf.get_input('plot15')
    cl1 = cf.get_parameter('capt_log1', 'string')
    cl2 = cf.get_parameter('capt_log2', 'string')
    cl3 = cf.get_parameter('capt_log3', 'string')
    cl4 = cf.get_parameter('capt_log4', 'string')
    cl5 = cf.get_parameter('capt_log5', 'string')
    cl6 = cf.get_parameter('capt_log6', 'string')
    cl7 = cf.get_parameter('capt_log7', 'string')
    cl8 = cf.get_parameter('capt_log8', 'string')
    cl9 = cf.get_parameter('capt_log9', 'string')
    cl10 = cf.get_parameter('capt_log10', 'string')
    cp1 = cf.get_parameter('capt_plot1', 'string')
    cp2 = cf.get_parameter('capt_plot2', 'string')
    cp3 = cf.get_parameter('capt_plot3', 'string')
    cp4 = cf.get_parameter('capt_plot4', 'string')
    cp5 = cf.get_parameter('capt_plot5', 'string')
    cp6 = cf.get_parameter('capt_plot6', 'string')
    cp7 = cf.get_parameter('capt_plot7', 'string')
    cp8 = cf.get_parameter('capt_plot8', 'string')
    cp9 = cf.get_parameter('capt_plot9', 'string')
    cp10 = cf.get_parameter('capt_plot10', 'string')
    cp11 = cf.get_parameter('capt_plot11', 'string')
    cp12 = cf.get_parameter('capt_plot12', 'string')
    cp13 = cf.get_parameter('capt_plot13', 'string')
    cp14 = cf.get_parameter('capt_plot14', 'string')
    cp15 = cf.get_parameter('capt_plot15', 'string')
    title = cf.get_parameter('sectionTitle', 'string')
    fragment = cf.get_output('fragment')

    latex = ''

    """use graphics package, and maybe also listing package for log. make title for each component log, and caption for each plot.
       \begin{figure}[h or h!] \caption{cp6} \includegraphics[scale=0.7]{p6} \end{figure}
    """
    latex = latex + '\\section*{' + title + '}\n'

    logfiles = [l1,l2,l3,l4,l5,l6,l7,l8,l9,l10]
    logcaptions = [cl1,cl2,cl3,cl4,cl5,cl6,cl7,cl8,cl9,cl10]
    plots = [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15]
    plotcaptions = [cp1,cp2,cp3,cp4,cp5,cp6,cp7,cp8,cp9,cp10,cp11,cp12,cp13,cp14,cp15]

    for i in xrange(len(logfiles)):
        latex = addLog(latex, logfiles[i], logcaptions[i])
    for i in xrange(len(plots)):
        latex = addPlot(latex, plots[i], plotcaptions[i])


    outf = open(fragment, 'w')
    outf.write(latex)
    outf.close()

    return 0

component_skeleton.main.main(execute)
