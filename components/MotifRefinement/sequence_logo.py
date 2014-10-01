import os

def generate_sequence_logo(fname, resfile):
    """
    Mylogo script that we used initially is bad. Weblogo isn't much better either...
    Logos in PDF format have damaged XREF table and can not be used.
    So: First create EPS vector graphics, then convert it into PNG for website and PDF for PDF-report.
    """

    Rscript_path = "/import/bc2/home/nimwegen/GROUP/local/bin/Rscript"
    mylogo_script = "mylogo.R"

    cmd = ' '.join([
            Rscript_path,
            mylogo_script,
            fname,
            'ps',
            resfile.replace('.pdf', '')])

    print "Sequence Logo: "
    print cmd 
    os.system(cmd)

    cmd = ' '.join([
            'convert -flatten -rotate 90',
            resfile.replace('.pdf', '.ps'),
            resfile])
    os.system(cmd)

    cmd = ' '.join([
            'convert -flatten -rotate 90',
            resfile.replace('.pdf', '.ps'),
            resfile.replace('.pdf', '.png'), 
            ])
    os.system(cmd)

    cmd = ' '.join([
        'rm',
        resfile.replace('.pdf', '.ps')
        ])
    os.system(cmd)
    return 0
