import os

def generate_sequence_logo(fname, resfile, weblogo_path, desc=' ', errorbars='NO'):
    cmd = ' '.join([
        weblogo_path,
        '--sequence-type dna',
        '--size large',
        '--color-scheme classic',
        '--errorbars %s' % errorbars,
        '--fineprint %s' % desc,
        '< %s >' % fname,
        resfile
        ])
    print "Sequence Logo: "
    print cmd 
    os.system(cmd)
    return 0
