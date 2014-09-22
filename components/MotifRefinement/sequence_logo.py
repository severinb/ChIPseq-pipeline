import os

def generate_sequence_logo(fname, resfile, weblogo_path, desc=' ', errorbars='NO'):
    cmd = ' '.join([
        weblogo_path,
        '--sequence-type dna',
        '--size large',
        '--color-scheme classic',
        '--errorbars %s' % errorbars,
        '--fineprint \'%s \'' % desc,
        '--format png', 
        '< %s >' % fname,
        resfile.replace('.pdf', '.png')
        ])
    # print "Sequence Logo: "
    # print cmd 
    os.system(cmd)
    cmd = ' '.join([
        'convert',
        resfile.replace('.pdf', '.png'),
        resfile, 
        ])
    os.system(cmd)
    return 0
