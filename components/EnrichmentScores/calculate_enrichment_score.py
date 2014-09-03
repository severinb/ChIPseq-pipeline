from motevo_stuff import *
from fitting_beta import *
from Bio import SeqIO
from enrichment_score import calculate_enrichment_scores
from concatenate_motifs import concatenate


def arguments():
    import argparse
    parser = argparse.ArgumentParser(description='Fits beta and prior, and also run calculate the enrichment scores')
    parser.add_argument('-w', '--wm',
                    action="store", dest="WM", type=str, nargs='+'
                    )    
    parser.add_argument('-t', '--trainseq',
                    action="store", dest="trainSeq", type=str
                    )
    parser.add_argument('-s', '--testseq',
                        action="store", dest="testSeq", type=str
                        )
    parser.add_argument('-o', '--outdir',
                        action="store", dest="outdir", type=str
                        )
    parser.add_argument('-g', '--genome',
                        action="store", dest="GENOME", type=str
                        )    
    results = parser.parse_args()
    return results


def cleanup(infiles):
    """
    To delete the files in the scratch directory
    """
    cmd = 'rm -f %s' % (' '.join(infiles))
    os.system(cmd)
    return 0    


def fittingParameters(WM, trainingPool, trainingLength, outdir, genome):
    siteFilename, priorFilename, paramFilename, motifName = run_motevo(WM, trainingPool, outdir, genome, priorFile=None, minposterior=0.0001)
    prior = extract_priors(priorFilename)
    beta = fit_beta(siteFilename, outdir, WM, trainingLength)
    cleanup([siteFilename, paramFilename])
    return {'prior': prior, 'beta':beta}, priorFilename


def calculateEnrichmetScores(WM, testPool, testLength, params, outdir, genome, priorFile):
    siteFilename, priorFilename, paramFilename, \
                  WM = run_motevo(WM, testPool, outdir, \
                                  genome, priorFile=priorFile, \
                                  minposterior=0.0001)
    motifName = os.path.basename(WM)
    scores = calculate_enrichment_scores(siteFilename, params['beta'], testLength, os.path.join(outdir, '%s.enrichment_score' % motifName))
    cleanup([siteFilename, priorFilename, paramFilename, priorFile])    
    return scores, motifName


def lengthOfWM(WMfiles):
    minLength = np.inf
    for WMfile in WMfiles:
        length = 0
        with open(WMfile) as inf:
            for line in inf:
                if re.search('^\d+\s+[(\.)0-9]+\s+[(\.)0-9]+\s+', line):
                    length += 1
        if length < minLength:
            minLength = length
    return minLength    


def lengthOfSequences(trainingPool, testPool, WMfile):
    wmLength = lengthOfWM(WMfile)
    trainingLength, testLength = {}, {}
    with open(trainingPool) as inf:
        for record in SeqIO.parse(inf, 'fasta'):
            trainingLength[re.sub('^>', '', record.id)] = (len(record.seq) - wmLength)*2
    with open(testPool) as inf:
        for record in SeqIO.parse(inf, 'fasta'):
            testLength[re.sub('^>', '', record.id)] = (len(record.seq) - wmLength)*2
    return trainingLength, testLength


def main():
    args = arguments()
    trainingLength, testLength = lengthOfSequences(args.trainSeq, args.testSeq, args.WM)
    params, priorFile = fittingParameters(args.WM, args.trainSeq, trainingLength, \
                                          args.outdir, args.GENOME)
    enrichmentScores, motifName = calculateEnrichmetScores(args.WM, \
                                                args.testSeq, testLength, params, \
                                                args.outdir, args.GENOME, priorFile)
    resFilename = os.path.join(args.outdir, motifName + '.results')
    with open(resFilename, 'w') as outf:
        outf.write('\t'.join([
            '%&%'.join(args.WM),
            str(enrichmentScores['mean']),
            str(enrichmentScores['std']),
            str(params['beta']),
            str(params['prior'])
            ]) + '\n')
    
    
    
if __name__ == '__main__':
    main()
