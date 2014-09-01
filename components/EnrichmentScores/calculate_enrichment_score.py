from motevo_stuff import *
from fitting_beta import *
from enrichment_score import calculate_enrichment_scores 


def arguments():
    import argparse
    parser = argparse.ArgumentParser(description='Fits beta and prior, and also run calculate the enrichment scores')
    parser.add_argument('-w', '--wm',
                    action="store", dest="WM", type=str
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


def fittingParameters(WM, trainingPool, outdir, genome):
    (siteFilename, priorFilename) = run_motevo(WM, trainingPool, outdir, genome)
    prior = extract_priors(priorFilename)
    beta = fit_beta(siteFilename, outdir, WM)
    os.system('rm %s %s' % (siteFilename, priorFilename))
    return {'prior': prior, 'beta':beta}


def calculateEnrichmetScores(WM, testPool, params, outdir, genome):
    siteFilename, priorFilename = run_motevo(WM, testPool, outdir, \
                                             genome, prior=params['prior'], \
                                             minposterior=0.0)
    motifName = os.path.basename(WM)
    enrichmentFile = os.path.join(outdir, motifName + '.enrichment_score')
    scores = calculate_enrichment_scores(siteFilename, params['beta'], enrichmentFile)
    return scores


def main():
    args = arguments()
    params = fittingParameters(args.WM, args.trainSeq, args.outdir, args.GENOME)
    enrichmentScores = calculateEnrichmetScores(args.WM, \
                                                args.testSeq, params, \
                                                args.outdir, args.GENOME)
    motifName = os.path.basename(args.WM)
    resFilename = os.path.join(args.outdir, motifName + '.results')
    with open(resFilename, 'w') as outf:
        outf.write('\t'.join([
            motifName,
            str(enrichmentScores['mean']),
            str(enrichmentScores['std']),
            str(params['beta']),
            str(params['prior'])
            ]) + '\n')
    
    
    
if __name__ == '__main__':
    main()
