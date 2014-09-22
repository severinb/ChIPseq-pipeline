filterAdaptors.pl	preprocesses raw reads. removes adaptors and filters sequences
importSample.pl		imports a sample to the repository
prepareMicroarray.pl	prepare microarray sample for use with importsSample.pl or createSubsample.pl
createSubsample.pl	create a 'virtual' sample that is a subset of an existing sample and uses its alignments
removeSample.pl		removes a sample from the repository
submitAnnotationTask.pl	submits an annotation task for an imported sample
annotationScheduler.pl	inspects the queue and starts annoation tasks if necessary (should not be called manually)
annotateSample.pl	manually enforce annotation of a sample(bypass the queue)
extractData.pl		extract weighted quantification data
frag2wig.pl             create a wiggle file for upload to UCSC genome browser from a fragment report created by extractData.pl
frag2bed.pl             create a bed file for use with peak finders (e.g. macs)
frag2error.pl		get details on errors in read alignments
frag2sga.pl             create a sga file for use with chipseq tools from Bucher lab
createErrorPlots.pl	create plots (error profile, fraction aligned, number of hits) from a fragment report created by extractData.pl
createFullReport.pl	count reads per annotation type using waterfall-principle
fullReport2bars.R	draw a bar chart from a full report
checkGenomes.pl		verify naming conventions (expected by extractData.pl) for available genomes
checkAnnot.pl		verify consistency of .fa files and familiesAndPriorities.txt file for available annotation databases
