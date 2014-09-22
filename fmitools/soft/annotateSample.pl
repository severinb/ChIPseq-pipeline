#!/usr/bin/env perl

use strict;
use warnings;
use FindBin '$RealBin';
use lib "$RealBin";
use File::stat;
use Time::localtime;
use DeepSeqTools;

# exit(0); annotation worked
# exit(1); annotation did not work

# absolute path to the repository
my $repDir             = "$RealBin/..";
#my $samplesDir         = "$repDir/samples";
my $seqBaseFileName    = "seqs";
my $groupDbDir         = getConfig('DbDir');
my $genomesDir         = "$groupDbDir/genomes";
my $seqAnnotDir        = "$groupDbDir/seqAnnotDB";
my $alnSoftDir         = "$repDir/soft/aln";
my $auxPrefix          = "AUX";
my $reposQueueFile     = "annotationQueue.tab";
my $errorLogDirName    = "errorLogs";
my $fullReportsDirName = "fullReports";
my $genomeDbDir        = "genomes";
my $annotDbDir         = "seqAnnotDB";
#my $tempRootDir        = "$repDir/tmp";
my $deleteTempTime     = getConfig('tmpDelDays');    # in days

if (@ARGV != 4)
{
    print STDERR "\nUSAGE: ./annotateSample.pl sampleID annotationOptions\n\n";
    print STDERR "   example: ./annotateSample.pl GSE1232 hg18-hgV01-aln1-AUX_1_3_8\n\n";
    exit(0);
}

# check if the annotation database contains all the files necessary to start annotation
# if yes submit, otherwise reject

# 1.11.13 modified by severin. Added samplesDir as argument
# 18.09.14 modified by severin. Added thread_number as argument. Passed to aln2.3.pl in all of its instances.
my ($sampleId, $annotOptions, $thread_number, $samplesDir) = ($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3]);
my $childSampleId = $sampleId;

my $tempRootDir = "$samplesDir/tmp";

# check if $sampleId is a virtual sample. if yes, replace it by the mother sampleId
if (-l "$samplesDir/$sampleId/mappings")
{

    # mappings directory is a symbolic link, follow link

    my $fileStatOut = `file $samplesDir/$sampleId/mappings`;

    if ($fileStatOut =~ /symbolic link.+\.\.\/(\S+)\/mappings/)
    {
        my $motherSampleId = $1;

        # replacing the id of the current sample by the mother sample
        $sampleId = $motherSampleId;
    }
    else
    {
        die "error: file operation does not work in order to resolve symbolic links\n";
    }
}

my ($genomeVersion, $annotVersion, $alnVersion, $auxString, $auxFileListRef) = parseAnnotationOptions($annotOptions);
my @auxFileList = @{$auxFileListRef};

my $genomeMappingDir = "$samplesDir/$sampleId/mappings/$genomeDbDir\_$genomeVersion\-$alnVersion";
my $annotMappingDir  = "$samplesDir/$sampleId/mappings/$annotDbDir\_$annotVersion\-$alnVersion";
my $auxMappingDir    = "$samplesDir/$sampleId/mappings/$annotDbDir\_$auxPrefix\-$alnVersion";

if (!(-d "$samplesDir/$sampleId"))
{
    die "error: sample $sampleId does not exist in $samplesDir\n";
}
if (!(-d "$genomesDir/$genomeVersion"))
{
    die "error: $genomeVersion does not exist in $genomesDir\n";
}
if (!(-d "$seqAnnotDir/$annotVersion"))
{
    die "error: $annotVersion does not exist in $seqAnnotDir\n";
}
if (!(-f "$alnSoftDir/$alnVersion.pl"))
{
    die "error: $alnVersion does not exist in $alnSoftDir\n";
}

# read the filenames (only id) of all the auxiliary sequences and store them in a hash
opendir(AUXDIR, "$seqAnnotDir/$auxPrefix/")
  or die("Cannot open directory $seqAnnotDir/aux/");
my @auxDbfiles = readdir(AUXDIR);
my %auxDbHash;
foreach (@auxDbfiles)
{
    if ($_ =~ /^$auxPrefix(\d+)\_/) { $auxDbHash{$1} = $_; }
}

# check if they exist
for my $auxFile (@auxFileList)
{
    if (!(defined($auxDbHash{$auxFile})))
    {
        die "error: auxiliary file $auxPrefix\_$auxFile does not exist in $seqAnnotDir/$auxPrefix/\n";
    }
}

# all the directories with the database files exist to perform the alignments

# genome mappings if not present yet
if (system("mkdir $genomeMappingDir 2> /dev/null") == 0)
{
    my $genTmpDir = "$tempRootDir/$sampleId\-$genomeDbDir\_$genomeVersion\-$alnVersion";

    # check if temporary dir exists, if yes delete everything in it
    if (system("mkdir $genTmpDir 2> /dev/null") != 0)
    {
        if( system("rm $genTmpDir/* -r 2> /dev/null") != 0)
        {
			die "Error in creating or removing $genTmpDir. Check permissions!\n";
		}
    }

    my $ret = system(
        "$alnSoftDir/$alnVersion.pl $genomesDir/$genomeVersion $samplesDir/$sampleId/$seqBaseFileName.tab g $genTmpDir $genomeMappingDir $thread_number 2>>$samplesDir/$sampleId/$errorLogDirName/$annotOptions.txt"
    );

	# piotr: I uncommented this
    system("rm $genTmpDir/ -r 2> /dev/null");

    if ($ret != 0)
    {
        system("rm -rf $genomeMappingDir");    # delete the whole directory containing the alignments
        die "error during genome mapping\n";
    }

}

# annotation mappings if not present yet
if (system("mkdir $annotMappingDir 2> /dev/null") == 0)
{
    my $annotTmpDir = "$tempRootDir/$sampleId\-$annotDbDir\_$annotVersion\-$alnVersion";

    # check if temporary dir exists, if yes delete everything in it
    if (system("mkdir $annotTmpDir 2> /dev/null") != 0)
    {
        system("rm $annotTmpDir/* -r 2> /dev/null");
    }

    my $ret = system(
        "$alnSoftDir/$alnVersion.pl $seqAnnotDir/$annotVersion $samplesDir/$sampleId/$seqBaseFileName.tab a $annotTmpDir $annotMappingDir $thread_number 2>>$samplesDir/$sampleId/$errorLogDirName/$annotOptions.txt"
    );

	# piotr: I uncommented this.
    system("rm $annotTmpDir/ -r 2> /dev/null");

    if ($ret != 0)
    {
        system("rm -rf $annotMappingDir");    # delete the whole directory containing the alignments
        die "error during annotation mapping\n";
    }

}

# aux mappings
system("mkdir $auxMappingDir 2> /dev/null");

# read the filenames (only id) of all the auxiliary sequences and store them in a hash
opendir(AUXDIRM, $auxMappingDir) or die("Cannot open directory $auxMappingDir");
my @auxMappingfiles = readdir(AUXDIRM);
my %auxMapHash;
foreach (@auxMappingfiles)
{
    if ($_ =~ /^$auxPrefix(\d+)\_/) { $auxMapHash{$1} = $_; }
}

# copy the relevant aux files into a temporary directory and execute the alignments
my $auxTmpDir            = "$tempRootDir/$sampleId\-$annotDbDir\_$auxPrefix\-$alnVersion";
my $auxTmpDirCopiedFiles = "$tempRootDir/$sampleId\-$annotDbDir\_$auxPrefix\-$alnVersion-CP";
if (system("mkdir $auxTmpDir 2> /dev/null") != 0)
{
    system("rm $auxTmpDir/* -r 2> /dev/null");
}
if (system("mkdir $auxTmpDirCopiedFiles 2> /dev/null") != 0)
{
    system("rm $auxTmpDirCopiedFiles/* -r 2> /dev/null");
}

my $auxFilesToMapTo = 0;
for (my $i = 0 ; $i < @auxFileList ; $i++)
{
    if (!(defined $auxMapHash{$auxFileList[$i]}))
    {
        system("cp $seqAnnotDir/$auxPrefix/$auxDbHash{$auxFileList[$i]} $auxTmpDirCopiedFiles");
        $auxFilesToMapTo++;
    }
}

# call mapping prog
if ($auxFilesToMapTo > 0)
{
    my $ret = system(
        "$alnSoftDir/$alnVersion.pl $auxTmpDirCopiedFiles $samplesDir/$sampleId/$seqBaseFileName.tab a $auxTmpDir $auxMappingDir $thread_number -r 2>>$samplesDir/$sampleId/$errorLogDirName/$annotOptions.txt"
    );
    if ($ret != 0) { die "error during aux file mapping\n"; }
}
else
{
    system("rm $auxTmpDir -r 2> /dev/null");
    system("rm $auxTmpDirCopiedFiles -r 2> /dev/null");
}

# create full report if not present

if (!(-f "$samplesDir/$childSampleId/$fullReportsDirName/$annotOptions.tab"))
{
    my $ret = system(
        "$repDir/soft/createFullReport.pl $childSampleId $annotOptions $samplesDir > $samplesDir/$childSampleId/$fullReportsDirName/$annotOptions.tab"
    );
    if ($ret != 0) { die "error when creating full report\n"; }
}

# piotr: just delete the tmp dir right now. The FMI team wanted to have it for 15 days
# available, but I see no reason.

# piotr: this didn't work anyway.
### go through the temp directory and remove everything that is older than $deleteTempTime days
#opendir(IMD, $tempRootDir) || die("error: Cannot open temp directory");
#my @tempDirEntries = readdir(IMD);

#my $currentTimeSecondsEpoch = time();
#foreach (@tempDirEntries)
#{
#    if (/\S+\-\S\S\S\S\S\S\S+/)
#    {    # make sure its not . or ..
#        my $dirTimeSecondsEpoch = stat("$tempRootDir/$_")->mtime;
#        my $dirAgeDays = ($currentTimeSecondsEpoch - $dirTimeSecondsEpoch) / (60 * 60 * 24);
#
#        if ($dirAgeDays > $deleteTempTime)
#        {
#
#            # check if the temp dir contains the keyword
#            system("rm $tempRootDir/$_ -r 2> /dev/null");
#        }
#
#    }
#
#}

sub parseAnnotationOptions
{
    my ($annotOptions) = @_;

    if ($annotOptions =~ /^([^\-]+)\-([^\-]+)\-([^\-]+)(\-$auxPrefix\_(\S+)){0,1}$/)
    {
        my ($genomeVersion, $annotVersion, $alnVersion, $auxString, $auxFileString) = ($1, $2, $3, $4, $5);
        if (!(defined $5)) { $auxFileString = ""; }
        if (!(defined $4)) { $auxString     = ""; }
        my @auxFileList = split("\_", $auxFileString);

        return ($genomeVersion, $annotVersion, $alnVersion, $auxString, \@auxFileList);

    }
    else
    {
        die "error: annotation options $annotOptions are in a wrong format\n";
    }
}