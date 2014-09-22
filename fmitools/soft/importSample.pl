#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use FindBin;
use lib "$FindBin::RealBin";
use DeepSeqTools;

# reads one file of unique sequences with counts and inserts the data into deepSeqRepos
# input can be either fasta file (default) or table (-t tab)
#
# fasta:
#
# >id1     800
# ACGTTTATTAG
# >id2     23
# AGCCTTATTAGTT
#
#
# table
#
# ACGTTTATTAG    800   id1(optional)
# AGCCTTATTAGTT  23    id2(optional)

# return values:
# 0 - went well
# 1 - sample exists. Severin will continue.
# 255 - die was called for some other reason (failed)

# absolute path to the repository
my $repDir             = "$Bin/..";
#my $samplesDir         = "$repDir/samples";
my $seqBaseFileName    = "seqs";
my $descrFileName      = "description.txt";
my $fullReportsDirName = "fullReports";
my $errorLogDirName    = "errorLogs";
my $mappingsDirName    = "mappings";

my $fileFormat     = "fasta";    # or tab
my $keepOldSeqId   = 0;
my $filterDisabled = 0;
my $results = GetOptions(
                         'f:s' => \$fileFormat,
                         'i'   => \$keepOldSeqId,
                         'v'   => \$filterDisabled
                        );

if (@ARGV != 4)
{
    print STDERR "\nUSAGE: importSample.pl sequencefile sampleID \"Sample Description\" [options]\n\n";
    print STDERR "   -f [fasta,tab]     sequence file format, default=fasta\n";
    print STDERR "   -v                 disable filter\n";
    print STDERR "   -i                 keep the orgnial id of the imported sequences\n\n";
    print STDERR "\n";
    exit(0);
}

my $sampleId    = $ARGV[1];
my $sampleDescr = $ARGV[2];

# 1.11.13 modified by severin. Added samplesDir as argument
my $samplesDir  = $ARGV[3];

print "loading sequences...\n";

# load the sequencing data
my $seqs;
if ($fileFormat eq "fasta")
{
    $seqs = loadOneLineFasta($ARGV[0]);
}
elsif ($fileFormat eq "tab")
{
    $seqs = loadOneLineTable($ARGV[0]);
}
else { die "unknown -f option $fileFormat\n"; }

# do some checks

my $nrSeqs                        = @{$seqs};
my $seqsWithNonRegularNucleotides = 0;

foreach (@{$seqs})
{

    # replace non ACGTN characters by N
    if (!($_->[1] =~ /^[ACGTN]+$/))
    {
        $_->[1] =~ tr/ACGTN/N/c;    # replace all non ACGTNs by Ns
        $seqsWithNonRegularNucleotides++;
    }
}

# filter sequences
print "filtering sequences...\n";
my @seqsFiltered;
my $totalSeqCounts = 0;
foreach (@{$seqs})
{
    my $filCode = 0;
    if (!($filterDisabled)) { $filCode = filter($_->[1]); }

    if ($filCode == 0)
    {
        push(@seqsFiltered, $_);
        $totalSeqCounts += $_->[2];
    }
}
my $nrSeqsFiltered = @seqsFiltered;

print "\n";

if ($sampleId =~ /\//) { die "error: sample identifier $sampleId must not contain a slash(/)\n\n"; }

if ($seqsWithNonRegularNucleotides != 0)
{
    print "warning: $seqsWithNonRegularNucleotides out of $nrSeqs sequences carried non ACGTN characters and were converted to Ns\n\n";
}

print "sample ID:      \t$sampleId\n";
print "description     \t$sampleDescr\n";
print "unique sequences\t$nrSeqs\n";
print "unique sq after filter\t$nrSeqsFiltered\n";
print "total counts    \t$totalSeqCounts\n";
print "\n";

# copy sample into repository
my $newSampleDir = "$samplesDir/$sampleId";

#my $ret=system("mkdir $newSampleDir 2> /dev/null");
mkdir($newSampleDir);
#my $ret = 


# replace the system given sequence id by the old one if given the option -i
if ($keepOldSeqId == 1)
{

    # check if the old id is unique
    my %oldIdHash;
    for (my $i = 0 ; $i < @{$seqs} ; $i++)
    {
        if (defined $oldIdHash{$seqs->[$i][3]})
        {
            die "error: option -i is used but the original sequence id $seqs->[$i][3] is not unique\n";
        }
        else
        {
            $oldIdHash{$seqs->[$i][3]} = 1;
        }
    }

    # replace the automatic id by the original one
    for (my $i = 0 ; $i < @{$seqs} ; $i++)
    {
        $seqs->[$i][0] = $seqs->[$i][3];
    }
}


if (! -e "$newSampleDir/$seqBaseFileName.tab" )  # piotr: $ret before
{

    # copy data into the new sample directory

	# piotr: I removed the unfiltered part, because it was redundant. Instead I
	# I am now creating symlinks. According to Michael:
	# Hi Piotr,
	# For sequence files that were processed with filterAdaptors.pl, the files
	# will be identical. The simplest solution to reclaim disk space for such
	# samples is to delete the seqs_unfiltered.tab file and create a symbolic
	# link with that name pointing to seqs.tab.
	# Regards,
	# Michael
	
    ## save table file for unfiltered
    #open(OUTFILE, ">" . "$newSampleDir/$seqBaseFileName\_unfiltered.tab")
    #  or die "File system error. cannot write to $newSampleDir/$seqBaseFileName\_unfiltered.tab\n";
    #foreach (@{$seqs})
    #{
    #    print OUTFILE join("\t", @{$_}), "\n";
    #}
    #close(OUTFILE);

    # save table file
    open(OUTFILE, ">" . "$newSampleDir/$seqBaseFileName.tab") or die "File system error. cannot write to $newSampleDir/$seqBaseFileName.tab\n";
    foreach (@seqsFiltered)
    {
        print OUTFILE join("\t", @{$_}), "\n";
    }
    close(OUTFILE);

	# piotr:
    system("ln -s $seqBaseFileName.tab $newSampleDir/$seqBaseFileName\_unfiltered.tab");
    # it is intentionally a link without a full path in case we want to move the repository somewhere in the future.

    # save description line in the description.txt
    open(OUTFILE, ">" . "$newSampleDir/$descrFileName") or die "File system error. cannot write to $newSampleDir/$descrFileName\n";
    print OUTFILE "$sampleDescr\n";
    close(OUTFILE);

    # create some directories
    my $ret1 = mkdir("$newSampleDir/$fullReportsDirName");
    my $ret2 = mkdir("$newSampleDir/$errorLogDirName");
    my $ret3 = mkdir("$newSampleDir/$mappingsDirName");

    print "The new sample $sampleId was successfully imported\n\n";

}
else
{

    #die "error: cannot import sample because $sampleId already exists\n\n";
    #die "error: cannot import sample: $!\n\n";
    warn "error: The sample was already imported: file $newSampleDir/$descrFileName exists!\n\n";
    exit 1;
}

sub loadOneLineFasta
{
    my ($filename) = @_;

    my @seqs;
    my $c = 1;

    my %h;
    open(F, $filename) or die "Cannot open $filename\n";
    while (<F>)
    {
        if (/^\>(\S+)\s+(\d+)/)
        {
            my ($oldId, $counts) = ($1, $2);
            $_ = <F>;
            chomp($_);
            my $seq = $_;

            $seq =~ s/\s//g;            # remove spaces
            $seq =~ tr/[a-z]/[A-Z]/;    # capitalize
            $seq =~ tr/U/T/;            # convert U to T

            if (!($seq =~ /^[A-Z]+$/)) { die "error: input sequence $seq contains non regular characters\n"; }
            if   (defined $h{$seq}) { die "error: sequence $seq is not unique\n"; }
            else                    { $h{$seq} = 1; }

            push(@seqs, ["sq$c", $seq, $counts, $oldId]);
            $c++;
        }
        else
        {
            die "error: input file $filename is not a one line per sequence fasta file with counts\n";
        }
    }
    return \@seqs;
}

sub loadOneLineTable
{
    my ($filename) = @_;

    my @seqs;
    my $c = 1;

    my %h;
    open(F, $filename) or die "Cannot open $filename\n";
    while (<F>)
    {
        chomp();
        my @l = split("\t");
        if (/^\#/) { next; }
        if (/^(\S+)\s+(\d+)/)
        {
            my ($seq, $counts) = ($1, $2);
            my $oldId = "sq$c";

            $seq =~ s/\s//g;            # remove spaces
            $seq =~ tr/[a-z]/[A-Z]/;    # capitalize
            $seq =~ tr/U/T/;            # convert U to T

            if (@l == 3) { $oldId = $l[2]; }

            if (!($seq =~ /^[A-Z]+$/)) { die "error: input sequence $seq contains unknown characters\n"; }
            if   (defined $h{$seq}) { die "error: sequence $seq is not unique\n"; }
            else                    { $h{$seq} = 1; }

            push(@seqs, ["sq$c", $seq, $counts, $oldId]);
            $c++;
        }
        else
        {
            die "error: input file $filename does not contain one sequence per line with counts.\n";
        }
    }
    return \@seqs;
}
