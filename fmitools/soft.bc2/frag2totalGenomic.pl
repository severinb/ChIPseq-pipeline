#!/usr/bin/env perl

# this script is intended to be piped into for genomic fragment report
# to get the total mapped read count. It does not suffer from the floor(count/#mappings)
# done in extractBed

# can process multiple samples at a time, but they must be
# coming one after another, and not mixed lines.


use strict;
use warnings;
use Getopt::Long;
use List::Util qw(min max);

my $quiet         = 1;

if (scalar(@ARGV) != 1)
{
    die(  "usage:\n\t$0 fragFile\n\n" );
}

# open IO
my ($FIN);

if ($ARGV[0] =~ m/\.gz$/)
{
    open($FIN, "gunzip -dc $ARGV[0]|")
      or die "error decompressing $ARGV[0]: $!\n";
}
elsif ($ARGV[0] =~ m/\.bz2$/)
{
    open($FIN, "bunzip2 -dc $ARGV[0]|")
      or die "error decompressing $ARGV[0]: $!\n";
}
else
{
    open($FIN, "<$ARGV[0]") or die "error reading $ARGV[0]: $!\n"; # can be '-'
}

# store one sequence at a time
my $currentSample      = "";
my $nSample            = 0;
#my $currentSampleCount = 0;
my %sampleCounts;
my %seqCounts;
while (my $line = <$FIN>)
{
    if ($line =~ m/^\#/) { next; }
    chomp $line;
    my @f = split("\t", $line);

    # format: $f[0] : read id
    #           [1] : target id
    #           [2] : target strand
    #           [3] : target start
    #           [4] : target end
    #           [5] : number of errors
    #           [6] : sequence alignment string
    #           [7] : target alignment string
    #           [8] : sampleId
    #           [9] : read count
    #          [10] : read inverse weight

    if ($f[8] ne $currentSample)
    {
        if ($nSample > 0)
        {
            $sampleCounts{$currentSample} = getSampleCounts($currentSample);
            %seqCounts = ();
        }
        $nSample++;
        $currentSample      = $f[8];
        #$currentSampleCount = 0;
        print STDERR "\tparsing $currentSample\n" if ($quiet == 0);
    }
    $seqCounts{$f[0]} = $f[9];
    #$currentSampleCount += $f[9];
}
$sampleCounts{$currentSample} = getSampleCounts($currentSample);
foreach my $sample (keys %sampleCounts)
{
	print $sample, "\t", $sampleCounts{$sample}, "\n";
}
close $FIN;

sub getSampleCounts
{
	my ($s) = @_;
	my $ret = 0;
	foreach my $v (values %seqCounts)
	{
		$ret += $v;
	}
	return $ret;
}
