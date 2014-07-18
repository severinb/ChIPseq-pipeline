#!/usr/bin/env perl

# file coded by Piotr based on frag2bed.pl from the original FMI toolbox.
# It outputs the ,,bedweight'' file which stores in the score column the
# number or reads / number of mappings

use strict;
use warnings;
use Getopt::Long;

#use FindBin '$RealBin';

# global variables
my $header = 0;
my $maxHits  = 99999;
my $ignCnts  = 0;
my $noUCSC     = 0;
my $verbose    = 0;

# digest command line
my $res = GetOptions(
                     "m|max|maxhits=i"  => \$maxHits,
                     "h|header" => 		\$header,
                     "i|ignorecounts"   => \$ignCnts,
                     "noUCSC"      		=> \$noUCSC,
                     "v|verbose"          => \$verbose,
                    );

if (scalar(@ARGV) != 1)
{
    die(  "usage:\n\t$0 fragFile|- [>bedFile] [options]\n\n"
        . "\t  -m int   : maximum number of hits per tag [$maxHits]\n"
        . "\t  -h       : print header with a sample name [no header]\n"
        . "\t  -i       : ignore sequence counts (all counts are set to one)\n"
        . "\t  -noUCSC  : do not produce UCSC/BEDtoos-compatible file (numbering -p- DNA bonds)\n"
        . "\t             instead number the nucleotides in a 1-based coordinates\n"
        . "\t  -v       : verbose\n"
        . "\n\n");
}
my $startshift = -1;
if ($noUCSC)
{
    $startshift = 0;
}

# open IO
my ($FIN, $FOUT);
if ($ARGV[0] eq '-')
{
    $FIN = \*STDIN;
}
elsif (not -r $ARGV[0])
{
    die "error: could not read input file $ARGV[0]: $!\n";
}
elsif ($ARGV[0] =~ m/\.gz$/)
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
    open($FIN, "<$ARGV[0]") or die "error reading $ARGV[0]: $!\n";
}
$FOUT = \*STDOUT;

# parse fragment report
print STDERR "creating bed file...\n" if ($verbose);
my $sep = "\t";
my $currentSample = "";
while (my $line = <$FIN>)
{
    if ($line =~ m/^#/) { next; }
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
    if ($f[10] <= $maxHits)
    {

        # output header line
        if ($f[8] ne $currentSample)
        {
            if($currentSample)
            {
				warn "Warning: multiple samples will be contained in a single bedweight file $currentSample, $f[8]\n";
			}
            $currentSample = $f[8];
            print $FOUT "track name='$currentSample'\n" if $header;
        }

        # output $n bed lines
        #  $n := round(read_count / read_weight)
        my $n = ($ignCnts ? 1 : $f[9]) / $f[10];
        my $bedline = join($sep, $f[1], $f[3] + $startshift, $f[4], $f[0], $n, $f[2]) . "\n";
        print $FOUT $bedline;
    }
}
close $FOUT;
print STDERR "done\n" if $verbose;

# clean up and leave
exit(0);

