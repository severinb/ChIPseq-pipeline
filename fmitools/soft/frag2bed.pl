#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
#use FindBin '$RealBin';

# global variables
my $noheader  = 0;
my $maxHits   = 99;
my $ignCnts   = 0;
my $tabSep    = 0;
my $UCSC      = 0;
my $quiet     = 0;


# digest command line
my $res = GetOptions("m|max|maxhits=i" => \$maxHits,
		     "s|suppressheader" => \$noheader,
		     "i|ignorecounts" => \$ignCnts,
		     "t|tab" => \$tabSep,
		     "U|UCSC" => \$UCSC,
		     "q|quiet" => \$quiet,
                    );


if(scalar(@ARGV) != 1) {
    die ("usage:\n\t$0 fragFile|- [>bedFile] [options]\n\n" .
         "\t  -m int   : maximum number of hits per tag [$maxHits]\n" .
         "\t  -s       : suppress bed track header [standard bed header]\n" .
	 "\t  -i       : ignore sequence counts (all counts are set to one)\n" .
	 "\t  -t       : use tab (instead of space) as field separator\n" .
	 "\t  -U       : produce UCSC-compatible bed file (zero-based start coordinate)\n" .
	 "\t  -q       : quiet (don't report progress)\n" .
	 "\n\n");
}
my $startshift = 0;
if($UCSC) {
    $startshift = -1;
}


# open IO
my ($FIN, $FOUT);
if($ARGV[0] eq '-') {
    $FIN = \*STDIN;
} elsif(not -r $ARGV[0]) {
    die "error: could not read input file $ARGV[0]: $!\n";
} elsif($ARGV[0] =~ m/\.gz$/) {
    open($FIN, "gunzip -dc $ARGV[0]|") or die "error decompressing $ARGV[0]: $!\n";
} elsif($ARGV[0] =~ m/\.bz2$/) {
    open($FIN, "bunzip2 -dc $ARGV[0]|") or die "error decompressing $ARGV[0]: $!\n";
} else {
    open($FIN, "<$ARGV[0]") or die "error reading $ARGV[0]: $!\n";
}
$FOUT = \*STDOUT;


# parse fragment report
print STDERR "creating bed file...\n" unless($quiet);
my $sep = $tabSep ? "\t" : " ";
my $currentSample = "";
while (my $line = <$FIN>) {
    if($line =~ m/^#/) { next; }
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
    if ($f[10]<=$maxHits) {
	# output header line
	if($f[8] ne $currentSample) {
	    $currentSample = $f[8];
	    unless($noheader) { print $FOUT "track name='$currentSample'\n"; }
	}

	# output $n bed lines
	#  $n := round(read_count / read_weight)
	my $n = int(($ignCnts ? 1 : $f[9])/$f[10] + 0.5);
	my $bedline = join($sep,$f[1],$f[3]+$startshift,$f[4],$f[0],$n,$f[2])."\n";
	foreach (1..$n) { print $FOUT $bedline; }
    }
}
close $FOUT;
print STDERR "done\n" unless($quiet);

# clean up and leave
exit(0);



