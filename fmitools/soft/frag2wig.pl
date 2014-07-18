#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use FindBin '$RealBin';
use lib "$RealBin";
#use DeepSeqTools;
#use Inline C => "DATA",
#    "NAME" => "deepSeqTools",
    #"DIRECTORY" => "$FindBin::RealBin/_Inline_".(hostname())."/",
#    "DIRECTORY" => "$FindBin::RealBin/_Inline_"."BC2"."/",
    #"FORCE_BUILD" => 1,
#    "PRINT_INFO" => 0,
 #   "BUILD_NOISY" => 0,
 #   ;

#use FindBin '$RealBin';
#use Sys::Hostname;
use Inline ("CPP" => "DATA",
	    "NAME" => "frag2wig_coverage",
    	"DIRECTORY" => "$FindBin::RealBin/_Inline_"."BC2"."/",
#	    "DIRECTORY" => "$FindBin::RealBin/_Inline_".(hostname())."/",
#	    "CLEAN_AFTER_BUILD" => 1,
#	    "OPTIMIZE" => "-g -O3",
#	    #"PRINT_INFO" => 1,
#	    #"BUILD_NOISY" => 1,
#	    #"FORCE_BUILD" => 1,
#	    #"GLOBAL_LOAD" => 1,
    );

# global variables
my $name      = "";
my $desc      = "weighted coverage";
my $window    = 100;
my $fragLen   = 0;
my $compress  = 0;
my $maxHits   = 99;
my $seqfilter = "";
my $fixedStep = 0;
my $repTrack  = 0;
my $ignCnts   = 0;
my $pool      = 0;
my $quiet     = 0;
my $log2Out   = 0;

# digest command line
my $res = GetOptions("n|name=s" => \$name,
		     "d|desc=s" => \$desc,
		     "w|win|window=i" => \$window,
		     "l|len|length=i" => \$fragLen,
                     "c|compress" => \$compress,
                     "m|max|maxhits=i" => \$maxHits,
                     "L|log2" => \$log2Out,
		     "s|seqfilter=s" => \$seqfilter,
		     "f|fixed|fixedstep" => \$fixedStep,
		     "t|repeattrack" => \$repTrack,
		     "i|ignorecounts" => \$ignCnts,
		     "p|poolsamples" => \$pool,
                     "q|quiet" => \$quiet,
                    );


if(scalar(@ARGV) != 1) {
    die ("usage:\n\t$0 fragFile|- [>wigFile] [options]\n\n" .
	 "assumes that fragment report is sorted by sample and targetId\n\n" .
         "\t  -n str   : track name suffix (added to sample name) [$name]\n" .
         "\t  -d str   : track description suffix (added to sample name) [$desc]\n" .
         "\t  -w int   : window size [$window]\n" .
         "\t  -l int   : fragment length [$fragLen]\n" .
	 "\t             (reads will be shifted by len/2 towards their 3'-end)\n" .
         "\t  -c       : compress output with bzip2 [uncompressed]\n" .
         "\t  -m int   : maximum number of hits [$maxHits]\n" .
	 "\t  -L       : report coverage in log2(x+1) [linear]\n" .
	 "\t  -i       : ignore sequence counts (all counts are set to one)\n" .
	 "\t  -s str   : sequnce filter (ids to keep) [$seqfilter]\n" .
	 "\t  -f       : create fixed step wiggle output [variable step]\n" .
	 "\t  -t       : repeat track header line for each target sequence [single header per sample]\n" .
	 "\t  -p       : pool samples [output different samples as separate tracks]\n" .
	 "\t  -q       : quiet (no progress report on STDERR)\n\n");
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
if($compress) {
    open($FOUT, "| bzip2 -zc") or die "error opening output compressor: $!\n";
} else {
    $FOUT = \*STDOUT;
}


# store one sequence at a time
print STDERR "creating wiggle file...\n" if($quiet==0);
my $shift = sprintf("%.0f",$fragLen/2); # int($fragLen/2 + 0.5) would not work for negative fragment lengthes
my $currentId = "";
my $currentSample = "";
my $nId = 0;
my $nSample = 0;
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
	# new sequence or sample?
	if($f[1] ne $currentId or ($pool==0 and $f[8] ne $currentSample)) {
	    # output former sequence?
	    if($currentId ne ""  and
	       ($seqfilter eq "" or index(uc($seqfilter), uc($currentId))>-1) ) {
		$nId++;
		if($repTrack or $nId==1) {
		    if($pool==0) {
			print $FOUT trackheader("$currentSample$name",
						"$currentSample$desc",
						$nSample);
		    } else {
			print $FOUT trackheader("pooled_samples$name",
						"pooled_samples$desc",
						$nSample);
		    }
		}
		my $len = getLen();
		if ($fixedStep) {
		    # fixedStep
		    print $FOUT "fixedStep chrom=$currentId start=1 step=$window span=$window\n";
		    print $FOUT sprintf("%.2f\n", getAvg(1, $window,$log2Out));
		    for(my $pos=1+$window; $pos<=($len-$window); $pos+=$window) {
			print $FOUT sprintf("%.2f\n", getAvg($pos, $window,$log2Out));
		    }
		} else {
		    # variableStep (suppress zeros)
		    print $FOUT "variableStep chrom=$currentId span=$window\n";
		    print $FOUT "1 ".sprintf("%.2f\n", getAvg(1, $window,$log2Out));
		    for(my $pos=1+$window; $pos<=($len-$window); $pos+=$window) {
			my $val = getAvg($pos, $window,$log2Out);
			if ($val > 0) {
			    print $FOUT "$pos ".sprintf("%.2f\n", $val);
			}
		    }
		}
	    }
	    clear();

	    if($pool==0 and $f[8] ne $currentSample) {
		$nId = 0;
		$nSample++;
	    }

	    $currentId = $f[1];
	    $currentSample = $f[8];
	    print STDERR "\tparsing $currentId\n" if($quiet==0);
	}

	add(($ignCnts ? 1 : $f[9])/$f[10], $f[3], $f[4], ($f[2] eq '+' ? $shift : -$shift));
    }
}
# output last one
if($currentId ne ""  and
   ($seqfilter eq "" or index(uc($seqfilter), uc($currentId))>-1) ) {
    $nId++;
    if($repTrack or $nId==1) {
	if($pool==0) {
	    print $FOUT trackheader("$currentSample$name",
				    "$currentSample$desc",
				    $nSample);
	} else {
	    print $FOUT trackheader("pooled_samples$name",
				    "pooled_samples$desc",
				    $nSample);
	}
    }
    my $len = getLen();
    if ($fixedStep) {
	# fixedStep
	print $FOUT "fixedStep chrom=$currentId start=1 step=$window span=$window\n";
	print $FOUT sprintf("%.2f\n", getAvg(1, $window,$log2Out));
	for(my $pos=1+$window; $pos<=($len-$window); $pos+=$window) {
	    print $FOUT sprintf("%.2f\n", getAvg($pos, $window,$log2Out));
	}
    } else {
	# variableStep (suppress zeros)
	print $FOUT "variableStep chrom=$currentId span=$window\n";
	print $FOUT "1 ".sprintf("%.2f\n", getAvg(1, $window,$log2Out));
	for(my $pos=1+$window; $pos<=($len-$window); $pos+=$window) {
	    my $val = getAvg($pos, $window,$log2Out);
	    if ($val > 0) {
		print $FOUT "$pos ".sprintf("%.2f\n", $val);
	    }
	}
    }
}
clear();
close $FOUT;
print STDERR "done\n" if($quiet==0);

# clean up and leave
exit(0);


sub trackheader {
    my ($nm, $dsc, $i) = @_;
    my @cols = ('27,158,119','217,95,2','117,112,179','231,41,138',
		'102,166,30','230,171,2','166,118,29','102,102,102');
    my $col = $cols[($i-1) % scalar(@cols)];
    return ("track type=wiggle_0 name='$nm' description='$dsc' visibility=full " .
	    "color=$col altColor=$col priority=100 autoscale=off " .
	    "gridDefault=on maxHeightPixels=128:128:11 graphType=bar " .
	    "yLineMark=0.0 yLineOnOff=off windowingFunction=maximum " .
	    "smoothingWindow=off\n");
}


__DATA__
__CPP__
#include <iostream>
#include <vector>
using namespace std;


//# globals
vector<float> cover;

//# functions
void add(float wgt, int start, int end, int shift) {
    int i;
    start += shift;
    if(start<1) {
	start = 1;
    }
    end += shift;
    if(end<start) {
	end = start;
    }
    if(end > cover.size()) {
	cover.resize(end, 0.0);
    }
    for(i=start-1; i<end; i++)
	cover[i] += wgt;
}

int getLen() {
    return cover.size();
}

float get(int pos) {
    return cover[pos-1];
}

float getAvg(int pos, int win, int log2Out) {
    //# average of window centered around pos
    float sum=0.0;
    float avg;
    int i, st, en;
    st = (pos-1) - win/2;
    en = st + win - 1;
    if(st < 0)
	st = 0;
    if(en > cover.size()-1)
	en = cover.size()-1;
    for(i=st; i<=en; i++)
	sum += cover[i];

    avg = sum/(en-st+1);

    //# compute log2(x+1)
    if(log2Out > 0){avg=log(avg+1)/log(2);} 

    return avg;
}

void clear() {
    cover.clear();
}

