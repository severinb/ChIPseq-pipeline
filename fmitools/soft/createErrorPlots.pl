#!/usr/bin/env perl
#
# createErrorPlots.pl (takes about 4-5 minutes per average sample; ~10% longer than extractData.pl alone)
# usage example:
# time ./extractData.pl Dm_H3K27me3_1.5pM,Dm_H3K27me3_2.5pM,Dm_H3K27me3_3.5pM dm3-dmV01-aln1 genome -f | /work/gbioinfo/deepSeqRepos/soft/createErrorPlots.pl -o Dm_H3K27me3 -


use strict;
use warnings;
use Getopt::Long;
use FindBin '$RealBin';
use Sys::Hostname;
use Inline ("C" => "DATA",
	    "NAME" => "createErrorPlots",
#	    "DIRECTORY" => "$FindBin::RealBin/_Inline_".(hostname())."/",
	    "DIRECTORY" => "$FindBin::RealBin/_Inline_"."BC2"."/",
	    "CLEAN_AFTER_BUILD" => 1,
	    "OPTIMIZE" => "-g -O3",
	    #"PRINT_INFO" => 1,
	    #"BUILD_NOISY" => 1,
	    #"FORCE_BUILD" => 1,
	    #"GLOBAL_LOAD" => 1,
	    );


# global variables
#my $sampleDir  = "$RealBin/../samples";
my @hitBins    = (0,1,5,10,50,100,"many"); # don't modify, it's hardcoded below
my $repFname   = "hitReport.tab";
my $seqsFname  = "seqs_unfiltered.tab";
my $maxSamples = 1024;
my $readLen    = 100;
my $maxErr     = 3;
my $ignoreCnt  = undef;
my $outprefix  = sprintf("error_plots_%d-%02d-%02d", (localtime)[5]+1900, (localtime)[4]+1, (localtime)[3]);
my $quiet      = 0;


# digest command line

# 1.11.13 modified by severin: Added sampleDir as argument 
my ($dd, $sampleDir) = @ARGV;
print STDERR "args:\n";
print STDERR "$dd\n";
print STDERR "$sampleDir\n";

my $res = GetOptions('l|len|length=i' => \$readLen,
		     'm|max|maxerr=i' => \$maxErr,
		     'i|ignorecounts' => \$ignoreCnt,
		     'o|out|outprefix=s' => \$outprefix,
		     'q|quiet' => \$quiet,
                    );
if(@ARGV != 2  or  ($ARGV[0] ne "-"  and  ! -r $ARGV[0])) {
    print STDERR "\nUSAGE: createErrorPlots.pl [options] -|fragmentReportFile\n\n";
    print STDERR "  -l int  maximum read length [$readLen]\n";
    print STDERR "          warning: all reads are assumed to be shorter or equal to 'l'\n";
    print STDERR "  -m int  maximum number of errors in alignment [$maxErr]\n";
    print STDERR "  -i      ignore sequence counts\n";
    print STDERR "  -q      quiet (don't report progress)\n";
    print STDERR "  -o str  prefix of output files [$outprefix]\n";
    print STDERR "          output files are:\n";
    print STDERR "           $outprefix.txt (numerical data for error plots)\n";
    print STDERR "           $outprefix.R   (R script to draw error plots)\n\n";
    exit(0);
}


# open IO
my ($FIN, $FOUTdat, $FOUTdraw);
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
open($FOUTdat,  ">$outprefix.txt") or die "error writing to $outprefix.txt: $!\n";
open($FOUTdraw, ">$outprefix.R")   or die "error writing to $outprefix.R: $!\n";
print $FOUTdat "#output of $0, generated on ".scalar(localtime)."\n";


# init further global variables
init_err($readLen, $maxErr+1, $maxSamples);
my @err = map {$_."err"} 0..$maxErr;
my @samples;
my $currentSample = "";
my $mapdir = "";
my $i = -1; # $currentSample in @samples
my @nHits;
for(1..$maxSamples) { push @nHits, [(0) x scalar(@hitBins)]; }
my @totUniqueSeqs = (0) x $maxSamples;
my @totCount = (0) x $maxSamples;
my (%cnt, %iw, %ambig);

# added by severin to count maximum read length:
my $max_readlen = 0;

# parse input (fragment report)
while (my $line = <$FIN>) {
    if($line =~ m/^#/) {
	if ($line =~ m/^#seqId;([^-]+)-[^-]+-(\S+?)(-AUX\S+)?(;[^\t]+)?\t/) {
	    $mapdir = "genomes_" . $1 . "-" . $2;
	} else {
	    die "error: could not get mapping directory from $line\n";
	}
	next;
    }
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


    my $readlen_i = $f[4] - $f[3];
    if ($readlen_i > $max_readlen) {
    	$max_readlen = $readlen_i;
    }

   if($currentSample ne $f[8]) {
	# new sample --> get counts and weights
	$currentSample = $f[8];
	push @samples, $currentSample;
	$i++;
	print STDERR "sample $currentSample\n" unless($quiet);

	my $seqFile = "$sampleDir/$currentSample/$seqsFname";
	my $repFile = "$sampleDir/$currentSample/mappings/$mapdir/$repFname";

	# store sequence counts
	print STDERR "\treading sequence counts..." unless($quiet);
	%cnt = ();
	($totUniqueSeqs[$i], $totCount[$i]) = countTab($seqFile, \%cnt, $ignoreCnt);
	print STDERR "done (n=$totUniqueSeqs[$i], N=$totCount[$i])\n" unless($quiet);

	# store inverse sequence weights and count no. of hits
	print STDERR "\treading sequence weights..." unless($quiet);
	%iw = ();
	%ambig = ();
	my ($id, $err, $n);
	open(IN, "<$repFile") or die;
	while(<IN>) {
	    chomp;
	    ($id, $err, $n) = split(/\t/);

	    if($n==-1) {
		$nHits[$i][6] += $cnt{$id};
		$ambig{$id} = 1;
	    } else {
		$iw{$id} = $n;

		if($n==1) {
		    $nHits[$i][1] += $cnt{$id};
		} elsif($n<=5) {
		    $nHits[$i][2] += $cnt{$id};
		} elsif($n<=10) {
		    $nHits[$i][3] += $cnt{$id};
		} elsif($n<=50) {
		    $nHits[$i][4] += $cnt{$id};
		} elsif($n<=100) {
		    $nHits[$i][5] += $cnt{$id};
		}
	    }
	}
	close IN;
	foreach $id (grep {!exists($iw{$_}) and !exists($ambig{$_})} keys %cnt) {
	    $nHits[$i][0] += $cnt{$id};
	}
	print STDERR "done\n" unless($quiet);

	print STDERR "\tparsing alignments...\n" unless($quiet);
   }

    # upper case the sequence strings (aln3: t --> T, only count sequencing errors, ignore t = C --> T mismatches)
    count_err(uc($f[6]), uc($f[7]), $i, $cnt{$f[0]}/$iw{$f[0]}); # $f[7] is already revcom'ed for $f[2] eq '-'
}

$readLen = $max_readlen + 1;

print "max read length:\n";
print $readLen;
print "\n";

# generate output
print STDERR "generating data output ($outprefix.txt)..." unless($quiet);
print $FOUTdat "HITS\n";
print $FOUTdat join("\t", "sample", @hitBins)."\n";
foreach $i (0..$#samples) {
    print $FOUTdat join("\t", defined($ignoreCnt) ? "$samples[$i] (ignore cnts)" : $samples[$i], @{$nHits[$i]})."\n";
}
print $FOUTdat "ERRORS\n";
print $FOUTdat join("\t", "sample", map {"err$_"} 1..$readLen)."\n";
foreach $i (0..$#samples) {
    print $FOUTdat join("\t", defined($ignoreCnt) ? "$samples[$i] (ignore cnts)" : $samples[$i], map {sprintf("%.3f", $_)} get_err($i, $readLen))."\n";
}
print $FOUTdat "ALIGNED\n";
print $FOUTdat join("\t", "sample", "err", map {"cycle$_"} 1..$readLen)."\n";
foreach $i (0..$#samples) {
    my $cum = 0; # 0+1+2 errors, the same for all cycles
    for(my $err=0; $err<scalar(@err); $err++) {
	my @cnt = map {get_cov_n($i, $_-1, $err)} 1..$readLen;
	$cum += $cnt[0]; # the same for all cycles (add numbers for cycle 1)
	print $FOUTdat join("\t", defined($ignoreCnt) ? "$samples[$i] (ignore cnts)" : $samples[$i], $err[$err], map {sprintf("%.3f",100*$_/$totCount[$i])} @cnt)."\n";
    }
    print $FOUTdat join("\t", defined($ignoreCnt) ? "$samples[$i] (ignore cnts)" : $samples[$i], "ambiguous", map {sprintf("%.3f",100*($nHits[$i][6])/$totCount[$i])} 1..$readLen)."\n";
    print $FOUTdat join("\t", defined($ignoreCnt) ? "$samples[$i] (ignore cnts)" : $samples[$i], "unmapped",  map {sprintf("%.3f",100*($nHits[$i][0])/$totCount[$i])} 1..$readLen)."\n";
}
close $FOUTdat;
print STDERR "done\n" unless($quiet);

print STDERR "generating drawing output ($outprefix.R)..." unless($quiet);
my $nSamples = scalar(@samples);
my $cycleHeader = join(",", (0) x $readLen);
my $Rcode = <<EOD;
fname    <- \"$outprefix.txt\"

# get line offsets
H.fh  <- pipe(paste('grep -n ^HITS', fname))
H.off <- scan(H.fh, what='', nmax=1, quiet=TRUE)
H.off <- as.integer(unlist(strsplit(H.off,':'))[1])
close(H.fh)

E.fh  <- pipe(paste('grep -n ^ERRORS', fname))
E.off <- scan(E.fh, what='', nmax=1, quiet=TRUE)
E.off <- as.integer(unlist(strsplit(E.off,':'))[1])
close(E.fh)

C.fh  <- pipe(paste('grep -n ^ALIGNED', fname))
C.off <- scan(C.fh, what='', nmax=1, quiet=TRUE)
C.off <- as.integer(unlist(strsplit(C.off,':'))[1])
close(C.fh)

# read data
H.format <- list('',0,0,0,0,0,0,0)
E.format <- list('',$cycleHeader)
C.format <- list('','',$cycleHeader)

H <- as.data.frame(scan(file=fname, what=H.format, nmax=E.off-H.off-2,
                        sep='\t', skip=H.off+1, multi.line=FALSE, quiet=TRUE))
names(H) <- c('sample','0','1','5','10','50','100','many')

E <- as.data.frame(scan(file=fname, what=E.format, nmax=C.off-E.off-2,
                        sep='\t', skip=E.off+1, multi.line=FALSE, quiet=TRUE))
names(E) <- c('sample',sprintf('err%d',1:$readLen))

C <- as.data.frame(scan(file=fname, what=C.format, sep='\t',
                        skip=C.off+1, multi.line=FALSE, quiet=TRUE))
names(C) <- c('sample','mm',sprintf('cycle%d',1:$readLen))

library(RColorBrewer)

# repeatedness
pdf('$outprefix-repeatedness.pdf', height=max($nSamples+1.7,4.2), width=12, pointsize=16)
cols <- c(brewer.pal(3,'Set1')[1], brewer.pal(5,'Greens')[5:1], brewer.pal(3,'Greys')[2])
par(las=1, mar=c(5,4+7,4+4,1)+.1, xpd=TRUE)
mp <- barplot(t(H[,2:ncol(H)]), horiz=TRUE, beside=FALSE, col=cols, border=NA, ylim=c(0,nrow(H)+1),
              main='', xlab='No. of reads', ylab='', names.arg=H[,'sample'])
legend(x=max(rowSums(H[,2:ncol(H)]))/2, y=mp[length(mp)]+ifelse(length(mp)==1,1.5,1.5*diff(mp[1:2])), xjust=.5, yjust=.5, bty='n',
       title='No. of hits:', fill=cols, ncol=ncol(H)-1, legend=colnames(H)[2:ncol(H)])
mp <- barplot(t(H[,2:ncol(H)]/rowSums(H[,2:ncol(H)])*100), horiz=TRUE, beside=FALSE, col=cols, border=NA, ylim=c(0,nrow(H)+1),
              main='', xlab='Percent of reads', ylab='', names.arg=H[,'sample'])
legend(x=50, y=mp[length(mp)]+ifelse(length(mp)==1,1.5,1.5*diff(mp[1:2])), xjust=.5, yjust=.5, bty='n',
       title='No. of hits:', fill=cols, ncol=ncol(H)-1, legend=colnames(H)[2:ncol(H)])
dev.off()

# error profile by cycle
pdf('$outprefix-error_profiles.pdf', height=8, width=$readLen/6+2, pointsize=16)
cols <- brewer.pal(nrow(E), name='Set2')
nms <- as.character(E[,'sample'])
plot(1:$readLen, E[1,2:ncol(E)], type='l', lwd=3, col=cols[1],
     ylim=c(0,max(E[,2:ncol(E)], na.rm=TRUE)+.5), main='Error Profile',
     xlab='Cycle (read position)', ylab='Sequencing errors (%)')
if(nrow(E)>1) { for(i in 2:nrow(E)) { lines(1:$readLen, E[i,2:ncol(E)], lwd=3, col=cols[i]) } }
abline(h=0, lty=2, col='gray')
abline(v=c(12,25), lty=3, col='red')
legend(x='topleft', bty='n', ncol=1, lwd=3, cex=0.8, col=cols, legend=nms)
dev.off()

# aligned by error and cycle
pdf('$outprefix-fraction_aligned.pdf', height=8, width=8, pointsize=16, onefile=TRUE)
cols <- brewer.pal(nrow(C)/length(unique(C[,'sample'])), 'Set1')
for(sample in levels(C[,'sample'])) {
  par(mar=c(5, 4, 4, 2) + 0.1)
  inx <- C[,'sample'] == sample
  dat <- as.matrix(C[inx,3:ncol(C)])
  layout(matrix(1:2, ncol=1), widths=1, heights=c(1,.1))
  mp <- barplot(dat, space=0, border=NA, beside=FALSE, col=cols, names.arg=1:$readLen,
                xlab='Read position', ylab='Percent of all reads', main=sample)
  for(y in seq(20,80,by=20)) {
    inc <- diff(mp[1:2]); l <- length(mp)
    segments(mp[1]-0.5*inc, y, mp[l]+0.5*inc, y, col='black', lty=3)
  }
  par(mar=c(0,0,0,0))
  plot(0,0,type='n', axes=FALSE)
  legend(x='top', bty='n', fill=cols, ncol=ceiling(nrow(C)/length(unique(C[,'sample']))/2), legend=C[inx,'mm'])
}
dev.off()

EOD

print $FOUTdraw $Rcode;
close $FOUTdraw;
print STDERR "done\n" unless($quiet);

# clean up and leave
free_all();
print STDERR "create error plots by calling:\n\tRscript --vanilla $outprefix.R\n\n" unless($quiet);
exit(0);


# subs
sub countTab {
    my ($fname, $cnt, $ignore) = @_;
    my ($totS, $totC) = (0, 0);
    open(IN, "<$fname") or die "error reading $fname: $!\n";
    my ($id, $seq, $count);
    if(defined $ignore) {
	while (<IN>) {
	    chomp();
	    ($id, $seq, $count) = split("\t");
	    $cnt->{$id} = 1;
	    $totS++;
	    $totC++;
	}
    } else {
	while (<IN>) {
	    chomp();
	    ($id, $seq, $count) = split("\t");
	    $cnt->{$id} = $count;
	    $totS++;
	    $totC += $count;
	}
    }
    close IN;
    return ($totS, $totC);
}



__DATA__
__C__
#include <stdio.h>
#include <string.h>

int max_err = 0;
int len = 0;
int nSamples = 0;
double**  myerr = NULL;   //# myerr[smpl][cycle] = count
double**  mytot = NULL;   //# mytot[smpl][cycle] = count
double*** mycov = NULL;   //# mycov[smpl][cycle][err] = count

void init_err(int l, int err, int nl) {
    int i, j, k;
    max_err = err;
    len = l;
    nSamples = nl;
    myerr = (double**) malloc(nl * sizeof(double*));
    mytot = (double**) malloc(nl * sizeof(double*));
    mycov = (double***)malloc(nl * sizeof(double**));
    for(i=0; i<nl; i++) {
        myerr[i] = (double*) calloc(l, sizeof(double));
        mytot[i] = (double*) calloc(l, sizeof(double));
        mycov[i] = (double**)malloc(l* sizeof(double*));
        for(j=0; j<len; j++) {
            mycov[i][j] = (double*)calloc(max_err+1, sizeof(double));
        }
    }
}

void free_all() {
    int i, j;
    for(i=0; i<nSamples; i++) {
        free(myerr[i]);
        free(mytot[i]);
        for(j=0; j<len; j++) {
            free(mycov[i][j]);
        }
        free(mycov[i]);
    }
    free(myerr);
    free(mytot);
    free(mycov);
}

void get_err(int smpl, int readlen) {
    int i;

    printf("%d\n", smpl);
    printf("%d\n", len);
    printf("%d\n", readlen);
    printf("%d\n", myerr);
    Inline_Stack_Vars;

    /* push back values to STACK */
    Inline_Stack_Reset;
    if(smpl < nSamples) {
        for(i=0; i<readlen; i++) {
            Inline_Stack_Push( newSVnv(100.0*myerr[smpl][i]/mytot[smpl][i]) ); //# see perldoc perlapi
        }
    }
    Inline_Stack_Done;
}

void get_cov(int smpl, int err) {
    int i;
    Inline_Stack_Vars;

    /* push back values to STACK */
    Inline_Stack_Reset;
    if(smpl < nSamples  &&  err <= max_err) {
        for(i=0; i<len; i++) {
            Inline_Stack_Push( newSVnv(100.0*mycov[smpl][i][err]/mytot[smpl][i]) ); //# see perldoc perlapi
        }
    }
    Inline_Stack_Done;
}

double get_cov_n(int smpl, int cycle, int err) {
    int i;

    if(smpl < nSamples && cycle < len && err <= max_err) {
	return mycov[smpl][cycle][err];
    } else {
	return -1.0;
    }
}

int count_err(char* s1, char* s2, int smpl, double wgt) {
    int i, err=0, l=strlen(s1), g=0;
    l = l>len ? len : l;
    for(i=0; i<l; i++) {
	if(s1[i] == '-') { g++; }
        mytot[smpl][i-g] += wgt;
        if(s1[i] != '-' && s1[i] != s2[i]) {
            //# err = err < max_err ? err+1 : err;
	    err++;  //# mustn't become greater than max_err!
            myerr[smpl][i-g] += wgt;
        }
        mycov[smpl][i-g][err] += wgt;
    }
    return err;
}

int count_err_revcom(char* s1, char* s2, int smpl, double wgt) {
    int i, j, err=0;
    int l=strlen(s1);
    l = l>len ? len : l;
    for(i=0, j=l-1; i<l; i++, j--) {
        mytot[smpl][i] += wgt;
        if((s1[i] == 'A' && s2[j] != 'T') ||
           (s1[i] == 'C' && s2[j] != 'G') ||
           (s1[i] == 'G' && s2[j] != 'C') ||
           (s1[i] == 'T' && s2[j] != 'A')) {
            //# err = err < max_err ? err+1 : err;
	    err++;  //# mustn't become greater than max_err!
            myerr[smpl][i] += wgt;
        }
        mycov[smpl][i][err] += wgt;
    }
    return err;
}



__END__
