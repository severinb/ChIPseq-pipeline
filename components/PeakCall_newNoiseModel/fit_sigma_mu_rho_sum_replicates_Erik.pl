#!/usr/bin/env perl
use strict;
use warnings;

#just one input file with damned counts
my $infile = $ARGV[0]; # output saeeds read counter
my $DIR = $ARGV[1]; # output directory
my $plotfile = $ARGV[2]; #file for histogram plot
#my $rawCutoff = int($ARGV[3]);
#NOTE: removed window lengths since we didn't use them


## load counts, apply cut-off on reads, and combine to 'input' to fitting
## use output from script: combine.pl
my $nbrfg = 0; # nbr of fg samples
my $nbrbg = 0; # nbr of fg samples

open(F,$infile) || die "cannot read input file\n";
while(<F>) {
    chomp($_);
    my @s = split(/\s+/,$_);
    for(my $i=0;$i<@s;$i++) {
	if($s[$i] =~ /fg/) { $nbrfg++; }
	if($s[$i] =~ /bg/) { $nbrbg++; }
    }
    last;
}
close F;

print STDERR "Number of samples:\n";
print STDERR "\t-fg: $nbrfg\n";
print STDERR "\t-bg: $nbrbg\n\n";

## Read in the data and store arrays needed to do the EM
my $rawCutoff = 2; # tags cut-off.
my $twopi = 3.141592658979 * 2.0;
my $oneovertwopi = 1.0/$twopi;


#array of log-frequency differences
my @sd = ();
#Information on the segments
my @xval = ();
#contribution to the variance from the Poisson sampling
my @sam = ();
#For the minimal and maximal log-ratio fg/bg
my $mindd = 1000;
my $maxdd = -1000;

my $totfg = 0;
my $totbg = 0;
my $L = 0; #S: Number of lines or windows
#read in count totals:
print STDERR "Reading $infile\n";
open(F,$infile) || die "cannot read input file\n";
while(<F>){
    #comment lines contain # and are skipped
    if( $_ =~ /\#/ ){ 
	next; 
    }
    chomp($_);
    my @s = split(/\s+/,$_);
    for(my$i=0;$i<$nbrfg;$i++) { 
	$totfg += $s[$i+4]; 
    }
    for(my$i=0;$i<$nbrbg;$i++) { 
	$totbg += $s[$i+4+$nbrfg]; 
    }
    ++$L; #S: increments L by one.
}
close(F);

my $pseudo = 0.5;
#needed when calculating difference in the log-frequencies
my $ldiftot = log(($totbg+$pseudo*$L)/($totfg+$pseudo*$L));
#to calculate pooled frequency
my $totdenom = 1.0/($totfg+$totbg+$pseudo*$L);
print STDERR "total fg $totfg total bg $totbg tot num windows $L\n";

my $passed = 0;
open(F,$infile) || die "cannot read input file\n";
##S: Get for each window log(f)-log(b) and 1/nf + 1/nb
while(<F>){
    if( $_ =~ /\#/ ) { next; }

    chomp($_);
    my @s = split(/\s+/,$_);

    my $f = 0;#raw pooled counts for fg and bg
    my $b = 0;
    for(my$i=0;$i<$nbrfg;$i++){
	$f += $s[$i+4];
    }
    for(my$i=0;$i<$nbrbg;$i++){ 
	$b += $s[$i+4+$nbrfg];
    }
    #if either below cut-off, go to next window
    if($f<$rawCutoff || $b<$rawCutoff){
	next; 
    } 
    $passed++;

    #difference in log frequencies (note, only 1 log to calculate)
    my $dd = log(($f+$pseudo)/($b+$pseudo))+$ldiftot;
    push @sd , $dd;
    if($dd > $maxdd){ 
	$maxdd = $dd; 
    }
    if($dd < $mindd){ 
	$mindd = $dd; 
    }

    # overall frequencies estimate
    my $pooled = ($f+$b+$pseudo)*$totdenom; 
    #the expected raw counts are nf = pooled*totfg, nb = pooled*$totbg 
    #Poisson variance is 1/nf + 1/nb
    my $samp = 1.0/$totfg;
    $samp += 1.0/$totbg;
    $samp /= $pooled;
    push @sam, $samp;

    # some shit for scatters and output
    push @xval, $_;
}
close(F);


print STDERR "Total number of windows (passed): $L ($passed)\n\n";

my $range = $maxdd-$mindd;
my $W = 1.0/$range;
my $num = @sd;
print STDERR "range: [$mindd,$maxdd]\tnum: $num\tW = $W\n\n";

my $mean = 0;
my $var = 0;
for(my$i=0;$i<$num;++$i){
    #$sd[$i] -= $lrat;
    $mean += $sd[$i]; #S: compute population mean and variance
    $var += $sd[$i]*$sd[$i];
}
$mean /= $num;
$var /= $num;
$var -= $mean*$mean;

print STDERR "Some statistics:\n";
print STDERR "\t- mean difference: ".sprintf("%.3f",$mean)."\n";
print STDERR "\t- standard-dev: ".sprintf("%.3f",sqrt($var))." on $num windows\n\n";

##mu is actually not a parameter that occurs in the noise model (only sigma and rho are there). But I think we allow here to additionally fit mu. But why and how exactly?
##Shifting mu is allowing to shift the whole noise distribution towards left or right. Thus if data is biased to either foreground or background (right or left) due to systematic bias we can bring the distribution back to 0.
##I think without shifting mu, we assume that the distribution o
##S: Here EM starts:
print STDERR "##### Fitting ###### \n";
my $sqmax = 1.25*$var;
my $sqmin = 0.001;
my $sq = undef; # used to be two sigma squared; now its sigma^2 * (1/#fg_samples^2 + 1/#bg_samples^2) #S: variable for sigma
my $pi = 0.005; #S: pi is the variable for rho
my $mu = $mean; #S: take population mean as initial value for mu
my $prefac = 2;  #was 1/nbrfg + 1/nbrfg, but always nbrfg=nbrbg=1 because we sum counts of replicates. bug: earlier it was 1/$nbrfg**2 + 1/$nbrbg**2

while (abs(2.0 * ($sqmax - $sqmin) / ($sqmax + $sqmin)) > 0.01)
{
    #go halfway
    $sq = 0.5 * ($sqmax + $sqmin); #sq is the Gaussian noise component that we want to fit. It is 2*sigma**2

    #use previous pi unless zero or out of range
    if($pi<= 0.000001 || $pi >= 1){
	$pi = 0.005;
    }
	
    #start with muvalue set to mean from above
    $mu = $mean;
    my $dev = 1.0;
    my $dsig = undef;
    my $munew;
    while ($dev > 0.01)
    {
	#calculate new pi and mu
	my $tot  = 0;
	$dsig = 0;
	my $totLogLikelihood = 0;
	#calculate mu for this pi and sigma
	$munew = 0;
	my $totv = 0;
	for (my $i = 0 ; $i < $num ; ++$i){
	    my $Ldif = $pi * $W; #S: uniform distribution rho times W=1/(max(log(nf/nb)) - min(log(nf/nb)))
	    my $v = 1.0 / ($sq + $sam[$i]); # usedo to be: 1/(2sigma^2 + 1/n_f + 1/n_b); sq=2sigma^2; sam= 1/n_f + 1/n_b #S: v is the whole denominator (2sigma**2+1/nf +1/nb). It's one over variance.
	    my $dd = $sd[$i];
	    my $dsq = ($dd-$mu)*($dd-$mu); #S: this is the numerator
	    my $Lsame = (1.0 - $pi) * exp(-0.5 * $dsq * $v) * sqrt($v * $oneovertwopi);
	    my $perc = $Ldif / ($Lsame + $Ldif);
	    $tot += $perc;
	    $dsig += 0.5 * $v * ($dsq * $v - 1) * (1.0-$perc) * $prefac; #derivative of likelihood with respect to sq
	    # the numerator looks good. This is d(lsame)/d(sq)
	    # the denominator comes from d(log(L))/d(sq) = d(lsame)/d(sq) / (L)
	    $totLogLikelihood += log($Lsame + $Ldif);
	    $munew += $dd*$v*(1.0-$perc);
	    $totv += $v*(1.0-$perc);
	}
	my $oldpi = $pi;
	$pi    = $tot / $num; #update equation where tot is the probability of the data being different (Ldif) and num is the number of data points.
	$munew /= $totv; #update equation for mu. both update equations are the derivatives of the log-likelihood with respect to the parameters.

	if ($pi < 0.000001){
	    $pi = 0;
	}
	if ($oldpi > 0 || $pi > 0){
	    $dev = 0.5 * abs($pi - $oldpi) / ($pi + $oldpi); #
	}
	else{
	    $dev = 0;
	}
	if($mu != 0 || $munew != 0){
	    $dev += 0.5 * abs($mu-$munew)/abs($mu+$munew);
	}
	print STDERR "$sqmax,$sqmin, 2sigma^2 $sq totLogLikelihood $totLogLikelihood dsig $dsig pi(rho) $pi munew $munew mu $mu dev $dev\n"; 
	$mu = $munew;
    }
    print STDERR "=============================================================\n";
    if ($dsig < 0) { 
	$sqmax = $sq; 
    } 
    else { 
	$sqmin = $sq; 
    } 
    print STDERR "\n\n";
}
my $sigma = sqrt($sq/$prefac); #prefac = 2
my $rho = $pi; 


print STDERR "Estimates:\n";
print STDERR "\t-sigma: $sigma\n";
print STDERR "\t-rho: $rho\n";
print STDERR "\t-mu: $mu\n\n";
open W, ">".$DIR."/stats_get_zvals";
print W "Fitted parameters:\n";
print W "\t-sigma: $sigma\n";
print W "\t-rho: $rho\n";
print W "\t-mu: $mu\n\n";
close W;

#now print out a distribution of z-values
my @count;
my @countfg;
my @countbg;

for(my $i=0;$i<=200;++$i){
    $count[$i] = 0;
    $countfg[$i] = 0;
    $countbg[$i] = 0;
}

my $score_col = -99;

open(G,">".$DIR."/outzdist");
open(H,">".$DIR."/outzvals");
my $totcount = 0;
for(my $i=0;$i<$num;++$i){
    my $Ldif = $pi * $W;
    my $v = 1.0 / ($sq + $sam[$i]);
    my $dd = $sd[$i];
    my $dsq = ($dd-$mu)*($dd-$mu);
    my $Lsame = (1.0 - $pi) * exp(-0.5 * $dsq * $v) * sqrt($v * $oneovertwopi);
    my $perc = $Ldif / ($Lsame + $Ldif);
    my $z = ($dd-$mu)*sqrt($v);

    my $string = $xval[$i]."\t".$z."\t".$perc;
    my @abc = split(/\s+/,$string);
    $score_col = @abc;
    $score_col--;

    print H $xval[$i], "\t", $z, "\t", $perc,"\n";
    #turn into a bin;
    my $bin = int(10*($z+10));
    if($bin > 200){
	$bin = 200;
    }
    elsif($bin < 0){
	$bin = 0;
    }
    ++$count[$bin];
    $countfg[$bin] += $perc;
    $countbg[$bin] += 1.0-$perc;
    ++$totcount;
}
print G "#fitted: sigma=$sigma rho=$rho mu=$mu range=$range windows=$totcount fg-reads $totfg bg-reads $totbg", "\n";
for(my $i=0;$i<=200;++$i){
    my $z = -10+0.1*$i;
    my $perc = $count[$i]/$totcount;
    my $percfg = $countfg[$i]/$totcount;
    my $percbg = $countbg[$i]/$totcount;
    print G $z, "\t", $perc,"\t", $percfg, "\t", $percbg,"\n";
}
close(G);
close(H);

## prepare R plot
open(P,">".$DIR."/hist.R");
#print P 'pdf("'.$plotfile.'")'."\n";
print P 'library(Hmisc)'."\n";
print P 'hx <- rnorm(10000000,m=0,sd=1)'."\n";
print P 'd = read.csv(file="'.$DIR.'/outzvals",sep="\t")'."\n";
print P 'x=hist(d[,'.$score_col.'],xlim=c(-12,12),col="red",xlab="z-value",breaks=500,plot=F);'."\n";

print P 'pdf("'.$plotfile.'")'."\n";
print P 'plot(x$mids,x$density,type="s",log="y",xlim=c(-12,12),xlab="Window Z-scores",ylab="Density",main="Histogram of Window Z-scores")'."\n";
print P 'lines(density(hx), type="l", lwd=2, lty=1,col="red")'."\n";
print P 'n=5'."\n";
print P 'minor.tick(nx=n, tick.ratio=1)'."\n";
print P 'dev.off()'."\n";

close P;

system('R CMD BATCH '.$DIR.'/hist.R '.$DIR.'/hist.Rout');
