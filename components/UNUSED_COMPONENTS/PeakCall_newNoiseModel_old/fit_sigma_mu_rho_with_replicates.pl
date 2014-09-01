#!/import/bc2/soft/app/perl/5.10.1/Linux/bin/perl -w
use strict;
use warnings;

#just one input file with damned counts
my $infile = $ARGV[0]; # output saeeds read counter
#if you know the optimal 2*sigma^2 then you can give it on the command line and it will produce the output
#otherwise, give a negative number
my $insq = -1;

my $win_fg = $ARGV[1]; # window length foreground
my $win_bg = $ARGV[2]; #               background
my $DIR = $ARGV[3]; # output directory

# get the window lengths used for the sliding window (bg and fg)
my $ratio = $win_fg/$win_bg;
print STDERR "Window lengths used:\n";
print STDERR "\t-fg: $win_fg\n";
print STDERR "\t-bg: $win_bg\n";
print STDERR "\t-ratio: $ratio\n\n";

## load counts, apply cut-off on reads, and combine to 'input' to fitting
## use output from script: combine.pl
my $nbrfg = 0; # nbr of fg samples
my $nbrbg = 0; # nbr of fg samples

open(F,$infile) || die "cannot read input file\n";
while(<F>) {

    chomp($_);
    my @s = split(/\t/,$_);

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

## start EM shit
my $rawCutoff = 2; # tags cut-off 
my $twopi = 3.141592658979 * 2.0;
my $oneovertwopi = 1.0/$twopi;
my $pc = 0.0; # pseudo-count. Only really needed when cut-offs are zero

#log-differences
my @sd = ();
#scatter information on the segments
my @xval = ();
#sum inverse raw counts 1/n + 1/m
my @sam = ();
#For the minimal and maximal log-ratio fg/bg
my $mindd = 1000;
my $maxdd = -1000;

my @totfg = (); for(my$i=0;$i<$nbrfg;$i++) { push(@totfg,0); }
my @totbg = (); for(my$i=0;$i<$nbrbg;$i++) { push(@totbg,0); }


my $L = 0;
#read in counts and totals:
print STDERR "Reading $infile\n";
open(F,$infile) || die "cannot read input file\n";
while(<F>){
    
    if( $_ =~ /\#/ ) { next; }

    chomp($_);
    my @s = split(/\s+/,$_);

    for(my$i=4;$i<(4+$nbrfg);$i++) { $totfg[$i-4] += $s[$i]; }
    for(my$i=($nbrfg+4);$i<(4+$nbrfg+$nbrbg);$i++) { $totbg[$i-4-$nbrfg] += $s[$i]; }

    ++$L;
}
close(F);

my @pseudofg = (); 
for(my$i=0;$i<$nbrfg;$i++) { push(@pseudofg,$totfg[$i]/2/$L); }

my @pseudobg = (); 
for(my$i=0;$i<$nbrbg;$i++) { push(@pseudobg,$totbg[$i]/2/$L); }

print STDERR "'Total' number of counts and pseudo counts\n";
for(my$i=0;$i<$nbrfg;$i++) { print STDERR "\t-fg $i: $totfg[$i] ($pseudofg[$i])\n"; }
for(my$i=0;$i<$nbrbg;$i++) { print STDERR "\t-bg $i: $totbg[$i] ($pseudobg[$i])\n"; }
print STDERR "\n";


my $totfg = 0;
my $totbg = 0;

my $pcfg = $totfg/(2*$L);
my $pcbg = $totbg/(2*$L);

my $totwin = 0;
my $passed = 0;
open(F,$infile) || die "cannot read input file\n";
while(<F>){

    if( $_ =~ /\#/ ) { next; }

    chomp($_);
    my @s = split(/\s+/,$_);

    # if one window in the background and one window in the foreground make the cut-off -> include window
    my $makecut_fg = 0;
    my $makecut_bg = 0;
    for(my$i=4;$i<(4+$nbrfg);$i++) {  if( $s[$i]>=$rawCutoff ) { $makecut_fg++; } }
    for(my$i=($nbrfg+4);$i<(4+$nbrfg+$nbrbg);$i++) { if( $s[$i]>=$rawCutoff ) { $makecut_bg++; } }

    $totwin++;
    if( $makecut_fg==0 || $makecut_bg==0 ) { next; }
    $passed++;
    
    # account for the different window length: do we really need to do that? 
    # the way we calculated the tags mapping to the window and the total number does already account for that, doesn't it?
    my $f = 0;
    my $b = 0;
    for(my$i=4;$i<(4+$nbrfg);$i++) { $f += log( ($s[$i]+$pseudofg[$i-4]) ); }
    for(my$i=($nbrfg+4);$i<(4+$nbrfg+$nbrbg);$i++) { $b += log( ($s[$i]+$pseudobg[$i-4-$nbrfg]) ); }
    $f /= $nbrfg;
    $b /= $nbrbg;
    my $dd = $f-$b;
    if($dd > $maxdd){ $maxdd = $dd; }
    if($dd < $mindd){ $mindd = $dd; }

    $f = 0;
    $b = 0;
    for(my$i=4;$i<(4+$nbrfg);$i++) { $f += log( ($s[$i]+$pseudofg[$i-4])/$totfg[$i-4] ); }
    for(my$i=($nbrfg+4);$i<(4+$nbrfg+$nbrbg);$i++) { $b += log( ($s[$i]+$pseudobg[$i-4-$nbrfg])/$totbg[$i-4-$nbrfg] ); }
    $f /= $nbrfg;
    $b /= $nbrbg;
    $dd = $f-$b;
    push @sd , $dd;

    # (sum_i 1/fg_i)/r^2 + (sum_j 1/bg_j)/k^2
    my $fv = 0;
    my $bv = 0;
    my $pc = 0.5;

    for(my$i=4;$i<(4+$nbrfg);$i++)                 { $fv +=  1/($s[$i]+$pc); }
    for(my$i=($nbrfg+4);$i<(4+$nbrfg+$nbrbg);$i++) { $bv +=  1/($s[$i]+$pc); }
    $fv /= $nbrfg**2;
    $bv /= $nbrbg**2;
    my $samp = $fv+$bv;
    push @sam, $samp;

    # some shit for scatters and output
    push @xval, $_;
}
close(F);


print STDERR "Total number of windows (passed): $totwin ($passed)\n\n";

my $range = $maxdd-$mindd;
my $W = 1.0/$range;
my $num = @sd;
print STDERR "range: [$mindd,$maxdd]\tnum: $num\tW = $W\n\n";

my $mean = 0;
my $var = 0;
for(my$i=0;$i<$num;++$i){
    #$sd[$i] -= $lrat;
    $mean += $sd[$i];
    $var += $sd[$i]*$sd[$i];
}
$mean /= $num;
$var /= $num;
$var -= $mean*$mean;

print STDERR "Some statistics:\n";
print STDERR "\t- mean difference: ".sprintf("%.3f",$mean)."\n";
print STDERR "\t- standard-dev: ".sprintf("%.3f",sqrt($var))." on $num windows\n\n";



print STDERR "##### Fitting ###### \n";
my $sqmax = 1.25*$var;
my $sqmin = 0.001;
my $sq = undef; # used to be two sigma squared; now its sigma^2 * (1/#fg_samples^2 + 1/#bg_samples^2)
my $pi = 0.005;
my $mu = $mean;
my $prefac = 1/$nbrfg**2 + 1/$nbrbg**2;

while (abs(2.0 * ($sqmax - $sqmin) / ($sqmax + $sqmin)) > 0.01)
{
    #go halfway
    $sq = 0.5 * ($sqmax + $sqmin);

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
	    my $Ldif = $pi * $W;
	    my $v = 1.0 / ($sq + $sam[$i]); # usedo to be: 1/(2sigma^2 + 1/n_f + 1/n_b); sq=2sigma^2; sam= 1/n_f + 1/n_b
	    my $dd = $sd[$i];
	    my $dsq = ($dd-$mu)*($dd-$mu);
	    my $Lsame = (1.0 - $pi) * exp(-0.5 * $dsq * $v) * sqrt($v * $oneovertwopi);
	    my $perc = $Ldif / ($Lsame + $Ldif);
	    $tot += $perc;
	    $dsig += 0.5 * $v * ($dsq * $v - 1) * (1.0-$perc) * $prefac;
	    # the numerator looks good. This is d(lsame)/d(sq)
	    # the denominator comes from d(log(L))/d(sq) = d(lsame)/d(sq) / (L)
	    $totLogLikelihood += log($Lsame + $Ldif);
	    $munew += $dd*$v*(1.0-$perc);
	    $totv += $v*(1.0-$perc);
	}
	my $oldpi = $pi;
	$pi    = $tot / $num;
	$munew /= $totv;

	if ($pi < 0.000001)
	{
	    $pi = 0;
	}
	if ($oldpi > 0 || $pi > 0)
	{
	    $dev = 0.5 * abs($pi - $oldpi) / ($pi + $oldpi);
	}
	else
	{
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
my $sigma = sqrt($sq/$prefac);
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
print G "#fitted: sigma=$sigma rho=$rho mu=$mu range=$range windows=$totcount fgreads=$totfg bgreads=$totbg\n";
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

print P 'library(Hmisc)'."\n";
print P 'hx <- rnorm(10000000,m=0,sd=1)'."\n";
print P 'd = read.csv(file="'.$DIR.'/outzvals",sep="\t")'."\n";
print P 'x=hist(d[,'.$score_col.'],xlim=c(-12,12),col="red",xlab="z-value",breaks=500,plot=F);'."\n";

print P 'jpeg("'.$DIR.'/hist_zvals.jpg")'."\n";
print P 'plot(x$mids,x$density,type="s",log="y",xlim=c(-12,12),xlab="Z-values",ylab="Density",main="Rev Cum Z-values")'."\n";
print P 'lines(density(hx), type="l", lwd=2, lty=1,col="red")'."\n";
print P 'n=5'."\n";
print P 'minor.tick(nx=n, tick.ratio=1)'."\n";
print P 'dev.off()'."\n";

close P;

system('R CMD BATCH '.$DIR.'hist.R');
