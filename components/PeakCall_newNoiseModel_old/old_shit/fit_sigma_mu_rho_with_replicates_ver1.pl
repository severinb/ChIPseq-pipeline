#!/import/bc2/soft/app/perl/5.10.1/Linux/bin/perl -w
#use strict;
#use fileOperations;
#use List::Util qw(first max min sum);
#use Math::Random qw(:all);

#just one input file with damned counts
$infile = $ARGV[0]; # output saeeds read counter
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

#read in counts and totals:
print STDERR "Reading $infile\n";
open(F,$infile) || die "cannot read input file\n";
#log-differences
my @sd = ();
#scatter information on the segments
my @xval = ();
#sum inverse raw counts 1/n + 1/m
my @sam = ();
#For the minimal and maximal log-ratio fg/bg
my $mindd = 1000;
my $maxdd = -1000;

$L = 0;
while(<F>){
    chomp($_);
    @s = split(/\s+/,$_);
    $totfg += $s[4];
    $totbg += $s[5];
    ++$L;
}
close(F);
$pcfg = $totfg/(2*$L);
$pcbg = $totbg/(2*$L);

open(F,$infile) || die "cannot read input file\n";
while(<F>){

    chomp($_);
    my @s = split(/\t/,$_);

    if($s[4] >= $rawCutoff && $s[5] >= $rawCutoff){

	# account for the different window length: ! only change the fraction (f/F), not the raw counts!
	my $dd = log(($s[4]+$pcfg)/($s[5]*$ratio+$pcbg));
	if($dd > $maxdd){
	    $maxdd = $dd;
	}
	if($dd < $mindd){
	    $mindd = $dd;
	}

	push @sd , $dd;
	$xx = $s[0] . "\t" . $s[1] . "\t" . $s[2] . "\t" . $s[3] . "\t" . $s[4] . "\t" . $s[5];
	push @xval, $xx;
	# 1/n + 1/m
	my $samp = (1.0/($s[4]+$pc)) + (1.0 /($s[5]+$pc));
	push @sam, $samp;
    }
}
close(F);
$totfg /= 2;
$totbg /= 2;
$totratio = $totfg/$totbg;
$lrat = log($totratio);
print STDERR "totcount fg $totfg and bg $totbg\n";
my $range = $maxdd-$mindd;
print STDERR "maximum $maxdd minimum $mindd difference ", $range, " compared to ", log(2000.), " that we used before\n";
#change to set the range to maximal difference
my $W = 1.0/$range;

#normalize the log-differences and calculate mean and var
$num = @sd;
print STDERR "number of lines that made the cut-off: $num\n";


$mean = 0;
$var = 0;
for($i=0;$i<$num;++$i){
    $sd[$i] -= $lrat;
    $mean += $sd[$i];
    $var += $sd[$i]*$sd[$i];
}
$mean /= $num;
$var /= $num;
$var -= $mean*$mean;
print STDERR "mean difference $mean standard-dev ", sqrt($var), " on $num windows\n";

print STDERR "Fitting\n";
my $sqmax = 1.25*$var;
my $sqmin = 0.001;
my $sq = undef; # two sigma squared
my $pi = 0.005;
my $mu = $mean;

if($insq>0){
    $sqmax = 1.02*$insq;
    $sqmin = 0.98*$insq;
}
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
	$totv = 0;
	for (my $i = 0 ; $i < $num ; ++$i){
	    my $Ldif = $pi * $W;
	    my $v = 1.0 / ($sq + $sam[$i]);
	    my $dd = $sd[$i];
	    my $dsq = ($dd-$mu)*($dd-$mu);
	    my $Lsame = (1.0 - $pi) * exp(-0.5 * $dsq * $v) * sqrt($v * $oneovertwopi);
	    my $perc = $Ldif / ($Lsame + $Ldif);
	    $tot += $perc;
	    $dsig += 0.5 * $v * ($dsq * $v - 1) * (1.0-$perc);
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
my $sigma = sqrt(0.5*$sq);
my $rho = $pi; 


open W, ">".$DIR."/fitted_params";
print W join("\t",'sigma:',$sigma)."\n";
print W join("\t",'rho:',$rho)."\n";
print W join("\t",'mu:',$mu)."\n";
print W join("\t",'totfg:',$totfg)."\n";
print W join("\t",'totbg:',$totbg)."\n";
close W;
#print join("\t","RESULT","sigma=$sigma","rho=$rho","mu=$mu")."\n";

#now print out a distribution of z-values
for($i=0;$i<=200;++$i){
    $count[$i] = 0;
    $countfg[$i] = 0;
    $countbg[$i] = 0;
}
open(G,">".$DIR."/outzdist");
open(H,">".$DIR."/outzvals");
$totcount = 0;
for($i=0;$i<$num;++$i){
    my $Ldif = $pi * $W;
    my $v = 1.0 / ($sq + $sam[$i]);
    my $dd = $sd[$i];
    my $dsq = ($dd-$mu)*($dd-$mu);
    my $Lsame = (1.0 - $pi) * exp(-0.5 * $dsq * $v) * sqrt($v * $oneovertwopi);
    my $perc = $Ldif / ($Lsame + $Ldif);
    $z = ($dd-$mu)*sqrt($v);
    print H $xval[$i], "\t", $z, "\t", $perc,"\n";
    #turn into a bin;
    $bin = int(10*($z+10));
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
for($i=0;$i<=200;++$i){
    $z = -10+0.1*$i;
    $perc = $count[$i]/$totcount;
    $percfg = $countfg[$i]/$totcount;
    $percbg = $countbg[$i]/$totcount;
    print G $z, "\t", $perc,"\t", $percfg, "\t", $percbg,"\n";
}
close(G);
close(H);


#system('R CMD BATCH hist.R');
