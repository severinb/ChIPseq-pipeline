use warnings;
use strict;

# INPUT

my $FILE_COUNTS = $ARGV[0]; # file containing (summed) background counts: output.readcounts.bg/result (last column is used!)
my $DIR = $ARGV[1];
my $plotfile = $ARGV[2];
my %hist;
my %histcum;
my %histscum;
my %histscum_log;


my $REM = 10**(-4);
my $MAX = 4000;
my $SCALE = 5;

my %windows;

my $counter = 0;
print STDERR "Creating histogram of background counts\n";
open F, $FILE_COUNTS; # file with background counts (outzvals)
while(<F>) {

    chomp($_);
    my @s = split(/\t/,$_);

    my $count = sprintf("%.0f",$s[-1]); # last column has to contain counts!
    
    for(my$i=$count;$i>=0;$i--) { $histcum{$i}++; }
    $hist{$count}++;

    #$windows{$counter} = $s[-1];
    #$counter++;
}
close F;

# file non-existing bins with 0
for(my$i=0;$i<$MAX;$i++) { if( !exists($hist{$i}) ) { $hist{$i} = 0.0; } }
for(my$i=0;$i<$MAX;$i++) { if( !exists($histcum{$i}) ) { $histcum{$i} = 0.0; } }


print STDERR "Creating smoothed histogram of background counts\n";
for(my$i=0;$i<$MAX;$i++) {

    for(my$j=0;$j<$MAX;$j++) {

	my $dist = abs($j-$i);
	my $w = exp(-$dist/$SCALE);

	$histscum{$j} += $w*$histcum{$i};
    }
}

for(my$i=0;$i<$MAX;$i++) { 
    if($histscum{$i}>0) { $histscum_log{$i} = log( $histscum{$i} ); }
    if($histscum{$i}==0) { $MAX=$i-1; }
}

open W, ">".$DIR."/hist.summedbgcounts";
for(my$i=0;$i<$MAX;$i++) { print W join("\t",$i,$hist{$i},$histscum_log{$i},$histscum{$i},$histcum{$i})."\n";  }
close W;




### analysis 


my $score = \%histscum_log;
my $BOUND = 25;
my %distance;
my %info;

my @diffs = ();

for(my$i=$BOUND;$i<($MAX-$BOUND);$i++) {

    my $left = 0;
    my $right = 0;

    my $c = 0;
    for(my $j=($i-$BOUND);$j<$i;$j++) {
	$left += $score->{$j};
	$c++;
    }
    $left /= $c;
    
    $c = 0;
    for(my $j=$i;$j<($BOUND+$i);$j++) {
	$right += $score->{$j};
	$c++;
    }
    $right /= $c;

    my $dist = ($left-$right) / ($left+$right);
    my $slope = -($left-$right)/$BOUND;
    my $btmp = ($left+$right)/2;
    my $b = $btmp - $slope*$i;

    $distance{$dist} = join("\t",$i,$slope,$b);

#    print $distance{$dist}."\n";
    push(@diffs,$distance{$dist});
   
#    print join("\t",'bg= '.$i,'dist= '.$dist,'slope= '.$slope,'b= '.$b,'l= '.$left,'r= '.$right)."\n";
}

my $final_cut = undef;
my $bgcut = 0;
my $slope = 0;
my $bee = 0;
for(my$i=0;$i<(@diffs+1);$i++) {
    
    my ($in,$slopen,$bn) = split(/\t/,$diffs[$i]);
    my ($inp1,$slopenp1,$bnp1) = split(/\t/,$diffs[$i+1]);

    if( ($slopen-$slopenp1)<0 ) {
	$bgcut = $in;
	$slope = $slopen;
	$bee = $bn;
	last;
    }
}


sub g_of_x {

    my $s = $_[0];
    my $b = $_[1];
    my $x = $_[2];

    return ($slope*$x+$b);
}

# mean & var around bgcut;
my $mean = 0;
my $absmean = 0;
my $tmp = 0;


my $AROUND = 50;
if($AROUND>$bgcut) { $AROUND = ($bgcut-1); }
$counter = 0;
for(my$i=-$AROUND+$bgcut;$i<($bgcut+$AROUND);$i++) {

    $mean += $histscum_log{$i}-&g_of_x($slope,$bee,$i);
    $absmean += abs( $histscum_log{$i}-&g_of_x($slope,$bee,$i) );
    $tmp += ($histscum_log{$i}-&g_of_x($slope,$bee,$i))**2;
    $counter++;
}

$mean /= $counter;
$absmean /= $counter;
my $var = sqrt( (1.0*$tmp)/($counter*1.0) - $mean**2 );

print STDERR "selected: $bgcut\n";
print STDERR "mean: $mean\tabsmean: $absmean\n";
print STDERR "sd: $var\n";



open W, ">".$DIR."/stats_bg_cutoff";
for(my$i=$bgcut;$i<$MAX;$i++) {

    my $delta = $histscum_log{$i}-&g_of_x($slope,$bee,$i);

#    print $delta."\n";

    if( abs($delta) > (10*$absmean)  ) {
	$i+=5;
	$final_cut = $i;
	print W $final_cut."\n";
	last;
    }
}

my $total = $histcum{0};
my $atcut = $histcum{$final_cut};

my $perc = ($atcut*1.0)/($total*1.0)*100;

print W join(' ','removed: ',$perc.'%','(',$total,$atcut,')')."\n";
if( !( $perc>(10**(-5)) && $perc<(10**(-3)) ) ) {
    
    print W "Warning: something might be wrong with the selected background cut-off!\n";
}


close W;


# create output file

open R, ">".$DIR.'/make_plot.R';

print R 'pdf("'.$plotfile.'")'."\n";
print R 'd = read.csv(file="'.$DIR.'/hist.summedbgcounts",sep="\t")'."\n";
print R 'mymax=1000+'.$final_cut."\n";
print R 'plot(d[,1],d[,3],t="l",lwd=3,xlab="summed background counts",ylab="density [log]", xlim=c(0,mymax),ylim=c(4,20))'."\n";
print R 'y = '.$slope.' * d[,1] + '.$bee.''."\n";
print R 'lines(d[,1],y)'."\n";
print R 'bbb='.$final_cut."\n";
print R 'lines( c(bbb,bbb), c(-100,1000), lwd=3, col="red" )'."\n";
print R 'dev.off()'."\n";

close R;

system('R CMD BATCH '.$DIR.'/make_plot.R '.$DIR.'/make_plot.Rout');
