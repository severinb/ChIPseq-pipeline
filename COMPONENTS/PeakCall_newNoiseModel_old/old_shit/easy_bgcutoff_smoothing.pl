use warnings;
use strict;

# INPUT

my $FILE_COUNTS = $ARGV[0]; # file containing (summed) background counts: output.readcounts.bg/result (last column is used!)
my $DIR = $ARGV[1];
my %hist;
my %histcum;
my %histscum;
my %histscum_log;

my $MAX = 4000;
my $SCALE = 5;

print STDERR "Creating histogram of background counts\n";
open F, $FILE_COUNTS; # file with background counts (outzvals)
while(<F>) {

    chomp($_);
    my @s = split(/\t/,$_);

    my $count = sprintf("%.0f",$s[-1]); # last column has to contain counts!
    
    for(my$i=$count;$i>=0;$i--) { $histcum{$i}++; }
    $hist{$count}++;
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


open W, ">hist.summedbgcounts";
for(my$i=0;$i<$MAX;$i++) { print W join("\t",$i,$hist{$i},$histscum_log{$i},$histscum{$i},$histcum{$i})."\n";  }
close W;


### analysis 


my $score = \%histscum_log;
my $BOUND = 25;
my %distance;
my %info;

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

    #print $distance{$dist}."\n";
    
#    print join("\t",'bg= '.$i,'dist= '.$dist,'slope= '.$slope,'b= '.$b,'l= '.$left,'r= '.$right)."\n";
}

my $bgcut = undef;
my $slope = undef;
my $bee = undef;

for my $dist (sort {$a <=> $b} keys %distance) {

    ($bgcut,$slope,$bee) = split(/\t/,$distance{$dist});
}

print STDERR "Some statistics:\n";
print STDERR "\t-read at max: $bgcut\n";
print STDERR "\t-slope at max: $slope\n";
print STDERR "\t-b at max: $bee\n";
print STDERR "\t-->   g(x) = $slope * x + $bee\n\n";


sub g_of_x {

    my $s = $_[0];
    my $b = $_[1];
    my $x = $_[2];

    return ($slope*$x+$b);
}

my $mean = 0;
my $tmp = 0;
my $var = undef;
my $cnt = 0;
for(my$i=$bgcut-$BOUND;$i<($bgcut+$BOUND);$i++) {

    my $delta = $histscum_log{$i}-&g_of_x($slope,$bee,$i);

    $mean += $delta;
    $tmp += $delta**2;

    $cnt++;
    #print join("\t",$histscum_log{$i},&g_of_x($slope,$bee,$i))."\n";
}
$mean /= $cnt;
$tmp /= $cnt;
$var = $tmp-$mean**2;
my $sigma = sqrt($var);

print STDERR "Fitting: mean: $mean\tsigma: $sigma\n";

# now get the fucking bg cut-off
my $final_cut = undef;

open W, ">".$DIR."/stats_bg_cutoff";
for(my$i=$bgcut;$i<$MAX;$i++) {

    my $delta = $histscum_log{$i}-&g_of_x($slope,$bee,$i);
    
    if( $delta > 1 ) {
	$i+=5;
	$final_cut = $i;
	print W $final_cut."\n";
	last;
    }
}
close W;

print STDERR "Suggested/Used bg cut-off: $final_cut\n\n";

# create output file

open R, ">".$DIR.'/make_plot.R';

print R 'jpeg("'.$DIR.'/bgcut.jpg")'."\n";
print R 'd = read.csv(file="'.$DIR.'/hist.summedbgcounts",sep="\t")'."\n";
print R 'mymax=300+'.$final_cut."\n";
print R 'plot(d[,1],d[,3],t="l",lwd=3,xlab="summed background counts",ylab="density [log]", xlim=c(0,mymax))'."\n";
print R 'y = '.$slope.' * d[,1] + '.$bee.''."\n";
print R 'lines(d[,1],y)'."\n";
print R 'bbb='.$final_cut."\n";
print R 'lines( c(bbb,bbb), c(-100,1000), lwd=3, col="red" )'."\n";
print R 'dev.off()'."\n";

close R;

system('R CMD BATCH '.$DIR.'/make_plot.R');
