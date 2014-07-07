use warnings;
use strict;

my $T = 5000;
my %hist;
for(my$i=0;$i<$T;$i++) { $hist{$i} = 0; }


# histogram
print STDERR "Creating histogram of background counts\n";
open F, $ARGV[0]; # file with background counts (outzvals)
while(<F>) {

    chomp($_);
    my @s = split(/\t/,$_);

    my $count = sprintf("%.0f",$s[$ARGV[1]]); # which column (starting from 0) contain the background counts (5)

#    for(my $i=0;$i<=$count;$i++) { $hist{$i}++; }
    $hist{$count}++;

}
close F;


my $LIM = 20;
if( defined($ARGV[2]) ) {
    $LIM = $ARGV[2];
}
my %d;

print STDERR "Getting distances ...\n";
for(my$i=$LIM;$i<($T-$LIM);$i++) { 


    my $left = 0;
    my $right = 0;

    for(my $j=$i-$LIM;$j<$i;$j++) { $left += $hist{$j}; }
    for(my $j=$i;$j<($i+$LIM);$j++) { $right += $hist{$j}; }

    if( ($left+$right)<100 ) { next; }

    if( ($left+$right)==0 ) { next; }


    my $diff = ($left-$right) / ($left+$right);

    $d{$diff} = $i;
    #print join("\t",$i,$diff)."\n";

}

my $bgcut;
for my $diff ( sort {$a <=> $b} keys %d ) {

    
    $bgcut = $d{$diff};
}

# usually, the bg cut-off estimated above is too small
$bgcut += 10;

my $tot = 0;
my $removed = 0;

open W, ">hist.bgcounts";
for(my$i=0;$i<100000;$i++) {
    
    if( !exists($hist{$i}) ) { 
	print W join("\t",$i,0.0)."\n";
	next; 
    }
    else { print W join("\t",$i,$hist{$i})."\n"; }

    $tot += $hist{$i};
    if( $i>=$bgcut ) {
	$removed += $hist{$i};
    }

}
close W;
my $ratio = sprintf("%.2f",$removed/$tot*100);

open W, ">estimated_bg_cutoff";
print W join("\t","bg cut-off: ",$bgcut)."\n";
print W join("\t","removed windows:",$ratio,"% ($removed)")."\n";
close W;

