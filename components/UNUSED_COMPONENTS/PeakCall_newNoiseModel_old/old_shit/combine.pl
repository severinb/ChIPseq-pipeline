use warnings;
use strict;

my $samples_fg = 0;
my $samples_bg = 0;

my $STEP = $ARGV[3]; # step size
my $SIZEFG = $ARGV[4]; # win sizefg
my $SIZEBG = $ARGV[5]; # win sizefg
my $PATH = $ARGV[6]; # output directory
my $CUTBG = $ARGV[7]; # negative number: read from file: stats_bg_cutoff, otherwise use specified number as bg cut-off


if( !defined($ARGV[7]) || $ARGV[7]<0 ) {

    print STDERR "Loading from file: stats_bg_cutoff\n";

    open F, $PATH.'/stats_bg_cutoff';
    while(<F>) {
	
	chomp($_);
	$CUTBG = $_;
	last;
    }
    close F;
}
else {
    print STDERR "Using specified cut-off: $CUTBG\n";
}

    
my %one;
open F, $ARGV[0]; # results fg
while(<F>) {

    chomp($_);
    my @s = split(/\t/,$_);
    my $chr = $s[0];
    my $middle = $s[3];

    $samples_fg = (@s-5);

    my $count = $s[4];
    for(my$i=5;$i<(@s-1);$i++) {

	$count .= "\t".$s[$i];
    }
    $one{$chr}{$middle} = $count;

}
close F;


my %two;
open F, $ARGV[1]; # results bg
while(<F>) {

    chomp($_);
    my @s = split(/\t/,$_);
    my $chr = $s[0];
    my $middle = $s[3];
    my $count = $s[4];

    $samples_bg = (@s-5);
    
    for(my$i=5;$i<(@s-1);$i++) {

	$count .= "\t".$s[$i];
    }
    $two{$chr}{$middle} = $count;
}
close F;

my %chrs;

open F, $ARGV[2]; # file with chrom information
while(<F>) {

    chomp($_);
    my ($chr,$length) = split(/\t/,$_);

    $chrs{$chr} = $length;

}
close F;



print STDERR "Number of samples:\n";
print STDERR "\t-fg: $samples_fg\n";
print STDERR "\t-bg: $samples_bg\n\n";

open W, ">".$PATH.'/all_counts';
print W join("\t",'#chr','start','stop','middle');
for(my$i=0;$i<$samples_fg;$i++) { print W "\t".'fg_'.$i; }
for(my$i=0;$i<$samples_bg;$i++) { print W "\t".'bg_'.$i; }
print W "\n";

my $total = 0;
my $rejected = 0;

for my $chr (keys %chrs) {

    my $length = $chrs{$chr};

    for(my $p=0;$p<($length-$SIZEFG);$p+=$STEP) {

	my $middle = int( $p+$SIZEFG/2 );

	if( ($middle < ($SIZEBG/2)) || ($middle > ($chrs{$chr}-$SIZEBG/2)) ) { next; }

	my $countsfg = "0"; for(my$i=1;$i<$samples_fg;$i++) { $countsfg .= "\t0"; }
	my $countsbg = "0"; for(my$i=1;$i<$samples_bg;$i++) { $countsbg .= "\t0"; }

	if( exists( $one{$chr}{$middle} ) ) { $countsfg = $one{$chr}{$middle}; }
	if( exists( $two{$chr}{$middle} ) ) { $countsbg = $two{$chr}{$middle}; }

	my @f = split(/\t/,$countsfg);
	my @b = split(/\t/,$countsbg);

	my ($zerof,$zerob) = (1,1);

	my $summedbg = 0;
	for(my$i=0;$i<@f;$i++) { if($f[$i]>0) { $zerof=0; } }
	for(my$i=0;$i<@b;$i++) { if($b[$i]>0) { $zerob=0; } $summedbg += $b[$i]; }

	if( $zerof==1 && $zerob==1 ) { next; }

	$total++;
	if( $summedbg > $CUTBG ) { 
	    $rejected++;
	    next; 
	}
	
	print W join("\t",$chr,$p,$p+$SIZEFG,$middle,$countsfg,$countsbg)."\n";
	
    }
}
close W;

my $per = sprintf("%.2f",$rejected/$total*100.);
open S, ">".$PATH.'/stats_combine_windows';
print S "cut-off on summed bg windows: $CUTBG\n";
print S "total number of windows: $total\n";
print S "rejected number of windows: $rejected ($per %)\n";
close S;
