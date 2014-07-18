#!/import/bc2/soft/app/perl/5.10.1/Linux/bin/perl -w
use strict;
use Getopt::Long;
# If you modify this script, it is important not to use locale.

# depends on:
# /import/bc2/home/nimwegen/GROUP/DeepSeqPipeline.Piotr/soft/pioSortBed9
# /import/bc2/home/nimwegen/GROUP/DeepSeqPipeline.Piotr/soft/bedmap

#### USAGE ####
# ./binReads.pl -b bgfile_1.bedweight[.gz] ... -b bfile_n -f fgfile_1 ... -f fgfile_n -l chromLenFile --fgwin=500 --shift=250 --bgwin=2000 --outfile=result
# optionally: --sorted --fmi-repo=...

# writes tmp files to the current directory and deletes them in the end.
# because it calls pioSortBed9, it needs to be run on a high memory node.

# alternative binning using bedtools.
my @fgFs;
my @bgFs;
my $chrLenF = "";
my $fgwin;
my $bgwin;
my $fgShift;
my $outF = '';
my $fmiRepo = '/import/bc2/home/nimwegen/GROUP/DeepSeqPipeline.Piotr';
my $res = GetOptions(
                     "b=s" => \@bgFs,
                     "f=s" => \@fgFs,
                     "l=s" => \$chrLenF,
                     "fgwin=i" => \$fgwin,
                     "shift=i" => \$fgShift,
                     "bgwin=i" => \$bgwin,
                     "outfile=s" => \$outF,
                     "fmi-repo=s" => \$fmiRepo,	#optional
                    );
die "Options missing\n" if ! $chrLenF || (! $fgwin && ! $bgwin) || ! $fgShift || ! $outF || (! @fgFs && ! @fgFs);

# prepare bins for the fg and bg centered at the same positions:
my %lenOchr;
open(my $cH, $chrLenF) or die;
while(<$cH>)
{
	chomp:
	my($chr, $len) = split(/\s+/);
	$lenOchr{$chr} = $len;
}
close $cH;
my $i = 0;

open(my $ofH, '>', "fgWindows.bed") or die;
open(my $obH, '>', "bgWindows.bed") or die;
foreach my $chr (sort keys %lenOchr)
{
	for(my $beg=0; $beg+$fgwin <= $lenOchr{$chr}; $beg+=$fgShift)
	{
		my $center = $beg+int($fgwin/2);
		print $ofH join("\t", $chr, $beg, $beg+$fgwin, "win.$chr.$center", 0, '+'), "\n";
		print $obH join("\t", $chr, $center-int($bgwin/2) < 0 ? 0 : $center-int($bgwin/2), $center+int($bgwin/2), "win.$chr.$center", 0, '+'), "\n";
	}
}
close $ofH;
close $obH;

# run counting:
foreach my $f (@fgFs)
{
	# sort, truncate to the 5' end and stack reads so that bedtools run much faster on the sorted data:
	system("zcat '$f' | $fmiRepo/soft/pioSortBed9 -s5 --collapse --input-file - > fgCollapsed.$i");
	system("$fmiRepo/soft/bedmap --faster --indicator --sum fgWindows.bed fgCollapsed.$i > fgBinned.$i");
	$i++;
}
$i = 0;
foreach my $f (@bgFs)
{
	system("zcat -f '$f' |  $fmiRepo/soft/pioSortBed9 -s5 --collapse --input-file - > bgCollapsed.$i");
	system("$fmiRepo/soft/bedmap --faster --indicator --sum bgWindows.bed bgCollapsed.$i > bgBinned.$i");
	$i++;
}
open(STDOUT, '>', $outF) or die $!;

# join them now together:
my @fHs;
$i = 0;
foreach my $f (@fgFs)
{
	open(my $h, "fgBinned.$i") or die $!;
	push @fHs, \$h;
	$i++;
}
$i = 0;
foreach my $f (@bgFs)
{
	open(my $h, "bgBinned.$i") or die $!;
	push @fHs, \$h;
	$i++;
}
open(F, 'fgWindows.bed') or die;

print "#chr\tstart\tstop\tmiddle";
foreach(my $i=0; $i<@fgFs; $i++) { print "\t", "fg_$i" }
foreach(my $i=0; $i<@bgFs; $i++) { print "\t", "bg_$i" }
print "\n";

for(my $i=0; 1; $i++)
{
	my @counts = ();
	my $flagNonZero = 0;
	# parse the window part:
	$_ = <F>; last if ! $_; chomp;
	my ($chr, $beg, $end, $id, $count, $str) = split(/\t/);
	$id =~ /win\.\S+\.(\d+)/ or die;
	my $lineBeg = join("\t", $chr, $beg, $end, $1);
	
	# parse all the count files:
	for(my $j=0; $j<@fHs; $j++)
	{
		my $h = ${$fHs[$j]};
		$_ = <$h>;
		chomp;
		my ($flag, $count) = split(/\|/);
		$count = $flag ? $count : '0';
		push @counts, $count;
		if ( $count  )
		{
			$flagNonZero = 1;
		}
	}
	if($flagNonZero)
	{
		print join ("\t", $lineBeg, @counts), "\n";
	}
}

close F;
foreach my $h (@fHs)
{
	close $$h;
}

close STDOUT;

# clean after it is finished:
`rm -f bgBinned.* bgCollapsed.*`;
`rm -f fgBinned.* fgCollapsed.*`;
`rm -f fgWindows.bed bgWindows.bed`;


