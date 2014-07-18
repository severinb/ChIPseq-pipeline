#!/import/bc2/soft/app/perl/5.10.1/Linux/bin/perl -w
use strict;
use Getopt::Long;
# If you modify this script, it is important not to use locale.

# depends on:
# /import/bc2/home/nimwegen/GROUP/DeepSeqPipeline.Piotr/soft/pioSortBed9
# /import/bc2/home/nimwegen/GROUP/DeepSeqPipeline.Piotr/soft.bc2/countReads.pl
# /import/bc2/home/nimwegen/GROUP/DeepSeqPipeline.Piotr/soft/bedtools (through countReads.pl)

#### USAGE ####
# ./binReads.pl -b bgfile_1.bedweight[.gz] ... -b bfile_n -f fgfile_1 ... -f fgfile_n -l chromLenFile --fgwin=500 --shift=250 --bgwin=2000 --outfile=result
# optionally: --sorted --fmi-repo=... --bedtools-dir=...

# I would use --sorted because bedtools are very inefficient in memory management. And there is no big speedup when you not use --sorted
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
my $bedToolsPath = '/import/bc2/home/nimwegen/GROUP/DeepSeqPipeline.Piotr/soft';
my $sorted = "";  # by default don't use the sweep algorithm, enable to save memory when you know that your files are sorted
my $res = GetOptions(
                     "b=s" => \@bgFs,
                     "f=s" => \@fgFs,
                     "l=s" => \$chrLenF,
                     "fgwin=i" => \$fgwin,
                     "shift=i" => \$fgShift,
                     "bgwin=i" => \$bgwin,
                     "outfile=s" => \$outF,
                     "fmi-repo=s" => \$fmiRepo,	#optional
                     "bedtools-dir=s" => \$bedToolsPath,	#optional
                     "sorted" => \$sorted	# optional, use the alternative bedtools algorithm which scans files linearly
                    );
die "Options missing\n" if ! $chrLenF || ! $fgwin || ! $bgwin || ! $fgShift || ! $outF || ! @fgFs || ! @fgFs;
$sorted = $sorted ? '--sorted' : '';
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
	system("less '$f' | $fmiRepo/soft/pioSortBed9 -s5 --collapse --input-file - > fgCollapsed.$i");
	system("$fmiRepo/soft.bc2/countReads.pl", '--reads', "fgCollapsed.$i",
		'--regions', 'fgWindows.bed', $sorted, '-5', '--output', "fgBinned.$i",
		'--bedtools-dir', $bedToolsPath); #  '--no-zeros',
	$i++;
}
$i = 0;
foreach my $f (@bgFs)
{
	system("less '$f' |  $fmiRepo/soft/pioSortBed9 -s5 --collapse --input-file - > bgCollapsed.$i");
	system("$fmiRepo/soft.bc2/countReads.pl", '--reads', "bgCollapsed.$i",
		'--regions', 'bgWindows.bed', $sorted, '-5', '--output', "bgBinned.$i",
		'--bedtools-dir', $bedToolsPath);
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

print "#chr\tstart\tstop\tmiddle";
foreach(my $i=0; $i<@fgFs; $i++) { print "\t", "fg_$i" }
foreach(my $i=0; $i<@bgFs; $i++) { print "\t", "bg_$i" }
print "\n";
loopTop:
for(my $i=0; 1; $i++)
{
	my @counts = ();
	my $lineBeg = "";
	my $flagNonZero = 0;
	for(my $j=0; $j<@fHs; $j++)
	{
		my $h = ${$fHs[$j]};
		$_ = <$h> || last loopTop;
		chomp;
		my ($chr, $beg, $end, $id, $count, $str) = split(/\t/);
		if(!$j)
		{
			$id =~ /win\.\S+\.(\d+)/ or die;
			$lineBeg = join("\t", $chr, $beg, $end, $1);
		}
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

foreach my $h (@fHs)
{
	close $$h;
}

close STDOUT;

# clean after it is finished:
#`rm -f bgBinned.* bgCollapsed.*`;
#`rm -f fgBinned.* fgCollapsed.*`;
#`rm -f fgWindows.bed bgWindows.bed`;


