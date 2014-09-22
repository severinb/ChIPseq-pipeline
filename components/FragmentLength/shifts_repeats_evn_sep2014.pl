#!/usr/bin/perl

#We will assume .bed file format

#Background chromatin sample
#/import/bc2/home/nimwegen/GROUP/ChIP_Sequencing_Data/ChIP_seq/BG.chromatin/oneline/MmES_Input_chromatin.mm9-mmV02-aln2.m100.bed
#K27 WT files
#/import/bc2/home/nimwegen/GROUP/ChIP_Sequencing_Data/ChIP_seq/ESC_NP_TN/H3K27me3/onefiles/
# The ab file is the concatenation of the a and b files. Same position can thus occur twice in the ab file

#K27 KO files
#/import/bc2/home/nimwegen/GROUP/ChIP_Sequencing_Data/ChIP_seq/ESC_NP_TN/H3K27me3_RESTko/
#Use the files with unique.bed

#File with mouse chromosome IDs
#/import/bc2/home/nimwegen/nimwegen/PhilChIPseq/all_chrom_ids

#Repeat annotation
#/import/bc2/home/nimwegen/nimwegen/PhilChIPseq/RepeatAnno
#files end with .out

#mouse CAGE promoters /import/bc2/home/nimwegen/GROUP/GNP/mm9_29.09.2008/regions

#$infile = "/import/bc2/home/nimwegen/GROUP/ChIP_Sequencing_Data/ChIP_seq/BG.chromatin/oneline/MmES_Input_chromatin.mm9-mmV02-aln2.m100.bed";
#$infile = "/import/bc2/home/nimwegen/GROUP/ChIP_Sequencing_Data/ChIP_seq/ESC_NP_TN/H3K27me3/onefiles/MmES_K27me3_b.mm9-mmV01-aln1.m100.bed";

my $infile = shift(@ARGV);
my $outputfilename = shift(@ARGV);
my $multi = shift(@ARGV);
my $repeatPath = shift(@ARGV);
print STDERR "$infile";
my @v = split(/\//, $infile);
my $file = pop @v;
my ($sample) = split(/\./, $file);
my $maxLen=600;


print STDERR "sample is $sample\n";

#$chromfile="/import/bc2/home/nimwegen/GROUP/hseq_pipeline/general/repeatMask/".$genome."/all_chrom_ids";
$chromfile=$repeatPath."/all_chrom_ids";

#get all chromosome IDs
open(F, $chromfile)
  || die "cannot open chromosome file $chromfile\n";
while (<F>)
{
    if ($_ =~ /^(chr\S+)\s*/)
    {
        $allID{$1} = 1;
  #      print STDOUT "$1\n"
    }
}
close(F);

#histogram correlation function
for ($i = 0 ; $i <= $maxLen ; ++$i)
{
    $count[$i] = 0;
}

#total number of mappers on plus and minus strand
$totplus     = 0;
$totmin      = 0;
$tot_revcomp = 0;

#go over all chromosomes
foreach $chrom (keys(%allID))
{

    print STDERR "doing chromosome $chrom\n";

    my @starts = ();
    my @ends   = ();
    $numrepeats=0;
    $numgood=0;
    #$repeatfile =
     #   "/import/bc2/home/nimwegen/GROUP/hseq_pipeline/general/repeatMask/" . $genome . "/" 
     # . $chrom
     # . ".fa.out";
    $repeatfile = $repeatPath."/".$chrom.".fa.out";
    open(F, $repeatfile) || die "cannot open file with repeat $repeatfile\n";

    #print STDERR "reading repeats\n";
    while (<F>)
    {
        if ($_ =~ /$chrom/)
        {
            chomp($_);

            #remove initial spaces
            $_ =~ s/^\s*//;
            @s = split(/\s+/, $_);
            push @starts, $s[5];
            push @ends,   $s[6];
        }
    }
    close(F);

    open(F, $infile) || die "cannot open $infile\n";
    my %readplus;
    my %readmin;

    $curnext   = 1;
    $numrepeat = @ends;

    $last = -999;
    my @lines = ();
    while (<F>)
    {

        #check this line is the right chromosome
        if ($_ =~ /^$chrom\s+/)
        {
            chomp($_);
            push @lines, $_;
        }
    }
    close(F);
    #sort should not be needed if the bedweight file is sorted
    #@ol = sort { sortbystart($a, $b) } @lines;
    $numlines = @lines; #@ol;
    for ($i = 0 ; $i < $numlines ; ++$i)
    {
        @s = split(/\s+/, @lines[$i]); #$ol[$i]);
        if ($s[1] < $last)
        {
            die "error at line $i $last $s[1] $_\n";
        }
        $last = $s[1];

        #position is $s[1];
        #check if position before first or after last
        if ($s[1] <= $starts[0] || $s[1] > $ends[$numrepeat - 1])
        {
            $inrepeat = 0;
	    $numgood+=1
        }
        else
        {
            while ($starts[$curnext] <= $s[1] && $curnext < ($numrepeat - 1))
            {
                ++$curnext;
            }
            if ($ends[$curnext - 1] >= $s[1])
            {
                $inrepeat = 1;
		$numrepeats=$numrepeats+1;
#		print "omitted due to repeat\n";
            }
            else
            {
                $inrepeat = 0;
		$numgood=$numgood+1;
            }
        }

        #NOTE: This assumes we are not doing multi
        if (!$inrepeat) {
            if ($s[5] eq "+") {
                $pos = $s[1];
                $readplus{$pos} = 1;
            }
            else {
                $pos = $s[2];
                $readmin{$pos} = 1;
            }
        }
    }

    #now count relative occurrences
    foreach $pos (keys(%readplus))
    {
        for ($relpos = 0 ; $relpos <= $maxLen ; ++$relpos)
        {
            $abspos = $pos + $relpos;
            if (defined($readmin{$abspos}))
            {
                ++$count[$relpos];
            }
        }
    }

  print "number reads: $numlines\n";
  print "number repeats:  $numrepeats,  num not $numgood\n";

}

#print "totals were $totplus $totmin\n";
#$outfile = $sample . ".reg_cor";
$outfile=$outputfilename;
open(G, ">", $outfile) or die;
for ($j = 0 ; $j <= $maxLen ; ++$j)
{
    $frac = $count[$j];
    print G $j, "\t", $frac, "\n";
}
close(G);

sub sortbystart
{
    my ($aa, $bb) = @_;
    my @sa = split(/\s+/, $aa);
    my @sb = split(/\s+/, $bb);
    return $sa[1] <=> $sb[1];
}

sub inrepeat
{
    my $pos    = shift(@_);
    my @cs     = @_;
    my $len    = @cs;
    my $bottom = 0;
    my $top    = $len - 1;

    #print "got position $pos\n";
    #print "bottom $bottom value ", $cs[$bottom], " top $top value ", $cs[$top], "\n";

    if ($cs[$bottom] > $pos)
    {
        return -1;
    }
    elsif ($cs[$top] < $pos)
    {
        return $top;
    }

    while ($top - $bottom > 1)
    {
        my $middle = int(($top + $bottom) / 2);
        if ($cs[$middle] > $pos)
        {
            $top = $middle;
        }
        elsif ($cs[$middle] <= $pos)
        {
            $bottom = $middle;
        }
        else
        {
            die "bottom $bottom top $top middle $middle pos $pos values",
              $cs[$bottom], " ", $cs[$top], " ", $cs[$middle], "\n";
        }

        #print "bottom $bottom value ", $cs[$bottom], " top $top value ", $cs[$top], "\n";
    }
    return $bottom;
}

#Fitting peak
#find optimum starting at 60
#Maybe this has an impact on the dynamic range

