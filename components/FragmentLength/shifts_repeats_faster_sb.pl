#!/usr/bin/perl

#We will assume .bed file format
# To fasten up this code:
# First put all positions from the input file (unique, sorted bedweight file) in one array. In a second array store the strand of each position.
# Go over all positions of the first array. If it is on the positive strand go maxLen forward and store the distances to the negative strand positions.

my $infile = shift(@ARGV);
my $outputfilename = shift(@ARGV);
my $repeatPath = shift(@ARGV);
print STDERR "$infile";
my @v = split(/\//, $infile);
my $file = pop @v;
my ($sample) = split(/\./, $file);
my $maxLen=600;

print STDERR "sample is $sample\n";

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
    my @read_positions;
    my @read_strands;

    $curnext   = 1;
    $numrepeat = @ends;

    my $seen_it = 0;
    $last = -999;
    my @lines = ();
    while (<F>)
    {
        #put before/after (seen_it) flag. Since file is sorted we can break the file loop once we passed the chromosome
        #check this line is the right chromosome

        if ($_ =~ /^$chrom\s+/)
        {
            $seen_it = 1;
            chomp($_);
            push @lines, $_;
        }

        elsif ($seen_it) {
            last;
        }
    }
    close(F);

    #sort should not be needed if the bedweight file is sorted
    #@ol = sort { sortbystart($a, $b) } @lines;

    $numlines = @lines; #@ol;
    for ($i = 0 ; $i < $numlines ; ++$i)
    {
        @s = split(/\s+/, @lines[$i]); #$ol[$i]);
        # Note: input reads are sorted by 5' end. Also use 5' end to check if it is within a repeat.
        if ($s[5] eq "+") {
            $read_5prime = $s[1];
        }
        else {
            $read_5prime = $s[2];
        }

        if ($read_5prime < $last)
        {
            die "error at line $i $last $s[1] $_\n";
        }
        $last = $read_5prime;

        #position is $s[1];
        #check if position before first or after last
        if ($read_5prime <= $starts[0] || $read_5prime > $ends[$numrepeat - 1])
        {
            $inrepeat = 0;
	    $numgood+=1
        }
        else
        {
            while ($starts[$curnext] <= $read_5prime && $curnext < ($numrepeat - 1))
            {
                ++$curnext;
            }
            if ($ends[$curnext - 1] >= $read_5prime)
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
            push @read_positions, $read_5prime;
            push @read_strands, $s[5];
        }

    }


    #now count relative occurrences
    $numreads = @read_positions;
    for ($i = 0 ; $i < $numreads ; ++$i){

        if ($read_strands[$i] eq "+")
            { $plus_pos = $read_positions[$i];
              $j = 1;
              while (1) {

                  if ($i + $j >= $numreads) {
                      last;
                  }

                  $cur_pos = $read_positions[$i+$j];
                  $rel_pos = $cur_pos - $plus_pos;

                  if ($rel_pos > $maxLen) {
                      last;
                  }

                  if ($read_strands[$i+$j] eq "-") {
                      ++$count[$rel_pos]
                  }

                  ++$j;
            }
        }
    }
    ##
  print "number inside repeats: $numrepeats\tnumber outside repeats $numgood\n";

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

