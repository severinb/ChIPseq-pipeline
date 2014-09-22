#!/usr/bin/env perl

use strict;
use warnings;
use FindBin '$Bin';
use Getopt::Long;

# this is alignment script with updated Bowtie version
# and with the sorting of RAL files optimized by Piotr
# based on original aln2.pl, by FMI

my $mthreads      = getConfig('Threads');                    # number of threads to use by oligomap and soap
my $m             = 100;                                     # maximum number of genomic hits per sequence
my $bowtieOptions = "-f -p $mthreads -v 2 -a -B 1 --quiet --best --strata"
  ;    # allow up to 2 mismatches, ignore quality scores, find all hits up to m. coord start 1
my $seqBins = [14, 20, 24];    # seqence lengths for 0,1 and 2 mismatches

my $dontPrintReport = 0;
my $results = GetOptions('r' => \$dontPrintReport);

if (@ARGV != 5)
{
    print STDERR "\nUSAGE: ./aln2.pl targetDir queryFile typeOfTarget tempDir outDir\n\n";
    print STDERR "   tempDir:  g: genome, a: annotation\n";
    print STDERR "   -r     :  dont print report\n\n";

    exit(0);
}

my $bowtieIndexSubdirName = "bowtieIndex";

my $repDir    = "$Bin/../..";
my $bowtieDir = "$repDir/soft/bowtie-0.12.7";

my ($targetDir, $queryFile, $typeOfTarget, $tempDir, $outDir) = ($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3], $ARGV[4]);

copyReadsIntoFastaFile($queryFile, "$tempDir/seqs.fa");

# map to genome
if ($typeOfTarget eq "g")
{

	
    # check if bowtie index exists. if not create one
    if (system("mkdir $targetDir/$bowtieIndexSubdirName 2> /dev/null") == 0)
    {

        # concatenate all chromosome files (.fa) into one big temporary file
        system("cat $targetDir/\*.fa > $targetDir/$bowtieIndexSubdirName/catChrs.fa");

        # build bowtie index
        my $ret1 = system(
            "$bowtieDir/bowtie-build $targetDir/$bowtieIndexSubdirName/catChrs.fa $targetDir/$bowtieIndexSubdirName/$bowtieIndexSubdirName > /dev/null"
        );
        system("rm $targetDir/$bowtieIndexSubdirName/catChrs.fa");
        if ($ret1 != 0)
        {
            system("rm \-rf $targetDir/$bowtieIndexSubdirName");    # delete index directoy if bowtie failed
            die "bowtie error when building the index\n";
        }
    }

    # map the reads to the genome
    my $ret2 = system(
        "$bowtieDir/bowtie $bowtieOptions -m $m --max $tempDir/seqsOverHits.fa $targetDir/$bowtieIndexSubdirName/$bowtieIndexSubdirName $tempDir/seqs.fa > $tempDir/seqs.bal"
    );
    if ($ret2 != 0) { die "bowtie error aligning the reads\n"; }

    # check if the file seqsOverHits.fa was created, otherwise create empty one (needed later)
    if (!(-f "$tempDir/seqsOverHits.fa"))
    {
        system("touch $tempDir/seqsOverHits.fa");
    }

    my $nrHitsPerSeq = splitAlignmentFilesIntoChromosomesAndConvertToRal("$tempDir/seqs.bal", $tempDir, $seqBins);

    # create hit report
    my $overHitSeqs = loadOneLineFasta("$tempDir/seqsOverHits.fa");
    open(F, ">$outDir/hitReport.tab")
      or die "Cannot open $outDir/hitReport.tab\n";
    for my $seqErrorId (sort keys %{$nrHitsPerSeq})
    {
        $seqErrorId =~ /^(\S+)\:(\d+)/;
        my ($seqId, $errors) = ($1, $2);
        print F "$seqId\t$errors\t$nrHitsPerSeq->{$seqErrorId}\n";
    }
    for my $seqId (sort keys %{$overHitSeqs}) { print F "$seqId\t0\t-1\n"; }
    close(F);

    # sort alignments for each chromosome alignment file and write to final destination folder
    opendir(RD, $tempDir) or die("Cannot open directory $tempDir");
    my @unsortedRalDir = readdir(RD);
    foreach my $fileName (@unsortedRalDir)
    {
        if ($fileName =~ /\.ral$/)
        {
            sortRAL("$tempDir/$fileName", "$outDir/$fileName");

            #my $ret3=system("sort $tempDir/$fileName -k 4 -n -S 4000000 -t\'".chr(9)."\' > $outDir/$fileName");
            #if($ret3 != 0){die "sort command does not work properly on that system\n";}
        }
    }

    # map to annotation
}
elsif ($typeOfTarget eq "a")
{

    # read the annotation files and cat them with a modified header (+filename, allows to split later)
    opendir(ANDIR, "$targetDir") or die("Cannot open directory $targetDir\n");
    my @annDbfiles = readdir(ANDIR);
    my @allAnnotFileNames;
    foreach (@annDbfiles)
    {
        if ($_ =~ /(\S+)\.fa$/)
        {
            my ($fileStem) = ($1);
            push(@allAnnotFileNames, $fileStem);
        }
    }

    # check if there are annotation files at all
    if (@allAnnotFileNames == 0)
    {

        #create empty hitreport if needed
        if (!($dontPrintReport)) { system("touch $outDir/hitReport.tab"); }
    }
    else
    {

        # check if bowtie index exists. if not create one
        if (system("mkdir $targetDir/$bowtieIndexSubdirName 2> /dev/null") == 0)
        {

            open(FOUT, ">$targetDir/$bowtieIndexSubdirName/catSeqs.fa")
              or die "cannot open $targetDir/$bowtieIndexSubdirName/catSeqs.fa\n";
            foreach (@annDbfiles)
            {
                if ($_ =~ /(\S+)\.fa$/)
                {
                    my ($fileStem) = ($1);

                    open(FIN, "$targetDir/$fileStem.fa")
                      or die "Cannot open $targetDir/$fileStem.fa\n";
                    while (<FIN>)
                    {
                        if (/^\s+$/) { next; }
                        if (/^\>(\S+)/)
                        {
                            print FOUT "\>$fileStem\/$1\n";
                        }
                        else
                        {

                            # convert U to T. bowtie will not find hits to target seqs with U
                            $_ =~ tr/Uu/Tt/;
                            print FOUT $_;
                        }
                    }
                    close(FIN);
                }
            }
            close(FOUT);

            # build bowtie index
            my $ret1 = system(
                "$bowtieDir/bowtie-build $targetDir/$bowtieIndexSubdirName/catSeqs.fa $targetDir/$bowtieIndexSubdirName/$bowtieIndexSubdirName > /dev/null"
            );
            system("rm $targetDir/$bowtieIndexSubdirName/catSeqs.fa");
            if ($ret1 != 0) { die "bowtie error when building the index\n"; }
        }

        # map the reads to annotation
        my $ret2 = system(
            "$bowtieDir/bowtie $bowtieOptions  $targetDir/$bowtieIndexSubdirName/$bowtieIndexSubdirName $tempDir/seqs.fa > $tempDir/seqs.bal"
        );
        if ($ret2 != 0) { die "bowtie error aligning the reads\n"; }

        foreach (@allAnnotFileNames)
        {
            system("touch $tempDir/$_\.ral");
        }    # create empty files

        my $maxHitsPerTarget =
          splitAlignmentFilesIntoAnnotFilesAndConvertToRal("$tempDir/seqs.bal", $tempDir, $seqBins);

        # sort alignments for each annotation file and write to final destination folder
        for my $annotFileName (@allAnnotFileNames)
        {
            my $ret = fastSortRAL("$tempDir/$annotFileName\.ral", "$outDir/$annotFileName\.ral");
            if($ret)
            {
				die "Sorting of file $tempDir/$annotFileName\.ral failed\n";
			}

            #my $ret3=system("sort $tempDir/$annotFileName\.ral  -k 2,2 -k 4,4n -S 4000000 -t\'".chr(9)."\' > $outDir/$annotFileName\.ral");
            #if($ret3 != 0){die "sort command does not work properly on that system\n";}
        }

        #create hitreport if needed
        if (!($dontPrintReport))
        {
            open(F, ">$outDir/hitReport.tab")
              or die "Cannot open $outDir/hitReport.tab\n";
            for my $seqErrorId (sort keys %{$maxHitsPerTarget})
            {
                $seqErrorId =~ /^(\S+)\:(\d+)/;
                my ($seqId, $errors) = ($1, $2);
                print F "$seqId\t$errors\t$maxHitsPerTarget->{$seqErrorId}\n";
            }
            close(F);
        }
    }

}
else { die "error: unknown typeOfTarget parameter $typeOfTarget\n"; }

# reads the seqs.tab file and copies the sequences into a new file in fasta format
sub copyReadsIntoFastaFile
{
    my ($seqsTabFile, $outFastaFile) = @_;

    open(FO, ">$outFastaFile") or die "Cannot open $outFastaFile\n";
    open(FI, $seqsTabFile)     or die "Cannot open $seqsTabFile \n";
    while (<FI>)
    {
        chomp();
        if (/^(\S+)\s+(\S+)\s+(\S+)/)
        {
            my ($id, $seq, $counts) = ($1, $2, $3);
            print FO ">$id\n$seq\n";
        }
        else
        {
            die "error: input file $seqsTabFile does not contain one sequence per line with counts.\n";
        }
    }
    close(FI);
    close(FO);
}

# reads the bowtie alignment output and splits the alignments according to chromosomes
# it also returns stats about the number of hits for every sequence
sub splitAlignmentFilesIntoChromosomesAndConvertToRal
{
    my ($balFile, $outDir, $seqBins) = @_;

    my %nrHitsPerSeq;
    my %outFilesOpened;    # contains all file handles opened at a given time
    open(F, $balFile) or die "Cannot open $balFile\n";
    while (<F>)
    {
        my @l = split("\t");
        if (/^(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t\S+\t\S+\t(\S*)$/)
        {
            my ($seqId, $strand, $chr, $from, $readSeq, $mmInfo) = ($1, $2, $3, $4, $5, $6);

            # convert to ral and print
            my $to         = $from + length($readSeq) - 1;
            my @mmInfoList = split("\,", $mmInfo);
            my $errors     = @mmInfoList;

            if (length($readSeq) < $seqBins->[$errors])
            {
                next;
            }    # check if read is long enough to allow the errors that it has

            if (!(defined $outFilesOpened{$chr}))
            {    # open new filehandle if new chromosome
                open(my $out, ">$outDir/$chr\.ral")
                  or die "Couldt read $outDir\n";
                $outFilesOpened{$chr} = $out;
            }

            # parse the mismatch field
            my @mmInfoListDetail;
            foreach (@mmInfoList)
            {
                if (/^(\d+)\:(\S)\>(\S)$/)
                {
                    my ($mmPos, $targetNuc, $queryNuc) = ($1, $2, $3);
                    push(@mmInfoListDetail, [$mmPos, $targetNuc, $queryNuc]);
                }
                else
                {
                    die "fatal error: cannot parse bowtie mismatch string\n";
                }
            }

            # reconstruct the target seq based on the error string of bowtie output
            my $targetSeq = $readSeq;
            if ($strand eq "+")
            {
                foreach (@mmInfoListDetail)
                {
                    if (substr($readSeq, $_->[0], 1) ne $_->[2])
                    {
                        die "fatal error 23423. bowtie error string does not match query seq\n";
                    }
                    substr($targetSeq, $_->[0], 1) = $_->[1];
                }
            }
            else
            {
                foreach (@mmInfoListDetail)
                {
                    if (substr($readSeq, length($readSeq) - 1 - $_->[0], 1) ne $_->[2])
                    {
                        die "fatal error 23423. bowtie error string does not match query seq\n";
                    }
                    substr($targetSeq, length($readSeq) - 1 - $_->[0], 1) = $_->[1];
                }

                # reverse complement query and target. ral is from query perspective
                $readSeq =~ tr/[ACGTUNRYWSMKBDHVacgtnrywsmkbdhv]/[TGCAANYRWSKMVDHBtgcanyrwskmvdhb]/;
                $readSeq = reverse($readSeq);
                $targetSeq =~ tr/[ACGTUNRYWSMKBDHVacgtnrywsmkbdhv]/[TGCAANYRWSKMVDHBtgcanyrwskmvdhb]/;
                $targetSeq = reverse($targetSeq);
            }

            my $fhl = $outFilesOpened{$chr};

            print $fhl "$seqId\t$chr\t$strand\t$from\t$to\t$errors\t$readSeq\t$targetSeq\n";

            # count hits per sequence
            my $seqErrorId = "$seqId\:$errors";
            if (defined $nrHitsPerSeq{$seqErrorId})
            {
                $nrHitsPerSeq{$seqErrorId}++;
            }
            else { $nrHitsPerSeq{$seqErrorId} = 1; }

        }
        else { die "error: bowtie output is not interpretable\n" }
    }

    foreach (sort keys %outFilesOpened)
    {
        close $outFilesOpened{$_};
    }    # close all file handles

    # test bowtie output: for each quey, the number of errors must be identical
    my %nrTestHash;
    foreach (keys %nrHitsPerSeq) { /^(\S+)\:\d+$/; $nrTestHash{$1} = 1; }
    my $a = keys %nrTestHash;
    my $b = keys %nrHitsPerSeq;
    if ($a != $b) { die "fatal bowtie output interpretation error\n \%nrTestHash has $a keys, \%nrHitsPerSeq nas $b keys\n"; }

    return \%nrHitsPerSeq;
}

# reads the bowtie alignment output and splits the alignments according to chromosomes
# it also returns stats about the number of hits for every sequence
sub splitAlignmentFilesIntoAnnotFilesAndConvertToRal
{
    my ($balFile, $outDir, $seqBins) = @_;

    my %nrHitsPerSeq;
    my %outFilesOpened;    # contains all file handles opened at a given time
    open(F, $balFile) or die "Cannot open $balFile\n";
    while (<F>)
    {
        my @l = split("\t");
        if (/^(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t\S+\t\S+\t(\S*)$/)
        {
            my ($seqId, $strand, $compositeAnnotId, $from, $readSeq, $mmInfo) = ($1, $2, $3, $4, $5, $6);

            my $annotFileId = "";
            my $targetSeqId = "";
            if ($compositeAnnotId =~ /^([^\/]+)\/(\S+)/)
            {
                ($annotFileId, $targetSeqId) = ($1, $2);
            }
            else { die "fatal error 43534\n"; }

            # convert to ral
            my $to         = $from + length($readSeq) - 1;
            my @mmInfoList = split("\,", $mmInfo);
            my $errors     = @mmInfoList;

            if (length($readSeq) < $seqBins->[$errors])
            {
                next;
            }    # check if read is long enough to allow the errors that it has

            if (!(defined $outFilesOpened{$annotFileId}))
            {    # open new filehandle if new annotation file
                open(my $out, ">$outDir/$annotFileId\.ral")
                  or die "Couldt read $outDir\n";
                $outFilesOpened{$annotFileId} = $out;
            }

            # parse the mismatch field
            my @mmInfoListDetail;
            foreach (@mmInfoList)
            {
                if (/^(\d+)\:(\S)\>(\S)$/)
                {
                    my ($mmPos, $targetNuc, $queryNuc) = ($1, $2, $3);
                    push(@mmInfoListDetail, [$mmPos, $targetNuc, $queryNuc]);
                }
                else
                {
                    die "fatal error: cannot parse bowtie mismatch string\n";
                }
            }

            # reconstruct the target seq based on the error string of bowtie output
            my $targetSeq = $readSeq;
            if ($strand eq "+")
            {
                foreach (@mmInfoListDetail)
                {
                    if (substr($readSeq, $_->[0], 1) ne $_->[2])
                    {
                        die "fatal error 23423. bowtie error string does not match query seq\n";
                    }
                    substr($targetSeq, $_->[0], 1) = $_->[1];
                }
            }
            else
            {
                foreach (@mmInfoListDetail)
                {
                    if (substr($readSeq, length($readSeq) - 1 - $_->[0], 1) ne $_->[2])
                    {
                        die "fatal error 23423. bowtie error string does not match query seq\n";
                    }
                    substr($targetSeq, length($readSeq) - 1 - $_->[0], 1) = $_->[1];
                }

                # reverse complement query and target. ral is from query perspective
                $readSeq =~ tr/[ACGTUNRYWSMKBDHVacgtnrywsmkbdhv]/[TGCAANYRWSKMVDHBtgcanyrwskmvdhb]/;
                $readSeq = reverse($readSeq);
                $targetSeq =~ tr/[ACGTUNRYWSMKBDHVacgtnrywsmkbdhv]/[TGCAANYRWSKMVDHBtgcanyrwskmvdhb]/;
                $targetSeq = reverse($targetSeq);
            }

            my $fhl = $outFilesOpened{$annotFileId};

            print $fhl "$seqId\t$targetSeqId\t$strand\t$from\t$to\t$errors\t$readSeq\t$targetSeq\n";

            # count hits per sequence
            my $seqErrorId = "$seqId\:$errors";
            if (defined $nrHitsPerSeq{$seqErrorId}{$compositeAnnotId})
            {
                $nrHitsPerSeq{$seqErrorId}{$compositeAnnotId}++;
            }
            else { $nrHitsPerSeq{$seqErrorId}{$compositeAnnotId} = 1; }

        }
        else { die "error: bowtie output is not interpretable\n" }
    }

    foreach (sort keys %outFilesOpened)
    {
        close $outFilesOpened{$_};
    }    # close all file handles

    my %maxHitsPerTarget;
    for my $seqErrorId (sort keys %nrHitsPerSeq)
    {
        my $maxHits = 0;
        for my $targetId (sort keys %{$nrHitsPerSeq{$seqErrorId}})
        {
            if ($nrHitsPerSeq{$seqErrorId}{$targetId} > $maxHits)
            {
                $maxHits = $nrHitsPerSeq{$seqErrorId}{$targetId};
            }
        }
        $maxHitsPerTarget{$seqErrorId} = $maxHits;
    }

    return \%maxHitsPerTarget;
}

sub loadOneLineTable
{
    my ($filename) = @_;

    my %seqs;

    open(F, $filename) or die "Cannot open $filename\n";
    while (<F>)
    {
        chomp();
        if (/^(\S+)\s+(\S+)\s+(\S+)/)
        {
            my ($id, $seq, $counts) = ($1, $2, $3);

            $seqs{$id} = $seq;
        }
        else
        {
            die "error: input file $filename does not contain one sequence per line with counts.\n";
        }
    }
    return \%seqs;
}

sub loadOneLineFasta
{
    my ($filename) = @_;

    my %seqs;

    open(F, $filename) or die "Cannot open $filename\n";
    while (<F>)
    {
        if (/^\>(\S+)/)
        {
            my $id = $1;
            $_ = <F>;
            chomp($_);
            my $seq = $_;

            $seqs{$id} = $seq;
        }
        else
        {
            die "error: input file $filename is not a one line per sequence fasta file\n";
        }
    }
    return \%seqs;
}

sub loadRAL
{
    my ($filename) = @_;

    my @RAL;
    my @fromCoords;

    open(F, $filename) or die "Cannot open $filename\n";
    while (<F>)
    {
        if (/\S+/)
        {
            my @l = split("\t");

            push(@RAL,        $_);
            push(@fromCoords, $l[3]);
        }
    }
    return (\@RAL, \@fromCoords);
}

sub fastSortRAL
{
	my ($inputFile, $outputFile) = @_;
	my $ret = system("$repDir/soft/pioSortBed6 --ral '$inputFile' > '$outputFile'");
	return $ret;
}

# sorts a ral file by the targetId first and then alignment start coordinate
sub sortRAL
{
    my ($inputFile, $outputFile) = @_;

    open(FI, $inputFile) or die "Cannot open $inputFile\n";
    my %lH;
    my %lC;
    while (<FI>)
    {
        chomp();
        my @l = split("\t");

        if (defined $lH{$l[1]})
        {
            push(@{$lH{$l[1]}}, $_);
            push(@{$lC{$l[1]}}, int($l[3]));
        }
        else
        {
            $lH{$l[1]} = [$_];
            $lC{$l[1]} = [int($l[3])];
        }
    }
    close(FI);

    open(FOUT, ">$outputFile") or die "cannot write to $outputFile\n";

    for my $targetId (sort keys %lH)
    {
        my $lines             = $lH{$targetId};
        my $startCoords       = $lC{$targetId};
        my $startCoordsElemsQ = @{$startCoords} - 1;

        my @list_order =
          sort { $startCoords->[$a] <=> $startCoords->[$b] } 0 .. $startCoordsElemsQ;

        foreach (@list_order)
        {
            print FOUT $lines->[$_], "\n";
        }
    }
    close(FOUT);
}

sub getConfig
{
    my @keys = @_;

    my %conf;
    my $configFile = "$FindBin::RealBin/../../samples/config.tab";
    open(CONF, "<$configFile") or die "ERROR reading $configFile: $!\n";
    while (my $l = <CONF>)
    {
        chomp $l;
        my @f = split(/\t/, $l);
        $conf{$f[0]} = $f[1];
    }
    close CONF;

    if (@keys == 0)
    {
        return \%conf;
    }
    elsif (@keys == 1)
    {
        return $conf{$keys[0]};
    }
    elsif (wantarray)
    {
        return (map { $conf{$_} } @keys);
    }
    else
    {
        die "ERROR: DeepSeqTools::getConfig called with multiple keys in scalar context\n";
    }
}
