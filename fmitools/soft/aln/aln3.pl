#!/usr/bin/env perl

use strict;
use warnings;
use FindBin '$Bin';
use Getopt::Long;

my $mthreads=getConfig('Threads'); # number of threads to use by oligomap and soap
my $m=100; # maximum number of genomic hits per sequence
my $bowtieOptions="-f -p $mthreads --best --strata -v 3 -B 1 --norc --quiet"; # allow up to 3 mismatches, ignore quality scores. coord start 1
my $seqBins=[24,28,33,38]; # seqence lengths for 0,1,2 and 3 mismatches


# THIS ALN VERSION CONVERTS ALL C'S TO T'S IN THE GENOME AND THE READS BEFORE ALIGNING

my $dontPrintReport=0;
my$results=GetOptions('r'  => \$dontPrintReport
);

if(@ARGV != 5){
    print STDERR "\nUSAGE: ./aln3.pl targetDir queryFile typeOfTarget tempDir outDir\n\n";
    print STDERR "   tempDir:  g: genome, a: annotation\n";
    print STDERR "   -r     :  dont print report\n\n";

    exit(0);
}

my $bowtieIndexSubdirName="bowtieIndexCtoT";

my $repDir="$Bin/../..";
my $bowtieDir="$repDir/soft/bowtie-0.10.0.1";


my($targetDir,$queryFile,$typeOfTarget,$tempDir,$outDir)=($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4]);

# convert the reads and save in a fasta file
copyReadsIntoFastaFile($queryFile,"$tempDir/seqs_O.fa",0);
copyReadsIntoFastaFile($queryFile,"$tempDir/seqs_R.fa",1);


# map to genome
if($typeOfTarget eq "g"){

    # check if bowtie index exists. if not create one
    if(system("mkdir $targetDir/$bowtieIndexSubdirName 2> /dev/null")==0){
	
	# concatenate both strand C to T converted chromosome files (.fa) into one big temporary file
	convertTtoCforDir($targetDir,"$targetDir/$bowtieIndexSubdirName/catChrs_O.fa",0);
	convertTtoCforDir($targetDir,"$targetDir/$bowtieIndexSubdirName/catChrs_C.fa",1);

	# build bowtie index for plus strand and complement (not rc)
	my $ret1=system("$bowtieDir/bowtie-build $targetDir/$bowtieIndexSubdirName/catChrs_O.fa $targetDir/$bowtieIndexSubdirName/$bowtieIndexSubdirName\_O > /dev/null");
	system("rm $targetDir/$bowtieIndexSubdirName/catChrs_O.fa");
	my $ret2=system("$bowtieDir/bowtie-build $targetDir/$bowtieIndexSubdirName/catChrs_C.fa $targetDir/$bowtieIndexSubdirName/$bowtieIndexSubdirName\_C > /dev/null");
	system("rm $targetDir/$bowtieIndexSubdirName/catChrs_C.fa");

	if(($ret1 != 0) or ($ret2 != 0)){
	    system("rm \-rf $targetDir/$bowtieIndexSubdirName"); # delete index directoy if bowtie failed
	    die "bowtie error when building the index\n";
	}
    }

    my $mp=$m+1;

    # map the reads to the converted plus strand genome and the complement (keep also overmapper hits with a limit of $m+1
    my $ret2=system("$bowtieDir/bowtie $bowtieOptions -k $mp $targetDir/$bowtieIndexSubdirName/$bowtieIndexSubdirName\_O $tempDir/seqs_O.fa > $tempDir/seqs_O.bal");
    if($ret2 != 0){die "bowtie error aligning the reads\n";}

    my $ret3=system("$bowtieDir/bowtie $bowtieOptions -k $mp $targetDir/$bowtieIndexSubdirName/$bowtieIndexSubdirName\_C $tempDir/seqs_R.fa > $tempDir/seqs_R.bal");
    if($ret3 != 0){die "bowtie error aligning the reads\n";}

 
    my $nrHitsPerSeq=splitAlignmentFilesAndConvertToRal("$tempDir/seqs",$tempDir,$seqBins,$m,$typeOfTarget,$targetDir);

    # sort alignments for each chromosome alignment file and write to final destination folder
    opendir(RD, $tempDir) or die("Cannot open directory $tempDir");
    my @unsortedRalDir= readdir(RD);
    foreach my $fileName (@unsortedRalDir){
	if ($fileName=~/\.ral$/){
	    sortRAL("$tempDir/$fileName","$outDir/$fileName");
	}
    }

    open(F,">$outDir/hitReport.tab") or die "Cannot open $outDir/hitReport.tab\n";
    for my $seqErrorId (sort {$nrHitsPerSeq->{$a} <=> $nrHitsPerSeq->{$b}} keys %{$nrHitsPerSeq}){
	$seqErrorId=~/^(\S+)\:(\d+)/;
	my($seqId,$errors)=($1,$2);
	my $hits=$nrHitsPerSeq->{$seqErrorId};
	
	if($hits<$m+1){
	    print F "$seqId\t$errors\t$hits\n";
	}else{
	    print F "$seqId\t0\t-1\n";
	}
    }
    close(F);

   

# map to annotation
}elsif($typeOfTarget eq "a"){

    # read the annotation files and cat them with a modified header (+filename, allows to split later)
    opendir(ANDIR, "$targetDir") or die("Cannot open directory $targetDir\n");
    my @annDbfiles= readdir(ANDIR);
    my @allAnnotFileNames;
    foreach (@annDbfiles){if ($_=~/(\S+)\.fa$/){my($fileStem)=($1);push(@allAnnotFileNames,$fileStem);}}

    # check if there are annotation files at all
    if(@allAnnotFileNames==0){
	#create empty hitreport if needed
	if(!($dontPrintReport)){system("touch $outDir/hitReport.tab");}
    }else{

	# check if bowtie index exists. if not create one
	if(system("mkdir $targetDir/$bowtieIndexSubdirName 2> /dev/null")==0){
	    
	    # concatenate both strand C to T converted chromosome files (.fa) into one big temporary file
	    convertTtoCforDir($targetDir,"$targetDir/$bowtieIndexSubdirName/catChrs_O.fa",0);
	    convertTtoCforDir($targetDir,"$targetDir/$bowtieIndexSubdirName/catChrs_C.fa",1);
	    
	    # build bowtie index for plus strand and complement (not rc)
	    my $ret1=system("$bowtieDir/bowtie-build $targetDir/$bowtieIndexSubdirName/catChrs_O.fa $targetDir/$bowtieIndexSubdirName/$bowtieIndexSubdirName\_O > /dev/null");
	    system("rm $targetDir/$bowtieIndexSubdirName/catChrs_O.fa");
	    my $ret2=system("$bowtieDir/bowtie-build $targetDir/$bowtieIndexSubdirName/catChrs_C.fa $targetDir/$bowtieIndexSubdirName/$bowtieIndexSubdirName\_C > /dev/null");
	    system("rm $targetDir/$bowtieIndexSubdirName/catChrs_C.fa");
	    
	    if(($ret1 != 0) or ($ret2 != 0)){
		system("rm \-rf $targetDir/$bowtieIndexSubdirName"); # delete index directoy if bowtie failed
		die "bowtie error when building the index\n";
	    }
	}
	
	# map the reads to the converted plus strand genome and the complement
	my $ret2=system("$bowtieDir/bowtie $bowtieOptions -a $targetDir/$bowtieIndexSubdirName/$bowtieIndexSubdirName\_O $tempDir/seqs_O.fa > $tempDir/seqs_O.bal");
	if($ret2 != 0){die "bowtie error aligning the reads\n";}
	
	my $ret3=system("$bowtieDir/bowtie $bowtieOptions -a $targetDir/$bowtieIndexSubdirName/$bowtieIndexSubdirName\_C $tempDir/seqs_R.fa > $tempDir/seqs_R.bal");
	if($ret3 != 0){die "bowtie error aligning the reads\n";}
	
	foreach (@allAnnotFileNames){system("touch $tempDir/$_\.ral");} # create empty files
	
	my $maxHitsPerTarget=splitAlignmentFilesAndConvertToRal("$tempDir/seqs",$tempDir,$seqBins,$m,$typeOfTarget,$targetDir);
	
	# sort alignments for each chromosome alignment file and write to final destination folder
	opendir(RD, $tempDir) or die("Cannot open directory $tempDir");
	my @unsortedRalDir= readdir(RD);
	foreach my $fileName (@unsortedRalDir){
	    if ($fileName=~/\.ral$/){
		sortRAL("$tempDir/$fileName","$outDir/$fileName");
	    }
	}
	
	#create hitreport if needed
	if(!($dontPrintReport)){
	    open(F,">$outDir/hitReport.tab") or die "Cannot open $outDir/hitReport.tab\n";
	    for my $seqErrorId (sort keys %{$maxHitsPerTarget}){
		$seqErrorId=~/^(\S+)\:(\d+)/;
		my($seqId,$errors)=($1,$2);
		print F "$seqId\t$errors\t$maxHitsPerTarget->{$seqErrorId}\n";
	    }
	    close(F);
	}
    }

}else{die "error: unknown typeOfTarget parameter $typeOfTarget\n";}




# reads the seqs.tab file, converts all C's to T's and copies the sequences into a new file in fasta format
sub copyReadsIntoFastaFile{
    my ($seqsTabFile,$outFastaFile,$reverse)=@_;
 
    open(FO,">$outFastaFile") or die "Cannot open $outFastaFile\n";
    open(FI,$seqsTabFile) or die "Cannot open $seqsTabFile \n";
      while(<FI>){
	chomp();
	if(/^(\S+)\s+(\S+)\s+(\S+)/){
	    my($id,$seq,$counts)=($1,$2,$3);
	    my $seqConv=$seq;
	    $seqConv=~tr/cC/tT/;

	    if($reverse==1){$seqConv=reverse($seqConv);}

	    print FO ">$seq\:$id\n$seqConv\n";
	}else{
	    die "error: input file $seqsTabFile does not contain one sequence per line with counts.\n";
	}
    }
    close(FI);
    close(FO);
}

# reads the bowtie alignment output and splits the alignments according to chromosomes
# it also returns stats about the number of hits for every sequence
sub splitAlignmentFilesAndConvertToRal{
    my ($balFile,$outDir,$seqBins,$m,$typeOfTarget,$targetDir)=@_;
 
    my %nrHitsPerSeqRaw;
    my %outFilesOpened; # contains all file handles opened at a given time
    

    for my $cstrand(("O","R")){
	open(F,"$balFile\_$cstrand.bal") or die "Cannot open $balFile\n";

	while(<F>){
	    my @l=split("\t");
	    if(/^(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t\S+\t\S+\t(\S*)$/){
		my($seqIdComp,$strand,$compositeAnnotId,$from,$readSeqCtoT,$mmInfo)=($1,$2,$3,$4,$5,$6);

		my $annotFileId="";
		my $targetSeqId="";
		if($compositeAnnotId=~/^([^\/]+)\/(\S+)/){($annotFileId,$targetSeqId)=($1,$2);}else{die "fatal error 43534\n";}

		my $readSeq="";
		my $seqId="";
		if($seqIdComp=~/^([^\:]+)\:(\S+)$/){($readSeq,$seqId)=($1,$2);}else{die "fatal error 43535\n";}

		# convert to ral and print
		my $to=$from+length($readSeq)-1;
		my @mmInfoList=split("\,",$mmInfo);
		my $errors=@mmInfoList;
		
		if(length($readSeq)<$seqBins->[$errors]){next;} # check if read is long enough to allow the errors that it has
		
		if(!(defined $outFilesOpened{$annotFileId})){ # open new filehandle if new chromosome
		    open(my $out, ">$outDir/$annotFileId\.pral" )   or die "Couldt read $outDir\n";
		    $outFilesOpened{$annotFileId}=$out;
		}
		
		# set the correct strand based on plus or minus strand of converted genome
		my $rStrand="+";
		if($cstrand eq "R"){$rStrand="-";}
		
		my $fhl=$outFilesOpened{$annotFileId};
		
		# replace C's in read by lowercase t
		$readSeq=~tr/C/t/;

		print $fhl "$seqId\t$targetSeqId\t$rStrand\t$from\t$to\t$errors\t$readSeq\n";
		
		if($typeOfTarget eq "g"){$compositeAnnotId="g";} # only one "genome" sequence

		# count hits per sequence
		if(defined $nrHitsPerSeqRaw{$seqId}{$errors}{$compositeAnnotId}){$nrHitsPerSeqRaw{$seqId}{$errors}{$compositeAnnotId}++;}else{$nrHitsPerSeqRaw{$seqId}{$errors}{$compositeAnnotId}=1;}

	    }else{die "error: bowtie output is not interpretable\n"};
	}
    }
    foreach(sort keys %outFilesOpened){close $outFilesOpened{$_};} # close all file handles

    # deterime minimum number of errors for each sequence (necessary due to genome strand separations).
    # check if there are not more than $m hits (in genomic alignment case)
    my %nrHitsPerSeq;
    for my $seqId(keys %nrHitsPerSeqRaw){
	my @errorList=sort keys %{$nrHitsPerSeqRaw{$seqId}};
	my $minErrors=$errorList[0];

	my $maxHits=0;
	my $maxHitsTaget="";
	for my $compId (keys %{$nrHitsPerSeqRaw{$seqId}{$minErrors}}){
	    if($nrHitsPerSeqRaw{$seqId}{$minErrors}{$compId} > $maxHits){
		$maxHits=$nrHitsPerSeqRaw{$seqId}{$minErrors}{$compId};
		$maxHitsTaget=$compId;
	    }
	}
	my $seqErrorId="$seqId\:$minErrors";
	
	$nrHitsPerSeq{$seqErrorId}=$maxHits;
    }

    # walk through all the files and extract target sequence from unconverted genome
    opendir(DIR, $outDir) or die("Cannot open directory $outDir\n");
    my @filesInDir= readdir(DIR);
    my @pralFileNames;
    foreach (@filesInDir){if ($_=~/(\S+)\.pral$/){push(@pralFileNames,$1);}}

    for my $fileStem(@pralFileNames){
	# read fasta file with the original sequences
	my $tagetSeqs=loadFastaFileIntoHash("$targetDir/$fileStem\.fa");

	open(FI,"$outDir/$fileStem\.pral") or die "Cannot open $outDir/$fileStem\.pral\n";
	open(FO,">$outDir/$fileStem\.ral") or die "Cannot open $outDir/fileStem\.ral\n";
	while(<FI>){
	    chomp();
	    my @l=split("\t");
	    my $subSeq=substr($tagetSeqs->{$l[1]},$l[3]-1,$l[4]-$l[3]+1);
	    if($l[2] eq "-"){ # reverse complement
		$subSeq =~ tr/[ACGTUNRYWSMKBDHVacgtnrywsmkbdhv]/[TGCAANYRWSKMVDHBtgcanyrwskmvdhb]/; 
		$subSeq = reverse($subSeq);
	    }
	    # replace C's in target by lowercase t
	    $subSeq=~tr/C/t/;

	    my $seqErrorId="$l[0]\:$l[5]";
	    if(defined $nrHitsPerSeq{$seqErrorId}){
		# print always except in case of genomic with more than $m hits
		if(!(($typeOfTarget eq "g") and ($nrHitsPerSeq{$seqErrorId} > $m))){
		    print FO "$_\t$subSeq\n";
		}
	    }
	}
	close(FI);
	close(FO);
    }

    return \%nrHitsPerSeq;
}



sub loadOneLineTable{
    my($filename)=@_;

    my %seqs;
 
    open(F,$filename) or die "Cannot open $filename\n";
    while(<F>){
	chomp();
	if(/^(\S+)\s+(\S+)\s+(\S+)/){
	    my($id,$seq,$counts)=($1,$2,$3);

	    $seqs{$id}=$seq;
	}else{
	    die "error: input file $filename does not contain one sequence per line with counts.\n";
	}
    }
    return \%seqs;
}


sub loadOneLineFasta{
    my($filename)=@_;

    my %seqs;

    open(F,$filename) or die "Cannot open $filename\n";
    while(<F>){
	if(/^\>(\S+)/){
	    my $id=$1;
	    $_=<F>;
	    chomp($_);
	    my $seq=$_;
	    
	    $seqs{$id}=$seq;
	}else{
	    die "error: input file $filename is not a one line per sequence fasta file\n";
	}
    }
    return \%seqs;
}


sub loadRAL{
    my($filename)=@_;

    my @RAL;
    my @fromCoords;

    open(F,$filename) or die "Cannot open $filename\n";
    while(<F>){
	if(/\S+/){
	    my@l=split("\t");

	    push(@RAL,$_);
	    push(@fromCoords,$l[3]);
	}
    }
    return (\@RAL,\@fromCoords);
}

# sorts a ral file by the targetId first and then alignment start coordinate
sub sortRAL{
    my($inputFile,$outputFile)=@_;

    open(FI,$inputFile) or die "Cannot open $inputFile\n";
    my %lH;
    my %lC;
    while(<FI>){
	chomp();
	my @l=split("\t");

	if(defined $lH{$l[1]}){
	    push(@{$lH{$l[1]}},$_);
	    push(@{$lC{$l[1]}},int($l[3]));
	}else{
	    $lH{$l[1]}=[$_];
	    $lC{$l[1]}=[int($l[3])];
	}
    }
    close(FI);

    open(FOUT,">$outputFile") or die "cannot write to $outputFile\n";

    for my $targetId(sort keys %lH){
	my $lines=$lH{$targetId};
	my $startCoords=$lC{$targetId};
	my $startCoordsElemsQ=@{$startCoords}-1;

	my @list_order =  sort { $startCoords->[$a] <=> $startCoords->[$b] } 0 .. $startCoordsElemsQ;

	foreach(@list_order){
	    print FOUT $lines->[$_],"\n";
	}
    }
    close(FOUT);
}



sub getConfig {
    my @keys = @_;

    my %conf;
    my $configFile = "$FindBin::RealBin/../../samples/config.tab";
    open(CONF, "<$configFile") or die "ERROR reading $configFile: $!\n";
    while(my $l = <CONF>) {
	chomp $l;
	my @f = split(/\t/, $l);
	$conf{$f[0]} = $f[1];
    }
    close CONF;

    if(@keys == 0) {
	return \%conf;
    } elsif(@keys == 1) {
	return $conf{$keys[0]};
    } elsif(wantarray) {
	return (map {$conf{$_}} @keys);
    } else {
	die "ERROR: DeepSeqTools::getConfig called with multiple keys in scalar context\n";
    }
}


# read all fasta files in directory and create a concatenated
# file with all of them in a T to C converted form performing
# the conversion either on plus strand or on complement (not reverse complement)
sub convertTtoCforDir{
    my($indir,$outFile,$complement)=@_;

    # find all .fa files in the input directory
    opendir(INDIR, $indir) or die("Cannot open directory $indir\n");
    my @filesInDir= readdir(INDIR);
    my @fastaFileNames;
    foreach (@filesInDir){if ($_=~/(\S+)\.fa$/){push(@fastaFileNames,$1);}}


    open(FOUT,">$outFile") or die "cannot open $outFile\n";
    for my $fileStem(@fastaFileNames){
	
	open(FIN,"$targetDir/$fileStem.fa") or die "Cannot open $targetDir/$fileStem.fa\n";
	while(<FIN>){
	    if(/^\s+$/){next;}
	    if(/^\>(\S+)/){
		print FOUT "\>$fileStem\/$1\n";
	    }else{
		# convert U to T. bowtie will not find hits to target seqs with U
		$_=~tr/Uu/Tt/;
		
		# make the complement (not reverse complement, does not require coordinates from end)
		if($complement==1){
		   $_ =~ tr/[ACGTUNRYWSMKBDHVacgtnrywsmkbdhv]/[TGCAANYRWSKMVDHBtgcanyrwskmvdhb]/; 
		}

		# convert C's to T's
		$_ =~tr/cC/tT/;

		print FOUT "$_";
	    }
	}
	close(FIN);

    }
}



sub loadFastaFileIntoHash{
    my ($filename)=@_;
    open(FHANDLE,$filename) or die "Cannot open $filename \n";
  
    my $id="";
    my $tempSeq="";
    my %seqs;

    while(<FHANDLE>){
	chomp();
	if((/^\>(\S+)/)){
	    if($id ne ""){
		if(!(defined $seqs{$id})){
		    $seqs{$id}=$tempSeq;
		}else{die "error loading from $filename containing non uniqe ids: $id\n";}
	    }

	    $id=$1;
	    $tempSeq="";
	}else{
	    tr/[a-z]/[A-Z]/;
	    tr/U/T/;
	    tr/ACGTN/N/c; # replace all non ACGTNs by Ns

	    $tempSeq.=$_;
	}

    }
    close(FHANDLE);

    if($id ne ""){
	if(!(defined $seqs{$id})){
	    $seqs{$id}=$tempSeq;
	}else{die "error loading from $filename containing non uniqe ids: $id\n";}
    }

    return \%seqs;
}
