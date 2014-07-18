#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use FindBin '$RealBin';
use lib "$RealBin";
use DeepSeqTools;


# exit(0); worked
# exit(1); did not work


# absolute path to the repository
my $repDir="$RealBin/..";
#my $samplesDir="$repDir/samples";
my $seqBaseFileName="seqs";
my $groupDbDir=getConfig('DbDir');
my $genomesDir="$groupDbDir/genomes";
my $seqAnnotDir="$groupDbDir/seqAnnotDB";
my $alnSoftDir="$repDir/soft/aln";
my $auxPrefix="AUX";
my $genomeDbDir="genomes";
my $annotDbDir="seqAnnotDB";
my $famsBaseFileName="familiesAndPriorities.txt";
my $extractDataProg="soft/extractData.pl";
my $noAnnotTag="none";
my $annotTagUnkNoGenomicHit="unknown_unmapped";
my $annotTagUnkTooManyGenomicHits="unknown_overmapped";
my $annotTagUnkMappedNonRepeat="unknown_non-repeat";
my $annotTagUnkMappedRepeat="unknown_repeat";



my $ignoreCounts=0;
my $fragmentReport=0;
my $newPriorityFile="";
my$results=GetOptions('c'  =>  \$ignoreCounts,
		      'f'  =>  \$fragmentReport,
		      'p:s'  =>  \$newPriorityFile
);


if(@ARGV != 3){
    print STDERR "\nUSAGE: ./createFullReport.pl sampleIDs/filename annotationOptions [options]\n\n";
    print STDERR "   -c                 ignore read counts\n\n";
    print STDERR "   -f                 fragment report (default matrix report)\n\n";
    print STDERR "   -p priorityFile    replace default annotation priorities (first col of familiesAndPriorities.txt)\n\n";
    print STDERR "   example: ./createFullReport.pl \"GSM261957,GSM261958\" hg18-hgV01-aln1-AUX_8_3_5 -c\n\n";
    print STDERR "\n";
    exit(0);
}

# 1.11.13 modified by severin. Added samplesDir as argument.
my($sampleIdString,$annotOptions, $samplesDir)=($ARGV[0],$ARGV[1],$ARGV[2]);

my @sampleIds;
if(-f $sampleIdString){ # load ids from file
    open(F,$sampleIdString) or die "Cannot open $sampleIdString\n";
    while(<F>){if(/^(\S+)/){push(@sampleIds,$1);}}
    close(F);
}else{
    @sampleIds=split(",",$sampleIdString);
}

my($genomeVersion,$annotVersion,$alnVersion,$auxString,$auxFileListRef)=parseAnnotationOptions($annotOptions);
my @auxFileList=@{$auxFileListRef};
my $famsAndPriorities="$seqAnnotDir/$annotVersion/$famsBaseFileName";
my $rmskFile="$genomesDir/$genomeVersion/rmsk/allChrs_rmsk.tab";
#if(-f $rmskFile){print STDERR "found repeat annotation\n";}else{print STDERR "no repeat annotation found, continue without\n";}

# read the new priority file if available
my @newPriorityList=();
if($newPriorityFile ne ""){
    open(F,$newPriorityFile) or die "Cannot open $newPriorityFile\n";
    while(<F>){if(/(^\S+)/){push(@newPriorityList,$1);}}
    close(F);
}

if($fragmentReport){ # write header in case of fragment report
    print "#$annotOptions\ttargetId\ttargetStrand\ttargetStart\ttargetEnd\tmm\tseqAlnString\ttargetAlnString\tsampleId\tseqCount\tseqInvWgt\n";
}

# check if for all samples all the mapping files are available
my %auxMappingHash;
for my $sampleId(@sampleIds){
    my $genomeMappingDir="$samplesDir/$sampleId/mappings/$genomeDbDir\_$genomeVersion\-$alnVersion";
    my $annotMappingDir="$samplesDir/$sampleId/mappings/$annotDbDir\_$annotVersion\-$alnVersion";
    my $auxMappingDir="$samplesDir/$sampleId/mappings/$annotDbDir\_$auxPrefix\-$alnVersion";
    
    if(!(-d "$samplesDir/$sampleId")){die "error: sample $sampleId does not exist in $samplesDir\n";}
    if(!(-d "$genomeMappingDir")){die "error: $genomeMappingDir is missing\n";}
    if(!(-d "$annotMappingDir")){die "error: $annotMappingDir is missing\n";}
    
    # read the filenames (only id) of all the auxiliary mapping files that are available (have been mapped)
    opendir(AUXDIR, "$auxMappingDir") or die("Cannot open directory $auxMappingDir");
    my @auxMappingfiles= readdir(AUXDIR);
    foreach (@auxMappingfiles){if ($_=~/^$auxPrefix(\d+)\_/){$auxMappingHash{$sampleId}{$1}=$_;}}
    
    # check if the selected aux mappings really exist
    for my $auxFile(@auxFileList){
	if(!(defined($auxMappingHash{$sampleId}{$auxFile}))){die "error: auxiliary mapping file $auxPrefix\_$auxFile does not exist in $auxMappingDir\n";}
    }
}


# read the hit reports as well as the aux alignments to determine the minimum number of errors for each read

my %annotCountHash;
my @uniqueFamsOrderedOutsideLoop;
my $withRepeatAnnot=0;
for my $sampleId(@sampleIds){
    my $genomeMappingDir="$samplesDir/$sampleId/mappings/$genomeDbDir\_$genomeVersion\-$alnVersion";
    my $annotMappingDir="$samplesDir/$sampleId/mappings/$annotDbDir\_$annotVersion\-$alnVersion";
    my $auxMappingDir="$samplesDir/$sampleId/mappings/$annotDbDir\_$auxPrefix\-$alnVersion";
   
    my ($alnFilesToTraverse,$uniqueFamsOrdered)=collectAnnotationFiles($famsAndPriorities,$annotMappingDir,$auxMappingDir,\@auxFileList,$auxMappingHash{$sampleId},\@newPriorityList);

    my $genomeHitReport=readHitReport("$genomeMappingDir/hitReport.tab");

    my $annotHitReport;
    my $auxHitReport;
    if($newPriorityFile eq ""){ # default priority files
	$annotHitReport=readHitReport("$annotMappingDir/hitReport.tab");
	$auxHitReport=createAuxHitReport($auxMappingHash{$sampleId},\@auxFileList,$auxMappingDir);
    }else{ 
        # in this case maybe some annotation files are completely excluded. get the files that have to be traversed from $alnFilesToTraverse
	# and pass the absolute filepaths provided to createAuxHitReport (3. variable is "").
	my @annotFileListTemp;
	my %annotMappingHashTemp;
	for (my $i=0;$i<@{$alnFilesToTraverse};$i++){
	    push(@annotFileListTemp,$i);
	    $annotMappingHashTemp{$i}=$alnFilesToTraverse->[$i][1];
	}
	$auxHitReport=createAuxHitReport(\%annotMappingHashTemp,\@annotFileListTemp,"");
    }

    my $finalHitReport=intersectHitReports($genomeHitReport,$annotHitReport,$auxHitReport,$noAnnotTag);

    my ($seqs,$seqsReads,$nrUnmappedReads)=readSeqs("$samplesDir/$sampleId/$seqBaseFileName.tab",$finalHitReport);
    
    annotateReads($alnFilesToTraverse,$finalHitReport,$seqs);

    foreach(@{$uniqueFamsOrdered}){$annotCountHash{$_}=0;} # initialize with zero counts

    for my $seqId(keys %{$finalHitReport}){
	my($mismatches,$hits,$annot)=@{$finalHitReport->{$seqId}};
	if($annot eq $noAnnotTag){
	    if($hits==-1){
		$finalHitReport->{$seqId}[2]=$annotTagUnkTooManyGenomicHits;
	    }else{
		$finalHitReport->{$seqId}[2]=$annotTagUnkMappedNonRepeat;
	    }
	}
    }

    # check if there exists a file containing repeat annotation
    
    if(-f "$rmskFile"){
	open(F,"$repDir/$extractDataProg \"$sampleId\,\" \"$annotOptions\" $rmskFile $samplesDir -f|") or die "Fatal extractdata pipe error.\n"; 
	while(<F>){
	    chomp();
	    if(/^\#/){next;}
	    my @l=split("\t");
	    my $seqId=$l[0];
	    if(!defined $finalHitReport->{$seqId}){die "fatal inconsistency error\n";}
	    
	    if($finalHitReport->{$seqId}[2] eq $annotTagUnkMappedNonRepeat){
		$finalHitReport->{$seqId}[2]=$annotTagUnkMappedRepeat;
	    }
	}
	close(F);
	$withRepeatAnnot=1;
    }

    if($fragmentReport){
	for my $seqId(sort keys %{$finalHitReport}){
	    print "$seqId\t$finalHitReport->{$seqId}[2]\t\.\t\.\t\.\t\.\t$seqsReads->{$seqId}\t\.\t$sampleId\t$seqs->{$seqId}\t\.\n";
	}

	# walk through seqs.tab file to print all the sequences that have no mappings
	open(F,"$samplesDir/$sampleId/$seqBaseFileName.tab") or die "Cannot open $samplesDir/$sampleId/$seqBaseFileName.tab\n";
	while(<F>){
	    chomp();
	    my @l=split("\t");
	    if(!(defined $finalHitReport->{$l[0]})){
		print "$l[0]\t$annotTagUnkNoGenomicHit\t\.\t\.\t\.\t\.\t$l[1]\t\.\t$sampleId\t$l[2]\t\.\n";
	    }
	}
    	close(F);

    }

    @uniqueFamsOrderedOutsideLoop=(@{$uniqueFamsOrdered},$annotTagUnkNoGenomicHit,$annotTagUnkTooManyGenomicHits,$annotTagUnkMappedNonRepeat,$annotTagUnkMappedRepeat);

    for my $seqId(keys %{$finalHitReport}){
	my($mismatches,$hits,$annot)=@{$finalHitReport->{$seqId}};

	if($ignoreCounts){
	    $annotCountHash{$sampleId}{$annot}++;
	}else{
	    $annotCountHash{$sampleId}{$annot}+=$seqs->{$seqId};
	}
    }
    
    $annotCountHash{$sampleId}{$annotTagUnkNoGenomicHit}=$nrUnmappedReads;

}

if(!($fragmentReport)){
    print "#$annotOptions\_";
    if(-f $rmskFile){print "rmsk";}else{print "normsk";}

    for my $sampleId(@sampleIds){print "\t$sampleId";}
    print "\n";

    for my $fam(@uniqueFamsOrderedOutsideLoop){
	print "$fam";
	for my $sampleId(@sampleIds){
	    if(defined $annotCountHash{$sampleId}{$fam}){
		print "\t$annotCountHash{$sampleId}{$fam}";
	    }else{
		print "\t0";
	    }
	}
	print "\n";
    }
}

sub parseAnnotationOptions{
    my($annotOptions)=@_;
    
    if($annotOptions=~/^([^\-]+)\-([^\-]+)\-([^\-]+)(\-$auxPrefix\_(\S+)){0,1}$/){
	my($genomeVersion,$annotVersion,$alnVersion,$auxString,$auxFileString)=($1,$2,$3,$4,$5);
	if(!(defined $5)){$auxFileString="";}
	if(!(defined $4)){$auxString="";}
	my @auxFileList=split("\_",$auxFileString);

	return ($genomeVersion,$annotVersion,$alnVersion,$auxString,\@auxFileList);

    }else{die "error: annotation options $annotOptions are in a wrong format\n";}
}

sub readSeqs{
    my($filename,$finalHitReport)=@_;

    my %seqsHash;
    my %seqsReads;
    my $nrUnmappedReads=0;
    open(F,$filename) or die "Cannot open $filename\n";
    while(<F>){
	/^(\S+)\t(\S+)\t(\S+)/;

	# check if seq has mapped at all, otherwise dont include
	if(defined $finalHitReport->{$1}){
	    $seqsHash{$1}=$3;
	    $seqsReads{$1}=$2;
	}else{
	    $nrUnmappedReads+=$3;
	}
    }
    close(F);
   
    # remove all seqs from $finalHitReport that are not in the sample
    # (in case of a virtual sample
    foreach(keys %{$finalHitReport}){
	if(!(defined $seqsHash{$_})){
	    delete($finalHitReport->{$_});
	}
    }


    return (\%seqsHash,\%seqsReads,$nrUnmappedReads);
}

sub readHitReport{
    my($filename)=@_;

    my %hitReportHash;
    open(F,$filename) or die "Cannot open $filename\n";
    while(<F>){
	chomp();
	my @l=split("\t");

	$hitReportHash{$l[0]}=[int($l[1]),int($l[2])];
    }
    close(F);

    return \%hitReportHash;
}

sub createAuxHitReport{
    my($auxMappingHash,$auxFileList,$auxMappingDir)=@_;

    my %minNumberOfMismatches;
    for my $auxNr(@{$auxFileList}){
	my $auxFileName=$auxMappingHash->{$auxNr};
	my $auxFileNameTotal="$auxMappingDir/$auxFileName";

	# read the alignment file for the current AUX
	open(F,$auxFileNameTotal) or die "Cannot open $auxFileNameTotal\n";

	while(<F>){
	    chomp();
	    my @l=split("\t");
	    my($seqId,$errors)=($l[0],$l[5]);


	    if(defined $minNumberOfMismatches{$seqId}){
		if($errors < $minNumberOfMismatches{$seqId}){
		    $minNumberOfMismatches{$seqId}=$errors;
		}
	    }else{
		$minNumberOfMismatches{$seqId}=$errors;
	    }
	}
	close(F);
    }

    my %hitsPerTranscript;

    for my $auxNr(@{$auxFileList}){
	my $auxFileName=$auxMappingHash->{$auxNr};
	my $auxFileNameTotal="$auxMappingDir/$auxFileName";

	# read the alignment file for the current AUX
	open(F,$auxFileNameTotal) or die "Cannot open $auxFileNameTotal\n";
	while(<F>){
	    chomp();
	    my @l=split("\t");
	    my($seqId,$target,$errors)=($l[0],$l[1],$l[5]);

	    # only take into account alignments that have the least numbers of errors
	    if($minNumberOfMismatches{$seqId}==$errors){
		my $tKey="AUX$auxNr\_$target";
		if(defined $hitsPerTranscript{$seqId}{$tKey}){
		    $hitsPerTranscript{$seqId}{$tKey}++;
		}else{
		    $hitsPerTranscript{$seqId}{$tKey}=1;
		}
	    }
	}
	close(F);
    }

    my %hitReport;
    for my $seqId(sort keys %minNumberOfMismatches){

	my $maxNrOfHitsPerSeq=0;
	for my $tKey(sort keys %{$hitsPerTranscript{$seqId}}){
	    if($hitsPerTranscript{$seqId}{$tKey}>$maxNrOfHitsPerSeq){
		$maxNrOfHitsPerSeq=$hitsPerTranscript{$seqId}{$tKey};
	    }
	}
	$hitReport{$seqId}=[$minNumberOfMismatches{$seqId},$maxNrOfHitsPerSeq];
    }

    return \%hitReport;
}


sub intersectHitReports{
    my($genomeHitReport,$annotHitReport,$auxHitReport,$noAnnotTag)=@_;

    my %finalHitReport;
    # copy the $genomeHitReport into %finalHitReport
    for my $seqId(keys %{$genomeHitReport}){
	my($mismatches,$hits)=@{$genomeHitReport->{$seqId}};
	$finalHitReport{$seqId}=[$mismatches,$hits,$noAnnotTag];
    }
    $genomeHitReport = {};

    for my $seqId(keys %{$annotHitReport}){
	my($mismatches,$hits)=@{$annotHitReport->{$seqId}};

	if(defined $finalHitReport{$seqId}){
	    if($finalHitReport{$seqId}->[0]>$mismatches){
		$finalHitReport{$seqId}=[$mismatches,$hits,$noAnnotTag];
	    }
	}else{
	    $finalHitReport{$seqId}=[$mismatches,$hits,$noAnnotTag];
	}
    }
    $annotHitReport = {};

    for my $seqId(keys %{$auxHitReport}){
	my($mismatches,$hits)=@{$auxHitReport->{$seqId}};

	if(defined $finalHitReport{$seqId}){
	    if($finalHitReport{$seqId}->[0]>$mismatches){
		$finalHitReport{$seqId}=[$mismatches,$hits,$noAnnotTag];
	    }
	}else{
	    $finalHitReport{$seqId}=[$mismatches,$hits,$noAnnotTag];
	}
    }

    $auxHitReport = {};

    return \%finalHitReport;
}

sub collectAnnotationFiles{
    my($famsAndPriorities,$annotMappingDir,$auxMappingDir,$auxFileList,$auxMappingHash,$newPriorityList)=@_;

    open(F,$famsAndPriorities) or die "Cannot open $famsAndPriorities\n";
    my @alnFilesToTraverse;
    while(<F>){
	chomp();
	my @l=split("\t");
	my($fam,$class,$fileName)=@l;
	
	# replace suffix
	$fileName=~s/\.fa$/\.ral/g;
	my $filenamePath="$annotMappingDir/$fileName";

	push(@alnFilesToTraverse,[$fam,$filenamePath]);

    }
    close(F);

    # add aux
    for my $auxNr(@{$auxFileList}){
	my $auxFileName=$auxMappingHash->{$auxNr};
	my $filenamePath="$auxMappingDir/$auxFileName";
	my $fam="AUX$auxNr";

	push(@alnFilesToTraverse,[$fam,$filenamePath]);
    }

    # if a new priority list exists, change the order of the entries in @alnFilesToTraverse
    if(@{$newPriorityList}>0){
	my %alnFilesToTraverseHash;
	for my $aftt(@alnFilesToTraverse){
	    my($fam,$filenamePath)=@{$aftt};
	    if(defined $alnFilesToTraverseHash{$fam}){
		push(@{$alnFilesToTraverseHash{$fam}},$filenamePath);
	    }else{
		$alnFilesToTraverseHash{$fam}=[$filenamePath];
	    }
	}
	my @alnFilesToTraverseReordered;
	foreach my $ptFam(@{$newPriorityList}){
	    if(defined $alnFilesToTraverseHash{$ptFam}){
		foreach(@{$alnFilesToTraverseHash{$ptFam}}){
		    push(@alnFilesToTraverseReordered,[$ptFam,$_]);
		}
	    }else{ # list all available annotation terms to correct the input
		foreach(sort keys %alnFilesToTraverseHash){print "$_\n";}
		die "\nThe term $ptFam in the new priority file is not known. see above for available terms:\n";
	    }
 	}
	@alnFilesToTraverse=@alnFilesToTraverseReordered;
    }

    my %uniqueFams;
    my @uniqueFamsOrdered;
    # populate the entries in @alnFilesToTraverse to @uniqueFamsOrdered
    for my $aftt(@alnFilesToTraverse){
	my$fam=$aftt->[0];
	if(!defined $uniqueFams{$fam}){
	    $uniqueFams{$fam}=1;
	    push(@uniqueFamsOrdered,"$fam\+");
	    push(@uniqueFamsOrdered,"$fam\-");
	}
    }

    return (\@alnFilesToTraverse,\@uniqueFamsOrdered);
}

# does not directly return anything, it updates the annotation field of finalHitReport
sub annotateReads{
    my($alnFilesToTraverse,$finalHitReport,$seqs)=@_;

     # go through the aln files in backwards order (respect the priorities principle)
    for (my $i=@{$alnFilesToTraverse}-1;$i>=0;$i--){
	my($fam,$filenamePath)=@{$alnFilesToTraverse->[$i]};

	open(F,$filenamePath) or die "Cannot open $filenamePath\n";
	while(<F>){
	    chomp();
	    my @l=split("\t");
	    my($seqId,$strand,$mismatchesAln)=($l[0],$l[2],$l[5]);
	    if(!defined $seqs->{$seqId}){next;} # do not consider seq in case of virtual sample

	    if(!defined $finalHitReport->{$seqId}){die "fatal error. sequence id is unknown (not allowed to happen)";}
	    
	    my($mismatches,$hits,$annot)=@{$finalHitReport->{$seqId}};

	    if($mismatchesAln == $mismatches){
		#$finalHitReport->{$seqId}->[2]="$fam$strand";

		if($strand eq "+"){ # priority of plus is higher
		    $finalHitReport->{$seqId}->[2]="$fam$strand";
		}else{ # check if the sequence was not annotated with plus
		    if(!($finalHitReport->{$seqId}->[2] eq "$fam\+")){
			$finalHitReport->{$seqId}->[2]="$fam$strand";
		    }

		}

	    }
	}
	close(F);
    }

}
