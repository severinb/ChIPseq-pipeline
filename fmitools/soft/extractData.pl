#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use FindBin '$RealBin';
use lib "$RealBin";
use DeepSeqTools;
#use Time::HiRes qw/gettimeofday tv_interval/;

### TESTERS
# ./extractData.pl GSM261957,GSM261958,GSM261959 mm8-mmV01-aln1 genome | less
### TODO

# extract weighted quantification data for given samples, annotation options and type

#my $t0 = [gettimeofday];
# global variables
my $repDir="$RealBin/..";
#my $samplesDir="$repDir/samples";
my $seqBaseFileName="seqs";
my $DbDir=getConfig('DbDir');
my $genomesDir="$DbDir/genomes";
my $seqAnnotDir="$DbDir/seqAnnotDB";
my $alnSoftDir="$repDir/soft/aln";
my $annotFile="familiesAndPriorities.txt";
my $hitFile="hitReport.tab";
my $seqFile="seqs.tab";
my $auxPrefix="AUX";

my %keepStrand=('p' => '+', 'm' => '-');

my $opt_help = 0;
my $opt_frag = 0;
my $opt_check = 0;
my $opt_over = 1;
my $opt_strand = undef;
my $opt_regionref = "genome";
my $opt_regionabsolute = 0;
my $opt_norm = 0;
my $opt_maxhits = 99;
my $opt_ignorecounts = 0;
my $opt_fraglen = 0;
my $results = GetOptions('h|help' => \$opt_help,
			 'f|frag|fragment' => \$opt_frag,
			 'c|check' => \$opt_check,
			 'o|over|overlap=i' => \$opt_over,
			 'r|ref|regionref=s' => \$opt_regionref,
			 's|strand|regionstrand:s' => \$opt_strand,
			 'a|absolute' => \$opt_regionabsolute,
			 'n|norm|normalize' => \$opt_norm,
			 'm|max|maxhits=i' => \$opt_maxhits,
			 'i|ignorecounts' => \$opt_ignorecounts,
			 "l|len|length=i" => \$opt_fraglen,
			);

#1.11.13, modified by severin: added samplesDir as argument
my ($sampleIdsString, $alnBasis, $type, $samplesDir) = @ARGV;
print STDERR "$sampleIdsString, $alnBasis, $type, $samplesDir";
my @params = (defined $alnBasis ? $alnBasis : "", "maxhits=$opt_maxhits");
if($opt_norm==1) { push @params, "norm"; }
if($opt_ignorecounts==1) { push @params, "ignorecounts"; }
if($opt_fraglen!=0) { push @params, "fraglen=$opt_fraglen"; }
if(defined $opt_strand) {
    if($opt_strand eq '') {
	$opt_strand = 'p';
    } elsif($opt_strand ne 'p' and $opt_strand ne 'm') {
	die "error: unknown strand mode '$opt_strand'; must be one of: p, m\n";
    }
    push @params, "respectStrand=$opt_strand";
}
if(defined $type) { push @params, "type=$type"; }
if(defined $type and -r $type) { push @params, "regionRef=$opt_regionref", "regionMinOverlap=$opt_over"; }
if(defined $type and -r $type and $opt_frag) { push @params, "regionCoords2Target=$opt_regionabsolute"; }
my $params = join(";", @params);

if(@ARGV != 4 or $opt_help){
    print STDERR "\nUSAGE: ./extractData.pl [options] samples genome-annotDB-alignmentAlgorithm[-auxAnnotation] type\n\n";
    print STDERR "    -h      display this help and exit\n";
    print STDERR "    -c      only check availability of alignments and exit\n";
    print STDERR "\n  REPORT TYPE\n";
    print STDERR "    -f      fragment report (default: matrix report)\n";
    print STDERR "\n  GENERAL PARAMETERS (affect all report types)\n";
    print STDERR "    -s p|m  respect strand (default: ignore strand)\n";
    print STDERR "            depending on 'type' and the parameter value for -s, this will result in:\n";
    print STDERR "              for types 'genome', 'int' and 'annoType'\n";
    print STDERR "                ignore alignments that are not on the plus (p) or minus (m) strand\n";
    print STDERR "              for type 'regionFile'\n";
    print STDERR "                ignore overlaps that are not on the same (p) or opposite (m) strand relative to the region strand\n";
    print STDERR "            specifying '-s' without parameter is identical to '-s p'\n";
    print STDERR "    -n      normalize hit counts for target sequence length\n";
    print STDERR "    -m int  maxhits: ignore reads with more than 'int' hits (default: $opt_maxhits)\n";
    print STDERR "    -i      ignore sequence counts (all counts are set to one, default: use observed counts)\n";
    print STDERR "    -l int  fragment length (default: $opt_fraglen, reads will be shifted by len/2 towards their 3'-end)\n";
    print STDERR "            remark: target coordinates in fragment report will not correspond to target sequence when shifting\n";
    print STDERR "\n  REGION OVERLAP PARAMETERS (only apply if 'type' is a regions file)\n";
    print STDERR "    -o int  minimal overlap required to assign sequence to a region (default: $opt_over)\n";
    print STDERR "    -r str  string that specifies the reference for regions (default: $opt_regionref)\n";
    print STDERR "    -a      return absolute rather than relative-to-region coordinates in fragment report\n";
    print STDERR "\n";
    print STDERR "  'genome-annotDB-alignmentAlgorithm[-auxAnnotation]' select the alignment basis\n";
    print STDERR "    e.g. hg18-hgV01-aln1-AUX_3_4 will consider alignments versus 'hg18' genome assembly,\n";
    print STDERR "         'hgV01' annotation database and auxiliary sequences in AUX files '3' and '4'\n";
    print STDERR "         produced by alignment algorithm 'aln1'\n";
    print STDERR "    All components except auxiliary files are mandatory for defining the alignment basis\n";
    print STDERR "\n";
    print STDERR "  'type' defines what to extract and is one of the following:\n";
    print STDERR "     genome     : genmoic mappings\n";
    print STDERR "     int        : an integer, corresponding to auxiliary AUX_int\n";
    print STDERR "     regionFile : a tab-separated file with columns regionId, seqId, strand, start, end,\n";
    print STDERR "                  used to extract data that overlaps these regions (see -r above)\n";
    print STDERR "     annoType   : annotation type, as defined in annotDB/familiesAndPriorities.txt\n";
    print STDERR "\n";
    print STDERR "   examples: ./extractData.pl  GSM261957,GSM261958  mm8-mmV01-aln1-AUX_1_3_8  genome  -f\n";
    print STDERR "             ./extractData.pl  samplesFile  hg18-hgV01-aln1  mRNA_hg_refseq  -n\n";
    print STDERR "             ./extractData.pl  GSM261957    hg18-hgV01-aln1-AUX_3  3 -i\n";
    print STDERR "             ./extractData.pl  samplesFile  hg18-hgV01-aln1  genomicRegionsFile  -r mRNA_hg_refseq -a\n\n";
    exit(0);
}

# parse annotation options
my($genomeVersion,$annotVersion,$alnVersion,$auxString,$auxFileListRef)=parseAnnotationOptions($alnBasis);
my @auxFileList=@{$auxFileListRef};

# get available types
my %knownTypes;
my %validAnnot;
opendir(DIR,$seqAnnotDir) or die "error: could not open $seqAnnotDir: $!\n";
foreach my $annot (grep { /V\d+$/ } readdir DIR) {
    if (-r "$seqAnnotDir/$annot/$annotFile") {
        open(IN, "<$seqAnnotDir/$annot/$annotFile");
        while(<IN>) {
            chomp;
            my ($family, $tp, $fname) = split /\t/;
            $knownTypes{$tp}{$annot} = $fname;
            if($tp eq $type) { $validAnnot{$annot}++; }
        }
        close IN;
    }
}
closedir DIR;

# read the filenames (only id) of all the auxiliary sequences and store them in a hash
my %auxHash;
opendir(AUXDIR, "$seqAnnotDir/$auxPrefix/") or die("Cannot open directory $seqAnnotDir/$auxPrefix/");
foreach (readdir(AUXDIR)) {
    if ($_=~/^$auxPrefix(\d+)\_/) {
	$auxHash{$1} = $_;
    }
}
closedir AUXDIR;
# check if they exist
for my $auxFile(@auxFileList){
    if(!(exists($auxHash{$auxFile}))){die "error: auxiliary file $auxPrefix\_$auxFile does not exist in $seqAnnotDir/$auxPrefix/\n";}
}

# check 'type'
my $nbRegions = 0;
my $regDir;
my %reg;
if($type eq "genome") {
    if(!(-d "$genomesDir/$genomeVersion")) {
	die "error: $genomeVersion does not exist in $genomesDir\n";
    }

} elsif($type =~ m/^\d+$/) {
    if(!exists($auxHash{$type})) {
	die "error: auxiliary file $auxPrefix\_$type does not exist in $seqAnnotDir/$auxPrefix/\n";
    }
    if(scalar(grep {$_ == $type} @auxFileList) != 1) {
	die "error: must include $auxPrefix\_$type in annotation options, currently missing: $alnBasis\n";
    }

} elsif(-r $type) {
    # check region reference
    if($opt_regionref eq 'genome') {
	if(!(-d "$genomesDir/$genomeVersion")) {
	    die "error: $genomeVersion does not exist in $genomesDir\n";
	}
    } elsif($opt_regionref =~ m/^\d+$/) {
	if(!exists($auxHash{$opt_regionref})) {
	    die "error: auxiliary file $auxPrefix\_$opt_regionref does not exist in $seqAnnotDir/$auxPrefix/\n";
	}
    } else {
	if(!exists($knownTypes{$opt_regionref})) {
	    die "error: unkown region reference '$opt_regionref', should be one of:\n\t" .
		join("\n\t", "genome", "(integer)", 
		     grep {exists $knownTypes{$_}{$annotVersion}} sort keys %knownTypes) . "\n";
	}
	if(!exists($knownTypes{$opt_regionref}{$annotVersion})) {
	    die "error: annotation '$annotVersion' does not contain region reference '$opt_regionref', would be in one of:\n\t" .
		join("\n\t",sort keys %{$knownTypes{$opt_regionref}}) . "\n";
	}
    }

    # store regions
    open(IN, "<$type");
    while(<IN>) {
	if (/^#/) {
	    next;
	} else {
	    chomp;
	    my @fields = split /\t/; # id, seqId, strand, start, end
	    if(!exists $reg{$fields[1]}) {
		$reg{$fields[1]} = [[], [], [], []]; # id, strand, start, end
	    }
	    push @{$reg{$fields[1]}[0]}, $fields[0];
	    push @{$reg{$fields[1]}[1]}, $fields[2];
	    push @{$reg{$fields[1]}[2]}, $fields[3];
	    push @{$reg{$fields[1]}[3]}, $fields[4];
	    $nbRegions++;
	}
    }
    close IN;
    
    # sort regions by start
    foreach my $regSeqId (keys %reg) {
	my $starts = $reg{$regSeqId}[2];
	my @inx = sort {$starts->[$a] <=> $starts->[$b]} 0..scalar(@$starts)-1;
	for(my $i=0; $i<4; $i++) {
	    @{$reg{$regSeqId}[$i]} = @{$reg{$regSeqId}[$i]}[@inx];
	}
    }

    # get region directory
    $regDir = $opt_regionref eq 'genome' ? "genomes_$genomeVersion\-$alnVersion" :
              $opt_regionref =~ m/^\d+$/ ? "seqAnnotDB_AUX\-$alnVersion" :
                                           "seqAnnotDB_$annotVersion\-$alnVersion";

} else {
    if(!exists($knownTypes{$type})) {
        die "error: unkown type '$type', should be one of:\n\t" .
          join("\n\t", "genome", "(integer)", "(regionsFile)", grep {exists $knownTypes{$_}{$annotVersion}} sort keys %knownTypes) . "\n";
    }
    if(!exists($knownTypes{$type}{$annotVersion})) {
        die "error: annotation '$annotVersion' does not contain type '$type', would be in one of:\n\t" .
          join("\n\t",sort keys %{$knownTypes{$type}}) . "\n";
    }
}
if(!(-f "$alnSoftDir/$alnVersion.pl")) { die "error: $alnVersion does not exist in $alnSoftDir\n"; }

# check samples and get available alignments
my @sampleIds;
my %availableGenomes; # $availableGenomes{mm8}{aln1}{sampleId} = 1;
my %availableAnnot;   # $availableAnnot{mmV01}{aln1}{sampleId} = 1;
my %availableAUX;     # $availableAUX{{3}{aln1}{sampleId} = 1;
if (-f $sampleIdsString) {
    my %sampleIds;
    my $sampleNo = 0;
    open(IN, "<$sampleIdsString") or die "ERROR: could not read $sampleIdsString: $!\n";
    while(<IN>) {
        chomp;
        if ( /^(\S+)/ ) { $sampleIds{$1} = $sampleNo++; }
    }
    close IN;
    @sampleIds = sort {$sampleIds{$a}<=>$sampleIds{$b}} keys %sampleIds;
    if (scalar(@sampleIds)==0) { die "error: no known sample ids found in string '$sampleIdsString'\n"; }
} else {
    @sampleIds = split(",", $sampleIdsString);
}
foreach my $sampleId (@sampleIds) {
    if(!(-d "$samplesDir/$sampleId")) { die "error: sample '$sampleId' does not exist in $samplesDir\n"; }
    opendir(DIR,"$samplesDir/$sampleId/mappings") or die "error: could not open $samplesDir/$sampleId: $!\n";
    foreach my $alndir (grep { /^(genomes_|seqAnnotDB_)/ } readdir DIR) {
        my ($realm, $dbV, $alnV) = split(/[_-]/, $alndir);
        if ($realm eq 'genomes') {
            $availableGenomes{$dbV}{$alnV}{$sampleId}++;
        } elsif ($dbV eq 'AUX') {
	    opendir(AUXDIR, "$samplesDir/$sampleId/mappings/$alndir") or die("Cannot open directory $samplesDir/$sampleId/mappings/$alndir");
	    foreach (readdir(AUXDIR)) {
		if($_ =~ m/^$auxPrefix(\d+)\_.*\.ral$/) {
		    $availableAUX{$1}{$alnV}{$sampleId}++;
		}
	    }
	    closedir AUXDIR;
	} else {
            $availableAnnot{$dbV}{$alnV}{$sampleId}++;
        }
    }
    closedir DIR;
}

# check if all necessary alignments are available
my %missingAUX;
foreach my $auxFile (@auxFileList) {
    if (scalar(keys %{$availableAUX{$auxFile}{$alnVersion}}) != @sampleIds) {
	$missingAUX{$auxFile}++;
    }
}
if(
   !exists $availableGenomes{$genomeVersion} or
   !exists $availableGenomes{$genomeVersion}{$alnVersion} or
   scalar(keys %{$availableGenomes{$genomeVersion}{$alnVersion}})!=@sampleIds or
   !exists $availableAnnot{$annotVersion} or
   !exists $availableAnnot{$annotVersion}{$alnVersion} or
   scalar(keys %{$availableAnnot{$annotVersion}{$alnVersion}})!=@sampleIds or
   scalar(keys %missingAUX) > 0 or
   $opt_check == 1
  ) {
    if ($opt_check == 1) {
	print STDERR "\njust checking if necessary alignments are available for $alnBasis:\n";
    } else {
	print STDERR "\nnot all necessary alignments are available for $alnBasis:\n";
    }
    if(!exists $availableGenomes{$genomeVersion}) { $availableGenomes{$genomeVersion}{$alnVersion} = {}; }
    if(!exists $availableAnnot{$annotVersion}) { $availableAnnot{$annotVersion}{$alnVersion} = {}; }
    outputAvailableAlignments(\@sampleIds, \%availableGenomes, \%availableAnnot, \%availableAUX,
                              $genomeVersion, $annotVersion, $alnVersion, $auxFileListRef, $type, \%validAnnot);

    if(!exists $availableGenomes{$genomeVersion} or
       !exists $availableGenomes{$genomeVersion}{$alnVersion} or
       scalar(keys %{$availableGenomes{$genomeVersion}{$alnVersion}})!=@sampleIds or
       !exists $availableAnnot{$annotVersion} or
       !exists $availableAnnot{$annotVersion}{$alnVersion} or
       scalar(keys %{$availableAnnot{$annotVersion}{$alnVersion}})!=@sampleIds or
       scalar(keys %missingAUX) > 0) {
        exit(1);
    } else {
        exit(0);
    }
}

#my $t1 = [gettimeofday];
#print STDERR "part: ".tv_interval($t0, $t1)."s\n";
#
# all necessery alignments are available using $genomeVersion-$annotVersion-$alnVersion and @auxFileList
#

# define alignment directory
my $alnDir = ($type eq 'genome' ? "genomes_$genomeVersion\-$alnVersion" :
	      $type =~ m/^\d+$/ ? "seqAnnotDB_AUX\-$alnVersion" :
	      $nbRegions > 0    ? "$regDir" :
	      "seqAnnotDB_$annotVersion\-$alnVersion");

# collect list of alignment files to be parsed,
#   for matrix report: intialize count matrix to zero (%sum)
#   for normalize: get target sequences lengthes (%len) and $avgLen
my @alnFiles;
my %sum; # $sum{targetId} = [sample1Count, ...]
my %len; # $len{targetId} = seqLength_in_nucleotides
my $avgLen = 1;
if($type eq 'genome') {
    my @fastaFiles;
    opendir(DIR, "$genomesDir/$genomeVersion") or
	die "error: could not open $genomesDir/$genomeVersion: $!\n";
    foreach my $alnFile (grep {/\.fa$/} readdir DIR) {
	push @fastaFiles, "$genomesDir/$genomeVersion/$alnFile";
	$alnFile =~ s/\.fa$/\.ral/;
	push @alnFiles, $alnFile;
    }
    closedir DIR;
    if ($opt_norm==1  or  $opt_frag==0) {
	$avgLen = sumInitFromFastalib(\@fastaFiles, \%sum, \%len, scalar(@sampleIds));
    }

} elsif($type =~ m/^\d+$/) {
    my $alnFile = $auxHash{$type};
    if ($opt_norm==1  or  $opt_frag==0) {
	$avgLen = sumInitFromFastalib("$seqAnnotDir/AUX/$alnFile", \%sum, \%len, scalar(@sampleIds));
    }
    $alnFile =~ s/\.fa$/.ral/;
    push @alnFiles, $alnFile;

} elsif($nbRegions>0) {
    if ($opt_norm==1  or  $opt_frag==0) {
	$avgLen = sumInitFromRegions(\%reg, \%sum, \%len, scalar(@sampleIds));
    }
    if($opt_regionref eq 'genome') {
	opendir(DIR, "$genomesDir/$genomeVersion") or
	    die "error: could not open $genomesDir/$genomeVersion: $!\n";
	foreach my $alnFile (grep {/\.fa$/} readdir DIR) {
	    $alnFile =~ s/\.fa$/\.ral/;
	    push @alnFiles, $alnFile;
	}
	closedir DIR;
    } elsif($opt_regionref =~ m/^\d+$/) {
	my $alnFile = $auxHash{$opt_regionref};
	$alnFile =~ s/\.fa$/.ral/;
	push @alnFiles, $alnFile;
    } else {
	my $alnFile = $knownTypes{$opt_regionref}{$annotVersion};
	$alnFile =~ s/\.fa$/\.ral/;
	push @alnFiles, $alnFile;
    }

} else {
    my $alnFile = $knownTypes{$type}{$annotVersion};
    if ($opt_norm==1  or  $opt_frag==0) {
	$avgLen = sumInitFromFastalib("$seqAnnotDir/$annotVersion/$alnFile", \%sum, \%len, scalar(@sampleIds));
    }
    $alnFile =~ s/\.fa$/\.ral/;
    push @alnFiles, $alnFile;
}

# matrix report or fragment report?
my $shift = sprintf("%.0f",$opt_fraglen/2); # int($opt_fraglen/2 + 0.5) would not work for negative fragment lengthes
my $shiftside = 0;
my ($seqId, $targetId, $targetStrand, $targetStart, $targetEnd, $err, $seqAlnStr, $targetAlnStr, $normFact);
if($opt_frag == 1) {  # fragment report
    %sum = ();
    # parse and generate output
    print STDOUT join("\t", "#seqId;$params", "targetId", "targetStrand", "targetStart", "targetEnd",
                      "mm", "seqAlnString", "targetAlnString", "sampleId", "seqCount", "seqInvWgt")."\n";

    for(my $i=0; $i<@sampleIds; $i++) {
	# store inverse weights (number of hits)
	my %iw;
	storeInvWgts(\%iw, \@auxFileList, \%auxHash, $samplesDir, $sampleIds[$i], $alnVersion, $hitFile, $opt_maxhits); #  $iw{seqId} = [iw, err]

	# store sequence counts
	my %seqCounts;
	storeSeqCounts(\%seqCounts, "$samplesDir/$sampleIds[$i]/$seqFile", $opt_ignorecounts); # $seqCounts{seqId} = count

	if($nbRegions>0) {
	    my %alnString;
	    foreach my $alnFile (@alnFiles) {
		if((my $maxAln=countLines("$samplesDir/$sampleIds[$i]/mappings/$alnDir/$alnFile")) > 0 ) {
		    # store alignments of current chromosome
		    DeepSeqTools::_clearC();
		    DeepSeqTools::_initC($maxAln);
		    open(IN, "<$samplesDir/$sampleIds[$i]/mappings/$alnDir/$alnFile")
			or die "error: could not open $samplesDir/$sampleIds[$i]/mappings/$alnDir/$alnFile: $!\n";
		    while(<IN>) {
			chomp;
			($seqId, $targetId, $targetStrand, $targetStart, $targetEnd, $err, $seqAlnStr, $targetAlnStr) = split /\t/;
			if(exists($iw{$seqId}) and $iw{$seqId}[1] == $err and exists($seqCounts{$seqId})) {
			    # ignore if ambiguous hit, not a best hit, or if sequence not in virtual sample
			    $shiftside = ($targetStrand eq "+" ? 1 : ($targetStrand eq "-" ? -1 : 0));
			    DeepSeqTools::_register_regionC($seqId, $targetStrand, $targetStart+$shift*$shiftside, $targetEnd+$shift*$shiftside);
			    $alnString{"$seqId;$targetId;".($targetStart+$shift*$shiftside).";".($targetEnd+$shift*$shiftside).";$targetStrand"} = "$seqAlnStr\t$targetAlnStr";
			}
		    }
		    close IN;
		    DeepSeqTools::_sortOnEndC();

		    # find and output overlaps
		    my $regSeqId = substr($alnFile, 0, -4); # cut off .ral
		    if(exists $reg{$regSeqId}) {
			my ($regId, $regStrand, $regStart, $regEnd) = @{$reg{$regSeqId}};
			for(my $j=0; $j<scalar(@$regId); $j++) {
			    $normFact = $opt_norm==1 ? $len{$regId->[$j]}/$avgLen : 1;
			    my @regInfo = findOverlapsC($regSeqId, $regStart->[$j], $regEnd->[$j], defined($opt_strand) ? ($opt_strand eq 'm' ? oppositeStrand($regStrand->[$j]) : $regStrand->[$j]) : '*', $opt_over);
			    if($opt_regionabsolute) {
				foreach my $regInfo (@regInfo) {
				    print STDOUT join("\t", $regInfo->[0], $regInfo->[1],
						      $regInfo->[4], $regInfo->[2], $regInfo->[3],
						      $iw{$regInfo->[0]}[1], $alnString{join(";",@{$regInfo}[0..4])},
						      $sampleIds[$i], $seqCounts{$regInfo->[0]}, $iw{$regInfo->[0]}[0]*$normFact
					)."\n";
				}
			    } else {
				foreach my $regInfo (@regInfo) {
				    print STDOUT join("\t", $regInfo->[0], $regId->[$j],
						      relativeStrandStartEnd($regStrand->[$j], $regStart->[$j], $regEnd->[$j],
									     $regInfo->[4], $regInfo->[2], $regInfo->[3]),
						      $iw{$regInfo->[0]}[1], $alnString{join(";",@{$regInfo}[0..4])},
						      $sampleIds[$i], $seqCounts{$regInfo->[0]}, $iw{$regInfo->[0]}[0]*$normFact
					)."\n";
				}
			    }
			}
		    }
		}
	    }

	} else {
	    foreach my $alnFile (@alnFiles) {
		if(-e "$samplesDir/$sampleIds[$i]/mappings/$alnDir/$alnFile") {
		    open(IN, "<$samplesDir/$sampleIds[$i]/mappings/$alnDir/$alnFile") or
			die "error: could not open $samplesDir/$sampleIds[$i]/mappings/$alnDir/$alnFile: $!\n";
		    while(<IN>) {
			chomp;
			($seqId, $targetId, $targetStrand, $targetStart, $targetEnd, $err, $seqAlnStr, $targetAlnStr) = split /\t/;
			if(exists($iw{$seqId}) and $iw{$seqId}[1] == $err and exists($seqCounts{$seqId}) and (!defined($opt_strand) or $targetStrand eq $keepStrand{$opt_strand})) {
			    # ignore if ambiguous hit, not a best hit, if sequence not in virtual sample, or if on wrong strand with defined($opt_strand)
			    $shiftside = ($targetStrand eq "+" ? 1 : ($targetStrand eq "-" ? -1 : 0));
			    $normFact = $opt_norm==1 ? $len{$targetId}/$avgLen : 1;
			    print STDOUT join("\t", $seqId, $targetId, $targetStrand, $targetStart+$shift*$shiftside, $targetEnd+$shift*$shiftside,
					      $err, $seqAlnStr, $targetAlnStr, $sampleIds[$i], $seqCounts{$seqId}, $iw{$seqId}[0]*$normFact
				)."\n";
			}
		    }
		    close IN;
		}
	    }
	}
    }

} else { # matrix report
    for(my $i=0; $i<@sampleIds; $i++) {
	# store inverse weights (number of hits)
	my %iw;
	storeInvWgts(\%iw, \@auxFileList, \%auxHash, $samplesDir, $sampleIds[$i], $alnVersion, $hitFile, $opt_maxhits); #  $iw{seqId} = [iw, err]

	# store sequence counts
	my %seqCounts;
	storeSeqCounts(\%seqCounts, "$samplesDir/$sampleIds[$i]/$seqFile", $opt_ignorecounts); # $seqCounts{seqId} = count

	# parse and count
	if($nbRegions>0) {
	    foreach my $alnFile (@alnFiles) {
		if((my $maxAln=countLines("$samplesDir/$sampleIds[$i]/mappings/$alnDir/$alnFile")) > 0 ) {
		    # store alignments of current chromosome
		    DeepSeqTools::_clearC();
		    DeepSeqTools::_initC($maxAln);
		    open(IN, "<$samplesDir/$sampleIds[$i]/mappings/$alnDir/$alnFile")
			or die "error: could not open $samplesDir/$sampleIds[$i]/mappings/$alnDir/$alnFile: $!\n";
		    while(<IN>) {
			chomp;
			($seqId, $targetId, $targetStrand, $targetStart, $targetEnd, $err) = split /\t/;
			if(exists($iw{$seqId}) and $iw{$seqId}[1] == $err and exists($seqCounts{$seqId})) {
			    # ignore if ambiguous hit, not a best hit, or if sequence not in virtual sample
			    $shiftside = ($targetStrand eq "+" ? 1 : ($targetStrand eq "-" ? -1 : 0));
			    DeepSeqTools::_register_regionC($seqId, $targetStrand, $targetStart+$shift*$shiftside, $targetEnd+$shift*$shiftside);
			}
		    }
		    close IN;
		    DeepSeqTools::_sortOnEndC();

		    # find and count overlaps
		    my $regSeqId = substr($alnFile, 0, -4); # cut off .ral
		    if(exists $reg{$regSeqId}) {
			my ($regId, $regStrand, $regStart, $regEnd) = @{$reg{$regSeqId}};
			for(my $j=0; $j<scalar(@$regId); $j++) {
			    my @regInfo = findOverlapsC($regSeqId, $regStart->[$j], $regEnd->[$j], defined($opt_strand) ? ($opt_strand eq 'm' ? oppositeStrand($regStrand->[$j]) : $regStrand->[$j]) : '*', $opt_over);
			    foreach my $regInfo (@regInfo) {
				$sum{$regId->[$j]}[$i] += ( $seqCounts{$regInfo->[0]} / $iw{$regInfo->[0]}[0] );
			    }
			}
		    }
		}
	    }

	} else {
	    foreach my $alnFile (@alnFiles) {
		if(-e "$samplesDir/$sampleIds[$i]/mappings/$alnDir/$alnFile") {
		    open(IN, "<$samplesDir/$sampleIds[$i]/mappings/$alnDir/$alnFile") or
			die "error: could not open $samplesDir/$sampleIds[$i]/mappings/$alnDir/$alnFile: $!\n";
		    while(<IN>) {
			chomp;
			($seqId, $targetId, $targetStrand, $targetStart, $targetEnd, $err) = split /\t/;
			if(exists($iw{$seqId}) and $iw{$seqId}[1] == $err and exists($seqCounts{$seqId}) and (!defined($opt_strand) or $targetStrand eq $keepStrand{$opt_strand})) {
			    # ignore if ambiguous hit, not a best hit, if sequence not in virtual sample, or if on wrong strand with defined($opt_strand)
			    # no need to shift reads here for annot-matrix report (in contrast to region-fragment/annot-fragment/region-matrix reports)
			    $sum{$targetId}[$i] += ( $seqCounts{$seqId} / $iw{$seqId}[0] );
			}
		    }
		    close IN;
		}
	    }
	}
    }

    # generate output
    print STDOUT join("\t", "#targetId;$params", @sampleIds)."\n";
    if ($opt_norm==1) {
	foreach my $targetId (sort keys %sum) {
	    print STDOUT join("\t", $targetId, map {sprintf("%.5f",$sum{$targetId}[$_]/$len{$targetId}*$avgLen)} 0..$#sampleIds)."\n";
	}
    } else {
	foreach my $targetId (sort keys %sum) {
	    print STDOUT join("\t", $targetId, map {sprintf("%.5f",$_)} @{$sum{$targetId}})."\n";
	}
    }
}


exit(0);


sub countLines {
    my ($fname) = @_;
    my $n = 0;
    if(-r $fname) {
	open(FILE, "<$fname") or die "error: could not open $fname: $!\n";
	$n += tr/\n/\n/ while sysread(FILE, $_, 2 ** 20);
	close FILE;
    }
    return $n;
}

sub oppositeStrand {
    return ($_[0] eq '+' ? '-' : ($_[0] eq '-' ? '+' : '*'));
}

sub storeInvWgts {
    # obtain inverse weights as:
    # iw(readId) = max(nb_genomic_hits_at_i,
    #                  max_nb_single_annot_sequence_hits_at_i,
    #                  max_nb_single_aux_sequence_hits_at_i)
    #     where i is the minimal observed nb of alignment errors for readId
    #     (ignore any alignments with of readId with j > i errors)
    # remove $readId if $iw{readId} > $maxhits
    my ($iw, $auxFileListRef, $auxHashRef, $smplDir, $smplId, $alnVer, $hitF, $maxhits) = @_;
    my ($seqId, $err, $nbHits, $targetId);
    # $iw->{seqId} = [iw, err]
    my %ambig; # stores seqIds with nbHits == -1

    # aux (no hit reports, need to parse alignments)
    my %nbHits; # $nbHits{seqId}{targetId} = [iw, mm]
    foreach my $auxId (@$auxFileListRef) {
	my $alnFile = $auxHashRef->{$auxId};
	$alnFile =~ s/\.fa$/.ral/;
	if(-e "$smplDir/$smplId/mappings/seqAnnotDB_AUX\-$alnVer/$alnFile") {
	    open(IN, "<$smplDir/$smplId/mappings/seqAnnotDB_AUX\-$alnVer/$alnFile")
		or die "error: could not open $smplDir/$smplId/mappings/seqAnnotDB_AUX\-$alnVer/$alnFile: $!\n";
	    while(<IN>) {
		chomp;
		($seqId, $targetId, undef, undef, undef, $err) = split /\t/;
		if(!exists $nbHits{$seqId}{$targetId}) {
		    $nbHits{$seqId}{$targetId} = [1, $err];
		} elsif($nbHits{$seqId}{$targetId}[1] > $err) {
		    $nbHits{$seqId}{$targetId}[0] = 1;
		    $nbHits{$seqId}{$targetId}[1] = $err;
		} elsif($nbHits{$seqId}{$targetId}[1]==$err) {
		    $nbHits{$seqId}{$targetId}[0] += 1;
		}
	    }
	    close IN;
	}
    }
    foreach $seqId (keys %nbHits) {
	foreach my $targetId (keys %{$nbHits{$seqId}}) {
	    ($nbHits, $err) = @{$nbHits{$seqId}{$targetId}};
	    if(!exists $iw->{$seqId}) {
		$iw->{$seqId} = $nbHits{$seqId}{$targetId};
	    } elsif($iw->{$seqId}[1] > $err) {
		$iw->{$seqId}[0] = $nbHits;
		$iw->{$seqId}[1] = $err;
	    } elsif($iw->{$seqId}[1]==$err and $iw->{$seqId}[0]<$nbHits) {
		$iw->{$seqId}[0] = $nbHits;
	    }
	}
    }

    # genomic alignments
    open(IN, "<$smplDir/$smplId/mappings/genomes_$genomeVersion\-$alnVer/$hitF")
      or die "error: could not open $smplDir/$smplId/mappings/genomes_$genomeVersion\-$alnVer/$hitF: $!\n";
    while(<IN>) {
	chomp;
	($seqId, $err, $nbHits) = split /\t/;
	if($nbHits == -1) {
	    $ambig{$seqId} = 1;
	    if(exists $iw->{$seqId}) { delete $iw->{$seqId}; }
	} elsif(!exists $iw->{$seqId}) {
	    $iw->{$seqId} = [$nbHits, $err];
	} elsif($iw->{$seqId}[1] > $err) {
	    $iw->{$seqId}[0] = $nbHits;
	    $iw->{$seqId}[1] = $err;
	} elsif($iw->{$seqId}[1]==$err and $iw->{$seqId}[0]<$nbHits) {
	    $iw->{$seqId}[0] = $nbHits;
	}
    }
    close IN;

    # annotation alignments
    open(IN, "<$smplDir/$smplId/mappings/seqAnnotDB_$annotVersion\-$alnVer/$hitF")
      or die "error: could not open $smplDir/$smplId/mappings/seqAnnotDB_$annotVersion\-$alnVer/$hitF: $!\n";
    while(<IN>) {
	chomp;
	($seqId, $err, $nbHits) = split /\t/;
	if(!exists $ambig{$seqId}) {
	    if(!exists $iw->{$seqId}) {
		$iw->{$seqId} = [$nbHits, $err];
	    } elsif($iw->{$seqId}[1] > $err) {
		$iw->{$seqId}[0] = $nbHits;
		$iw->{$seqId}[1] = $err;
	    } elsif($iw->{$seqId}[1]==$err and $iw->{$seqId}[0]<$nbHits) {
		$iw->{$seqId}[0] = $nbHits;
	    }
	}
    }
    close IN;

    # filter by $maxhits
    foreach $seqId (keys %$iw) {
	if($iw->{$seqId}[0] > $maxhits) { delete $iw->{$seqId}; }
    }
}

sub storeSeqCounts {
    my ($seqCounts, $fname, $ignore) = @_;
    my ($seqId, $count) = ("a", 0);
    # $seqCounts->{seqId} = count
    open(IN, "<$fname") or die "error: could not open seqcount file $fname: $!\n";
    if($ignore) {
        while(<IN>) {
            chomp;
            ($seqId) = split /\t/;
            $seqCounts->{$seqId} = 1;
        }
    } else {
        while(<IN>) {
            chomp;
            ($seqId, undef, $count) = split /\t/;
            $seqCounts->{$seqId} = $count;
        }
    }
    close IN;
}

sub relativeStrandStartEnd { # return (strand,start,end) of 2 relative to 1
    my ($strand1, $start1, $end1, $strand2, $start2, $end2) = @_;
    if($strand1 eq '+' or $strand1 eq '*') {
	if($strand2 eq '+' or $strand2 eq '*') {
	    return ('+', $start2-$start1+1, $end2-$start1+1);
	} else {
	    return ('-', $start2-$start1+1, $end2-$start1+1);
	}
    } else {
	if($strand2 eq '+' or $strand2 eq '*') {
	    return ('-', $end1-$end2+1, $end1-$start2+1);
	} else {
	    return ('+', $end1-$end2+1, $end1-$start2+1);
	}
    }
}

sub relativeStrandStartEndReg { # return (strand,start,end) of target relative to region
    my ($rInfo, $tStrand, $tStart, $tEnd) = @_; #$rInfo = [regId, regSeqId, regStart, regEnd, regStrand]
    if($rInfo->[4] eq '+' or $rInfo->[4] eq '*') {
	if($tStrand eq '+' or $tStrand eq '*') {
	    return ('+', $tStart-$rInfo->[2]+1, $tEnd-$rInfo->[2]+1);
	} else {
	    return ('-', $tStart-$rInfo->[2]+1, $tEnd-$rInfo->[2]+1);
	}
    } else {
	if($tStrand eq '+' or $tStrand eq '*') {
	    return ('-', $rInfo->[3]-$tEnd+1, $rInfo->[3]-$tStart+1);
	} else {
	    return ('+', $rInfo->[3]-$tEnd+1, $rInfo->[3]-$tStart+1);
	}
    }
}

sub sumInitFromRegions {
    my ($regs, $sum, $len, $n) = @_;
    my %regIds;
    my $length;
    my $lensum = 0;
    foreach my $regSeqId (keys %$regs) {
	my $regIds    = $regs->{$regSeqId}[0];
	my $regStarts = $regs->{$regSeqId}[2];
	my $regEnds   = $regs->{$regSeqId}[3];
	for(my $i=0; $i<scalar(@$regIds); $i++) {
	    # for multiple identical ids, sum up lenghtes
	    if(defined $len->{$regIds->[$i]}) {
		$length = $regEnds->[$i] - $regStarts->[$i] + 1;
		$lensum += $length;
		$len->{$regIds->[$i]} += $length;
	    } else {
		$len->{$regIds->[$i]} = $regEnds->[$i] - $regStarts->[$i] + 1;
		$lensum += $len->{$regIds->[$i]};
	    }
	}
	foreach my $regId (@$regIds) {
	    $regIds{$regId}++;
	}
    }

    foreach my $regId (keys %regIds) {
	$sum->{$regId} = [(0) x $n];
    }

    return $lensum/scalar(keys %$len);
}

sub sumInitFromFastalib {
    my ($fname, $sum, $len, $n) = @_;
    my ($currId, $currLen, $lensum) = ("", 0, 0);
    my @fnames;
    if(ref($fname) eq "ARRAY") {
	@fnames = @$fname;
    } else {
	@fnames = ($fname);
    }
    foreach $fname (@fnames) {
	open(IN, "<$fname") or die "error: could not open $fname: $!\n";
	while(<IN>) {
	    if(/^>(\S+)/) {
		if ($currId ne "") {
		    $len->{$currId} = $currLen;
		    $lensum += $currLen;
		    $currLen = 0;
		}
		$currId = $1;
		$sum->{$currId} = [(0) x $n];
	    } else {
		$currLen += scalar($_ =~ tr/ACGTUNacgtun/ACGTUNacgtun/);
	    }
	}
	close IN;
	if ($currId ne "") {
	    $len->{$currId} = $currLen;
	    $lensum += $currLen;
	}
    }
    return $lensum/scalar(keys %$len);
}

sub parseAnnotationOptions{
    my($annotOptions)=@_;

    if($annotOptions=~/^([^\-]+)\-([^\-]+)\-([^\-]+)(\-$auxPrefix\_(\S+)){0,1}$/){
        my($genomeVersion,$annotVersion,$alnVersion,$auxString,$auxFileString)=($1,$2,$3,$4,$5);
        if(!(defined $5)){$auxFileString="";}
        if(!(defined $4)){$auxString="";}
        my @auxFileList=split("\_",$auxFileString);

        return ($genomeVersion,$annotVersion,$alnVersion,$auxString,\@auxFileList);

    }else{die "annotation options $annotOptions are in a wrong format\n";}

}

sub outputAvailableAlignments {
    my ($ids, $gens, $annots, $auxs, $genV, $annotV, $alnV, $auxList, $t, $validAnnots) = @_;
    my (@w, $i);
    my $gensMaxN = 0;
    my $annotsMaxN = 0;
    my $auxString = join("_", @$auxList ? "-AUX" : "", @$auxList);
    my ($gensMaxVer, $gensMaxAln, $annotsMaxVer, $annotsMaxAln) = ($genV, $alnV, $annotV, $alnV);
    push @w, (sort {$b<=>$a} map {length($_)} (@$ids, "sampleId"))[0];
    $i = 0;
    print STDERR padR("sampleId",$w[$i++]);
    foreach my $gen (sort keys %$gens) {
        foreach my $aln (sort keys %{$gens->{$gen}}) {
            $w[$i] = length("$gen-$aln") + 1;
            print STDERR padR("$gen-$aln",$w[$i++]);
        }
    }
    foreach my $annot (sort keys %$annots) {
        foreach my $aln (sort keys %{$annots->{$annot}}) {
            $w[$i] = length("$annot-$aln") + 1;
            print STDERR padR("$annot-$aln",$w[$i++]);
        }
    }
    if (@$auxList) {
	foreach my $aux (sort {$a<=>$b} keys %$auxs) {
	    foreach my $aln (sort keys %{$auxs->{$aux}}) {
		$w[$i] = length("AUX$aux-$aln") + 1;
		print STDERR padR("AUX$aux-$aln",$w[$i++]);
	    }
	}
    }
    print STDERR "\n";
    foreach my $sampleId (@$ids) {
        $i = 0;
        print STDERR padR($sampleId,$w[$i++]);
        foreach my $gen (sort keys %$gens) {
            foreach my $aln (sort keys %{$gens->{$gen}}) {
                print STDERR padR(exists $gens->{$gen}{$aln}{$sampleId} ? "1" : "0",$w[$i++]);
            }
        }
        foreach my $annot (sort keys %$annots) {
            foreach my $aln (sort keys %{$annots->{$annot}}) {
                print STDERR padR(exists $annots->{$annot}{$aln}{$sampleId} ? "1" : "0",$w[$i++]);
            }
        }
	if (@$auxList) {
	    foreach my $aux (sort {$a<=>$b} keys %$auxs) {
		foreach my $aln (sort keys %{$auxs->{$aux}}) {
		    print STDERR padR(exists $auxs->{$aux}{$aln}{$sampleId} ? "1" : "0",$w[$i++]);
		}
	    }
	}
	print STDERR "\n";
    }
    $i = 0;
    print STDERR padR("TOTAL",$w[$i++]);
    foreach my $gen (sort keys %$gens) {
        foreach my $aln (sort keys %{$gens->{$gen}}) {
            my $n = scalar(keys %{$gens->{$gen}{$aln}});
            print STDERR padR($n==@$ids ? "*$n" : $n, $w[$i++]);
            if($gensMaxN<$n) {
                $gensMaxN = $n;
                $gensMaxVer = $gen;
                $gensMaxAln = $aln;
            }
        }
    }
    foreach my $annot (sort keys %$annots) {
        foreach my $aln (sort keys %{$annots->{$annot}}) {
            my $n = scalar(keys %{$annots->{$annot}{$aln}});
            print STDERR padR($n==@$ids ? "*$n" : $n, $w[$i++]);
            if($annotsMaxN<$n and exists $validAnnots->{$annot}) {
                $annotsMaxN = $n;
                $annotsMaxVer = $annot;
                $annotsMaxAln = $aln;
            }
        }
    }
    if (@$auxList) {
	foreach my $aux (sort {$a<=>$b} keys %$auxs) {
	    foreach my $aln (sort keys %{$auxs->{$aux}}) {
		my $n = scalar(keys %{$auxs->{$aux}{$aln}});
		print STDERR padR($n==@$ids ? "*$n" : $n, $w[$i++]);
	    }
	}
    }
    print STDERR "\n\n";

    my %haveAllAux;
    foreach my $sampleId (@$ids) {
	if(scalar(grep {exists($auxs->{$_}) and exists($auxs->{$_}{$gensMaxAln}) and exists($auxs->{$_}{$gensMaxAln}{$sampleId})} @$auxList) == scalar(@$auxList)) {
	   $haveAllAux{$sampleId}++;
       }
    }
    if($gensMaxN==@$ids and $annotsMaxN==@$ids and $gensMaxAln eq $annotsMaxAln and scalar(keys %haveAllAux)==@$ids) {
        print STDERR "Samples are comparable using:\n";
        print STDERR "  $0 $sampleIdsString $gensMaxVer\-$annotsMaxVer\-$gensMaxAln$auxString $t\n\n";
    } else {
        print STDERR "Minimal required annotation for $gensMaxVer\-$annotsMaxVer\-$gensMaxAln$auxString:\n";
	my $n = 0;
        foreach my $sampleId (@$ids) {
            if(
	       !exists($gens->{$gensMaxVer}{$gensMaxAln}) or
               !exists($gens->{$gensMaxVer}{$gensMaxAln}{$sampleId}) or
               !exists($annots->{$annotsMaxVer}{$gensMaxAln}) or
               !exists($annots->{$annotsMaxVer}{$gensMaxAln}{$sampleId}) or
	       !exists($haveAllAux{$sampleId})
	       #or scalar(grep {exists($auxs->{$_}) and exists($auxs->{$_}{$gensMaxAln}) and exists($auxs->{$_}{$gensMaxAln}{$sampleId})} @$auxList) != scalar(@$auxList)
	      ) {
                print STDERR "  submitAnnotationTask.pl $sampleId $gensMaxVer\-$annotsMaxVer\-$gensMaxAln$auxString\n";
		$n++;
            }
        }
	if ($n==0) { print STDERR "  none\n"; }
        print STDERR "will allow comparison by:\n";
        print STDERR "  $0 ".join(",",@$ids)." $gensMaxVer\-$annotsMaxVer\-$gensMaxAln$auxString $t\n";
        print STDERR "\n";
    }
    if($gensMaxVer ne $genV or $gensMaxAln ne $alnV or $annotsMaxVer ne $annotV) {
        print STDERR "Required annotation for requested $genV\-$annotV\-$alnV$auxString:\n";
        foreach my $sampleId (@$ids) {
            if(!exists($gens->{$genV}) or
               !exists($gens->{$genV}{$alnV}) or
               !exists($gens->{$genV}{$alnV}{$sampleId}) or
               !exists($annots->{$annotV}) or
               !exists($annots->{$annotV}{$alnV}) or
               !exists($annots->{$annotV}{$alnV}{$sampleId})) {
                print STDERR "  submitAnnotationTask.pl $sampleId $genV\-$annotV\-$alnV$auxString\n";
            }
        }
        print STDERR "will allow comparison by:\n";
        print STDERR "  $0 ".join(",",@$ids)." $genV\-$annotV\-$alnV$auxString $t\n";
        print STDERR "\n";
    }
    print STDERR "\n";
}

sub padR { return sprintf('%*2$s', $_[0], $_[1] || 20); }

