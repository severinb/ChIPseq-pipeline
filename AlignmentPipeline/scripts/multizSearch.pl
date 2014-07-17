#!/import/bc2/soft/bin/perl5/perl -w

use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';

#use lib '/import/bc2/home/nimwegen/GROUP/AlignmentPipeline/scripts';
use lib '/import/wnz/home/crunch/AlignmentPipeline/scripts';
use DTable::DTableOps;

#assumes that we have the .maf files
#when we start from .axt formats we need to convert them with
#/import/bc2/home8/zavolan/GROUP/Perl/PerlScripts/dPerlScripts/axt_to_maf.pl indir anchorOrgChrLengths secondOrgChrLengths
#and to get the chromosome lengths use
#/import/bc2/home8/zavolan/GROUP/Perl/PerlScripts/dPerlScripts/chrLengthsInDir.pl chrdir >lenfile

my $outFormat='extMZ'; # or maf
my$results=GetOptions('f:s'   => \$outFormat);


if(@ARGV != 2) {
    die "usage: multizSearch multizDir regions.genCoords > regions.mz\n";
}


# read the table with the coordinates

my $lines=DTableOps::loadFileRef($ARGV[1],'');
my $tab=DTableOps::convertLinesToTableRef($lines,'\s+');

# split the entries according to the chromosomes

# create a set of all the used chromosomes
my %occChrs;
for my $i(0..$#{$tab}){
    $occChrs{$tab->[$i][0]}=1;
}


foreach (keys %occChrs){
    my @chrTab=DTableOps::selectRowsWhere($tab,0,$_);
    my ($chrTabSorted,$chrTabSortedInd)=DTableOps::sortTableForColumn(\@chrTab,1);

    #DTableOps::printTable($chrTabSorted);
    findRegionsInMZA($ARGV[0]."/$_.maf",$chrTabSorted,$outFormat);
}
#DTableOps::printTable($tab);


sub findRegionsInMZA{
    my($filename,$tab,$outFormat)=@_;

    my $fin;
    my $flag;
    my $ind=0;
    if (!open($fin, '<', $filename)) {
        print STDERR "ERROR!!! No file $filename found. Skipping $tab.\n";
        return 1;
    } 
    while(<$fin>) {
	# check if an msa starts
	if ( ($_=~/^a score\=(\S+)/) and ($#{$tab}>=$ind)){
	    my $score=$1;

	    # read the first line and check if coordinates match
	    $_ = <$fin>;
	    $_ =~ /^s\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\S+)/;
	    my ($name, $beg, $len, $mzOr, $chrLength, $seq) = ($1, $2, $3, $4, $5, $6);

	    # remove all hsps that are "behind" the actual line
	    while ( ( $#{$tab}>=$ind) and ($tab->[$ind][2]<$beg) ){
		$ind++;
	    }
	    # if all the entries are checked, then return
	    if ($ind>$#{$tab}){return;}

	    my $sitebeg=$tab->[$ind][1];
	    my $siteend=$tab->[$ind][2];
	    my $chr=$tab->[$ind][0];
	    my $or=$tab->[$ind][3];
	    my $hsp=$tab->[$ind][4];

	    
	    my $end = $beg+$len-1;

	    if(($sitebeg <= $end) and ( $beg <= $siteend)) {
		#overlap
		my @lines;
		while($_ =~ /\S/) {
		    push @lines, $_;
		    $_ = <$fin>;
		}
 
		# there can be more than one hsp in the alignment.
		my $localInd=$ind;
		while(($sitebeg <= $end) and ($beg <= $siteend) and ( $#{$tab}>=$localInd)){
		    extractMsaRegion($hsp,$chr,$or,$sitebeg,$siteend,\@lines,$score,$outFormat);
		    $localInd++;
		    if ($localInd<=$#{$tab}){
			$sitebeg=$tab->[$localInd][1];
			$siteend=$tab->[$localInd][2];
			$or=$tab->[$localInd][3];
			$chr=$tab->[$localInd][0];
			$hsp=$tab->[$localInd][4];
		    }
		}
	    }
	}
	
    }
    close($fin);

}

sub extractMsaRegion{
    my($hsp,$chr,$or,$sitebeg,$siteend,$linesRef,$score,$outFormat)=@_;
    my @lines=@{$linesRef};

    my (@names, @begins, @ends, @ors, @chrLengths, @seqs);
    my $maxl = 0;
    for(my $i = 0; $i < @lines; $i++) {
	if($lines[$i] =~ /^s\s+(\S+\s+)(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\S+)/) {
	    ($names[$i], $begins[$i], $ends[$i], $ors[$i], $chrLengths[$i],$seqs[$i]) = ($1, $2, $3, $4, $5 ,$6);
	    $begins[$i]++;
	    $ends[$i] = $begins[$i] + $ends[$i] - 1;
	    if(length($names[$i]) > $maxl) {
		$maxl = length($names[$i]);
	    }
	}
    }
    my @index;
    $index[0] = $begins[0];
    my ($starti, $endi) = (-1, -1);
    for($starti = 0; $index[0] < $sitebeg && $index[0] < $ends[0]; $starti++) {
	if(substr($seqs[0], $starti, 1) !~ /\-/) {
	    $index[0]++;
	}
    }
    #start site
    for($endi = $starti; $index[0] < $siteend && $index[0] < $ends[0]; $endi++) {
	if(substr($seqs[0], $endi, 1) !~ /\-/) {
	    $index[0]++;
	}
    }
    if($starti >= 0 && $endi >= 0) {

	# check if there are "---" (gaps) at the beginnning of the alignment in the reference organism->remove
	my $seq=substr($seqs[0], $starti, $endi-$starti+1);

	my $introGAPs=0;
	my $outroGAPs=0;
	if ($seq=~/^(\-+)/){
	    $introGAPs=length($1);
	}
	if ($seq=~/(\-+)$/){
	    $outroGAPs=length($1);
	}

	# check if the complete block is empty now
	if ((length($seq)!=$introGAPs) and (length($seq)!=$outroGAPs)){

	    if ($outFormat eq "extMZ"){

		print STDOUT "$hsp $chr $or $sitebeg $siteend\n";
		my (@sitebegins, @siteends);
		for(my $i = 0; $i < @begins; $i++) {
		    my ($first, $mid);
		    $first = substr($seqs[$i], 0, $starti);
		    $mid = substr($seqs[$i], $starti, $endi-$starti+1);
		    $first =~ s/\-//g;
		    $mid =~ s/\-//g;
		    $sitebegins[$i] = $begins[$i]+length($first);
		    $siteends[$i] = $begins[$i]+length($first)+length($mid)-1;
		    my $ps=substr($seqs[$i], $starti+$introGAPs, $endi-$starti+1-$introGAPs-$outroGAPs);
		    
		    # in multiz files coordinates on the neg strand and calculated from the end
		    # compute the right ones based on the total length of the chromosome
		    if ($ors[$i] eq "-"){
			my $tempEnd=$chrLengths[$i]-$sitebegins[$i]+1;
			my $tempBegin=$chrLengths[$i]-$siteends[$i]+1;
			$sitebegins[$i]=$tempBegin;
			$siteends[$i]=$tempEnd	
			}
		    
		    
		    # reverse complement the seqs if on - strand
		    if ($or eq "-"){
			$ps=revCompAln($ps);
			$ors[$i]=~tr/\+\-/\-\+/;
		    }
		    
		    print STDOUT pad($names[$i], $maxl), "\t",  $ps,
		    "\t$sitebegins[$i]\t$siteends[$i]\t$ors[$i]\n";
		}
		print STDOUT "\n";
	    }elsif($outFormat eq "maf"){

		print STDOUT "a score=0.0\n";
		my (@sitebegins, @siteends);
		for(my $i = 0; $i < @begins; $i++) {
		    my ($first, $mid);
		    $first = substr($seqs[$i], 0, $starti);
		    $mid = substr($seqs[$i], $starti, $endi-$starti+1);
		    $first =~ s/\-//g;
		    $mid =~ s/\-//g;
		    $sitebegins[$i] = $begins[$i]+length($first);
		    $siteends[$i] = $begins[$i]+length($first)+length($mid)-1;
		    my $ps=substr($seqs[$i], $starti+$introGAPs, $endi-$starti+1-$introGAPs-$outroGAPs);
		    
		    # in multiz files coordinates on the neg strand and calculated from the end
		    # compute the right ones based on the total length of the chromosome
		    if ($ors[$i] eq "-"){
			my $tempEnd=$chrLengths[$i]-$sitebegins[$i]+1;
			my $tempBegin=$chrLengths[$i]-$siteends[$i]+1;
			$sitebegins[$i]=$tempBegin;
			$siteends[$i]=$tempEnd	
			}
		    
		    
		    # reverse complement the seqs if on - strand
		    if ($or eq "-"){
			$ps=revCompAln($ps);
			$ors[$i]=~tr/\+\-/\-\+/;
		    }
		    
		    #print STDOUT "s pad($names[$i], $maxl), "\t",  $ps,
		    #"\t$sitebegins[$i]\t$siteends[$i]\t$ors[$i]\t$chrLengths[$i]\n";
		    
		    my $alnLen=$siteends[$i]-$sitebegins[$i]+1;

		    print STDOUT "s ",pad($names[$i], $maxl),pad($sitebegins[$i],10),pad($alnLen,6)," $ors[$i] ",pad($chrLengths[$i],10),"$ps\n";

		}
		print STDOUT "\n";



	    }else{die "unknown input argument $outFormat\n";}
	}
    }
}



sub pad {
    my ($s, $l) = @_;
    while(length($s) < $l) {
	$s .= " ";
    }
    return $s;
}

sub revCompAln{
    my ($seq)=@_;

    $seq=~tr/ACGTacgt\-/TGCAtgca\-/;
    $seq=reverse $seq;

}
