package DeepSeqTools;

use base 'Exporter';
our @EXPORT = qw/getConfig filter aln2str str2aln solexaFastqToQual regInitFromArrays regInitFromFile findOverlaps findOverlapsC/;

use strict;
use warnings;
use FindBin;
use Time::HiRes qw/gettimeofday tv_interval/;
use Sys::Hostname;
use Inline C => "DATA",
    "NAME" => "deepSeqTools",
    #"DIRECTORY" => "$FindBin::RealBin/_Inline_".(hostname())."/",
    "DIRECTORY" => "$FindBin::RealBin/_Inline_"."BC2"."/",
    #"FORCE_BUILD" => 1,
    "PRINT_INFO" => 0,
    "BUILD_NOISY" => 0,
    ;

# globals variables used by region-overlap-finder
my %reg; # $reg{seqId}{strand} = [ [start0, start1, ...], [end0, end1, ...], [id0, id1, ...], maxLen ]


=head1 SYNOPSIS

 use FindBin;
 use lib "$FindBin::RealBin";
 use DeepSeqTools;


 my $filtCode = filter("AGCTACGGACTAGCGACTAGC");


 my $s1 = "ACGTACGAGCTTAGC";

 my $a1 = "AC--GTACGAGCTTAGC";
 my $a2 = "ACCCGTtCGA-CTTAGC";

 my $alnStr    = aln2str($a1, $a2);
 my ($b1, $b2) = str2aln($s1, $alnStr);

 if ($a1 eq $b1  and  $a2 eq $b2) { print "this must be TRUE\n"; }

=head1 DESCRIPTION

This module exports function that are used by the script processing
deep sequencing data sets.

=cut

=head1 GENERAL FUNCTIONS

=head2 getConfig

 Title   : getConfig
 Usage   : my $configVal = getConfig($configKey);
 Function: reads deepseq pipeline config file and get's config values
 Returns : hashref of config key/value pairs (if called without arguments)
           list of values corresponding to keys (if called with arguments)
 Args    : list of config keys (optional)

=cut

sub getConfig {
    my @keys = @_;

    my %conf;
    my $configFile = "config.tab";
    if(! -e $configFile)
    {
		# fallback to the FMI config. And for being consistent with the previous runs of these scripts.
		$configFile = "$FindBin::RealBin/../samples/config.tab";
		print STDERR "Using the global config file: $configFile\n";
	}
    open(CONF, "<$configFile") or die "ERROR reading $configFile: $!\n";
    while(my $l = <CONF>) {
	chomp $l;
	my @f = split(/\s+/, $l);
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

=head2 filter

 Title   : filter
 Usage   : my $filtCode = filter("AGCTACGGACTAGCGACTAGC");
 Function: Checks if a sequence string passes filtering rules
 Returns : 0 (zero), if the string passes all rules
           1, if the string is too short (<14)
           2, if the string contains too many Ns (>2)
           3, if the string is low complexisty
 Args    : sequence string

=cut

sub filter {
    my $seq = shift;

    # length filter
    if (length($seq) < 14) {
	return 1;

    # non-base character filter
    } elsif (($seq =~ tr/N/N/) > 2) {
	return 2;

    # entropy filter
    } else {
	if (calcH($seq, length($seq))/0.9770337 < 0.5) {
	    return 3;
	} else {
	    return 0;
	}
    }
}

=head2 aln2str

 Title   : aln2str
 Usage   : my $alnStr = aln2str($a1, $a2);
 Function: Transform a sequence alignment (query and target
           alignment strings) into a condensed string descriptor,
           relative to the query sequence.
           The descriptor is a comma-separated list of number-character
           elements, describing the non-identical positions in the
           alignment, in the form (\d+)(D|I\w|M\w), where:
             $1 : is the (one-based) offset in the query sequence
                  as contained in the alignemnt
             $2 : is the type of alignment, with:
                    D\w : a deletion (gap) in the query aligned to
                          (\w) in the target
                    I   : an insertion in the query (gap in the target)
                    M\w : a mismatch with the query character aligned
                          to \w in the target
 Returns : alignment descriptor string
 Args    : query alignment string (reference for the descriptor)
           target alignment string

=cut

sub aln2str {
    my ($a1, $a2) = @_;
    my @alnStr;
    my $cumGap1 = 0;
    for (my $i=0; $i<length($a1); $i++) {

	if (substr($a1,$i,1) ne substr($a2,$i,1)) {
	    if (substr($a1,$i,1) eq '-') {
		push @alnStr, ($i+1-$cumGap1).'D'.substr($a2,$i,1);
		$cumGap1++;

	    } elsif (substr($a2,$i,1) eq '-') {
		push @alnStr, ($i+1-$cumGap1).'I';

	    } else {
		push @alnStr, ($i+1-$cumGap1).'M'.substr($a2,$i,1);
	    }
	}
    }
    return join(",", @alnStr);
}

=head2 str2aln

 Title   : str2aln
 Usage   : my ($b1, $b2) = str2aln($s1, $alnStr);
 Function: Inverse the transformation of aln2str(), i.e.
           from the original query sequence and the alignment descriptor
           string, recreate the alignment.
 Returns : list of two string (query and target alignment strings)
 Args    : query sequence string
           alignment descriptor string

=cut

sub str2aln {
    my ($s1, $alnStr) = @_;
    my ($a1, $a2) = ($s1, $s1);
    my $cumGap1 = 0;
    foreach my $desc (split(",", $alnStr)) {
	$desc =~ m/^(\d+)([MDI])(\w?)$/;

	if ($2 eq 'D') {
	    substr($a1, $1-1+$cumGap1, 0, '-');
	    substr($a2, $1-1+$cumGap1, 0, $3);
	    $cumGap1++;

	} elsif ($2 eq 'I') {
	    substr($a2, $1-1+$cumGap1, 1, '-');

	} else {
	    substr($a2, $1-1+$cumGap1, 1, $3);
	}
    }
    return ($a1, $a2);
}

=head2 solexaFastqToQual

 Title   : solexaFastqToQual
 Usage   : my @qValues = solexaFastqToQual($qString);
 Function: Transform Solexa/Illumina single-byte per value
           quality string into numerical quality values.
 Returns : array of numerical quality values
 Args    : quality string

=cut

sub solexaFastqToQual {
    my ($qstring) = @_;
    return (map { int( 10*log(1+10**((ord($_)-64)/10.0))/log(10) + 0.5) } split("",$qstring));
}

=head1 COORDINATE REGION FUNCTIONS

=head2 regInitFromArrays

 Title   : regInitFromArrays
 Usage   : regInitFromArrays(\@regSeqId, \@regStart, \@regEnd, \@regStrand, \@regId, $clear);
 Function: Initialize data structures for region-overlap-finder from arrays.
 Returns : nothing
 Args    : array ref to sequence ids
           array ref to starts
           array ref to ends
           array ref to strands
           array ref to and region ids
           [optional] clear_flag (if defined, existing regions will be cleared)


=cut

sub regInitFromArrays {
    my ($seqId, $start, $end, $strand, $regId, $clear) = @_;

    unless(scalar(@$seqId)==scalar(@$start) and scalar(@$seqId)==scalar(@$end) and
	   scalar(@$seqId)==scalar(@$strand) and  scalar(@$seqId)==scalar(@$regId)) {
	die "error [regInitFromArrays]: arrays are not all of equal lengths\n";
    }

    if (defined $clear) { %reg = (); }

    # split regions by sequence identifier
    # $reg{seqId}{strand} = [ [start0, start1, ...], [end0, end1, ...], [id0, id1, ...], maxLen ]
    my ($len, $lists);
    for(my $i=0; $i<scalar(@$seqId); $i++) {
	if(!exists $reg{$seqId->[$i]}) {
	    $reg{$seqId->[$i]} = { $strand->[$i] => [[], [], [], 0] };
	} elsif(!exists $reg{$seqId->[$i]}{$strand->[$i]}) {
	    $reg{$seqId->[$i]}{$strand->[$i]} = [[], [], [], 0];
	}
	$lists = $reg{$seqId->[$i]}{$strand->[$i]};
	push @{$lists->[0]}, $start->[$i];
	push @{$lists->[1]}, $end->[$i];
	push @{$lists->[2]}, $regId->[$i];
	$len = $end->[$i] - $start->[$i] + 1;
	if($lists->[3] < $len) {
	    $lists->[3] = $len;
	}
    }

    # sort regions (ascending on end)
    _sortOnEnd();
}

=head2 regInitFromFile

 Title   : regInitFromFile
 Usage   : regInitFromFile($fname, $clear);
 Function: Initialize data structures for region-overlap-finder from
           coordinate/ral/fragment-report file.
 Returns : nothing
 Args    : file name of coordinate file with tab-delimited columns
               regionId, sequenceId, strand, start, end, e.g.:
                 NM_175941_promoter  chr2L   +    6529    8029
                 NM_078715_promoter  chr2L   -   18083   19583
                 NM_164348_promoter  chr2L   -   20872   22372

               These first 5 columns are also found in .ral and fragment
               report files. Any further columns will be ignored.

               The .ral file has 3 additional columns
               (nrErrors, qString, tString), and a fragment report has 6
               (nrErrors, qString, tString, sampleId, seqCount, seqInvWgt).
             
           [optional] clear_flag (if defined, existing regions will be cleared)

=cut

sub regInitFromFile {
    my ($fname, $clear) = @_;

    if (defined $clear) { %reg = (); }

    # store regions by sequence identifier
    # $reg{seqId}{strand} = [ [start0, start1, ...], [end0, end1, ...], [id0, id1, ...], maxLen ]
    my ($regId, $seqId, $regStrand, $regStart, $regEnd, $len, $lists);
    open(IN, "<$fname") or die "error in regInitFromCoordFile: could not open $fname: $!\n";
    while(<IN>) {
	if (/^#/) {
	    next;
	} else {
	    chomp;
	    ($regId, $seqId, $regStrand, $regStart, $regEnd) = split /\t/;
	    if(!exists $reg{$seqId}) {
		$reg{$seqId} = { $regStrand => [[], [], [], 0] };
	    } elsif(!exists $reg{$seqId}{$regStrand}) {
		$reg{$seqId}{$regStrand} = [[], [], [], 0];
	    }
	    $lists = $reg{$seqId}{$regStrand};
	    push @{$lists->[0]}, $regStart;
	    push @{$lists->[1]}, $regEnd;
	    push @{$lists->[2]}, $regId;
	    $len = $regEnd - $regStart + 1;
	    if($lists->[3] < $len) {
		$lists->[3] = $len;
	    }
	}
    }
    close IN;

    # sort regions (ascending on end)
    _sortOnEnd();
}

=head2 findOverlaps

 Title   : findOverlaps
 Usage   : my @regions = findOverlaps($seqId, $start, $end, $strand);
 Function: Find regions overlapping with the specified region.
 Returns : array of arrayrefs of overlapping regions, of the
           form: [regId, regSeqId, regStart, regEnd, regStrand]
 Args    : sequence id
           region start
           region end
           [optional] region strand

=cut

my ($t0, $t1, $t_find, $t_walk);
sub findOverlaps {
    my ($id, $s, $e, $strand) = @_;

    if(!exists $reg{$id}) {
	# sequence $id has no registered regions
	return ();
    } else {
	my @regions;

	if (defined $strand) {
	    if(exists $reg{$id}{$strand}) {
		# search only one strand
		my $ar  = $reg{$id}{$strand};
		my $s2  = $ar->[0];
		my $e2  = $ar->[1];
		my $id2 = $ar->[2];
		my $maxLen = $ar->[3];
		my $n = scalar(@$s2);
		
		my $i = _find_largest_element_less_than_or_equal($s, $e2);
 
		if (defined $i  or  $e2->[0]-$maxLen+1 <= $e) {
		    $i = $i || 0;
		    while($i<$n  and  $e2->[$i]-$maxLen < $e) {
			if($s2->[$i]<=$e  and  $e2->[$i]>=$s) {
			    push @regions, [ $id2->[$i], $id, $s2->[$i], $e2->[$i], $strand ];
			}
			$i++;
		    }
		}
	    }
	} else {
	    # search both strands
	    foreach $strand (keys %{$reg{$id}}) {
		my $ar  = $reg{$id}{$strand};
		my $s2  = $ar->[0];
		my $e2  = $ar->[1];
		my $id2 = $ar->[2];
		my $maxLen = $ar->[3];
		my $n = scalar(@$s2);
		
		my $i = _find_largest_element_less_than_or_equal($s, $e2);
 
		if (defined $i  or  $e2->[0]-$maxLen+1 <= $e) {
		    $i = $i || 0;
		    while($i<$n  and  $e2->[$i]-$maxLen < $e) {
			if($s2->[$i]<=$e  and  $e2->[$i]>=$s) {
			    push @regions, [ $id2->[$i], $id, $s2->[$i], $e2->[$i], $strand ];
			}
			$i++;
		    }
		}
	    }
	}

	return reverse @regions;
    }
}

=head2 _register_region

 Title   : _register_region
 Usage   : [internal]
 Function: register a single region (push it to the end of existing regions)
           WARNING: use with care, as after this call, regions are not sorted
                    anymore. _sortOnEnd() will have to be called before
                    findOverlaps() can be used again.
 Returns : nothing
 Args    : regId, regSeqId, regStrand, regStart, regEnd

=cut

sub _register_region {
    # $_[0]  $_[1]  $_[2]      $_[3]     $_[4]
    # regId  seqId  regStrand  regStart  regEnd
    if(!exists $reg{$_[1]}) {
	$reg{$_[1]} = { $_[2] => [[], [], [], 0] };
    } elsif(!exists $reg{$_[1]}{$_[2]}) {
	$reg{$_[1]}{$_[2]} = [[], [], [], 0];
    }
    my $lists = $reg{$_[1]}{$_[2]};
    push @{$lists->[0]}, $_[3];
    push @{$lists->[1]}, $_[4];
    push @{$lists->[2]}, $_[0];
    my $len = $_[4] - $_[3] + 1;
    if($lists->[3] < $len) {
	$lists->[3] = $len;
    }
}

=head2 _register_region_sorted

 Title   : _register_region_sorted
 Usage   : [internal]
 Function: register a single region at its sorted position
           (push it to the appropriate place into array of existing regions)
 Returns : nothing
 Args    : regId, regSeqId, regStrand, regStart, regEnd

=cut

sub _register_region_sorted {
    # $_[0]  $_[1]  $_[2]      $_[3]     $_[4]
    # regId  seqId  regStrand  regStart  regEnd
    if(!exists $reg{$_[1]}) {
	$reg{$_[1]} = { $_[2] => [[], [], [], 0] };
	push @{$reg{$_[1]}{$_[2]}[0]}, $_[3];
	push @{$reg{$_[1]}{$_[2]}[1]}, $_[4];
	push @{$reg{$_[1]}{$_[2]}[2]}, $_[0];
	$reg{$_[1]}{$_[2]}[3] = $_[4] - $_[3] + 1;
    } elsif(!exists $reg{$_[1]}{$_[2]}) {
	$reg{$_[1]}{$_[2]} = [[], [], [], 0];
	push @{$reg{$_[1]}{$_[2]}[0]}, $_[3];
	push @{$reg{$_[1]}{$_[2]}[1]}, $_[4];
	push @{$reg{$_[1]}{$_[2]}[2]}, $_[0];
	$reg{$_[1]}{$_[2]}[3] = $_[4] - $_[3] + 1;
    } else {
	my $lists = $reg{$_[1]}{$_[2]};
	my $inx = _find_largest_element_less_than_or_equal($_[4], $lists->[1]);
	$inx = defined $inx ? $inx + 1 : 0; # WARNING: may be +n rather than +1, in case there are multiple identical!
	splice @{$lists->[0]}, $inx, 0, $_[3];
	splice @{$lists->[1]}, $inx, 0, $_[4];
	splice @{$lists->[2]}, $inx, 0, $_[0];
	my $len = $_[4] - $_[3] + 1;
	if($lists->[3] < $len) {
	    $lists->[3] = $len;
	}
    }
}

=head2 _clear

 Title   : _clear
 Usage   : [internal]
 Function: clear all registered regions
 Returns : nothing
 Args    : none

=cut

sub _clear {
    %reg = ();
}

=head2 _find_largest_element_less_than_or_equal

 Title   : _find_largest_element_less_than_or_equal
 Usage   : [internal]
 Function: Find largest element less than the specified value.
 Returns : index of largest element less than or equal to value
           undef if all elements are larger than value
 Args    : value to compare to
           reference to array of ascendingly sorted values

=cut

sub _find_largest_element_less_than_or_equal {
    my ($val, $aref) = @_;

    if ($aref->[0] > $val) {
	return undef;

    } else {
	my $low = 0;
	my $high = scalar(@$aref) - 1;
	my $mid;
	while($low <= $high) {
	    $mid = int(($low + $high) / 2); # prevent int overflow by: $mid = $low + (($high - $low) / 2);
	    if($aref->[$mid] > $val) {
		$high = $mid - 1;
	    } elsif($aref->[$mid] < $val) {
		$low = $mid + 1;
	    } else {
		while($mid > 0 and $aref->[$mid-1] == $val) {
		    $mid--; # find first of multiple identical
		}
		return $mid; # found
	    }
	}
	return $low-1; # not found
    }
}

=head2 _sortOnEnd

 Title   : _sortOnEnd
 Usage   : [internal]
 Function: Sort regions acendingly by end coordinate.
 Returns : nothing
 Args    : none

=cut

sub _sortOnEnd {
    foreach my $sId (keys %reg) {
	foreach my $strand (keys %{$reg{$sId}}) {
	    my $ends = $reg{$sId}{$strand}[1];
	    my @inx = sort {$ends->[$a] <=> $ends->[$b]} 0..(scalar(@$ends)-1);
	    @{$reg{$sId}{$strand}[0]} = @{$reg{$sId}{$strand}[0]}[@inx];
	    @{$reg{$sId}{$strand}[1]} = @{$reg{$sId}{$strand}[1]}[@inx];
	    @{$reg{$sId}{$strand}[2]} = @{$reg{$sId}{$strand}[2]}[@inx];
	}
    }
}

=head1 REGIION FUNCTIONS IN C

=head2 findOverlapsC

 Title   : findOverlapsC
 Usage   : my @regions = findOverlapsC($start, $end, $strand, $minOverlap);
 Function: Find regions overlapping with the specified region.
           Regions must be properly registered by _initC, _clearC, 
           _register_regionC and _sortOnEndC
 Returns : array of arrayrefs of overlapping regions, of the
           form: [regId, undef, regStart, regEnd, regStrand]
 Args    : region start
           region end
           region strand ('+', '-' or '*')
           minimum overlap

=cut

sub findOverlapsC {
    my @res = _findOverlapsC($_[1], $_[2], $_[3], $_[4]);
    my @res2;
    for(my $i=0; $i<@res; $i+=4) {
	push @res2, [ $res[$i], $_[0], $res[$i+1], $res[$i+2], $res[$i+3] ];
    }
    return reverse( @res2 );
}

1;

__DATA__
__C__
#include <stdlib.h>
#include <string.h>
#include <math.h>

double freq[16];

int base2int(char b) {
    switch(b) {
	case 'A': return 0;
	case 'C': return 1;
	case 'G': return 2;
	case 'T': return 3;
        default:  return -1;
    }
}

double calcH(char* seq, int l) {
    int i, n=0, d1, d2;
    double Hseq=0.0;

    for(i=0; i<16; i++) { freq[i] = 0.0; }

    for(i=0; i<l-1; i++) {
	d1 = base2int(seq[i]);
	d2 = base2int(seq[i+1]);

	if (d1!=-1 && d2!=-1) {
	    n++;
	    freq[d1+d2*4]++;
	}
    }

    for(i=0; i<16; i++) {
	if(freq[i]>0.0) {
	    freq[i] /= n;
	    Hseq += (freq[i] * log(freq[i])/log(16));
	}
    }

    return -Hseq;
}

int maxRegions, nRegions=0, maxLen=0;
char** regId;
char* regStrand;
int *regStart, *regEnd;
void _clearC() {
    int i;
    if(nRegions>0) {
        for(i=0; i<nRegions; i++) {
            free(regId[i]);
        }
        free(regId);
        free(regStrand);
        free(regStart);
        free(regEnd);
        maxRegions = nRegions = maxLen = 0;
    }
}

void _initC(int maxR) {
    regId = (char**)malloc(maxR * sizeof(char*));
    regStrand = (char*)malloc(maxR * sizeof(char));
    regStart = (int*)malloc(maxR * sizeof(int));
    regEnd = (int*)malloc(maxR * sizeof(int));
    maxRegions = maxR;
    nRegions = 0;
    maxLen = 0;
}

void _register_regionC(char* rId, char* rStrand, int rStart, int rEnd) {
    regId[nRegions] = (char*)malloc( (strlen(rId)+1) * sizeof(char) );
    strcpy(regId[nRegions], rId);
    regStrand[nRegions] = rStrand[0];
    regStart[nRegions] = rStart;
    regEnd[nRegions] = rEnd;

    if(rEnd - rStart + 1 > maxLen) {
        maxLen = rEnd - rStart + 1;
    }

    nRegions++;
}

int _compareEnds(const void * a, const void * b) {
  return ( regEnd[ *(int*)a ]  - regEnd[ *(int*)b ] );
}

void _sortOnEndC() {
    int i, *tmpStart, *tmpEnd, *order;
    char **tmpId, *tmpStrand;
    //# allocate new arrays
    tmpId = (char**)malloc(nRegions * sizeof(char*));
    tmpStrand = (char*)malloc(nRegions * sizeof(char));
    tmpStart = (int*)malloc(nRegions * sizeof(int));
    tmpEnd = (int*)malloc(nRegions * sizeof(int));
    //# find permutation
    order = (int*)malloc(nRegions * sizeof(int));
    for(i=0; i<nRegions; i++) { order[i] = i; }
    qsort(order, nRegions, sizeof(int), _compareEnds);
    //# permute
    for(i=0; i<nRegions; i++) {
        tmpId[i] = regId[ order[i] ];
	tmpStrand[i] = regStrand[ order[i] ];
	tmpStart[i] = regStart[ order[i] ];
	tmpEnd[i] = regEnd[ order[i] ];
    }
    //# clean up old and update to new
    free(order);
    free(regId); regId = tmpId;
    free(regStrand); regStrand = tmpStrand;
    free(regStart); regStart = tmpStart;
    free(regEnd); regEnd = tmpEnd;
}

//# Returns : index of largest element less than or equal to value
//#           -1 if all elements are larger than value
//# Args    : value to compare to
//#           reference to array of ascendingly sorted values
int _find_largest_element_less_than_or_equalC(int val, int* aref, int aSize) {
    int low, high, mid;

    if (aref[0] > val) {
	return -1;

    } else {
	low = 0;
	high = aSize - 1;
	while (low <= high) {
	    mid = low + ((high - low) / 2); //# prevents int overflows
	    if(aref[mid] > val) {
		high = mid - 1;
	    } else if(aref[mid] < val) {
		low = mid + 1;
	    } else {
		while(mid > 0 && aref[mid-1] == val) {
		    mid--; //# find first of multiple identical
		}
		return mid; //# found
	    }
	}
	return low-1; //# not found
    }
}



//# Returns : array of arrayrefs of overlapping regions, of the
//#           form: [regId, undef, regStart, regEnd, regStrand]
//# Args    : region start
//#           region end
//#           region strand ('+', '-' or '*')
//#           minimum overlap
void _findOverlapsC(int start, int end, char* str, int minOverlap) {
    int i, ov, n = 0;
    //# AV *ovInfo;
    char strand = str[0];
    Inline_Stack_Vars;
    Inline_Stack_Reset;

    if (nRegions > 0) {
	i = _find_largest_element_less_than_or_equalC(start, regEnd, nRegions);

	if (i != -1  ||  regEnd[0]-maxLen+1 <= end) {
	    i = (i != -1) ? i : 0;

	    if (strand != '*') {
		//# search only one strand
		while (i<nRegions  &&  regEnd[i]-maxLen < end) {
		    if (regStart[i]<=end  &&  regEnd[i]>=start  &&  regStrand[i]==strand) {
			ov = (end<regEnd[i] ? end : regEnd[i]) - (start>regStart[i] ? start : regStart[i]) + 1;
			if (ov >= minOverlap) {
			    //# ovInfo = (AV*) sv_2mortal( (SV*)newAV() ); //#see perldoc perlapi
			    //# av_push( ovInfo, sv_2mortal(newSVpv(regId[i], 0)) );     //# regId
			    //# av_push( ovInfo, &PL_sv_undef );                         //# undef (formerly targetSeqId)
			    //# av_push( ovInfo, sv_2mortal(newSViv(regStart[i])) );     //# regStart
			    //# av_push( ovInfo, sv_2mortal(newSViv(regEnd[i])) );       //# regEnd
			    //# av_push( ovInfo, sv_2mortal(newSVpvn(regStrand+i,1)) );  //# regStrand
			   
			    //# av_push( results, newRV((SV *)ovInfo) );
			    Inline_Stack_Push( sv_2mortal(newSVpv(regId[i], 0)) );     //# regId
			    Inline_Stack_Push( sv_2mortal(newSViv(regStart[i])) );     //# regStart
			    Inline_Stack_Push( sv_2mortal(newSViv(regEnd[i])) );       //# regEnd
			    Inline_Stack_Push( sv_2mortal(newSVpvn(regStrand+i,1)) );  //# regStrand
			    n++;
			}
		    }
		    i++;
		}

	    } else {
		//# search both strands
		while (i<nRegions  &&  regEnd[i]-maxLen < end) {
		    if (regStart[i]<=end  &&  regEnd[i]>=start) {
			ov = (end<regEnd[i] ? end : regEnd[i]) - (start>regStart[i] ? start : regStart[i]) + 1;
			if (ov >= minOverlap) {
			    //# ovInfo = (AV*) sv_2mortal( (SV*)newAV() ); //#see perldoc perlapi
			    //# av_push( ovInfo, sv_2mortal(newSVpv(regId[i], 0)) );     //# regId
			    //# av_push( ovInfo, &PL_sv_undef );                         //# undef (formerly targetSeqId)
			    //# av_push( ovInfo, sv_2mortal(newSViv(regStart[i])) );     //# regStart
			    //# av_push( ovInfo, sv_2mortal(newSViv(regEnd[i])) );       //# regEnd
			    //# av_push( ovInfo, sv_2mortal(newSVpvn(regStrand+i,1)) );  //# regStrand
			   
			    Inline_Stack_Push( sv_2mortal(newSVpv(regId[i], 0)) );     //# regId
			    Inline_Stack_Push( sv_2mortal(newSViv(regStart[i])) );     //# regStart
			    Inline_Stack_Push( sv_2mortal(newSViv(regEnd[i])) );       //# regEnd
			    Inline_Stack_Push( sv_2mortal(newSVpvn(regStrand+i,1)) );  //# regStrand
			    n++;
			    //# fprintf(stderr, "\tjust stored in C: [%s, %d, %d, %c]\n", regId[i], regStart[i], regEnd[i], regStrand[i]);
			}
		    }
		    i++;
		}
	    }
	}
	
    }
    Inline_Stack_Done;
    Inline_Stack_Return(4*n);
}
