#!/usr/bin/env perl
#
# filterAdaptors.pl
# filter out adaptor sequences from sequence reads
#
# Michael Stadler, started May 21, 2008
#
# HISTORY:
# March 9, 2009:  replaced step 1. (SW alignments) by mismatch-tolerant string matcher
# April 26, 2010: minimal removed adaptor length parameters
#
# TODO:

use strict;
use warnings;
use Getopt::Std;
use File::Temp;
use FindBin;
use lib "$FindBin::RealBin";
use DeepSeqTools;


### GLOBALS
my %options = ('h' => undef,                 # display help
	       'i' => "",                    # input sequences
	       'F' => "fasta",               # input format
	       '5' => "",                    # 5'-adaptor
	       '3' => "",                    # 3'-adaptor
	       'm' => 0,                     # minimal 5'-adaptor alignment length
	       'n' => 0,                     # minimal 3'-adpator alignment length
	       '1' => 7,                     # match length to tolerate 1 mismatch
	       '2' => 10,                    # match length ot lolerate 2 mismatches
	       'f' => 2,                     # max mismatches for full length adapter hits
	       'a' => undef,                 # allow indels as mismatches
	       'L' => undef,                 # log file
	       's' => undef,                 # short log file format
	       'v' => undef,                 # disable sequence filtering
	      );
my $usage=<<EOU;
$0 :
  unify sequences, filter out adaptor sequences, unify again and filter

The adaptor processing is done in three steps:
  1. full adaptor sequence matching with up to '-f' mismatches,
     with any sequence preceeding/trailing the 5'/3' adaptor, repectively.
  2. medium length pre-/suffix matching allowing 1 or 2 mismatches
     (starting with adaptor length minus 1, decreasing down to '-1')
  3. exact matching of short pre-/suffixes
     (starting with length '-1'-1, down to 1)

After adapters are removed, reads are unified, and filtered based on
length, non-base characters and entropy.

Results are written to STDOUT as a fasta library (identifier randomly
chosen amongst identical sequences, description is the sequence count).
Optionally, a log file is produced (-L).

usage: $0 [options] -i input_seqs.fa

with valid options:
    -h      : display help

  sequence parameters:
    -i file : input sequences                               [$options{i}]
    -F str  : input format (fasta|tab|fastq)                [$options{F}]
    -5 str  : 5'-adaptor sequence string                    [$options{5}]
    -3 str  : 3'-adaptor sequence string                    [$options{3}]
    -m int  : minimal 5'-adaptor alignment length           [$options{m}]
    -n int  : minimal 3'-adaptor alignment length           [$options{n}]
    -v      : disable sequence filtering                    [filter sequences]

  search paramters:
    -1 int  : match length to tolerate 1 mismatch           [$options{1}]
    -2 int  : match length ot lolerate 2 mismatches         [$options{2}]
    -f int  : max mismatches for full length adaptor hits   [$options{f}]
    -a      : allow indels as mismatches                    [no indels]

  log file:
    -L str  : log file with per read and summary info       [do not create a log file]
    -s      : produce short log file                        [produce long log file]


  comments:

  Formats supported for tab input (tab delimited, one sequence per line):
	    1. seq
	    2. id seq
	    3. id seq counts
	    4. seq counts

EOU


### MAIN
# digest command line
getopts('hi:F:5:3:m:n:1:2:f:aL:sv', \%options);
$options{F} = lc($options{F});
if ($options{F} !~ m/^(fasta|tab|fastq)$/) {
    print STDERR "ERROR: unknown input format '$options{F}'\n";
    $options{h} = 1;
}
$options{5} = uc($options{5});
if ($options{5} ne "" and $options{5} =~ m/[^ACGT]/) {
    print STDERR "ERROR: found non-base characters in 5'-adaptor: $options{5}\n";
    $options{h} = 1;
}
$options{3} = uc($options{3});
if ($options{3} ne "" and $options{3} =~ m/[^ACGT]/) {
    print STDERR "ERROR: found non-base characters in 3'-adaptor: $options{3}\n";
    $options{h} = 1;
}
if ($options{i} eq "") {
    print STDERR "ERROR: must specify input sequence file (-i)\n";
    $options{h} = 1;
}
if ($options{2} < $options{1}) {
    print STDERR "ERROR: -2 $options{2} must be greater than -1 $options{1}\n";
    $options{h} = 1;
}
if (defined $options{h}) { die $usage; }


# prepare
my $doReport = defined($options{s}) ? sub {} : \&report;
my %count = (nSeqs => count_seqs($options{i}, $options{F}),
	     nUniqueSeqs => 0,
	     rm5_exact => 0, rm5_mm => 0, rm5_full => 0,
	     rm3_exact => 0, rm3_mm => 0, rm3_full => 0,
	     filt_alnLen5 => 0, filt_alnLen3 => 0,
	     filt_len => 0, filt_N => 0, filt_H => 0,
	     passed => 0
	    );
my @filt = qw/passed filt_len filt_N filt_H/;
my (@adapt5, @adapt3, @var5, @var3, %full5, %full3);
my ($len5, $len3) = (length($options{5}), length($options{3}));
if ($len5 > 0) {
    %full5 = %{generate_variants($options{5}, 5, "", $options{f}, $options{a})};
    @adapt5 = ("", map {substr($options{5},-$_)} 1..$len5);
    while(scalar(@adapt5)<$options{1}) {
	push @adapt5, "";
    }
    for(my $i=$options{1}; $i<=$len5-1; $i++) {
	$var5[$i] = generate_variants(substr($options{5}, -$i), 5, substr($options{5}, -$i-1, 1), $i>=$options{2} ? 2 : 1, $options{a});
	#print STDERR scalar(keys %{$var5[$i]})." ".($i>=$options{2} ? 2 : 1)."mm-variants of ".substr($options{5}, -$i).":\n".join("\n", sort keys %{$var5[$i]})."\n";
    }
}
if ($len3 > 0) {
    %full3 = %{generate_variants($options{3}, 3, "", $options{f}, $options{a})};
    @adapt3 = ("", map {substr($options{3},0,$_)} 1..$len3);
    while(scalar(@adapt3)<$options{1}) {
	push @adapt3, "";
    }
    for(my $i=$options{1}; $i<=$len3-1; $i++) {
	$var3[$i] = generate_variants(substr($options{3}, 0, $i), 3, substr($options{3}, $i, 1), $i>=$options{2} ? 2 : 1, $options{a});
	#print STDERR scalar(keys %{$var3[$i]})." ".($i>=$options{2} ? 2 : 1)."mm-variants of ".substr($options{3}, 0, $i).":\n".join("\n", sort keys %{$var3[$i]})."\n";
    }
}
my $LOG;
# By Severin: Some modifications so that no paths or other more private information gets shown on the web report interface...
if (defined $options{L}) {
    open($LOG, ">$options{L}") or die "ERROR writing to log file $options{L}: $!\n";
    #print $LOG "# output of $0, started on ".scalar(localtime)."\n";
    #print $LOG "# input file          : $options{i} (n_unique=$count{nSeqs})\n";
    print $LOG "# 5'/3' adaptors      : $options{5} / $options{3}\n";
    print $LOG "# 5'/3' min aln. len. : $options{m} / $options{n}\n";
    print $LOG "# min len (1 mm)      : $options{1}\n";
    print $LOG "# min len (2 mm)      : $options{2}\n";
    print $LOG "# max mm (full)       : $options{f}\n";
    #print $LOG "# log file            : $options{L}\n";
    print $LOG "# sequence filter     : ".(defined($options{v}) ? "off" : "on")."\n";
}


# store sequences
if(defined $LOG) { print $LOG "# ".scalar(localtime).": READING SEQUENCES\n"; }
&$doReport("# id\tcomment");
my %seq; # format: $seq{sequence} = [id, count, newseq, rm5, rm3]
keys(%seq) = $count{nSeqs};
store_seqs($options{i}, $options{F}, \%seq);
# fix $count{nSeqs} (may be wrong as count_seqs() counts unique seqs and ignores seq counts)
$count{nSeqs} = 0;
$count{nUniqueSeqs} = scalar(keys %seq);
if(defined $LOG) { print $LOG "# ".scalar(localtime).": done reading sequences (n=$count{nUniqueSeqs})\n"; }


# removing adaptors
my %newseq; # cleaned seqs, format: $newseq{sequence} = [id, count]
keys(%newseq) = $count{nSeqs};
if(defined $LOG) { print $LOG "# ".scalar(localtime).": REMOVING ADAPTORS\n"; }
&$doReport("# ".join("\t", qw/id seq.before seq.after removed.5p type.5p removed.3p type.3p/));
my ($i, $seq, $newseq, $newseq2, $mmcount, $comment);

# 1. full adaptor matching
$count{nSeqs} = 0; # fix $count{nSeqs} (may be wrong as count_seqs() counts unique seqs and ignores seq counts)
foreach $seq (keys %seq) {
    $count{nSeqs} += $seq{$seq}[1];
    # ...match 5'-adaptor
    $mmcount = undef;
    if ($len5 > 0) {
	for($i=length($seq)-$len5+1; $i>=0; $i--) {
	#for($i=0; $i<length($seq)-$len5+2; $i++) {
	    # remark: loop runs over end of string, but substr() just truncates at end
	    if(exists($full5{substr($seq, $i, $len5-1)})) {
		# deletion hit?
		$mmcount = $full5{substr($seq, $i, $len5-1)};
		$newseq = substr($seq, $i+$len5-1);
		$comment = ($i+1)."-".($i+$len5-1)."\tfull $mmcount";
		last;
	    } elsif(exists($full5{substr($seq, $i, $len5)})) {
		# perfect/mismatch hit?
		$mmcount = $full5{substr($seq, $i, $len5)};
		$newseq = substr($seq, $i+$len5);
		$comment = ($i+1)."-".($i+$len5)."\tfull $mmcount";
		last;
	    } elsif(exists($full5{substr($seq, $i, $len5+1)})) {
		# insertion hit?
		$mmcount = $full5{substr($seq, $i, $len5+1)};
		$newseq = substr($seq, $i+$len5+1);
		$comment = ($i+1)."-".($i+$len5+1)."\tfull $mmcount";
		last;
	    }
	}
	if(defined $mmcount) {
	    $seq{$seq}[2] = $newseq;
	    $seq{$seq}[3] = $comment; 
	    $count{rm5_full} += $seq{$seq}[1];
	}
    }

    # ...match 3'-adaptor
    if ($len3 > 0) {
	if(!defined $mmcount) {
	    $newseq = $seq;
	} else {
	    $mmcount = undef;
	}
	for($i=0; $i<length($newseq)-$len3+2; $i++) {
	    # remark: loop runs over end of string, but substr() just truncates at end
	    if(exists($full3{substr($newseq, $i, $len3-1)})) {
		# deletion hit?
		$mmcount = $full3{substr($newseq, $i, $len3-1)};
		$newseq2 = substr($newseq, 0, $i);
		$comment = (length($seq)-length($newseq)+$i+1)."-".(length($seq)-length($newseq)+$i+$len3-1)."\tfull $mmcount";
		last;
	    } elsif(exists($full3{substr($newseq, $i, $len3)})) {
		# perfect/mismatch hit?
		$mmcount = $full3{substr($newseq, $i, $len3)};
		$newseq2 = substr($newseq, 0, $i);
		$comment = (length($seq)-length($newseq)+$i+1)."-".(length($seq)-length($newseq)+$i+$len3)."\tfull $mmcount";
		last;
	    } elsif(exists($full3{substr($newseq, $i, $len3+1)})) {
		# insertion hit?
		$mmcount = $full3{substr($newseq, $i, $len3+1)};
		$newseq2 = substr($newseq, 0, $i);
		$comment = (length($seq)-length($newseq)+$i+1)."-".(length($seq)-length($newseq)+$i+$len3+1)."\tfull $mmcount";
		last;
	    }
	}
	if(defined $mmcount) {
	    $seq{$seq}[2] = $newseq2;
	    $seq{$seq}[4] = $comment; 
	    $count{rm3_full} += $seq{$seq}[1];
	}
    }
}

#    5'-adaptor
if ($len5 > 0) {
    foreach $seq (keys %seq) {
	$newseq = $seq{$seq}[2];

	# 2. pre-/suffix matching (up to length $len5-1)
	if ($seq{$seq}[3] eq "\t") {
	    for($i=$len5-1; $i>=$options{1}; $i--) {
		if ( exists $var5[$i]{substr($newseq, 0, $i)} ) {
		    $newseq = substr($newseq, $i);
		    $count{rm5_mm} += $seq{$seq}[1];
		    $seq{$seq}[3] = "1-$i\tmismatch ".$var5[$i]{substr($newseq, 0, $i)};   # REMARK: throws warning if full match removed whole sequence
		    last;
		}
	    }

	    # 3. exact pre-/suffix matching (starting with length '-1'-1, down to 1)
	    if ($seq{$seq}[3] eq "\t") {
		for($i=$options{1}-1; $i>0; $i--) {
		    #if (index($newseq, $adapt5[$i]) == 0) {
		    if (substr($newseq, 0, $i) eq $adapt5[$i]) {
			$newseq = substr($newseq, $i);
			$count{rm5_exact} += $seq{$seq}[1];
			$seq{$seq}[3] = "1-$i\texact";
			last;
		    }
		}
	    }

	    $seq{$seq}[2] = $newseq;
	}
    }
}

#    3'-adaptor
if ($len3 > 0) {
    foreach $seq (keys %seq) {
	$newseq = $seq{$seq}[2];

	# 2. pre-/suffix matching (up to length $len3-1)
	if ($seq{$seq}[4] eq "\t") {
	    for($i=$len3-1; $i>=$options{1}; $i--) {
		if ( exists $var3[$i]{substr($newseq, -$i)} ) {
		    $count{rm3_mm} += $seq{$seq}[1];
		    $seq{$seq}[4] = (length($seq)-$i+1)."-".length($seq)."\tmismatch ".$var3[$i]{substr($newseq, -$i)};
		    $newseq = substr($newseq, 0, -$i);
		    last;
		}
	    }

	    # 3. exact pre-/suffix matching (starting with length '-1'-1, down to 1)
	    if ($seq{$seq}[4] eq "\t") {
		for($i=$options{1}-1, my $pos=length($newseq)-$options{1}+1; $i>0; $i--, $pos++) {
		    #if (index($newseq, $adapt3[$i], $pos) == $pos) {
		    if (substr($newseq, $pos, $i) eq $adapt3[$i]) {    # REMARK: throws warning for adaptors that are shorter than $options{1}
			$count{rm3_exact} += $seq{$seq}[1];
			$seq{$seq}[4] = (length($seq)-length($newseq)+$pos+1)."-".length($seq)."\texact";
			$newseq = substr($newseq, 0, $pos);
			last;
		    }
		}
	    }

	    $seq{$seq}[2] = $newseq;
	}
    }
}

# report on adaptor removal and condense all in %newseq
my (@filtAlnLen5, @filtAlnLen3); # remember filered sequence for filtering report
my ($nFiltAlnLen5, $nFiltAlnLen3) = (0, 0);
foreach $seq (keys %seq) {
    if ($seq{$seq}[3] ne "\t"  or  $seq{$seq}[4] ne "\t") {
	&$doReport(join("\t", $seq{$seq}[0], $seq, $seq{$seq}[2], $seq{$seq}[3], $seq{$seq}[4]));
    }

    # first check if minimal adaptor alignment lengthes were achieved
    if(adaptorAlignmentLength($seq{$seq}[3]) < $options{m}) {
	# ...5' adaptor alignment too short?
	$count{filt_alnLen5} += $seq{$seq}[1];
	$nFiltAlnLen5++;
	if(not defined($options{s})) { push @filtAlnLen5, $seq; }
	
    } elsif(adaptorAlignmentLength($seq{$seq}[4]) < $options{n}) {
	# ...3' adaptor alignment too short?
	$count{filt_alnLen3} += $seq{$seq}[1];
	$nFiltAlnLen3++;
	if(not defined($options{s})) { push @filtAlnLen3, $seq; }

    } else {
	if (exists $newseq{$seq{$seq}[2]}) {
	    $newseq{$seq{$seq}[2]}[1] += $seq{$seq}[1];
	} else {
	    $newseq{$seq{$seq}[2]} = $seq{$seq};
	}
	#delete $seq{$seq};
    }
}
$count{nUniqueSeqs2} = scalar(keys %newseq) + $nFiltAlnLen5 + $nFiltAlnLen3;
if(defined $LOG) { print $LOG "# ".scalar(localtime).": done removing adaptors (n=$count{nUniqueSeqs2})\n"; }


# filtering
if(defined $LOG) { print $LOG "# ".scalar(localtime).": FILTERING AND OUTPUT\n"; }
&$doReport(join("\t", "# id", "seq", "filter.reason"));
my $f;
if(not defined($options{s})) {
    foreach $seq (@filtAlnLen5) {  &$doReport(join("\t", $seq{$seq}[0], $seq{$seq}[2], "filt_alnLen5"));  }
    foreach $seq (@filtAlnLen3) {  &$doReport(join("\t", $seq{$seq}[0], $seq{$seq}[2], "filt_alnLen3"));  }
}
foreach $seq (keys %newseq) {
    $f = defined($options{v}) ? 0 : filter($seq);
    $count{ $filt[$f] } += $newseq{$seq}[1];
    if ($f == 0) {
	print STDOUT ">$newseq{$seq}[0] $newseq{$seq}[1]\n";
	print STDOUT "$seq\n";
    } else {
	&$doReport(join("\t", $newseq{$seq}[0], $seq, $filt[$f]));
    }
}
if(defined $LOG) { print $LOG "# ".scalar(localtime).": done FILTERING AND OUTPUT\n"; }


# make summary
if(not defined $LOG) { $LOG = \*STDERR; }
print $LOG "# SUMMARY:\n";
print $LOG "# Number of sequences:\n";
print $LOG sprintf("# % 38s : % 8d (%5.1f%%)\n", "input sequences",  $count{nSeqs}, 100*$count{nSeqs}/$count{nSeqs});
print $LOG sprintf("# % 38s : % 8d (%5.1f%%)\n", "unique sequences", $count{nUniqueSeqs}, 100*$count{nUniqueSeqs}/$count{nSeqs});
print $LOG "#\n";
print $LOG "# Adaptor removal:\n";
print $LOG sprintf("# % 38s : % 8d (%5.1f%%)\n", "full matches to 5'-adaptor",      $count{rm5_full},    100*$count{rm5_full}/$count{nSeqs});
print $LOG sprintf("# % 38s : % 8d (%5.1f%%)\n", "inexact suffix matches to 5'-adaptor", $count{rm5_mm},    100*$count{rm5_mm}/$count{nSeqs});
print $LOG sprintf("# % 38s : % 8d (%5.1f%%)\n", "exact suffix matches to 5'-adaptor",   $count{rm5_exact}, 100*$count{rm5_exact}/$count{nSeqs});
print $LOG sprintf("# % 38s : % 8d (%5.1f%%)\n", "full matches to 3'-adaptor",      $count{rm3_full},    100*$count{rm3_full}/$count{nSeqs});
print $LOG sprintf("# % 38s : % 8d (%5.1f%%)\n", "inexact prefix matches to 3'-adaptor", $count{rm3_mm},    100*$count{rm3_mm}/$count{nSeqs});
print $LOG sprintf("# % 38s : % 8d (%5.1f%%)\n", "exact prefix matches to 3'-adaptor",   $count{rm3_exact}, 100*$count{rm3_exact}/$count{nSeqs});
print $LOG "#\n";
print $LOG sprintf("# % 38s : % 8d (%5.1f%%)\n", "total matches to 5'-adaptor",   $count{rm5_exact}+$count{rm5_mm}+$count{rm5_full},
	       100*($count{rm5_exact}+$count{rm5_mm}+$count{rm5_full})/$count{nSeqs});
print $LOG sprintf("# % 38s : % 8d (%5.1f%%)\n", "total matches to 3'-adaptor",   $count{rm3_exact}+$count{rm3_mm}+$count{rm3_full},
	       100*($count{rm3_exact}+$count{rm3_mm}+$count{rm3_full})/$count{nSeqs});
print $LOG "#\n";
print $LOG sprintf("# % 38s : % 8d (%5.1f%%)\n", "unique sequences (w/o adaptors)", $count{nUniqueSeqs2}, 100*$count{nUniqueSeqs2}/$count{nSeqs});
print $LOG "#\n";
print $LOG "# Sequence filtering:\n";
print $LOG sprintf("# % 38s : % 8d (%5.1f%%)\n", "5'-adaptor alignment too short",  $count{filt_alnLen5}, 100*$count{filt_alnLen5}/$count{nSeqs});
print $LOG sprintf("# % 38s : % 8d (%5.1f%%)\n", "3'-adaptor alignment too short",  $count{filt_alnLen3}, 100*$count{filt_alnLen3}/$count{nSeqs});
print $LOG sprintf("# % 38s : % 8d (%5.1f%%)\n", "too short",                       $count{filt_len}, 100*$count{filt_len}/$count{nSeqs});
print $LOG sprintf("# % 38s : % 8d (%5.1f%%)\n", "too many non-base characters",    $count{filt_N},   100*$count{filt_N}/$count{nSeqs});
print $LOG sprintf("# % 38s : % 8d (%5.1f%%)\n", "low entropy",                     $count{filt_H},   100*$count{filt_H}/$count{nSeqs});
print $LOG "#\n";
print $LOG "# Final sequences:\n";
print $LOG sprintf("# % 38s : % 8d (%5.1f%%)\n", "passed", $count{passed}, 100*$count{passed}/$count{nSeqs});

if (defined $options{L}) { close $LOG; }
exit;






### SUBS
sub adaptorAlignmentLength {
    my $str = $_[0];
    my ($a, $b) = split("-", (split("\t", $str, 2))[0]);
    if(defined($b)) {
	return ($b-$a+1);
    } else {
	return 0;
    }
}
    
sub count_seqs {
    my ($fname, $format) = @_;
    my $cnt = 0;
    if ($format eq 'fasta') {
	$cnt = `grep -c '^>' $fname`;
	chomp $cnt;

    } elsif ($format eq 'tab') {
	#faster: open(FILE, "<$fname"); $cnt += tr/\n/\n/ while sysread(FILE, $_, 2 ** 20); close FILE;
	$cnt = `wc -l $fname`;
	if ($cnt =~ m/^\s*(\d+)/) {
	    $cnt = $1;
	} else {
	    $cnt = 0;
	}

    } elsif ($format eq 'fastq') {
	$cnt = `grep -c '^\@' $fname`;
	chomp $cnt;
    }

    return $cnt;
}

sub store_seqs {
    my ($fname, $format, $seq) = @_;

    my ($currSeq, $currId, $currDesc);
    open(SEQIN, "<$options{i}") or die "ERROR opening $options{i}: $!\n";

    if ($format eq 'fasta') {
	while (my $line = <SEQIN>) {
	    chomp $line;
	    if ($line =~ m/^>/o) {
		if(defined $currId) {
		    $currSeq =~ tr/acgtbdefhijklmnopqrsuvwxyzBDEFHIJKLMNOPQRSUVWXYZ.,/ACGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/;
		    if (exists $seq->{$currSeq}) {
			&$doReport("$currId\tidentical to $seq->{$currSeq}[0]");
			$seq->{$currSeq}[1] += $currDesc;
		    } else {
			&$doReport("$currId\t$currDesc");
			$seq->{$currSeq} = [$currId, $currDesc, $currSeq, "\t", "\t"];
		    }
		}
		($currId, $currDesc) = split(/\s/, substr($line, 1), 3);
		if (!defined($currDesc)  or
		    $currDesc !~ m/^\d+$/) { $currDesc = 1; }
		$currSeq = "";
	    } else {
		$currSeq .= $line;
	    }
	}
	if (defined $currSeq and $currSeq ne "") {
	    $currSeq =~ tr/acgtbdefhijklmnopqrsuvwxyzBDEFHIJKLMNOPQRSUVWXYZ.,/ACGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/;
	    if (exists $seq->{$currSeq}) {
		&$doReport("$currId\tidentical to $seq->{$currSeq}[0]");
		$seq->{$currSeq}[1] += $currDesc;
	    } else {
		&$doReport("$currId\t$currDesc");
		$seq->{$currSeq} = [$currId, $currDesc, $currSeq, "\t", "\t"];
	    }
	}

    } elsif ($format eq 'tab') {
	while (my $line = <SEQIN>) {
	    chomp $line;

	    # tab can be in 4 formats:
	    # 1. seq
	    # 2. id seq
	    # 3. id seq counts
	    # 4. seq counts

	    if($line=~/^#/){next;}
	    my($currId, $currSeq, $currDesc);
	    my @l=split("\t",$line);

	    # determine the format
	    if(@l==1){
		($currId, $currSeq, $currDesc)=($l[0],$l[0],1);
	    }elsif(@l==2){
		if($l[1]=~/^\d+$/){
		    ($currId, $currSeq, $currDesc)=($l[0],$l[0],$l[1]);
		}else{
		    ($currId, $currSeq, $currDesc)=($l[0],$l[1],1);
		}
	    }elsif(@l==3){
		($currId, $currSeq, $currDesc) =($l[0],$l[1],$l[2]);
	    }
	    
	    #($currId, $currSeq, $currDesc) = split(/\t/, $line, 4);
	    #if (! defined $currDesc) { $currDesc = 1; }
	    $currSeq =~ tr/acgtbdefhijklmnopqrsuvwxyzBDEFHIJKLMNOPQRSUVWXYZ.,/ACGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/;
	    if (exists $seq->{$currSeq}) {
		&$doReport("$currId\tidentical to $seq->{$currSeq}[0]");
		$seq->{$currSeq}[1] += $currDesc;
	    } else {
		&$doReport("$currId\t$currDesc");
		$seq->{$currSeq} = [$currId, $currDesc, $currSeq, "\t", "\t"];
	    }
	}

    } elsif ($format eq 'fastq') {
	$currDesc = 1; # fastq format does not support description
	while (my $line = <SEQIN>) {
	    if ($line =~ m/^\@(\S+)/) {
		$currId = $1;
		$line = <SEQIN>;
		chomp $line;
		$currSeq = $line; # sequence in fastq is always on a single line
		$currSeq =~ tr/acgtbdefhijklmnopqrsuvwxyzBDEFHIJKLMNOPQRSUVWXYZ.,/ACGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/;
		if (exists $seq->{$currSeq}) {
		    &$doReport("$currId\tidentical to $seq->{$currSeq}[0]");
		    $seq->{$currSeq}[1] += $currDesc;
		} else {
		    &$doReport("$currId\t$currDesc");
		    $seq->{$currSeq} = [$currId, $currDesc, $currSeq, "\t", "\t"];
		}
	    }
	}
    }
    close SEQIN;
}

sub generate_variants {
    my ($str, $side, $neighborBase, $mm, $doindels) = @_;
    my %hash = ($str => 0);
    my %not  = (A => [qw/C G T N/], C => [qw/A G T N/], G => [qw/A C T N/], T => [qw/A C G N/]);
    my @all  = qw/A C G T N/;
    my $var;
    for(my $m=0; $m<$mm; $m++) {
	foreach my $s (keys %hash) {
	    for(my $p=0; $p<length($s); $p++) {
		# mismatch variants
		foreach my $other (@{$not{substr($s,$p,1)}}) {
		    $var = substr($s,0,$p).$other.substr($s,$p+1);
		    if (not exists $hash{$var}) { $hash{$var} = $m+1; }
		}
		if(defined($doindels) and
		   $p>1 and $p<length($s)-2) { # no flanking indels!
		    # deletion variant
		    $var = ( $side==5 ?
			     $neighborBase.substr($s,0,$p).substr($s,$p+1) :
			     substr($s,0,$p).substr($s,$p+1).$neighborBase);
		    if (not exists $hash{$var}) { $hash{$var} = $m+1; }
		    # insertion variants
		    foreach my $other (@all) {
			$var = substr($s,0,$p).$other.substr($s,$p,length($s)-$p-1);
			if (not exists $hash{$var}) { $hash{$var} = $m+1; }
		    }
		}
	    }
	}
    }
    return \%hash;
}

sub report {
    my $msg = shift;
    if (defined $LOG) { print $LOG "$msg\n"; }
}
