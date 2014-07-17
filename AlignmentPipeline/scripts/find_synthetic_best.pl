#!/usr/bin/perl -w

use strict;

if (@ARGV != 1) {
  die "usage: find_syntenic_best.pl file\n";
}

my ($f) = @ARGV;
my %mapping;
my %regions;
open(my $fin, '<', $f) || die "cannot open $f\n";
while (<$fin>) {
  if ($_ =~ /^(ctc\S+)/) {
    my $name = $1;
    $_ = <$fin>;
    my($hchr, $hal, $hbeg, $hend, $hstrand) = split(/\s+/, $_);
    $_ = <$fin>;
    my ($chr, $al, $beg, $end, $strand) = split(/\s+/, $_);
    if ($hal =~ /^-+$/ || $al =~ /^-+$/) {
        print STDERR "hal $hal al $al Skipping...\n";
#        next;
    }
    push @{$regions{$name}}, "$chr $beg $end $strand";
    if (!defined($mapping{$chr}{$strand}{"$beg:$end"})) {
      $mapping{$chr}{$strand}{"$beg:$end"} = "$hchr $hbeg $hend $hstrand";
    } elsif ($mapping{$chr}{$strand}{"$beg:$end"} eq "$hchr $hbeg $hend $hstrand") {
      print STDERR "Warning! Segment \[$chr $beg $end $strand\] to ".
	"[$hchr $hbeg $hend $hstrand] met twice!\n";
    } else {
      print STDERR "Warning! Segment \[$chr $beg $end $strand\] aligned to ".
	$mapping{$chr}{$strand}{"$beg:$end"} . " and to $hchr $hbeg $hend $hstrand." .
	  "Leaving only the first alignment.\n";
    }
  }
}
close($fin);

foreach my $name (keys(%regions)) {
    print STDERR "region name $name\n";
  my ($hchr, $hbeg, $hend, $hstr);
  if ($name =~ /^ctc\-([^\-]+)\-([^\-]+)\-([^\-]+)\-(\S+)/) {
    ($hchr, $hbeg, $hend, $hstr) = ($1, $2, $3, $4);
  }
  my $len = ($hend-$hbeg+1);
  my %chrblocks = split_by_TU(@{$regions{$name}});
  my ($absbest, $abschr, $absbeg, $absend, $absstr);
  $absbest = 0;
  foreach my $chr (keys (%chrblocks)) {
    my @v = sort {sort_by_begin($a, $b)} @{$chrblocks{$chr}};
    print STDERR "chr $chr" .join("|",@v) ."\n";
    my $best= 0;
    my ($first, $last);
    for (my $i = 0; $i < @v; $i++) {
      for (my $j = $i; $j < @v; $j++) {
	my ($span, $coverage) = compute_span_coverage(@v[$i..$j]);
	print STDERR "span $span coverage $coverage len $len\n";
	if ($span <= 1.5 * $len || scalar(@v) == 1) {
	  if ($coverage > $best) {
	    $best = $coverage;
	    $first = $i;
	    $last = $j;
	  }
	}
      }
    }
    if ($best > $absbest) {
      print STDERR "asigning new best values\n";
      ($abschr, $absstr) = split(":", $chr);
      ($absbeg, $absend) = pick_boundaries($v[$first], $v[$last]);
      $absbest = $best;
    }
  }

  # get start and end of corresponding blocks in human genome
  print STDERR "abschr $abschr absstr $absstr\n";
  my ($hbeg1,$hend1,$hbeg2,$hend2)= (0,0,0,0);
  for my $key (keys %{$mapping{$abschr}{$absstr}}) {
    my ($start,$end) = split /:/, $key;
    my @data = split /\s+/,$mapping{$abschr}{$absstr}{$key};
    if ($start ==  $absbeg || $end == $absbeg) {
      ($hbeg1,$hend1) = ($data[1],$data[2]);
    }
    if ($start ==  $absend || $end == $absend) {
      ($hbeg2,$hend2) = ($data[1],$data[2]);
    }
    if ($hbeg1 && $hend2) {
      last;
    }
  }
  my @tmp = sort {$a <=> $b} ($hbeg1,$hbeg2,$hend1,$hend2);
  $hbeg = $tmp[0];
  $hend = $tmp[3];
  #find necessary extentions from block edges to edges of the segment

  $name  =~ /chr\S+-(\d+)-(\d+)-([\+-])/;
  my $segment_start = $1 < $2 ? $1 : $2;
  my $segment_end = $1 < $2 ? $2 : $1;
  my $segment_str = $3;
  if ($segment_start > $hbeg) {
    print STDERR "Multiz segment ($hbeg) starting before segment of interest ($segment_start). Skipping...\n";
    next;
  }
  if ($segment_end < $hend) {
    print STDERR "Multiz segment ($hend) ending after segment of interest ($segment_end). Skipping...\n";
    next;
  }
  my $extension_start = $hbeg - $segment_start;
  my $extension_end = $segment_end - $hend;
  print STDERR "human original start $segment_start\tend $segment_end\n";
  print STDERR "human myltiZ    start $hbeg\tend $hend\n";
  print STDERR "extension_start $extension_start extension_end $extension_end\n";
  # recalculate edge points for organism
  print STDERR "organism multiz start $absbeg\tend $absend\n";
  if ($segment_str eq $absstr) {
    $absbeg -= $extension_start;
    $absend += $extension_end;
  } else {
    $absbeg -= $extension_end;
    $absend += $extension_start;
  }
  if ($absbeg < 0) {
    $absbeg = 0;
  }
  print STDERR "organism extended start $absbeg\tend $absend\tstrand $absstr\n";
  # better make for end of chromosome too but what is the best way?

  # print out final values

  if ($absbest > 0) {
    print STDOUT "$name\t$abschr\t$absstr\t$absbeg\t$absend\n";
  }
}

sub split_by_TU {
  my @s = @_;
  my %r;
  for (my $i = 0; $i < @s; $i++) {
    my @ss = split(/\s+/, $s[$i]);
    push @{$r{"$ss[0]:$ss[3]"}}, $s[$i];
  }
  return %r;
}

sub sort_by_begin {
  my ($v1, $v2) = @_;
  my ($b1, $e1, $b2, $e2);
  if ($v1 =~ /\s(\d+)\s+(\d+)/) {
    ($b1, $e1) = ($1, $2);
  } else {
    die $v1;
  }
  if ($v2 =~ /\s(\d+)\s+(\d+)/) {
    ($b2, $e2) = ($1, $2);
  } else {
    die $v2;
  }
  return $b1 <=> $b2;
}

sub compute_span_coverage {
  my @vv = @_;
  my ($v1, $v2) = ($vv[0], $vv[$#vv]);
  my ($b1, $e1, $b2, $e2);
  my $totalcoverage = 0;
  my %b;
  for (my $i = 0; $i < @vv; $i++) {
    if ($vv[$i] =~ /\s(\d+)\s+(\d+)/) {
      my ($bb, $ee) = ($1, $2);
      for (my $j = $bb; $j <= $ee; $j++) {
	$b{$j} = 1;
      }
    } else {
      die $vv[$i];
    }
  }
  $totalcoverage = scalar(keys(%b));
  if ($v1 =~ /\s(\d+)\s+(\d+)/) {
    ($b1, $e1) = ($1, $2);
  } else {
    die $v1;
  }
  if ($v2 =~ /\s(\d+)\s+(\d+)/) {
    ($b2, $e2) = ($1, $2);
  } else {
    die $v2;
  }
  my @c = sort {$a<=>$b} ($b1, $e1, $b2, $e2);
  return ($c[$#c]-$c[0]+1, $totalcoverage);
}

sub pick_boundaries {
  my ($v1, $v2) = @_;
  my ($b1, $e1, $b2, $e2);
  if ($v1 =~ /\s(\d+)\s+(\d+)/) {
    ($b1, $e1) = ($1, $2);
  } else {
    die $v1;
  }
  if ($v2 =~ /\s(\d+)\s+(\d+)/) {
    ($b2, $e2) = ($1, $2);
  } else {
    die $v2;
  }
  my @c = sort {$a<=>$b} ($b1, $e1, $b2, $e2);
  return ($c[0], $c[$#c]);
}
