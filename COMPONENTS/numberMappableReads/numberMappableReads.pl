#!/usr/bin/perl

use strict;

my $sample = shift @ARGV;
my $out_file = shift @ARGV;
my $anno = shift @ARGV;
my $perlPATH = shift @ARGV;
my $rootD = shift @ARGV;

# necessery to exec the RIGHT perl in the system call
$ENV{PATH} = 'perlPATH:' . $ENV{PATH};

open(my $iH, "$rootD/soft/extractData.pl -f $sample $anno genome | $rootD/soft.bc2/frag2totalGenomic.pl - | ");
open(my $oH, ">", $out_file) or die;
while(<$iH>)
{
	print $oH $_;
}
close $oH;
close $iH;
