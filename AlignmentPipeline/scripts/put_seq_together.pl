#!/usr/bin/perl -w

use strict;

if(@ARGV != 3) {
    die "usage: put_seq_together.pl file outdir species_code\n";
}

my ($f, $d, $sp) = @ARGV;

open(F, $f) || die "cannot open $f\n";
my $namef = "";
my $nameo = "";
while(<F>) {
    my $line = $_;
    if($line =~ /\>ctc\-([^\|]+)\|([^\|]+)\|([^\|]+)\|([^\|]+)\|([^\|]+)/) {
	$namef = $1;
	$nameo = $sp . "_" . "$2" . "_" . "$4" . "_" . "$5" . "_" . "$3";
	if($namef =~ /^([^\-]+)\-([^\-]+)\-([^\-]+)\-(\S+)/) {
	    $namef = $1 . "_" . $2 . "_" . $3 . "_" . $4;
	}
	$line = ">$nameo\n";
    }
    else {
	$line =~ tr/[a-z]/[A-Z]/;
    }
    open(OF, ">>$d/$namef.fna") || die "cannot open $namef.fna\n";
    print OF $line;
    close(OF);
}
close(F);
