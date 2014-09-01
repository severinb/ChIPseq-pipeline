#!/import/bc2/soft/bin/perl5/perl -w

use strict;
use warnings;

my $dir = shift || "./";
my $out_file = shift || "final_out.aln";

# read contents of the alignment dir
opendir(DIR,$dir) || die "Can not open dir $dir\n";
my @list = grep{/\.aln$/} readdir(DIR);
closedir(DIR);

open(OUT,"> $out_file") || die "Can not open file $out_file\n";

for my $file (@list) {
    next if $file eq $out_file;
    open(IN,"$dir/$file") || die "Can not open file $dir/$file\n";
#    print STDERR "Opened file $dir/$file for reading\n";
    print OUT ">";
    while (<IN>) {
        my $string = $_;
        next if /\A \s*\n/x;
        print OUT $string;
    }
    close(IN);
#    print STDERR "Wrote to final_out and closed file $dir/$file\n";
}

close(OUT);
               
