use warnings;
use strict;


my $cut_bg = $ARGV[0];     # cut-off on background 
my $cut_score = $ARGV[1];  # cut-off on z-value/posterior
my $DIR = $ARGV[2]; # output directory
#my $in_dir = $ARGV[3]; #input idrectory (output directory from fit_sigma_mu_rho.pl)

my $string = 'awk \'{ if($6<='.$cut_bg.' && $7>'.$cut_score.') { print $0 } }\' outzvals | sort -g -k 7 -r > final_windows';

system($string);

open W, ">".$DIR."cut_offs_used";
print W join("\t",'bg cut-off:',$cut_bg)."\n";
print W join("\t",'score cut-off:',$cut_score)."\n";
close W;

