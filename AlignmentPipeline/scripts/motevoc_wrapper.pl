#!/import/bc2/home/nimwegen/pachko/local/bin/perl

use lib '/import/bc2/home/nimwegen/pachko/local/lib/';

use strict;
use warnings;
use File::Path;
use File::Temp;
use File::Copy;
use Sys::Hostname;
use Getopt::Long;
use Fcntl qw(:flock SEEK_END);

my $alignments_file = '../final_out.aln';
my $wm_dir          = '/import/bc2/home/nimwegen/GROUP/WMs/Mammals/CurrentWMs'; #'/import/bc2/home/nimwegen/GROUP/WMs/p1/P1s';
my $working_dir     = 'ResMotEvo';
my $motevo_path     = '/import/bc2/home/nimwegen/GROUP/software/motevo_ver1.02/bin/motevo'; # '/import/bc2/home/nimwegen/GROUP/MotEvoC/motevo';
my $DEBUG           = 1;

Getopt::Long::GetOptions('alignments=s'  => \$alignments_file,
                         'wm_dir=s'      => \$wm_dir,
                         'working_dir=s' => \$working_dir,
                         'motevo_path=s' => \$motevo_path,
                         'debug=i'       => \$DEBUG
                        );

my $taskid           = $ENV{"SGE_TASK_ID"};
my $hostname         = Sys::Hostname::hostname();

# wait a bit before start reading a file
my $delay            = ($taskid -1) % 50;
sleep $delay;

#read which WM we are talking about

opendir(my $wm_dir_h, $wm_dir) ||die "Can not open dir $wm_dir.\n";
my @wms =  sort grep {$_ !~ /^\.+/} readdir($wm_dir_h);
if (($taskid - 1) > scalar(@wms)) {
    print STDERR "No wm for task id $taskid\n";
    exit 0;
}
my $wm_name = $wms[$taskid -1];

closedir($wm_dir_h);

chomp($wm_name);
my $wm_file = "$wm_dir/$wm_name";

# count length of WM

my $wm_length = get_WM_length($wm_file);

# get wm id
my $wm_id = $wm_name;

if (-e "$working_dir/sites_$wm_id") {
    print STDERR "$hostname : $working_dir/sites_$wm_id file exists\n";
    exit 0;
}

# create scratch directory
my $username    = get_username();
my $motevo_scratch_dir = "/scratch/$username/motevo";
File::Path::make_path($motevo_scratch_dir) unless  (-e $motevo_scratch_dir);

my $scratch_dir = File::Temp::tempdir( 'ResMotEvo_XXXXX', DIR => $motevo_scratch_dir );
chmod 0777, $scratch_dir;

# make motevoc_param file for this matrix
my $paramfile = "params_prom_$wm_id";

open(F,"$working_dir/motevoc_params") || die "Can not open file $working_dir/motevoc_params\n";
open(G,"> $scratch_dir/$paramfile") || die "Can not open file $scratch_dir/$paramfile\n";

while(<F>){
    if(/sitefile/){
	print G "sitefile $scratch_dir/sites_$wm_id\n";
    }
    elsif(/priorfile/){
	print G "priorfile $scratch_dir/priors_$wm_id\n";
    }
    elsif (/otherwmlen/) {
        print G "otherwmlen $wm_length\n";
    }
    elsif (/UFEwmlen/) {
        print G "UFEwmlen $wm_length\n
";
    }
    else{
	print G $_;
    }
}
close(F);

close(G);


#run motevo
my $cmd = "$motevo_path \"$alignments_file\" \"$scratch_dir/$paramfile\" \"$wm_file\"";

# fix "bad" characters
#$cmd =~ s/([\[\]\{\}\(\)])/\\$1/g;
print STDERR "$hostname : $cmd\n" if $DEBUG;
system($cmd);

#remove all otherwm sites
my $sitefile  = "sites_" . $wm_id;
my $priorfile = "priors_" . $wm_id;
my $outfile   = "tmp_" . $wm_id;

print STDERR "sitefile $sitefile\n" if $DEBUG;

open(F,"$scratch_dir/$sitefile") || die "$hostname : Can not open file $scratch_dir/$sitefile\n";
open(G,"> $scratch_dir/$outfile") || die "$hostname : Can not open file $scratch_dir/$outfile\n";

my $insite = 0;
while(<F>){
    if($_ =~ /^\d+/){
	if(/otherwm/ || /UFEwm/){
	    $insite = 0;
	}
	else{
	    print G $_;
	    $insite = 1;
	}
    }
    elsif($insite){
	print G $_;
    }
}
close(F);
close(G);

# move results to working directory
move("$scratch_dir/$outfile", "$working_dir/$sitefile");
move("$scratch_dir/$priorfile", "$working_dir/$priorfile");

# remove scratch dir

File::Path::remove_tree($scratch_dir);

# echo to log
open(my $fout, ">>task_log") || die "motevo_wrapper: Can not open file task_log\n";
flock($fout, LOCK_EX) or die "Cannot lock task_log - $!\n";
# and, in case someone appended while we were waiting...
seek($fout, 0, SEEK_END) or die "Cannot seek - $!\n";
print $fout "finished task $taskid\n";
flock($fout, LOCK_UN) or die "Cannot unlock  - $!\n";
close($fout);

# get username

sub get_username {
    return id("u");
}

sub  id {
    my $type = shift;
    open(IN, "id -$type --name |") || die "Can not open pipe to id\n";
    my $out = <IN>;
    chomp $out;
    close(IN);

    return $out;
}

sub get_WM_length {

    my $wm_file = shift;

    open(WM, $wm_file) || die "Can not open file $wm_file\n";
    
    my @wm = <WM>;
    
    my @positions = grep {/\A \d+ /x} @wm;
    
    close(WM);

    return scalar(@positions);
}
