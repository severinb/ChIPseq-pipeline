#!/import/bc2/soft/bin/perl5/perl -w

use strict;
use warnings;
use Fcntl qw(:flock SEEK_END);

############################################################
#  programm starts the iprquery programm on the clusternode#
############################################################

my $taskid         = $ENV{"SGE_TASK_ID"};
my $method         = shift || "tcoffe"; #clustal or tcoffee
my $filedir        = shift || "OUTPUT"; #put directory with all files
my $outdir         = shift || "$method\_$filedir";
my $numfilesperrun = shift || 10;

my $filelist = "$filedir/filelist";

open(F, $filelist) || die "alignment_wrapperalignment_wrapper: Can not open file $filelist\n";
my @names = <F>;
close(F);
my $numlines = @names;
--$taskid;
$taskid *= $numfilesperrun;
my $end = $taskid + $numfilesperrun;
if ($end > $numlines) {
    $end = $numlines;
}

mkdir "$filedir/log" if(!-d "$filedir/log");

for (my $i=$taskid;$i<$end;++$i) {
    my $infile = $names[$i];
    $infile =~ s/\s*$//g;
    my $logfile = "$filedir/log/$infile.log";
    my $treefile = $infile;
    $treefile =~ s/\.fna/\.dnd/;
    my $outfile = $infile;
    $outfile =~ s/\.fna/\.tmp/;
    $infile = $filedir . "/" . $infile;
    $treefile = $filedir . "/" . $treefile;
    $outfile = $outdir . "/" . $outfile;
    # check if outfile is present
    # if so, then skip this file
    my $alignment_file = $outfile;
    $alignment_file =~ s/tmp$/aln/;
    logger($logfile, "$i:$infile:". $ENV{'HOSTNAME'} ." started\n\n");
    if (open(TMP,$alignment_file)) {
        close(TMP);
        print STDERR "Alignment file $alignment_file already exists. Skipping.\n";
        next;
    }
    if ($method eq "clustal") {
	system(join(" ",
                    "/import/bc2/soft/bin/clustalw",
                    "-infile=$infile",
#                    "-usetree=$treefile",
                    "-type=dna",
                    "-output=fasta",
                    "-outfile=$outfile",
                    "-case=UPPER",
                    "> /dev/null"
                   )
              );
    } elsif ($method eq "tcoffee") {
        my $cmd = join(' ',
                       'TMP_4_TCOFFEE=/scratch/'. $ENV{'USER'} .'/',
                       '/import/wnz/home/crunch/soft/T-COFFEE_distribution_Version_9.03.r1318/bin/t_coffee',
                       "$infile",
#                       "-usetree=$treefile",
                       '-type=dna',
                       "-outfile=$outfile",
                       '-output fasta_aln',
                       '-n_core=1',
                       ">> $logfile 2>&1"
                      );
        system($cmd);
    }
    logger($logfile, "\nfinished\n");
    #now fix the names in the outfile
    #read original names
    my (%original_names,@original_names_list);
    open(F,$infile) || die "alignment_wrapper: Can not open file $infile\n";
    while (<F>) {
	chomp($_);
	if (/^>(\S+)/) {
	    my $curname = $1;
	    my @data = split /_/, $curname;
	    $original_names{$data[0]} = $curname;
            push @original_names_list, $curname;
        }
    }
    close(F);

    my $finalout = $outfile;
    $finalout =~ s/\.tmp/\.aln/;
    $finalout = ">" . $finalout;

    # read names in outfile and put original names in final out
    open(F,$outfile) || die "alignment_wrapperalignment_wrapper: Can not open file $outfile\n";
    open(G,$finalout) || die "alignment_wrapperalignment_wrapper: Can not open file $finalout\n";
    my ($cur_name,%seqs);
    while (<F>) {
        chomp;
	if ($_ =~ /^>(\S+)/) {
	    my @data = split /_/, $1;
	    $cur_name = $original_names{$data[0]};
	    $seqs{$cur_name} = "";
	} else {
            my $string = $_;
            $string =~ tr/actgnx/ACTGNX/;
	    $seqs{$cur_name} .= $string;
	}
    }
    close(F);

    # find if first seq starts and ends with gaps
    my($gap_start_len,$gap_end_len);
    
    if ($seqs{$original_names_list[0]} =~ /^(-+)[actgnxACTGNX]/) {
        $gap_start_len = length($1);
    } else {
        $gap_start_len = 0;
    }
    ;
    if ($seqs{$original_names_list[0]} =~ /[actgnxACTGNX](-+)$/) {
        $gap_end_len = length($1);
    } else {
        $gap_end_len = 0;
    }

    for my $name (@original_names_list) {
        print G ">" , $name, "\n";
        my $sequence = substr $seqs{$name},$gap_start_len;
        $sequence = substr $sequence, 0, (length($sequence) - $gap_end_len);
        print G "$sequence\n";
    }
    close(G);
    system("rm $outfile");

    # remove .dnd file
    $outfile =~ s/\.tmp$/\.dnd/;
    my @data = split /\//, $outfile;
    $outfile = $data[-1];
    if (-e $outfile) {
        system("rm $outfile") == 0 || print STDERR "Failed to remove dnd file $outfile\n";
    }
}

sub logger {
    my $file = shift;
    my $message = shift;
    open(my $fout, ">>$file") || die "alignment_wrapper: Can not open file $file\n";
    flock($fout, LOCK_EX) or die "Cannot lock $file - $!\n";
    # and, in case someone appended while we were waiting...
    seek($fout, 0, SEEK_END) or die "Cannot seek - $!\n";
    print $fout $message;
    flock($fout, LOCK_UN) or die "Cannot unlock  - $!\n";
    close($fout);
}
