#!/import/bc2/home/nimwegen/pachko/local/bin/perl

use lib '/import/bc2/home/nimwegen/pachko/local/lib/';

use strict;
use warnings;
use YAML::Any;

# read private variables

my $config_file = '/import/bc2/home/nimwegen/GROUP/Pipeline2/global_pipeline.conf';
my $config      = YAML::Any::LoadFile($config_file);
my %REFSEQ = %{$config->{REFSEQ}};

# read user variables
$config_file = shift || "pipeline.conf";
$config      = YAML::Any::LoadFile($config_file);

my $root_organism   = $config->{ROOT_ORGANISM};
my $refseq_file     = $config->{REFSEQ_FILE};
my $regions_file    = $config->{REGIONS_FILE};
my $RANGE           = $config->{RANGE};
my $ALLOWED_OVERLAP = 0;


# check if we have data for specified root organism

if (!(grep {/^$root_organism$/} keys %REFSEQ)) {
    die "No data available for root organism $root_organism\n";
}

# read chromosome lengths

open(IN, '/import/bc2/data/databases/UCSC/hg18/hg18chrs.tab') || die "Can not open file /import/bc2/data/databases/UCSC/hg18/hg18chrs.tab\n";

my %chr_length;

while (<IN>) {
    my @data = split /\s+/, $_;
    $chr_length{$data[0]} = $data[1];
}

close(IN);


open(IN,$refseq_file) || die "Can not open file $refseq_file\n";
# refseqs should be given one per line
# make array with refseq id
my %refseqs;
while (<IN>) {
    chomp;
    if (/^\s*(\S+)/) {
        $refseqs{$1} = 0;
    }
}
close(IN);

#get coordinates


my %regions;
open(IN,$REFSEQ{$root_organism}) || die "Can not open file ". $REFSEQ{$root_organism} ."\n";

# file to write new id with coordinates
open(OUT,"> refseqs_id_with_coordinates") || die "Can not open file refseqs_id_with_coordinates\n";

while (<IN>) {
    chomp;
    if (/(\S+)\s+(\S+)\s+mRNA\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+ID=(\S+)/) {
        my $chr = $1;
        my ($start,$end) = sort {$a <=> $b} ($3,$4);
        my $strand = $6;
        my $refseq_id = $8;
        if (defined($refseqs{$refseq_id})) {
            print OUT $refseq_id .".". $refseqs{$refseq_id} ."\t". join("\t",
                                                                        $chr,
                                                                        $start,
                                                                        $end,
                                                                        $strand
                                                                       ) ."\n";
            my ($region_start, $region_end);
            if ($strand eq "+") {
                $region_start = $start - $RANGE->[1];
                $region_end = $start - $RANGE->[0];
            }
            elsif ($strand eq "-") {
                $region_start = $end + $RANGE->[0];
                $region_end = $end + $RANGE->[1];
            }
            else {
                die "Wrong strand for $refseq_id\n";
            }

            if (defined($regions{$chr}{$region_start})) {
                # check if strand the same
                if ($strand ne $regions{$chr}{$region_start}{'strand'}) {
                    die "Same start but different strand for refseq_id and ". $regions{$chr}{$region_start}{'id'} ."\n";
                }
                $regions{$chr}{$region_start}{'id'} .= $refseq_id .".". $refseqs{$refseq_id} ." ";
            }
            else {
                $regions{$chr}{$region_start}{'id'} = $refseq_id .".". $refseqs{$refseq_id} ." ";
                $regions{$chr}{$region_start}{'end'} = $region_end;
                $regions{$chr}{$region_start}{'strand'} = $strand;
            }

            ++$refseqs{$refseq_id};
        }
    }
}
close(OUT);

# check if all refseq ids were processed

for my $refseq_id (keys %refseqs) {
    if (!$refseqs{$refseq_id}) {
        print STDERR "No data found for refseq $refseq_id\n";
    }
}

# check for overlaps and merge overlaping regions

my %final_regions;

for my $chr (sort keys %regions) {
    for my $start (sort {$a <=> $b} keys %{$regions{$chr}}) {
        if (!defined($final_regions{$chr})) {
            $final_regions{$chr}{$start}{'id'} = $regions{$chr}{$start}{'id'};
            $final_regions{$chr}{$start}{'end'} = $regions{$chr}{$start}{'end'};
            $final_regions{$chr}{$start}{'strand'} = $regions{$chr}{$start}{'strand'};
        }
        else {
            my $overlap_flag = 0;
            for my $final_start (sort {$a <=> $b} keys %{$final_regions{$chr}}) {
                # skip if strand mismatch
                
                if ($final_regions{$chr}{$final_start}{'strand'} ne $regions{$chr}{$start}{'strand'}) {
                    next;
                }
                
                # check for overlap

                if ($final_regions{$chr}{$final_start}{end} >= ($start + $ALLOWED_OVERLAP)) {
                    if ($final_regions{$chr}{$final_start}{end} <= $regions{$chr}{$start}{'end'}) {
                        print STDERR "region $chr $start overlaps with final_region $chr $final_start. merging...\n";
                        $final_regions{$chr}{$final_start}{'end'} = $regions{$chr}{$start}{'end'};
                        $final_regions{$chr}{$final_start}{'id'} .= $regions{$chr}{$start}{'id'};
                        ++$overlap_flag;
                    }
                    else {
                        print STDERR "region $chr $start is inside final_region $chr $final_start. merging...\n";
                        $final_regions{$chr}{$final_start}{'id'} .= $regions{$chr}{$start}{'id'};
                        ++$overlap_flag;
                    }
                }
            }
            if (!$overlap_flag) {
                $final_regions{$chr}{$start}{'id'} = $regions{$chr}{$start}{'id'};
                $final_regions{$chr}{$start}{'end'} = $regions{$chr}{$start}{'end'};
                $final_regions{$chr}{$start}{'strand'} = $regions{$chr}{$start}{'strand'};
            }
        }
    }
}

# print regions

open(OUT,"> $regions_file") || die "Can not open file regions_file\n";

for my $chr (sort keys %final_regions) {
    for my $start (sort {$a <=> $b} keys %{$final_regions{$chr}}) {
        if (($final_regions{$chr}{$start}{'end'} - $start) > 2000) {
            print STDERR "Region is longer than 2000.Skipping...\t". join("\t",
                                                                          $chr,
                                                                          $start,
                                                                          $final_regions{$chr}{$start}{'end'},
                                                                          $final_regions{$chr}{$start}{'strand'},
                                                                          $final_regions{$chr}{$start}{'id'}
                                                                         ) ."\n";
            next;
        }
        print OUT join("\t",
                   $chr,
                   $start,
                   $final_regions{$chr}{$start}{'end'},
                   $final_regions{$chr}{$start}{'strand'},
                   $final_regions{$chr}{$start}{'id'}
                  ) ."\n";
    }
}

close(OUT);
