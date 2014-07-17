#!/import/bc2/soft/bin/perl5/perl -w

#directory with the fasta files
$d = shift(@ARGV);
$d =~ s/\/\s*$//;
#pattern for fasta files
$pat = shift(@ARGV);

opendir(D,$d);
@f = readdir(D);
closedir(D);
for($i=0;$i<@f;++$i){
    if($f[$i] =~ /$pat/){
	$infile = $d . "/" . $f[$i];
	$prefix = $f[$i];
	$prefix =~ s/\.fna\s*$//;
	$outfile = ">" . $d . "/" . $prefix . "\.dnd";
	open(F,$infile);
	#set all species as absent
	$pres{"hg18"} = "";
	$pres{"rheMac2"} = "";
	$pres{"canFam2"} = "";
	$pres{"bosTau2"} = "";
	$pres{"mm8"} = "";
	$pres{"monDom4"} = "";
        $pres{"equCab1"} = "";
	while(<F>){
	    if($_ =~ /\>([\S\d\_\-\+]+)\s*$/){
		$name = $1;
		@s = split(/\_/,$name);
		$spec = $s[0];
		print $spec, "\n";
		$pres{$spec} = $name;
	    }
	}
	close(F);
	if($pres{"hg18"} ne "" && $pres{"rheMac2"} ne ""){
            #	    die "no human present in input file $infile\n";
            $line = "((". $pres{"hg18"} .",". $pres{"rheMac2"} . ")";
	}
	elsif($pres{"rheMac2"} ne ""){
	    $line = "(" . $pres{"rheMac2"};
	}
        elsif ($pres{"hg18"} ne "") {
            $line = "(" . $pres{"hg18"};
        }
        else {
            $line = "(";
        }
	#dog, cow, horse
	if($pres{"canFam2"} ne "" && $pres{"bosTau2"} ne "" && $pres{"equCab1"} ne ""){
	    $line = "((" . $line . ",(" . $pres{"equCab1"} . "," . $pres{"bosTau2"} . ")),". $pres{"canFam2"} .")";
	}
        # dog & cow
        elsif ($pres{"canFam2"} ne "" && $pres{"bosTau2"}) {
          $line = "(". $line .",(". $pres{"canFam2"} .",". $pres{"bosTau2"} ."))";
        }
        # dog & horse
        elsif ($pres{"canFam2"} ne "" && $pres{"equCab1"}) {
          $line = "(". $line .",(". $pres{"canFam2"} .",". $pres{"equCab1"} ."))";
        }
        # cow & horse
        elsif ($pres{"bosTau2"} ne "" && $pres{"equCab1"}) {
          $line = "(". $line .",(". $pres{"bosTau2"} .",". $pres{"equCab1"} ."))";
        }
	#just dog
	elsif($pres{"canFam2"} ne ""){
	    $line = "(" . $line . "," . $pres{"canFam2"} . ")";
	}
	#just cow
	elsif($pres{"bosTau2"} ne ""){
	    $line = "(" . $line . "," . $pres{"bosTau2"} . ")";
	}
        #just horse
        elsif($pres{"equCab1"} ne ""){
            $line = "(" . $line . "," . $pres{"equCab1"} . ")";
          }
	#mouse
	if($pres{"mm8"} ne ""){
	    $line = "(" . $line . "," . $pres{"mm8"} . ")";
	}
	#opposum
	if($pres{"monDom4"} ne ""){
	    $line = "(" . $line . "," . $pres{"monDom4"} . ")";
	}
	#final bracket
	$line =~ s/^\(//;
        $line =~ s/\(,/\(/g; # fix "(," if any
	$line = $line . "\;";
	open(G,$outfile);
	print G $line, "\n";
	close(G);
    }
}


