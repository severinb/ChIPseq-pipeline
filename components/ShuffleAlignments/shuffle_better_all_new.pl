use warnings;
use strict;


## species
#my @specs = ('hg19','mm9','rheMac2','canFam2','equCab1','bosTau6','monDom4');
#my @nonrefs = ('mm9','rheMac2','canFam2','equCab1','bosTau6','monDom4');
#my $REF = 'mm9'; # = #specs - @nonsrefs
#my $REF = $ARGV[1]; # reference species

my @specs = ();
my @nonrefs = ();

my $REF;
for(my$i=1;$i<@ARGV;$i++) {

    if($i==1) {
	$REF = $ARGV[$i];
	push(@specs,$ARGV[$i]);
    }
    else {
	push(@specs,$ARGV[$i]);
	push(@nonrefs,$ARGV[$i]);
    }
}

#print STDERR "\n@specs\n";
#my $kk = scalar(@specs);
#print STDERR "$kk\n";
#exit;

## prepare pool with sequences, uses a shitload of memory
my %seqs;

print STDERR "Loading multi alignment ...";
open F, $ARGV[0]; # file with multi alignments
my @ma = <F>;
close F;

my %load_2_id;

my %patterns;
my $countload = 0;
for(my $i=0;$i<@ma;$i++) {

    if( $ma[$i] =~ /\>\>/ ) { # read multi alignments
	
	chomp($ma[$i]);
	$countload++;

	#print STDERR "==============new multi alignment: $ma[$i] ==== ( $countload ) ===============\n\n";

	# get id of region: >>mm9_chr5_24604500_24608500_+
	my @s = split(/\_/,$ma[$i]);
	my $id = join('_',$s[1],$s[2],$s[3]);

	my %tmpseq;

	$i++;  chomp($ma[$i]); $tmpseq{$REF} = $ma[$i];
	my $len = length($ma[$i]);

	for my $spec (@nonrefs) { $tmpseq{$spec}=''; }

	for(my$pos=0;$pos<$len;$pos++) {
	    for my $spec (@nonrefs) { $tmpseq{$spec} .= '-'; }
	}

	$i++;
	while( $ma[$i] !~ /\>\>/ && $i<(@ma-2)) {

	    if( $ma[$i] =~ /\>/ ) {
		
		chomp($ma[$i]);
		$ma[$i] =~ s/\>//g;
		my @s = split(/\_/,$ma[$i]);
	    
		$i++; chomp($ma[$i]);
		$tmpseq{$s[0]} = $ma[$i];
	    }
	    $i++;
	}

	#### idea #####
        #      %seqs{id}{spec}{pos} = char
        #      %pool{sel} = @(1,5,77,99, ...)
        ##############

	# get selection pattern, same nucl as reference, 
	for(my$pos=0;$pos<$len;$pos++) {

	    # get selection pattern
	    my $sel = '';
	    for my $spec (@specs) {
		my $char = substr($tmpseq{$spec},$pos,1);
		if( $char eq '-' ) { $sel .= '0'; }
		else { $sel .= '1'; }
	    }
#	    print "\t".$sel."\n";

	    # same nucl as reference
	    my $refchar = '';
	    my $selbase = '';
	    for my $spec (@specs) {

		if( $spec eq $REF ) {
		    
		    $refchar = substr($tmpseq{$spec},$pos,1);
		    next;
		}

		my $char = substr($tmpseq{$spec},$pos,1);

		if( $char eq $refchar ) { $selbase .= '1'; }
		else { $selbase .= '0'; }

	    }
#	    print "\t ".$selbase."\n";

	    # all bases
	    my $bases = ''; # A=0, C=1, G=2, T=3, -=4, N=5
	    for my $spec (@specs) {

		my $char = substr($tmpseq{$spec},$pos,1);

		if( $char eq 'A' || $char eq 'a' ) { $bases .= '0'; }
		if( $char eq 'C' || $char eq 'c' ) { $bases .= '1'; }
		if( $char eq 'G' || $char eq 'g' ) { $bases .= '2'; }
		if( $char eq 'T' || $char eq 't' ) { $bases .= '3'; }
		if( $char eq '-' ) { $bases .= '4'; }
		if( $char eq 'N' ) { $bases .= '5'; }
            }
	    if( length($bases)<scalar(@specs) ) { 
		print STDERR "Something is wrong:\t";
		print STDERR "$countload, $pos\n";
		exit;
	    }
#	    print "\t".$bases."\n";
#	    print join("\t",$pos,$sel,$selbase,$bases)."\n";
#	    exit;
#	    if( $pos>8 ) { exit; }

	    if( exists( $patterns{$sel}{$selbase}{$bases} ) ) { $patterns{$sel}{$selbase}{$bases}++; }
	    else { $patterns{$sel}{$selbase}{$bases}=1; }


	    

	}
	
	$i--;
    }


#    if( $countload > 10 ) { last; }

}
############################
# create the dummy regions #
############################
@ma = ();
print STDERR "Loading multi alignment ...";
open F, $ARGV[0]; # file with multi alignments
@ma = <F>;
close F;


$countload = 0;
for(my $i=0;$i<@ma;$i++) {

    if( $ma[$i] =~ /\>\>/ ) { # read multi alignments

	chomp($ma[$i]);
	$countload++;

	my $reffield = $ma[$i];
	$reffield =~ s/\>\>//g;

	#print STDERR "==============new multi alignment: $ma[$i] ==== ( $countload ) ===============\n\n";

	# get id of region: >>mm9_chr5_24604500_24608500_+
	my @s = split(/\_/,$ma[$i]);
	my $id = join('_',$s[1],$s[2],$s[3]);

	my %tmpseq;

	$i++;  chomp($ma[$i]); $tmpseq{$REF} = $ma[$i];
	my $len = length($ma[$i]);

	for my $spec (@nonrefs) { $tmpseq{$spec}=''; }

	for(my$pos=0;$pos<$len;$pos++) {
	    for my $spec (@nonrefs) { $tmpseq{$spec} .= '-'; }
	}

	$i++;
	while( $ma[$i] !~ /\>\>/ && $i<(@ma-2)) {

	    if( $ma[$i] =~ /\>/ ) {
		
		chomp($ma[$i]);
		$ma[$i] =~ s/\>//g;
		my @s = split(/\_/,$ma[$i]);
	    
		$i++; chomp($ma[$i]);
		$tmpseq{$s[0]} = $ma[$i];
	    }
	    $i++;
	}

#	my $pos = 3;
#	my $spec ='monDom4';
#	my $char = substr($tmpseq{$spec},$pos,1);
#	print join("\t",$char,$pos,$spec)."\n";
#	exit;

	my %dummyseq;

	# get selection pattern, same nucl as reference, 
	for(my$pos=0;$pos<$len;$pos++) {

	    # get selection pattern
	    my $sel = '';
	    for my $spec (@specs) {
		my $char = substr($tmpseq{$spec},$pos,1);
		if( $char eq '-' ) { $sel .= '0'; }
		else { $sel .= '1'; }
	    }
#	    print "\t".$sel."\n";

	    # same nucl as reference
	    my $refchar = '';
	    my $selbase = '';
	    for my $spec (@specs) {

		if( $spec eq $REF ) {
		    
		    $refchar = substr($tmpseq{$spec},$pos,1);
		    next;
		}

		my $char = substr($tmpseq{$spec},$pos,1);

		if( $char eq $refchar ) { $selbase .= '1'; }
		else { $selbase .= '0'; }

	    }
#	    print "\t ".$selbase."\n";

	    if( length($sel)!=scalar(@specs) ) { 
#                my $nwk=length(@specs);
#                print STDERR "$nwk\n";
#                print STDERR "$sel\n";
		print STDERR "Something is wrong with length of selection\n";
		exit;
	    }
	    if( length($selbase)!=scalar(@specs)-1 ) { 
                print STDERR "\n$selbase\n";
		print STDERR "Something is wrong with length of selection (base)\n";
		exit;
	    }

	    if( $pos== -1 ) { 
		print "pos = $pos\n"; 
		print "sel = $sel\n";
		print "selbase = $selbase\n";
		
		my $char = substr($tmpseq{$REF},$pos,1);
		print $char."\n";
		$char = substr($tmpseq{$nonrefs[-1]},$pos,1);
		print $char."\n";
		
		exit; 
	    }
	    # now, get a RANDOM column satisfying sel and selbase
	    my $N = 0;
	    for my $cols ( sort {$a cmp $b} keys %{$patterns{$sel}{$selbase}} ) {
		$N += $patterns{$sel}{$selbase}{$cols};
	    }

	    my $random = rand();
	    my ($min,$max) = (0,0);
	    for my $cols ( sort {$a cmp $b} keys %{$patterns{$sel}{$selbase}} ) {
		
		$max += $patterns{$sel}{$selbase}{$cols}/$N;

#		print join("\t",$min,'<=',$random,'<',$max,$cols)."\n";

		if( ($min<=$random) && ($max>$random) ) {

#		    print "At pos $pos, we selected: $cols\n";
		    for(my $k=0;$k<length($cols);$k++) {
			
			my $nbr = substr($cols,$k,1);
			my $letter = '?';
			if( $nbr==0 ) {  $letter='A'; }
			if( $nbr==1 ) {  $letter='C'; }
			if( $nbr==2 ) {  $letter='G'; }
			if( $nbr==3 ) {  $letter='T'; }
			if( $nbr==4 ) {  $letter='-'; }
			if( $nbr==5 ) {  $letter='N'; }

			if( $letter eq '?' ) {

			    print STDERR "Something is very wrong (?): $pos, $k, $specs[$k], $countload, $sel, $selbase, $cols\n";
			    exit;
			}

			$dummyseq{$pos}{$specs[$k]} = $letter;
			#{ print join("\t",'-->',$nbr,$letter,$specs[$k])."\n"; }
		    }
		    last;
		}

		$min += $patterns{$sel}{$selbase}{$cols}/$N;

	    }

#	    if( $pos>8 ) { last; }
	}

	# print alignment

	for my $spec (@specs) {


	    my $str_spec = '';
	    my $str_letters = '';


	    #if( $spec eq $REF ) { print ">"; }
	    if( $spec eq $REF ) { $str_spec .= ">"; }

	    #print ">$spec".'_'.$reffield."\n";
#	    $str_spec .= ">$spec".'_'.$reffield;
	    
	    my $ider = substr($reffield,length($REF)+1,length($reffield)-1-length($REF));

	    $str_spec .= ">$spec".'_'.$ider;

	    for(my$pos=0;$pos<$len;$pos++) {

		if( !exists($dummyseq{$pos}{$spec}) ) {

		    print STDERR "Something is very wrong: $pos, $spec, $countload\n";
		    exit;
		}

		#print $dummyseq{$pos}{$spec};
		$str_letters .= $dummyseq{$pos}{$spec};
	    }
	    
	    #print "\n";

	    if( $str_letters =~ /[A-Za-z]/ ) {

		print join("\n",$str_spec,$str_letters)."\n";
	    }
	    
	}
	
	
	$i--;
    }

#    if( $countload > 0 ) { last; }

}
