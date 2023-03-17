use strict;
use warnings;
use Carp;
my %pair;
my %GE;
my %agp;
my $anchors=$ARGV[0];
my $gff1=$ARGV[1];
my $gff2=$ARGV[2];
my $paf=$ARGV[3];
my $n=$ARGV[4];
open AN,"$anchors"or die "cannot open the file $anchors\n";
while(my $line=<AN>){
	chomp $line;
	next if($line=~/^#/);
	my $pr=(split /\t/,$line)[0];
	if($pr=~/"(\w+)-(\w+)"/){
		my $prA=$1;
		my $prB=$2;
		$pair{$prB}=$line;
	}

}
close AN;

open GFFA,"$gff1"or die "cannot open the file $gff1\n";
while(my $line=<GFFA>){
	chomp $line;
	my $id;
	my ($type,$attribute)=(split /\t/,$line)[2,8];
	if($type eq "gene"){
		if($attribute=~/ID=(\w+);/){
			$id=$1;
			$GE{$id}=$line;
		}
		
	}
}
close GFFA;

open GFFB,"$gff2"or die "cannot open the file $gff2\n";
while(my $line=<GFFB>){
        chomp $line;
        my $id;
        my ($type,$attribute)=(split /\t/,$line)[2,8];
        if($type eq "gene"){
                if($attribute=~/ID=(\w+);/){
                        $id=$1;
                        $GE{$id}=$line;
                }

        }
}
close GFFB;


open PAF,"$paf"or die "cannopt open the file $paf\n";
while(my $line=<PAF>){
	chomp $line;
	my ($chrB,$startB,$endB,$strand,$chrA,$startA,$endA,$cigar)=(split /\t/,$line)[0,2,3,4,5,7,8,-1];
	next unless($chrA eq $chrB);
	my $position="$startB,$endB";
#	chr1    1       128481886       1       W       h1tg000012l     1       128481886       -
	my $temp="$chrA\t$startA\t$endA\t-\tW\t$chrB\t$startB\t$endB\t$strand\t$cigar";
	$agp{$chrB}{$position}=$temp;
}
close PAF;


foreach my $prB(keys %pair){
	my $line=$pair{$prB};
	my $prA;
	if($line=~/"(\w+)-(\w+)"/){
		$prA=$1;
        }
	my $gfA=$GE{$prA};
	my $gfB=$GE{$prB};
	my $start;
	my $end;
	my($chrA,$startA,$endA,$strandA,$phaseA,$attributesA)=(split /\t/,$gfA)[0,3,4,6,7,8];
	my($chrB,$startB,$endB,$strandB,$phaseB,$attributesB)=(split /\t/,$gfB)[0,3,4,6,7,8];
	next unless($chrA eq $chrB);
	my $hash2 = $agp{$chrB};
	my $agp_line;
	my $cigar;
	foreach my $key (keys %$hash2){
		my ($s,$e)=(split /,/,$key)[0,1];
		if($startB >= $s && $endB <= $e){
			$agp_line=$agp{$chrB}{$key};
		}

	}
	
	unless($agp_line){
		carp "cannot find the $prB in paf file";
		next;
	}

	my $agp_strand=(split /\t/,$agp_line)[8];
	if($agp_strand eq "-"){
		$start=&agp_turn($endB,$agp_line);
		$end=&agp_turn($startB,$agp_line);
	}else{

		$start=&agp_turn($startB,$agp_line);
		$end=&agp_turn($endB,$agp_line);
	}
	
	if(abs($start-$startA)<$n or abs($end-$endA)<$n){
		print "$line\n";
	}

}

sub agp_turn{
#本子程序利用类似agp的方式将位置信息转移到对应的位置
#输入$in,agp
	my $pos=$_[0];
	my $agp=$_[1];
	my $out;
	#10 80  +
	#chr1    1       128481886       1       W       h1tg000012l     1       128481886       -	cigar
	my ($ref_chr,$ref_start,$ref_end,$qry_chr,$qry_start,$qry_end,$agp_strand,$cigar)=(split /\t/,$agp)[0,1,2,5,6,7,8,9];
	my $ref_length=abs($ref_end-$ref_start)+1;
	unless($cigar){
		$cigar="${ref_length}M";
	}

	if($pos >= $qry_start and $pos <= $qry_end){
        	if($agp_strand eq "-"){
        		$out=$qry_end-$pos+1;
			$out=&turn_cigar($out,$cigar);
			$out=$out+$ref_start-1;
        	}else{
			$out=$pos-$qry_start+1;
			$out=&turn_cigar($out,$cigar);
			$out=$out+$ref_start-1;
        	}
			return $out;
    	}else{
		my $warn="not in agp";
		carp "$pos not in agp $agp";
		return "$warn";
    	}
	

}

sub turn_cigar{
	my $pos=$_[0];
	my $cigar=$_[1];
	my $pointer_ref=0;
	my $pointer_qry=0;
	my $q_pos=$pos;
	while ($cigar =~ /(\d+)([MIDNSHP=X])/g){		
		my $length = $1;
		my $operation = $2;
		my $add;
		#print "$length\t$operation\n";
		next if($operation eq "P");
		if($operation eq "D" or $operation eq "N"){
			$add=0;
		}else{
			$add=$length;
		}
		#首先判断是否会溢出
		if(($pointer_qry+$add)<$pos){
			#没有溢出时
		#	print "$length$operation\t没溢出\n";
			if($operation eq "I" or $operation eq "S" or $operation eq "H"){
				$q_pos-=$length;
		#		print "-$q_pos\n";
				$pointer_qry+=$length;
			}elsif($operation eq "D" or $operation eq "N"){
				$q_pos+=$length;
		#		print "+$q_pos\n";
				$pointer_ref+=$length;
			}elsif($operation eq "M" or $operation eq "=" or $operation eq "X"){
				$pointer_ref+=$length;
				$pointer_qry+=$length;
			}else{
				print "unkown operation $length$operation\n";
				next;
			}
		}else{

			#出现溢出
		#	print "$length$operation\t溢出\n";
			if($operation eq "I" or $operation eq "S" or $operation eq "H"){
				$q_pos-=($pos-$pointer_qry);
				$q_pos++;
				return "$q_pos";
				last;
				
			}elsif($operation eq "D" or $operation eq "N"){
				return "$q_pos";
				last;
			}elsif($operation eq "M" or $operation eq "=" or $operation eq "X"){
				return "$q_pos";
				last;
			}
		}
	}

}
