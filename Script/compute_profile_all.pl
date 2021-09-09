#!/usr/bin/perl

use warnings;
my $gtf = shift; 
my $bam = shift;
my $strand = shift;
my $A_site_end = shift;
my $halfwindow = 100;
my %mid_plus;

open IN1,$gtf or die;
while(<IN1>){
my @a = split;
 if($a[2] eq 'start_codon'){
	if($a[6] eq '+'){
	 $mid_plus{$a[0]}{$a[4]} = 1;
	}
 }
}

close(IN1);


foreach $k1 (keys %mid_plus) {
	foreach $k2 ( keys %{$mid_plus{$k1}} ) {
        
 	 my $s_w = $k2 - $halfwindow;
 	 my $s_e = $k2 + $halfwindow;
 	 my @IN2="";

	 	if($strand eq 'pos_neg'){
  	  	 @IN2 = `samtools view -F 0x10 $bam "$k1:$s_w-$s_e"`;
	 	}else{
 	  	 @IN2 = `samtools view -f 0x10 $bam "$k1:$s_w-$s_e"`;
	 	}
	 
	 %pos = ();


		foreach (@IN2){
		 my @b = split;
		  
		 my $l = 0;
                 $b[5] =~ s/(\d+)[M]/$l+=$1/eg;
 		 my $length = $l;

		 if($A_site_end eq '3p'){
		  $b[3] = $b[3] + $length;
		 }
			if(($b[11] eq 'NH:i:1') and (($b[14] eq 'nM:i:1') | ($b[14] eq 'nM:i:0'))) { # Unique mapping read
		 		
				if(exists($pos{$length}{$b[3]})){
		  	 	 $pos{$length}{$b[3]}++;
				}
		    		else{
				$pos{$length}{$b[3]} = 1;
		     		}
			}
		}

		foreach $ke1 (keys %pos){
		 print("L:$ke1\t");
			for($i = $s_w; $i <= $s_e; $i++){
     				if(exists($pos{$ke1}{$i})){
      			 	 print("$pos{$ke1}{$i}\t");
      				}
      				else
      				{
      		 	 	 print("0\t");
      				}
			}
		print("\n");
		
		}
	}
}

