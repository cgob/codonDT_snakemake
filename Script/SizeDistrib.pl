#!/usr/bin/perl
use warnings;

my $f1= shift;
my %length;
open BAM,"samtools view $f1 |";
my %reads;

while(<BAM>){
 next if ($_ =~ m/^\@/);
 s/\n//;  
 s/\r//;
 my @a = split(/\t+/);
 my $l = 0;
 $a[5] =~ s/(\d+)[M]/$l+=$1/eg;  
 
 	if(($a[11] eq 'NH:i:1') and (not exists($reads{$a[0]}))){ 
		if(exists($length{$l})){
 		 $length{$l}++;
 		}
 		else{
 	 	 $length{$l} = 1;
 		}		 
	$reads{$a[0]}=1;
	}
}

foreach my $i (1..75){
 if(exists($length{$i})){
  print("$i\t$length{$i}\n");
 }else{
  print("$i\t0\n");
 }
}

close BAM;

