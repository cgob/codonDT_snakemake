#!/usr/bin/perl
use warnings;
 my $glob_strand;
 ############Function definition#######

 # Load Information from fasta CDS

 sub load_cds{
 my %master;
 my @list=@_;
 open CDS,$list[0] or die;
 local $/ = "\n>";
 my $sm ="";
 	
	if($list[1] eq 'pos_neg'){
	$sm=1;		
	$glob_strand="-f";
	}else{
	$sm=-1;
	$glob_strand="-F";
	}

        while(<CDS>){
          chomp;
          @a = split("\n",$_);
          @attr = split(/ /,$a[0]);
          $trans = $attr[0];
          $gene = (split(/\:/,$attr[3]))[1];
          $chr = (split(/\:/,$attr[2]))[2];
          $start = (split(/\:/,$attr[2]))[3];
          $stop = (split(/\:/,$attr[2]))[4];
          $strand = (split(/\:/,$attr[2]))[5];
          $comp = "$chr\:$start\-$stop";
          @a = @a[1 .. $#a];
          $seq = join('',@a);
	 $l2=length($seq);
                 if($attr[5] =~ /protein_coding/ and $l2 > 115){
                 $master{${gene}}{${trans}}{'seq'}=$seq;  
                 $master{${gene}}{${trans}}{'comp'}=$comp;
		 	if($strand == $sm){
		 	 $master{${gene}}{${trans}}{'strand'} = "-F";
			}else{
			 $master{${gene}}{${trans}}{'strand'} = "-f";
			}
                 }
        

}
 close CDS;
 return %master;
 }
 
 
############## Input files #############
 my $strand_mode = 'pos_neg';
 my $cds = shift;
 my %master = load_cds($cds, $strand_mode);
############# Output file #############
my %seq;
        foreach $k (sort keys %master){ # loop around the gene
	%seq=();
        	for $k2 (sort keys %{$master{$k}}){ # loop around the transcript
		 my @codon=unpack("(A3)*", $master{$k}{$k2}{'seq'});
			for my $i (0 .. ($#codon-39)) {
			my $stri= join("",@codon[$i..$i+39]);
			$seq{$k}{$stri}=1
			} 

                        }

for my $key1 ( keys %seq) {
    for my $key2 ( keys %{$seq{$key1}}) {
    	print("$key1\t$key2\n");
	}
}
}

