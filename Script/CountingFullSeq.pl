#!/usr/bin/perl
#use POSIX;
use warnings;
 ############Function definition#######

 # Load Information from fasta CDS

sub load_cds{
 my %master;
 my @list = @_;
 open CDS,$list[0] or die;
 local $/ = "\n>";
 my $sm ="";

        if($list[1] eq 'pos_neg'){
        $sm=1;
        }else{
        $sm=-1;
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
                 if($attr[5] =~ /protein_coding/ and $l2 > 120){
                # $seq = substr $seq, 45, -45;
                 $master{${gene}}{${trans}}{'seq'}=$seq;
                 $master{${gene}}{${trans}}{'comp'}=$comp;
		 $master{${gene}}{${trans}}{'strand_2'} = $strand;
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

 my $rfp = shift;
 my $l_1 = shift;
 my $l_2 = shift;
 my $strand_mode =shift;
 my $cds = shift; #'/scratch/cgobet/Genomes/Yeast/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa';
 my $file_out = shift;
 my %master = load_cds($cds, $strand_mode);

############# Output file #############
 open(my $fh, '>', $file_out) or die "Could not open file '$file_out' $!";

########## Call Reads via samtools ### 
 my %codon;
 my %read_tot=();
 my %codon_pos;
        foreach $k (sort keys %master){ # loop around the gene
	 %codon = ();
	 %codon_pos=();
	for $k2 (sort keys %{$master{$k}}){ # loop around the transcript
		print(""); 
		my @RE= `samtools view $master{$k}{$k2}{'strand'}  0x10 $rfp "$master{$k}{$k2}{'comp'}"`;
		foreach(@RE){ #loop around the reads
			my @read = split;
			 my $rseq = '';
			 my $length = '';
			 my $fivep=0;
				if($master{$k}{$k2}{'strand_2'} == -1){
			 	 my @pat = split(/(?<=\D)/, $read[5]);
                         	 my @pat_2 = reverse @pat;
                         	 my $lpat = join('',@pat_2);
                         	 $read[5] = $lpat;
                         	 $read[9] = reverse $read[9];
                         	 $read[9] =~ tr/NATCG/NTAGC/;
			 	}
			 	
			  my $lpi = $read[5];
			  my ($statut) = join('',$read[5] =~ /[A-Z]+/g);
			  			 
				  unless(exists($read_tot{$read[0]})){ # make sure we haven't seen the read globally
				 
				 	if($read[14] eq 'nM:i:0' and $read[11] eq 'NH:i:1'){
				 	#if($read[14] eq 'nM:i:0'){
					 my $l = 0;
		                         $lpi =~ s/(\d+)[M]/$l+=$1/eg;
		                         $length = $l;
		                          
						if($read[5] =~ /^(\d+)[S]/){
		                                 $rseq = substr $read[9], $1, $length;
						}else{
		                                 $rseq = substr $read[9], 0, $length;
						}
							
						if($master{$k}{$k2}{'seq'} =~ $rseq and $length > $l_1 and $length < $l_2){ # if sequence match fasta cds
						 #print("$master{$k}{$k2}{'seq'}\t$rseq\n");
						 my $posi = $-[0];
						 $read_tot{$read[0]}=1;	
					         my $posi2=0;	
							if(($posi + 1) % 3 == 0){ # Align the read in the right frame
							 $rseq = substr($rseq, 1);
							 $posi2=$posi+1;
							}

							if(($posi + 2) % 3 == 0){ # Align the read in the right frame
		                                         $rseq = substr($rseq, 2);
							 $posi2=$posi+2;
		                                        }
		                                        
					  	 my $ll = length($rseq);
					  	 my $rseq_2= substr($master{$k}{$k2}{'seq'}, $posi2-60, 120);
						 
							if(exists($codon{$rseq_2})){
						 	 $codon{$rseq_2}++;
						 	 $codon_pos{$rseq_2}{'length'}=$posi2;
							}else{
						 	 $codon{$rseq_2}=1;
						 	 $codon_pos{$rseq_2}{'trans'}=  $k2;	
						 	 $codon_pos{$rseq_2}{'length'}=  $posi2;	
							}
						}
					
				
				 	}
			}	
		}
	}
                foreach $w (keys %codon){
		 print $fh  "$k\t$w\t$codon{$w}\t$codon_pos{$w}{'length'}\t$codon_pos{$w}{'trans'}\n";	
		}		
                                        
 }


close $fh;

