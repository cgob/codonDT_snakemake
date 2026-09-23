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
          # optional shard filter: keep only genes on this shard's contigs
          next if defined($list[2]) && !exists($list[2]->{$chr});
          @a = @a[1 .. $#a];
          $seq = join('',@a);
         $l2=length($seq);
                 if($attr[5] =~ /protein_coding/ and $l2 > 120){
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
 my $strand_mode = shift;
 my $a_site_end = shift;
 my $cds = shift;
 my $file_apos = shift;
 my $file_out = shift;
 my $contig_file = shift;   # optional: file listing the contigs of this shard
 my $keep;
 if (defined $contig_file) {
   my %kc;
   open my $CF, '<', $contig_file or die "Cannot open $contig_file: $!";
   while (my $c = <$CF>) { chomp $c; $kc{$c} = 1 if length $c; }
   close $CF;
   $keep = \%kc;   # an empty shard file yields no genes, never "all genes"
 }
 my %master = load_cds($cds, $strand_mode, $keep);


############ Load A site position #####
my %As_pos;

open my $fh, '<', $file_apos or die "Cannot open: $!";

while (my $line = <$fh>) {
  chomp($line);
  my @array = split /\t/, $line;
  $As_pos{$array[0]}{'0'} = $array[1];
  $As_pos{$array[0]}{'1'} = $array[2];
  $As_pos{$array[0]}{'2'} = $array[3];
}
close $fh;

############# Output files #############

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
					 my $l = 0;
		                         $lpi =~ s/(\d+)[M]/$l+=$1/eg;
		                         $length = $l;
		                          
						if($read[5] =~ /^(\d+)[S]/){
		                                 $rseq = substr $read[9], $1, $length;
						}else{
		                                 $rseq = substr $read[9], 0, $length;
						}
							
						if($master{$k}{$k2}{'seq'} =~ $rseq and $length >= $l_1 and $length <= $l_2){ # if sequence match fasta cds
						 my $posi = $-[0];
						 $read_tot{$read[0]}=1;	
					         my $posi2=$posi;
						 # offsets are keyed by residue mod 3 of the anchored read end
						 my $frame = ($a_site_end eq '5p') ? $posi % 3 : ($posi + $length) % 3;

						 # no offset for this length: drop the read
						 next unless exists $As_pos{$length} and defined $As_pos{$length}{$frame};
							# A-site position; the 120 nt window starts 75 nt upstream of it
							my $a_pos = ($a_site_end eq '5p')
							          ? $posi2 + $As_pos{$length}{$frame}
							          : $posi2 + $l - $As_pos{$length}{$frame};
							my $shift_pos = $a_pos - 15 - 60;

					  	 my $rseq_2= substr($master{$k}{$k2}{'seq'}, $shift_pos, 120);
						 
							if(exists($codon{$rseq_2})){
						 	 $codon{$rseq_2}++;
						 	 $codon_pos{$rseq_2}{'pos'}=$a_pos;
							}else{
						 	 $codon{$rseq_2}=1;
						 	 $codon_pos{$rseq_2}{'trans'}=  $k2;	
						 	 $codon_pos{$rseq_2}{'pos'}=  $a_pos;	
							}
						}
					
				
				 	}
			}	
		}
	}
                foreach $w (keys %codon){
		 print $fh  "$k\t$w\t$codon{$w}\t$codon_pos{$w}{'pos'}\t$codon_pos{$w}{'trans'}\n";	
		}		
                                        
 }


close $fh;

