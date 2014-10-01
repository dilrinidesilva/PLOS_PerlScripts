#!/usr/bin/perl -w

## Script to determine SNPs in nonvariable and variable regions in CNEs using 
## an alignment of seven species

use strict;

#filter 1809 CNEs
my %filter =();
open F, "< 1809CNEs.txt"; # list of 1809 CNEs in analysis to be selected from the larger dataset
while(<F>){
    chomp;
    $filter{$_} =1;
}

close F;


my $outfile = ""; # output filename
open OUT, "> $outfile";

my $infile = ""; # file with tab separated CNE ID, genomic coordinates of CNE and SNPs 

open AF, "< $infile";
while(<AF>){
    chomp;
    my $line =$_;
    @_=split("\t",$_);
    my $cne=$_[8];
    if($filter{$cne}){
	my $cnestart=0; # genomic start coordinate of CNE from file
	my $snploc=0; # genomic coordinate of SNP from file
	my $state =""; # whether NVR or RVR
	my $pos=0; # position of SNP in the alignment string
	my $position = 0; # position for substr() command to extract bases from alignment 
	my %sequences=();
	my($humanseq, $macseq, $mouseseq, $chseq, $frogseq, $zseq, $fseq) = "";
	my $filename="$cne.aln"; # clustal alignment file
	$cnestart = $_[12]-1; 
	$snploc = $_[1]; 
	
	open ALN, "< $filename";
	
	while(<ALN>){
	    chomp;
	    if(!/^CLUSTAL/){
		if(/human/){
		    $humanseq .= substr($_, 16);
		}		
		if(/macaque/){
		    $macseq .= substr($_, 16);
		}
		if(/mouse/){
		    $mouseseq .= substr($_, 16);
		}
		if(/chicken/){
		    $chseq .= substr($_, 16);
		}
		if(/frog/){
		    $frogseq .= substr($_, 16);
		}
		if(/zfish/){
		    $zseq .= substr($_, 16);
		}
		if(/fugu/){
		    $fseq .= substr($_, 16);
		}
	    }
	    if(eof){
		%sequences = ("human",$humanseq,
			      "macaque",$macseq,
			      "mouse",$mouseseq,
			      "chicken",$chseq,
			      "frog", $frogseq,
			      "zfish", $zseq,
			      "fugu", $fseq);
		
		while($humanseq =~ /(A|T|G|C)/g){ # use the human sequence to get the position of the SNP in the alignment and then match to genomic coordinate
		    
		    $cnestart++; # iterate through each base in the human sequence of the alignment until the genomic coordinate of the SNP is reached 
		    
		    if($cnestart == $snploc){
			$pos=pos($humanseq);
		    
			$position = $pos-1; # offset by 1 for substr command 

			my($macbase,$mousebase, $chbase, $frogbase,$zbase, $fbase) ="";
 		       
                        #get the base at the SNP position for each species in the alignment (excluding human)
		
			$macbase = substr($sequences{macaque},$position,1);
			$mousebase = substr($sequences{mouse},$position,1);
			$chbase = substr($sequences{chicken},$position,1);
			$frogbase = substr($sequences{frog},$position,1);
			$zbase = substr($sequences{zfish},$position,1);
			$fbase = substr($sequences{fugu},$position,1);
	
			
			my $column = $macbase.$mousebase.$chbase.$frogbase.$zbase.$fbase;
			
			if($column =~ /A{6}|T{6}|C{6}|G{6}/){
		    
			    $state = "NVR";
		    
			}else{
			    $state = "RVR";
			}
		
			print OUT "$line\t$macbase\t$mousebase\t$chbase\t$frogbase\t$zbase\t$fbase\t$state\t$position\n";

		    }
		}
	    }
	}
	close ALN;
    }
}    
close AF;
close OUT;
