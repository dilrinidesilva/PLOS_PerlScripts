#!/usr/env/perl -w

## Script to count nonvariable and variable regions in CNEs using 
## an alignment of six species (i.e. excluding human)

use strict;

my(%filter,%sequences) =();

#filter 1809 CNEs

open F, "< 1809CNEs_list.txt"; # list of 1809 CNEs in analysis to be selected from the larger dataset
while(<F>){
    chomp;
    $filter{$_} =1;
}

close F;

my $outfile = ""; # specify output filename

open OUT, "> $outfile";

print OUT "CNEID\tAlnLength\tNVRcount\tRVRcount\n";



for my $cne (keys %filter){
    my$filename="";
 
    $filename="$cne.aln"; # clustal aligned file 
 
    my($nvrCount, $rvrCount) =0;
    open ALN, "< $filename";

    my($macseq, $mouseseq, $chseq, $frogseq, $zseq, $fseq) = ""; # variable to store the sequence for each species in the alignment
    my %sequences=(); # hash to store the sequences
    my $seqLength =0;

    while(<ALN>){
	chomp;
	if(!/^CLUSTAL/){
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
	        %sequences = ("macaque",$macseq,
			      "mouse",$mouseseq,
			      "chicken",$chseq,
			      "frog", $frogseq,
			      "zfish", $zseq,
			      "fugu", $fseq);
		    
		my $seqLength = length($sequences{macaque}); # only one line of the alignment is required to get the length
 
		my $position = 0; 

		while($position < $seqLength){ #iterate through each base in the alignment

		    my($macbase,$mousebase, $chbase, $frogbase,$zbase, $fbase) ="";
		       

		    $macbase = substr($sequences{macaque},$position,1);
		    $mousebase = substr($sequences{mouse},$position,1);
		    $chbase = substr($sequences{chicken},$position,1);
		    $frogbase = substr($sequences{frog},$position,1);
		    $zbase = substr($sequences{zfish},$position,1);
		    $fbase = substr($sequences{fugu},$position,1);
		    
		    my $column = $macbase.$mousebase.$chbase.$frogbase.$zbase.$fbase;

		    if($column =~ /A{6}|T{6}|C{6}|G{6}/){
                        
			$nvrCount++;			    
		    }

		    $position++;
		}

		$rvrCount = $seqLength-$nvrCount;

		print OUT "$cne\t$seqLength\t$nvrCount\t$rvrCount\n";
	}
    }

    close ALN;
}

close OUT;
