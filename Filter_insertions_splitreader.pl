#!/usr/bin/perl
use strict ;
use warnings ;

#############################################
#                                           #  
#               SPLITREADER                 #
#  filter insertions by positive coverage   #
#           Baduel et al. 2020              #
#                                           # 
#############################################

##Questions or comments to pbaduel(ar)bio.ens.psl.eu

## define subroutines
sub max {
    my ($max, @vars) = @_;
    for (@vars) {
        $max = $_ if $_ > $max;
    }
    return $max;
}
sub binary {
    my ($target, $direction, @array) = @_;
	my $size = $#array; 
    my ($left, $right) = (0, $size);

	if ($array[$left] <= $target && $array[$right] >= $target){
		while ($right - $left > 8){
			my $tmp = ($right+$left) / 2;
			if ($array[$tmp] > $target){
				$right = $tmp;
			} 
			elsif ($array[$tmp] < $target) {
				$left = $tmp;
			}
			else{ #on target by chance
				$right = $tmp;				
				$left = $tmp;
			}
		}
		if ($direction==1){ #finding first element above value
			for ($left .. $right){
				return $_ if $array[$_] >= $target;
			}
		}
		else{  #finding last element below value
			for (-$right .. -$left){ 
				return -$_ if $array[-$_] <= $target;
			}
		}
	}
	elsif ($array[$left] > $target ){
		return $left;}
	else{
		return $right;
		}
}


my $subsetname = $ARGV[0] ; # name of cohort
my $depth = $ARGV[1] ; # number of reads (split+discordant) required to call an insertion on 1st pass
my $project_dir = $ARGV[2]; # path to working directory
my $fam = $ARGV[3]; # TE family 
my $out_dir = $ARGV[4]; # path to output directory by TE family
my $TSDthresh = $ARGV[5]; # threshold above witch presence variants are rejected based on expected TSD length of TE family
my $noTSDthresh = $ARGV[6]; # threshold above witch presence variants are rejected based when no TSD defined for TE family

my @all_pops = @ARGV[7..$#ARGV]; # list of all the genomes analysed in the cohort 


# retrieve TSD length for all TE families
my %TSDlength;
open IN, "<$project_dir/TE_sequence/TE-information-all.txt" ;
while(<IN>){
	chomp $_ ;
	my @line = split(/\t/, $_) ; 
	$TSDlength{$line[0]}=$line[1];
}
close IN ; 


print STDERR "Filtering Insertions for ".$subsetname." at DP ".$depth."\n";


print STDOUT "number of individuals $#all_pops\n";

# find missing genomes because of BAM issues 
my %missing_inds;
foreach my $ind (0..$#all_pops){
	$missing_inds{$all_pops[$ind]}=0;
}
open IN, "<$project_dir/BAMs/$subsetname.missing_bams.txt" ;
while(<IN>){
	chomp $_ ;
	my @line = split(/\t/, $_) ; 
	if($#line>0){
	$line[0] =~ s/\s+$//;
	$missing_inds{$line[0]}=1;}
}
close IN ; 
print ".";

open IN, "<$project_dir/BAMs/$subsetname.truncated_bams.txt" ;
while(<IN>){
	chomp $_ ;
	my @line = split(/\t/, $_) ; 
	if($#line>0){
	$line[0]=~ s/\s+$//;
	$missing_inds{$line[0]}=1;}
}
close IN ; 
print ".";

my @pops_to_analyze;
foreach my $ind (0..$#all_pops){
	if ($missing_inds{$all_pops[$ind]}==1){
		print STDERR "$all_pops[$ind] missing", "\n";}
	else{
		push(@pops_to_analyze, $all_pops[$ind])
	}

}


# # retrieve all intervals where putative presence variant calls overlap for a given TE family
my %intersect_array_sort;
open IN, "<$out_dir/$fam.$subsetname-intersect.sort.bed" ;
while(<IN>){
	chomp $_ ;
	
	next if /^chrom/; #skip header

	my @line = split(/\t/, $_) ; 

	push(@{$intersect_array_sort{$line[0]}{'start'}}, $line[1]);
	push(@{$intersect_array_sort{$line[0]}{'stop'}}, $line[2]);
	push(@{$intersect_array_sort{$line[0]}{'indNB'}}, $line[3]);
}
close IN ; 

# # loop through all merged intervals to define minimal interval as most shared between genomes
my %intersect_array;
my %peak_array;
open IN, "<$out_dir/$fam.$subsetname-intersect.mrg.bed" ;
while(<IN>){
	chomp $_ ;
	
	next if /^chrom/; #skip header

	my @line = split(/\t/, $_) ; 
	my $pos_start = $line[1] ;
	my $pos_end = $line[2] ;
	push(@{$intersect_array{$line[0]}{'breadth'}}, $pos_end-$pos_start);
	push(@{$intersect_array{$line[0]}{'start'}}, $pos_start);
	push(@{$intersect_array{$line[0]}{'stop'}}, $pos_end);

	# # LOOK UP IN UNMERGED ARRAY FOR MOST COMMON POSITION
	my $intersect_inf = binary($pos_start, 1, @{$intersect_array_sort{$line[0]}{'start'}}); #first interval of intersect starting after or at pos_start
	my $intersect_sup = binary($pos_end, -1, @{$intersect_array_sort{$line[0]}{'stop'}}); # last interval of intersect ending before pos_end
	
	if ($intersect_sup < $intersect_inf){ # when first interval starting after or at pos_start ends after pos_end due to merging of neighboring insertions for intersect
		$intersect_sup=$intersect_inf;
	}

	# first define boundaries of most common position as highest peak
	my $max_intersect=0;
	my $max_start=${$intersect_array_sort{$line[0]}{'start'}}[$intersect_inf];
	my $max_stop=${$intersect_array_sort{$line[0]}{'stop'}}[$intersect_sup];
	my $peak_ind=$intersect_inf;
	foreach my $pos_ind ($intersect_inf..$intersect_sup){
		my $inter_start = ${$intersect_array_sort{$line[0]}{'start'}}[$pos_ind];
		my $inter_NB = ${$intersect_array_sort{$line[0]}{'indNB'}}[$pos_ind]; # nb of bed overlapping in region
		if ($inter_NB > $max_intersect ){
			$max_intersect=$inter_NB;
		 	$max_start= $inter_start;
			$max_stop = ${$intersect_array_sort{$line[0]}{'stop'}}[$pos_ind];
			$peak_ind=$pos_ind;
		}
	}

	# # refine will be done after taking into account split reads
	my $hard_start=0; # consider interval borders as flexible
	my $hard_stop=0; # consider interval borders as flexible

	# # add peak associated with merge interval
	push(@{$peak_array{$line[0]}{$pos_start}{'maxstart'}}, $max_start);
	push(@{$peak_array{$line[0]}{$pos_start}{'maxstop'}}, $max_stop);
	push(@{$peak_array{$line[0]}{$pos_start}{'indNB'}}, 0);

	push(@{$peak_array{$line[0]}{$pos_start}{'hardstart'}}, $hard_start);
	push(@{$peak_array{$line[0]}{$pos_start}{'hardstop'}}, $hard_stop);
	
}
close IN ; 


print ".";

## for each genome check that insertion boundaries defined by split reads do not conflict with maxstart and maxstop by more than 2bp
## if conflict between intervals then define new peak
foreach my $pop (@pops_to_analyze){
my $bedfile = "$out_dir/$fam.$pop-insertion-sites.bed";
	if (-f $bedfile){ #check if file exists
		open IN, "<$bedfile" ;
		while(<IN>){
			chomp $_ ;
			my @line = split(/\t/, $_) ; 
			my $pos_start = $line[1] ;
			my $pos_end = $line[2] ;
			my $TE_name = $line[3] ;

			my @sum_line;
			foreach my $index (5..6){  # total number of reads (split)
				$sum_line[$index] = 0;
				my @rnb = split(",", $line[$index]) ; 
				if ($#rnb>0){ # more than one read number (separate clusters of reads)
					foreach my $c (0..$#rnb){
						$sum_line[$index] += $rnb[$c] ;}
				}
				else{ # only one read number
				$sum_line[$index] += $line[$index] ;}
			}
			my $split_left = $sum_line[6]; ## SPLIT READS ON THE RIGHT DEFINE LEFT BOUNDARIES!!		 
			my $split_right = $sum_line[5]; ## SPLIT READS ON THE LEFT DEFINE RIGHT BOUNDARIES!!
			
			my $intersect_inf = binary($pos_start, -1, @{$intersect_array{$line[0]}{'start'}}); #last interval of merged intersect starting before or at pos_start
			my $intersect_sup = binary($pos_end, 1, @{$intersect_array{$line[0]}{'stop'}}); # first interval of merged intersect ending after or at pos_end
			
			if ($intersect_sup < $intersect_inf){ # when first interval starting after or at pos_start ends after pos_end due to merging of neighboring insertions for intersect
				$intersect_sup=$intersect_inf;
			}

			print STDERR "$pop\t$line[0]\t$pos_start\t$pos_end\t$split_left\t$split_right\t${$intersect_array{$line[0]}{'start'}}[$intersect_inf]\t${$intersect_array{$line[0]}{'stop'}}[$intersect_sup]\n";
			
			if ($intersect_sup == $intersect_inf) {
				my $max_ind = $intersect_inf;
				my $inter_start = ${$intersect_array{$line[0]}{'start'}}[$max_ind];
				# # look through peaks associated with merge interval
				my $disagreement = 1;
				my $peak_ind=0;
				while($disagreement>0 && $peak_ind<=$#{$peak_array{$line[0]}{$inter_start}{'maxstart'}}){
					# # check if pop interval agrees with max interval
					my $max_start = ${$peak_array{$line[0]}{$inter_start}{'maxstart'}}[$peak_ind];
					my $max_stop = ${$peak_array{$line[0]}{$inter_start}{'maxstop'}}[$peak_ind];
					my $hard_start = ${$peak_array{$line[0]}{$inter_start}{'hardstart'}}[$peak_ind];
					my $hard_stop = ${$peak_array{$line[0]}{$inter_start}{'hardstop'}}[$peak_ind];

					if($pos_start<=$max_stop && $pos_end>=$max_start){ #intervals overlap at least partially
						$disagreement = 0;
						# # compare max_start and stop to split-read-defined bounds
						if($split_left > 0 ){ # start defined by at least 1 split read
							if($pos_start<$max_start-1 || $pos_start>$max_start+1){ #starts disagree by 2bp or more
								if($hard_start>0){ #peak boundaries well defined
									$disagreement = 1;
								}
								else{ #peak boundaries still flexible > redefine start with split_start
									${$peak_array{$line[0]}{$inter_start}{'maxstart'}}[$peak_ind]=$pos_start;
								}
							}
							if($pos_start==$max_start){
								${$peak_array{$line[0]}{$inter_start}{'hardstart'}}[$peak_ind]+=$split_left;
							}
						}

						if($split_right > 0 ){ # stop defined by at least 1 split read
							if($pos_end<$max_stop-1 || $pos_end>$max_stop+1){ #stops disagree by 2bp or more
								if($hard_stop>0){ #peak boundaries well defined
									$disagreement = 1;
								}
								else{ #peak boundaries still flexible > redefine start with split_stop
									${$peak_array{$line[0]}{$inter_start}{'maxstop'}}[$peak_ind]=$pos_end;
								}
							}
							if($pos_end==$max_stop){
								${$peak_array{$line[0]}{$inter_start}{'hardstop'}}[$peak_ind]+=$split_right;
							}
						}

						
						# # no disagreement if one of the boundaries is not defined by split reads (then the other is not well defined by splitreader even if supported by split reads)
						if( $split_left ==0 || $split_right ==0 ){ 
							$disagreement = 0;								
						}

						# # no disagreement if one of the split boundaries match perfectly (+/-1) even if the other doesn't
						if( $split_left >0 && $hard_start>0 && ($pos_start >= $max_start-1 && $pos_start <= $max_start+1) ){ 
							$disagreement = 0;								
						}
						if( $split_right >0 && $hard_stop>0 && ($pos_end >= $max_stop-1 && $pos_end <= $max_stop+1) ){ 
							$disagreement = 0;					
						}
					}


					if($pos_end<$max_start && $hard_start==0){ #intervals don't overlap but merged boundaries are flexible
						$disagreement = 0;
						${$peak_array{$line[0]}{$inter_start}{'maxstart'}}[$peak_ind]=$pos_end;
					}

					if($pos_start>$max_stop && $hard_stop==0){ #intervals don't overlap but merged boundaries are flexible
						$disagreement = 0;
						${$peak_array{$line[0]}{$inter_start}{'maxstop'}}[$peak_ind]=$pos_start;
					}

					if($pos_end<$max_start && $pos_end>=$max_start-10 && $split_right==0 && $split_left==0){ #intervals don't overlap by much and both individual boundaries are flexible
						$disagreement = 0; ## do not create secondary peak
					}

					if($pos_start>$max_stop && $pos_start<=$max_stop+10 && $split_right==0  && $split_left==0){ #intervals don't overlap by much and both individual boundaries are flexible
						$disagreement = 0; ## do not create secondary peak
					}

					$peak_ind++;				
				}
				my $max_peak_ind = $peak_ind-1;

				if($disagreement > 0){
					# # no peak match > add new peak
					# # add peak associated with merge interval
					push(@{$peak_array{$line[0]}{$inter_start}{'maxstart'}}, $pos_start);
					push(@{$peak_array{$line[0]}{$inter_start}{'maxstop'}}, $pos_end);
					push(@{$peak_array{$line[0]}{$inter_start}{'indNB'}}, 1);

					push(@{$peak_array{$line[0]}{$inter_start}{'hardstart'}}, $split_left);
					push(@{$peak_array{$line[0]}{$inter_start}{'hardstop'}}, $split_right);

					$max_peak_ind++;
					}
				else{
					# Found peak overlap 
					${$peak_array{$line[0]}{$inter_start}{'indNB'}}[$max_peak_ind]+=1;
					${$peak_array{$line[0]}{$inter_start}{'hardstart'}}[$max_peak_ind]+=$split_left;
					${$peak_array{$line[0]}{$inter_start}{'hardstop'}}[$max_peak_ind]+=$split_right;
				}
			}
			else{
			print STDERR "No overlap for $pop over scaffold $line[0] from $pos_start to $pos_end\n";}
		}
		close IN ;
	}
	else{print STDERR "No putative insertions detected in $pop\n";}
	print ".";
}

## loop again through genomes to attribute them to the peaks just redefined
my %TE_sites ;
foreach my $pop (@pops_to_analyze){
my $bedfile = "$out_dir/$fam.$pop-insertion-sites.bed";
	if (-f $bedfile){ #check if file exists
		open IN, "<$bedfile" ;
		while(<IN>){
			chomp $_ ;
			my @line = split(/\t/, $_) ; 
			my $pos_start = $line[1] ;
			my $pos_end = $line[2] ;
			my $TE_name = $line[3] ;
			
			my @sum_line;
			foreach my $index (5..6){  # total number of reads (split)
				$sum_line[$index] = 0;
				my @rnb = split(",", $line[$index]) ; 
				if ($#rnb>0){ # more than one read number (separate clusters of reads)
					foreach my $c (0..$#rnb){
						$sum_line[$index] += $rnb[$c] ;}
				}
				else{ # only one read number
				$sum_line[$index] += $line[$index] ;}
			}
			my $split_left = $sum_line[6]; ## SPLIT READS ON THE RIGHT DEFINE LEFT BOUNDARIES!!		 
			my $split_right = $sum_line[5]; ## SPLIT READS ON THE LEFT DEFINE RIGHT BOUNDARIES!!
			
			my $intersect_inf = binary($pos_start, -1, @{$intersect_array{$line[0]}{'start'}}); #last interval of intersect starting before or at pos_start
			my $intersect_sup = binary($pos_end, 1, @{$intersect_array{$line[0]}{'stop'}}); # first interval of intersect ending after or at pos_end
			
			if ($intersect_sup < $intersect_inf){ # when first interval starting after or at pos_start ends after pos_end due to merging of neighboring insertions for intersect
				$intersect_sup=$intersect_inf;
			}

			if($intersect_sup == $intersect_inf){
				my $max_ind = $intersect_inf;
				my $inter_start = ${$intersect_array{$line[0]}{'start'}}[$max_ind];
				# # look through peaks associated with merge interval
				my $disagreement = 1;
				my $peak_ind=0;
				while($disagreement>0 && $peak_ind<=$#{$peak_array{$line[0]}{$inter_start}{'maxstart'}}){
					# # check if pop interval agrees with max interval
					my $max_start = ${$peak_array{$line[0]}{$inter_start}{'maxstart'}}[$peak_ind];
					my $max_stop = ${$peak_array{$line[0]}{$inter_start}{'maxstop'}}[$peak_ind];
					my $hard_start = ${$peak_array{$line[0]}{$inter_start}{'hardstart'}}[$peak_ind];
					my $hard_stop = ${$peak_array{$line[0]}{$inter_start}{'hardstop'}}[$peak_ind];

					if($pos_start<=$max_stop && $pos_end>=$max_start){ #intervals overlap at least partially
						$disagreement = 0;
						# # compare max_start and stop to split-read-defined bounds
						if($split_left >0 ){ # start defined by at least 1 split read
							if($pos_start<$max_start-1 || $pos_start>$max_start+1){ #starts disagree by 2bp or more
								if($hard_start>0){
									$disagreement = 1;
								}
							}
						}

						if($split_right>0){ # stop defined by at least 1 split read
							if($pos_end<$max_stop-1 || $pos_end>$max_stop+1){ #stops disagree by 2bp or more
								if($hard_stop>0){ #peak boundaries well defined
									$disagreement = 1;
								}
							}
						}

						

						# # no disagreement if one of the boundaries is not defined by split reads (then the other is not well defined by splitreader even if supported by split reads)
						if( $split_left ==0 || $split_right ==0 ){ 
							$disagreement = 0;								
						}
						# # no disagreement if one of the split boundaries match perfectly (+/-1) even if the other doesn't 
						if( $split_left >0 && $hard_start>0 && ($pos_start >= $max_start-1 && $pos_start <= $max_start+1) ){ 
							$disagreement = 0;								
						}
						if( $split_right >0 && $hard_stop>0 && ($pos_end >= $max_stop-1 && $pos_end <= $max_stop+1) ){ 
							$disagreement = 0;					
						}
					}
					

					if($pos_end<$max_start && $pos_end>=$max_start-10 && $split_right==0 && $split_left==0){ #intervals don't overlap by much and both individual boundaries are flexible
						$disagreement = 0; 
					}

					if($pos_start>$max_stop && $pos_start<=$max_stop+10 && $split_right==0  && $split_left==0){ #intervals don't overlap by much and both individual boundaries are flexible
						$disagreement = 0; 
					}

					
					$peak_ind++;
					
				}

				my $max_peak_ind = $peak_ind-1;

				if($disagreement == 0){# # found matching peak (always should based on previous loop)			
					my @sum_line;
					foreach my $index (5..8){  # total number of reads (split+discordant)
						$sum_line[$index] = 0;
						my @rnb = split(",", $line[$index]) ; 
						if ($#rnb>0){ # more than one read number (separate clusters of disc reads)
							foreach my $c (0..$#rnb){
								$sum_line[$index] += $rnb[$c] ;}
						}
						else{ # only one read number
						$sum_line[$index] += $line[$index] ;}
					}
					my $sum_left = $sum_line[5]+ $sum_line[7];		
					my $sum_right = $sum_line[6]+ $sum_line[8];
					my $sum = $sum_left + $sum_right;
						
					$TE_sites{$line[0]}{$max_ind}{$max_peak_ind}{$pop}{$TE_name}{'all'} = $sum;
					$TE_sites{$line[0]}{$max_ind}{$max_peak_ind}{$pop}{$TE_name}{'left'} = $sum_left;
					$TE_sites{$line[0]}{$max_ind}{$max_peak_ind}{$pop}{$TE_name}{'right'} = $sum_right;
					
				}
				else{
				print STDERR "No peak overlap for $pop over scaffold $line[0] from $pos_start ($split_left) to $pos_end ($split_right) vs ${$peak_array{$line[0]}{$inter_start}{'maxstart'}}[0] (${$peak_array{$line[0]}{$inter_start}{'hardstart'}}[0]) to ${$peak_array{$line[0]}{$inter_start}{'maxstop'}}[0] (${$peak_array{$line[0]}{$inter_start}{'hardstop'}}[0])\n";}
			}
			else{
			print STDERR "No overlap for $pop over scaffold $line[0] from $pos_start to $pos_end\n";}
		}
		close IN ;
	}
	else{print STDERR "No putative insertions detected in $pop\n";}
	print ".";
}


open OUT, ">$out_dir/$fam.$subsetname-insertions.filt3.DP$depth.bed" ;

print OUT "scaff \t start \t end \t TE_name \t overlap(NBsplit_start/NBsplit_stop) \t" ;
foreach my $ind (0..$#pops_to_analyze){
	print OUT $pops_to_analyze[$ind], "\t" ;
}
print OUT "\n" ;

foreach my $scaff (sort keys %TE_sites){
	foreach my $pos (sort{$a<=>$b} keys %{$TE_sites{$scaff}} ){
		my $inter_start=${$intersect_array{$scaff}{'start'}}[$pos];
		foreach my $peak (sort{$a<=>$b} keys %{$TE_sites{$scaff}{$pos}} ){
			my $peak_start=${$peak_array{$scaff}{$inter_start}{'maxstart'}}[$peak];
			my $peak_stop=${$peak_array{$scaff}{$inter_start}{'maxstop'}}[$peak];
			my $hard_start=${$peak_array{$scaff}{$inter_start}{'hardstart'}}[$peak];
			my $hard_stop=${$peak_array{$scaff}{$inter_start}{'hardstop'}}[$peak];
			my $peak_breadth=$peak_stop-$peak_start;
			# # ONLY KEEP ONE TE PER POS
			my $maxDP = 0 ;
			my $TENB = 0 ;
			my $popNB = 0 ;
			my $TE_name = 'NA';
			my %TE_hash; 
			foreach my $pop (keys %{$TE_sites{$scaff}{$pos}{$peak}} ){
				$popNB++;
				my $TEpopNB = 0 ;
				foreach my $TE (keys %{$TE_sites{$scaff}{$pos}{$peak}{$pop}} ){
					$TEpopNB ++ ; # nb of TEs with an insertion at this site in this pop
					$TE_hash{$TE}+=1;
					if($TE_sites{$scaff}{$pos}{$peak}{$pop}{$TE}{'all'}>$maxDP && $TE_sites{$scaff}{$pos}{$peak}{$pop}{$TE}{'left'}>=1 && $TE_sites{$scaff}{$pos}{$peak}{$pop}{$TE}{'right'}>=1){
					$maxDP = $TE_sites{$scaff}{$pos}{$peak}{$pop}{$TE}{'all'} ;
					$TE_name = $TE ; # only keeps TE with the highest DP and at least 1 read on each side
					}
				}
				if($TEpopNB>$TENB){$TENB=$TEpopNB;} #KEEP MAX NB OF TE PER SITE

			}

			my $print_bool=0;
			if($maxDP>=$depth){ #at least one individual with number of reads higher than $depth
				if(defined $TSDlength{$TE_name} && $hard_start>0  && $hard_stop>0 ){
					if($peak_breadth<=$TSDthresh+$TSDlength{$TE_name}){ # matching TSD length (+ $TSDthresh)
					$print_bool=1;
					}
				}
				else{
					if($peak_breadth<$noTSDthresh){ # 
					$print_bool=1;
					}
				}
			}

			if($print_bool>0){
				print OUT $scaff, "\t", $peak_start, "\t", $peak_stop,"\t", $TE_name, "\t", $popNB,"(", $hard_start,"/", $hard_stop,")", "\t"; 
				foreach my $ind (0..$#pops_to_analyze){
					if( exists $TE_sites{$scaff}{$pos}{$peak}{$pops_to_analyze[$ind]}{$TE_name} ){
					print OUT $TE_sites{$scaff}{$pos}{$peak}{$pops_to_analyze[$ind]}{$TE_name}{'all'}, "\t" ;
					}
					else{ print OUT "0 \t";}
				}
				print OUT "\n";
			}
		}
	}
}

close OUT ;

exit ;
