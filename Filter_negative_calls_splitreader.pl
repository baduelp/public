#!/usr/bin/perl

use strict ;
use warnings ;
use POSIX;
STDOUT->autoflush(1);

#############################################
#                                           #  
#               SPLITREADER                 #
#    Filter presence variant calls based    #
#  on negative coverage ratios (NC-filter)  #
#           Baduel et al. 2020              #
#                                           # 
#############################################

# # This script filters putative insertion sites based on the minimum negative coverage observed within each insertion site as well as 100bp up and 100bp downstream by BAM-readcount_wrapper.sh.
# # The BAMrc output files ($pop.$subsetname-insertion.[0/100up/100down].rc) are stored for each $pop in $workspace_dir/$subsetname/$pop/BAMrc/. 
# # If another directory architecture is used the paths on ll. 174, 189, and 204 need to be changed accordingly.
# # For further analysis four output files are generated :
# # $project_dir/$subsetname-insertions.$filtname.NC.DP$depth.bed contains all the NC information but is not filtered for the ratio of negative coverage (ratioNCfilt): for each putative insertion sites it shows for each sample the coverage supporting the non-reference insertion / the coverage supporting the reference absence / the coverage 100bp up / the coverage100bp down. 
# # $project_dir/$subsetname-insertions.ratioNC$filtname.bool.DP$depth.bed filters the putative insertions based on the ratioNCfilt with a 1 for samples where insertions have passed the filter, a 0 where they have not and the reference coverage is sufficient (>$DPrefmin reads) to be confident that the insertion is indeed absent, and - for the samples where the reference coverage is not sufficient to confirm the absence.
# # $project_dir/$subsetname-insertions.ratioNC$filtname.DP$depth.bed gives the coverage supporting the ratioNCfilt insertions
# # $project_dir/$subsetname-insertions.ratioNC$filtname.NConly.DP$depth.bed gives the reference coverage over ratioNCfilt insertion sites
## Questions or comments to pbaduel(ar)bio.ens.psl.eu

my $subsetname = $ARGV[0] ;# name of cohort
my $depth = $ARGV[1] ; # number of reads (split+discordant) used to call an insertion in Filter_insertion_splitreader.pl
my $project_dir = $ARGV[2]; # path to working directory where $subsetname-insertions.$filtname.DP$depth.bed is located
my $workspace_dir = $ARGV[3]; # path to input files by sample
my $filtname = $ARGV[4] ; # name of first pass filter (positive coverage) 

my @pops_to_analyze = @ARGV[5..$#ARGV]; # array of all samples in cohort

my $DPrefmin=5; # minimum negative coverage over insertion sites required to be confident that an insertion is indeed absent
my $DPrefmax=100; # maximum negative coverage over insertion sites to remove regions with aberrant sequencing depth (to be adjusted depending on library depth) 

print STDOUT "Calling Ref DP for $subsetname \n";


print STDOUT "number of individuals: ",$#pops_to_analyze+1,"\n";
my $ind_thresh=floor($#pops_to_analyze/10); # requires at least 10% of individuals covered
print STDOUT "number of individuals covered required by position: $ind_thresh\n";



print STDOUT "importing positions\t";
my %masked_pos ;
open IN, "<$project_dir/$subsetname-insertions.$filtname.DP$depth.bed" ;
while(<IN>){
	chomp $_ ;
	next if /^scaff/; # remove header line
	# Buffer array with values of every line
	my @line = split(/\t/, $_) ; 
	my $scaff = $line[0] ;
	my $pos_start = $line[1] ;
	my $pos_end = $line[2] ;

	$masked_pos{$scaff}{$pos_start}{$pos_end}=1;
}
close IN ;
print STDOUT "done\n";


my %TE_cov ;
my %ind_pop ;
open IN, "<$project_dir/$subsetname-insertions.$filtname.DP$depth.bed" ;


print STDOUT "finding sample index\t";
# FIND index of sample column IN header 
my @header = split(/\t/, <IN>) ;
foreach my $pop (@pops_to_analyze){
	# print STDOUT "$pop \t";
	foreach my $index (5..$#header) {
		chomp $header[$index] ;
		my $ind = $header[$index] ;		 
		#collect indices for each pop
		if($ind eq $pop){
			$ind_pop{$pop} = $index ; 
			# print STDOUT "$ind_pop{$pop}\n";
			}
	}
}
print STDOUT "done\n";

my %TE_hash ;
print STDOUT "importing DP\t";
while(<IN>){
	chomp $_ ;
	# Buffer array with values of every line
	my @line = split(/\t/, $_) ; 
	my $scaff = $line[0] ;
	my $pos_start = $line[1] ;
	my $pos_end = $line[2] ;
	my $TE_name = $line[3] ;
	my $TE_dups = $line[4] ;
	if ($masked_pos{$scaff}{$pos_start}{$pos_end}){
		# # # potentially several TEs can map to the same position
		my @cov_split = split /[()\/]+/, $TE_dups;
		$TE_cov{$scaff}{$pos_start}{$pos_end}{$TE_name}{'popnb'}=$cov_split[0];		
		$TE_cov{$scaff}{$pos_start}{$pos_end}{$TE_name}{'splitcov'}=$cov_split[1]+$cov_split[2];		
		foreach my  $pop (sort keys %ind_pop){
			$line[$ind_pop{$pop}] =~ s/\s+$//; #remove whitespace at end of value
			$TE_cov{$scaff}{$pos_start}{$pos_end}{$TE_name}{'DP'}{$pop}=$line[$ind_pop{$pop}];
			if($line[$ind_pop{$pop}]>0){
			$TE_hash{$scaff}{$pos_start}{$pos_end}{$pop}{$TE_name}+=1;
			}
		}
	}
}
close IN ;
print STDOUT "done\n";


my $max_cov;
my $max_pops;
my $max_TE;
print STDOUT "filtering duplicates\t";
foreach my $scf (sort keys %TE_hash){
	foreach my $posstart (sort{$a<=>$b} keys %{$TE_hash{$scf}} ){
		foreach my $posstop (sort{$a<=>$b} keys %{$TE_hash{$scf}{$posstart}} ){
			foreach my  $pop (sort keys %ind_pop){
				my $size = keys %{$TE_hash{$scf}{$posstart}{$posstop}{$pop}};
				if($size>1){
					print STDOUT "\n", $scf, "\t", $posstart, "\t", $posstop, "\t", $pop, "\t";
					# several TEs present at same location in sam pop
					$max_cov=0;
					$max_pops=0;
					foreach my $TE (keys %{$TE_hash{$scf}{$posstart}{$posstop}{$pop}} ){
						print STDOUT $TE, "\t", $TE_cov{$scf}{$posstart}{$posstop}{$TE}{'popnb'},"(", $TE_cov{$scf}{$posstart}{$posstop}{$TE}{'splitcov'},")", "\t";
						# find TE found in max nb of pops
						if ($TE_cov{$scf}{$posstart}{$posstop}{$TE}{'popnb'}>$max_pops){
							$max_pops=$TE_cov{$scf}{$posstart}{$posstop}{$TE}{'popnb'};
							$max_cov=$TE_cov{$scf}{$posstart}{$posstop}{$TE}{'splitcov'};	
							$max_TE=$TE;
						}
						# if tied find TE found with max nb of split reads
						elsif($TE_cov{$scf}{$posstart}{$posstop}{$TE}{'popnb'}==$max_pops){
							if ($TE_cov{$scf}{$posstart}{$posstop}{$TE}{'splitcov'}>$max_cov){
								$max_cov=$TE_cov{$scf}{$posstart}{$posstop}{$TE}{'splitcov'};	
								$max_TE=$TE;
							}
						}
					}

					# mark all the non-best as duplicates
					foreach my $TE (keys %{$TE_hash{$scf}{$posstart}{$posstop}{$pop}} ){
						if($TE ne $max_TE){
							$TE_cov{$scf}{$posstart}{$posstop}{$TE}{'dup'}+=1;
						}
						else{
							$TE_cov{$scf}{$posstart}{$posstop}{$TE}{'dup'}=0;
							print STDOUT $TE, "\t", $TE_cov{$scf}{$posstart}{$posstop}{$TE}{'popnb'},"(", $TE_cov{$scf}{$posstart}{$posstop}{$TE}{'splitcov'},")";
						}
					}


				}
				else{
					foreach my $TE (keys %{$TE_hash{$scf}{$posstart}{$posstop}{$pop}} ){
						$TE_cov{$scf}{$posstart}{$posstop}{$TE}{'dup'}=0;
					}
				}
			}
		}
	}
}
print STDOUT "done\n";


print STDOUT "importing NC cov\t";
my %ref_cov ;
foreach my  $pop (sort keys %ind_pop){
	print STDOUT "$pop\t";
	open IN2, "<$workspace_dir/$subsetname/$pop/BAMrc/$pop.$subsetname-insertion.0.rc" ;
	while(<IN2>){
		chomp $_ ;
		my @line = split(/\t/, $_) ;
		# coverage already processed and min cov output for each insertion site
		$ref_cov{$pop}{$line[0]}{$line[1]}{$line[2]}=$line[3]; # min coverage over insertion
	}
}
close IN2;
print STDOUT "done\n";

print STDOUT "importing NC cov 100up\t";
my %up_cov ;
foreach my  $pop (sort keys %ind_pop){
	print STDOUT "$pop\t";
	open IN3, "<$workspace_dir/$subsetname/$pop/BAMrc/$pop.$subsetname-insertion.100up.rc" ;
	while(<IN3>){
		chomp $_ ;
		my @line = split(/\t/, $_) ;
		# coverage already processed and min cov output for each insertion site
		$up_cov{$pop}{$line[0]}{$line[1]+100}{$line[2]+100}=$line[3]; # min coverage over insertion
	}
}
close IN3;
print STDOUT "done\n";

print STDOUT "importing NC cov 100down\t";
my %down_cov ;
foreach my  $pop (sort keys %ind_pop){
	print STDOUT "$pop\t";
	open IN4, "<$workspace_dir/$subsetname/$pop/BAMrc/$pop.$subsetname-insertion.100down.rc" ;
	while(<IN4>){
		chomp $_ ;
		my @line = split(/\t/, $_) ;
		# coverage already processed and min cov output for each insertion site
		$down_cov{$pop}{$line[0]}{$line[1]-100}{$line[2]-100}=$line[3]; # min coverage over insertion
	}
}
close IN4;
print STDOUT "done\n";

open OUT, ">$project_dir/$subsetname-insertions.$filtname.NC.DP$depth.bed" ;
open OUT3, ">$project_dir/$subsetname-insertions.ratioNC$filtname.NConly.DP$depth.bed" ;
open OUT4, ">$project_dir/$subsetname-insertions.ratioNC$filtname.DP$depth.bed" ;
open OUT5, ">$project_dir/$subsetname-insertions.ratioNC$filtname.bool.DP$depth.bed" ;

print OUT "scaff \t start \t end \t TE_name \t popNB(NBsplit) \t" ;
print OUT3 "scaff \t start \t end \t TE_name \t popNB(NBsplit) \t" ;
print OUT4 "scaff \t start \t end \t TE_name \t popNB(NBsplit) \t" ;
print OUT5 "scaff \t start \t end \t TE_name \t popNB(NBsplit) \t" ;
foreach my  $pop (sort keys %ind_pop){
	print OUT $pop, "\t" ;
	print OUT3 $pop, "\t" ;
	print OUT4 $pop, "\t" ;
	print OUT5 $pop, "\t" ;
}
print OUT "\n" ;
print OUT3 "\n" ;
print OUT4 "\n" ;
print OUT5 "\n" ;

my $output_ratio;
my $site_filter_bool;
my $site_filterNO_bool;
my $NC_filter_nb;
foreach my $scf (sort keys %TE_cov){
	foreach my $posstart (sort{$a<=>$b} keys %{$TE_cov{$scf}} ){
		foreach my $posstop (sort{$a<=>$b} keys %{$TE_cov{$scf}{$posstart}} ){
			foreach my $TE (keys %{$TE_cov{$scf}{$posstart}{$posstop}} ){
				if ($TE_cov{$scf}{$posstart}{$posstop}{$TE}{'dup'}<1){ # check if site is not a duplicate
					print OUT $scf, "\t", $posstart, "\t", $posstop, "\t", $TE, "\t", $TE_cov{$scf}{$posstart}{$posstop}{$TE}{'popnb'},"(", $TE_cov{$scf}{$posstart}{$posstop}{$TE}{'splitcov'},")", "\t";
					# find corresponding ref_cov
					$site_filter_bool=0;
					$site_filterNO_bool=0;
					$NC_filter_nb=0;
					foreach my  $pop (sort keys %ind_pop){
						if( defined $ref_cov{$pop}{$scf}{$posstart}{$posstop}){
							# # min cov already processed and output by bamrc
							$TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DPref'}{$pop}=$ref_cov{$pop}{$scf}{$posstart}{$posstop};
							$TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DPtot'}{$pop}=$TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DP'}{$pop}+$TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DPref'}{$pop};
							$TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DPup'}{$pop}=$up_cov{$pop}{$scf}{$posstart}{$posstop};
							$TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DPdown'}{$pop}=$down_cov{$pop}{$scf}{$posstart}{$posstop};
							
							print OUT "$TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DP'}{$pop}/$TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DPref'}{$pop}/$TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DPup'}{$pop}/$TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DPdown'}{$pop}\t" ;
							
							if($TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DPup'}{$pop}+$TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DPdown'}{$pop}>0){
								# if surrounding coverage is not null: outputs DPneg/mean(DPup, DPdown)
								$output_ratio=2*$TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DPref'}{$pop}/($TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DPup'}{$pop}+$TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DPdown'}{$pop});
								if($TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DP'}{$pop}>0){
									# positive call
									if($output_ratio<=0.5){
										# NOT FALSE POSITIVE 
										$site_filter_bool+=1;
									}
									else{
										$site_filterNO_bool+=1;
									}
									if($TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DPtot'}{$pop}<$DPrefmax){
										$NC_filter_nb+=1;
									}

								}
								else{
									# negative call
									if($TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DPtot'}{$pop}<$DPrefmax && $TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DPref'}{$pop}>=$DPrefmin){
										$NC_filter_nb+=1;
									}
									}
								
								}
						}
						else{
							print OUT "$TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DP'}{$pop}/-/-/-\t" ;
							print STDOUT "$scf $posstart $posstop $pop missing \n";
						}

					}
					print OUT "\n";

					# print site to OUT4 only if at least 1 TRUE POSITIVEs and at least 10% inds covered
					if($site_filter_bool>=1 && $NC_filter_nb>=$ind_thresh){
						print OUT3 $scf, "\t", $posstart, "\t", $posstop, "\t", $TE, "\t", $TE_cov{$scf}{$posstart}{$posstop}{$TE}{'popnb'},"(", $TE_cov{$scf}{$posstart}{$posstop}{$TE}{'splitcov'},")", "\t";
						print OUT4 $scf, "\t", $posstart, "\t", $posstop, "\t", $TE, "\t", $TE_cov{$scf}{$posstart}{$posstop}{$TE}{'popnb'},"(", $TE_cov{$scf}{$posstart}{$posstop}{$TE}{'splitcov'},")", "\t";
						print OUT5 $scf, "\t", $posstart, "\t", $posstop, "\t", $TE, "\t", $TE_cov{$scf}{$posstart}{$posstop}{$TE}{'popnb'},"(", $TE_cov{$scf}{$posstart}{$posstop}{$TE}{'splitcov'},")", "\t";
						foreach my  $pop (sort keys %ind_pop){
							if( defined $ref_cov{$pop}{$scf}{$posstart}{$posstop}){
								print OUT3 "$TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DPref'}{$pop}\t" ;
								print OUT4 "$TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DP'}{$pop}\t" ;
								if( $TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DP'}{$pop} > 0 ){ #carriers
									print OUT5 "1\t" ;
								}
								else{
									if( $TE_cov{$scf}{$posstart}{$posstop}{$TE}{'DPref'}{$pop} >= $DPrefmin ){ #NC covered non-carriers
									print OUT5 "0\t" ;
									}
									else{
									print OUT5 "-\t" ;
									}
								}
							}
							else{
							print OUT3 "-\t" ;
							print OUT4 "-\t" ;
							print OUT5 "-\t" ;

							}
						}
						print OUT3 "\n";
						print OUT4 "\n";
						print OUT5 "\n";
					}
				}
			}
		}
	}
}

close OUT ;
close OUT3 ;
close OUT4 ;
close OUT5 ;


exit ;
