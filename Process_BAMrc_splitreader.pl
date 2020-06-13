#!/usr/bin/perl
use strict ;
use warnings ;

#############################################
#                                           #  
#               SPLITREADER                 #
#  analyze and combine negative coverages   #
#           Baduel et al. 2020              #
#                                           # 
#############################################

##Questions or comments to pbaduel(ar)bio.ens.psl.eu

my $subsetname = $ARGV[0] ; # name of cohort
my $depth = $ARGV[1] ; # number of reads required to call an insertion on 1st pass
my $pop = $ARGV[2]; # name of genome
my $OutputDir = $ARGV[3];  #path/to/output
my $region = $ARGV[4]; #0 for insertion sites and 100up or 100down for 100bp up or down respectively

print STDERR "Calling Min DP for ".$pop."\n";

my %TE_cov ;
open IN, "<$OutputDir/$subsetname-insertion-sites.$region.bed" ;

while(<IN>){
	chomp $_ ;
	if ( $_ =~ m/start/ ) { 	
			next ; 
		}
	$_ =~ /^$/ and next;
	# Buffer array with values of every line
	my @line = split(/\t/, $_) ; 
	my $scaff = $line[0] ;
	my $pos_start = $line[1] ;
	my $pos_end = $line[2] ;

	$TE_cov{$scaff}{$pos_start}{$pos_end}+=1;
}
close IN ;


my %ref_cov ;
open IN2, "<$OutputDir/$pop.$subsetname-insertion.tmp.$region.rc" ;
while(<IN2>){
	# bamtools-readcount output
	chomp $_ ;
	my @line = split(/\t/, $_) ;
	my $refRC = 0 ;
	my $i = 5 ;
	while ($i <= 8 && $refRC == 0) {	
		my @readcount = split(/:/, $line[$i]) ; # coverage per allele
		if ($readcount[0] eq $line[2]){ # ONLY TAKE REF Allele
			$refRC = $readcount[1] ;
		}
		$i++; 
	}
	$ref_cov{$line[0]}{$line[1]}=$refRC;
}
close IN2;


open OUT, ">$OutputDir/$pop.$subsetname-insertion.$region.rc" ;

foreach my $scf (sort {$a cmp $b} keys %TE_cov){
	foreach my $posstart (sort{$a<=>$b} keys %{$TE_cov{$scf}} ){
		foreach my $posstop (sort{$a<=>$b} keys %{$TE_cov{$scf}{$posstart}} ){
			print OUT $scf, "\t", $posstart, "\t", $posstop, "\t";
			# find corresponding ref_cov
			# BAMtools-readcount output
			my $minDP=0;
			foreach my $pos (keys %{$ref_cov{$scf}}){
				if($pos>=$posstart && $pos<=$posstop){
					if($minDP==0){
					$minDP=$ref_cov{$scf}{$pos};}
					elsif($ref_cov{$scf}{$pos}<$minDP){
						$minDP=$ref_cov{$scf}{$pos};}
				}
			}
			print OUT "$minDP\n" ;
		}
	}
}

close OUT ;

exit ;
