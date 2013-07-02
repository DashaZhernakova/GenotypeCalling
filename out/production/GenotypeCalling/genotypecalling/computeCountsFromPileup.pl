#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $pileupFile = $ARGV[0];			# pileup file from samtools
my $countsFile = $ARGV[1];			# name for the output file
my $minNumReads = $ARGV[2];			# minimum number of reads
my $refAllele;
my $bases;
my $counts;
my $refAlleleCount;
my $nonRefAlleleCount;

my $line;
my @split;

open(FILE, "< $pileupFile") or die "Can't open file: $! for reading";
open(OUTFILE, "> $countsFile") or die "Can't open file: $! for writing";

print OUTFILE "Chromosome\tPosition\tReferenceAllele\tNrReportedReads\tNrCountedReads\tCountsRefAllele\tCountsNonrefAllele\tCountsA\tCountsC\tCountsG\tCountsT\n";

while(<FILE>) {
	$line = $_;
	chomp($line);
	@split = split("\t",$line);
		
	$bases = $split[4];
	$refAllele = uc($split[2]);
	
	# count alleles
	$counts = &countAlleles($refAllele, $bases);
	$nonRefAlleleCount = 0;
	
	foreach my $allele (keys %$counts) {
		if ($allele eq $refAllele) {
			$refAlleleCount = $counts->{$allele};
		}
		else {
			$nonRefAlleleCount += $counts->{$allele};
		}
	}
	
	
	# check number of reads for the position
	if (($counts->{"G"} + $counts->{"A"} + $counts->{"T"} + $counts->{"C"}) >= $minNumReads) {
	
		print OUTFILE join("\t", @split[0..3]);
		print OUTFILE "\t" . ($counts->{"G"} + $counts->{"A"} + $counts->{"T"} + $counts->{"C"});
		print OUTFILE "\t" . $refAlleleCount . "\t" . $nonRefAlleleCount; 
		
		foreach my $allele (sort keys %$counts) {
			print OUTFILE "\t" . $counts->{$allele}; 
		}
		
		print OUTFILE "\n";
	}
}

close FILE;
close OUTFILE;

sub countAlleles {
	my $refAllele = shift();
	my $bases = shift();
	my %alleleCounts;
	
	my @alleleA = $bases =~ /A|a/g;
	my @alleleC = $bases =~ /C|c/g;
	my @alleleT = $bases =~ /T|t/g;
	my @alleleG = $bases =~ /G|g/g;
	my @alleleRef = $bases =~ /\.|\,/g;
	
	if ($refAllele eq "A") {
		@alleleA = @alleleRef;
	}
	elsif ($refAllele eq "C") {
		@alleleC = @alleleRef;
	}
	elsif ($refAllele eq "T") {
		@alleleT = @alleleRef;
	}
	elsif ($refAllele eq "G") {
		@alleleG = @alleleRef;
	}
	
	$alleleCounts{"A"} = scalar @alleleA;
	$alleleCounts{"C"} = scalar @alleleC;
	$alleleCounts{"T"} = scalar @alleleT;
	$alleleCounts{"G"} = scalar @alleleG;
	
	return(\%alleleCounts);
}