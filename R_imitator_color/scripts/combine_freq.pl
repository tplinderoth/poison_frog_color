#!/usr/bin/perl

# combine_freq.pl

use warnings;
use strict;

die(qq/
combine_freq.pl <comma-delimited pop name prefixes> <all-pop MAF> <pop1 MAF> <pop2 MAF> ... <popn MAF>
\n/) if (!@ARGV || scalar @ARGV < 2);

my @names = split(',',$ARGV[0]);
die ("Number of names does not match number of MAF files: $!\n") if (scalar @names != (scalar @ARGV)-1);

open(my $metafh, '<', $ARGV[1]) or die("Unable to open all-pops MAF file $ARGV[0]\n");

my @site_order;
my %sites;

<$metafh>; #skip header

while (<$metafh>) {
	chomp;
	my @tok = split('\s+', $_);
	my $id = "${tok[0]}:${tok[1]}";
	push @site_order, $id;
	$sites{$id}{"${names[0]}_major"} = $tok[2];
	$sites{$id}{"${names[0]}_minor"} = $tok[3];
	$sites{$id}{"${names[0]}_maf"} = $tok[4];
	$sites{$id}{"${names[0]}_snp_pval"} = $tok[5];
	$sites{$id}{"${names[0]}_nind"} = $tok[6];
}
close $metafh;

for (my $i=2; $i <= $#ARGV; $i++) {
	my $pop = $names[$i-1];
	open(my $fh, '<', $ARGV[$i]) or die("Unable to open MAF file $ARGV[$i]: $!\n");
	<$fh>; #skip header
	while (<$fh>) {
		my @tok = split('\s+', $_);
		my $id = "${tok[0]}:${tok[1]}";
		die("$tok[0] $tok[1] not present in MAF file 1\n") if (!exists $sites{$id});
		$sites{$id}{"${pop}_maf"} = $tok[4];
		$sites{$id}{"${pop}_mafpval"} = $tok[5];
	}
	close $fh;
}

print STDOUT "chromo\tposition\t${names[0]}.major\t${names[0]}.minor\t${names[0]}.maf\t${names[0]}.snp.pval\t${names[0]}.nind";
for (my $i = 1; $i <= $#names; $i++) {
	print STDOUT "\t${names[$i]}.maf\t${names[$i]}.maf.pval";
}
print STDOUT "\n";

foreach (@site_order) {
	my @pos = split(':',$_);
	print STDOUT "$pos[0]\t$pos[1]\t", $sites{$_}{"${names[0]}_major"},"\t",$sites{$_}{"${names[0]}_minor"},"\t",$sites{$_}{"${names[0]}_maf"},"\t",$sites{$_}{"${names[0]}_snp_pval"},"\t",$sites{$_}{"${names[0]}_nind"};
	for (my $i=1; $i <= $#names; $i++) {
		if (exists $sites{$_}{"${names[$i]}_maf"}) {
			print STDOUT "\t",$sites{$_}{"${names[$i]}_maf"},"\t",$sites{$_}{"${names[$i]}_mafpval"};
		} else {
			print STDOUT "\tNA\tNA";
		}
	}
	print STDOUT "\n";
}

exit;
