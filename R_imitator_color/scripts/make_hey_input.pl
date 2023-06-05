#!/usr/bin/perl

# make_hey_input.pl

use warnings;
use strict;

die(qq/
make_hey_input.pl <hka file> <sequence length file> <ingroup code: 1 = imitator, 2 = band, 3 = stripe>
\n/) unless (@ARGV && scalar @ARGV == 3);

open(my $hkafh, '<', $ARGV[0]) or die("Unable to open HKA file $ARGV[0]: $!\n");
open(my $lenfh, '<', $ARGV[1]) or die("Unable to open sequence length file $ARGV[1]: $!\n");
my $s1code = $ARGV[2];
if ($s1code == 1) {
} elsif ($s1code == 2) {
} elsif ($s1code == 3) {
} else {
	die("Invalide ingroup code $s1code: $!\n");
}

my %len;
<$lenfh>; # skip header
while(<$lenfh>) {
	chomp;
	my @tok = split(/\s+/, $_);
	$len{$tok[0]} = $tok[1];
}
close $lenfh;

my $inheritance = 1;
my $nimi = 0;
if ($s1code == 1) {
	$nimi = 132; # number imitator haplotypes
} elsif ($s1code == 2 || $s1code == 3) {
	$nimi = 66; # number banded/striped haplotypes
} else {
	die("Invalide ingroup code $s1code: $!\n");
}

my $nmodel = 1; # number model species haplotypes (considering one outgroup sequence)

my %seen;
if ($s1code == 1) {
	print STDOUT "FULL.ID\tSHORT.ID\tINHERITANCE\tIMI.SEQ.LEN\tMODEL.SEQ.LEN\tSEQ.COMPARE.LEN\tIMI.N.SEQ\tMODEL.N.SEQ\tIMI.S\tMODEL.S\tD\n";
} elsif ($s1code == 2) {
	print STDOUT "FULL.ID\tSHORT.ID\tINHERITANCE\tBAND.SEQ.LEN\tMODEL.SEQ.LEN\tSEQ.COMPARE.LEN\tBAND.N.SEQ\tMODEL.N.SEQ\tBAND.S\tMODEL.S\tD\n";
} elsif ($s1code == 3) {
	print STDOUT "FULL.ID\tSHORT.ID\tINHERITANCE\tSTRIPE.SEQ.LEN\tMODEL.SEQ.LEN\tSEQ.COMPARE.LEN\tSTRIPE.N.SEQ\tMODEL.N.SEQ\tSTRIPE.S\tMODEL.S\tD\n";
} else {
	die("Invalid ingroup code $s1code: $!\n");
}
<$hkafh>; # skip header
while(<$hkafh>) {
	chomp;
	my @larr = split('\t',$_);
	my $id = $larr[0];
	my ($s1, $d);
	if ($s1code == 1) {
		$s1 = $larr[1];
		$d = $larr[2];
	} elsif ($s1code == 2) {
		$s1 = $larr[13];
		$d = $larr[14];
	} elsif ($s1code == 3) {
		$s1 = $larr[25];
		$d = $larr[26];
	} else {
		die("Invalid ingroup code $s1: $!\n");
	}

	my $s2 = 0;

	my $namearr = "";
	$namearr = $1 if $id =~ /_combined_(\S+)$/;
	if ($namearr) {
		my @tok = split('_', $namearr);
		if ($tok[0] =~ /\D/ && $tok[0] ne "NA") {
			if (exists $tok[0]) {
				$seen{$tok[0]}++;
			} else {
				$seen{$tok[0]} = 1;
			}
			my $idshort = scalar(@tok == 3) ? "${tok[0]}_${seen{$tok[0]}}" : $tok[0];
			print STDOUT "${id}\t${idshort}\t${inheritance}\t$len{$id}\t$len{$id}\t$len{$id}\t$nimi\t$nmodel\t$s1\t$s2\t$d\n" if (exists $len{$id});
		}
	}
}

exit;
