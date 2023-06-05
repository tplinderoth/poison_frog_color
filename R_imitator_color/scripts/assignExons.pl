#!/usr/bin/perl

# assignExons.pl <mafs file> <exon bed file>

use warnings;
use strict;

die(qq/
assignExons.pl <mafs file> <exon bed file>
\n/) if (!@ARGV || scalar @ARGV < 2);

open(my $maffh, '<', $ARGV[0]) or die("Unable to open MAFs file $ARGV[0]: $!\n");
open(my $bedfh, '<', $ARGV[1]) or die("Unable to open bed file $ARGV[1]: $!\n");

my %exons;
while (<$bedfh>) {
	chomp;
	my ($contig, $start, $stop) = ($1, $2, $3) if $_ =~ /^(.+)_(\d+)_(\d+)\s+/;
	push @{$exons{$contig}{start}}, $start;
	push @{$exons{$contig}{stop}}, $stop;
}
close $bedfh;

$" = "\t";
chomp(my @header = split("\t",<$maffh>));
print STDOUT "@header\n";

while (<$maffh>) {
	chomp;
	my @tok = split("\t", $_);
	my $found = 0;
	my $ex = "";
	if (exists $exons{$tok[0]}) {
		for (my $i = 0; $i <= $#{$exons{$tok[0]}{start}}; $i++) {
			my $a = ${$exons{$tok[0]}{start}}[$i];
			my $b = ${$exons{$tok[0]}{stop}}[$i];
			if ($tok[1] >= $a && $tok[1] <= $b) {
				$ex = "${tok[0]}_${a}_${b}";
				$found = 1;
				last;
			}
		}
	} else {
		die("$tok[0] not in BED file: $!\n");
	}
	die("Unable to find exon for $tok[0] $tok[1]: $!\n") if (!$found);
	$tok[0] = $ex;
	print STDOUT "@tok\n";
}

close $maffh;

exit;
