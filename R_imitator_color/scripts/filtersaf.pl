#!/usr/bin/perl

# filterSaf.pl

use warnings;
use strict;
use IO::Zlib;

die(qq/
filterSaf.pl [Nminor] [region file] [saf TSV file]

Discards variable sites with fewer than Nminor minor alleles. Fixed categories are retained.
E.g. Nminor=6 would discard sites with 1-5 minor alleles. 
\n/) if (!@ARGV || scalar @ARGV < 3);

my $nminor = $ARGV[0];
die("Nminor should be at least 2\n") if $nminor < 2;

open(my $rf, '<', $ARGV[1]) or die("Couldn't open region file $ARGV[1]: $!\n");
my %contigs;
my @contig_order;
while(<$rf>) {
	chomp;
	@{$contigs{$_}} = ();
	push @contig_order, $_;
}

my $saf = $ARGV[2];

open(my $saffh, '<', $saf) or die("Unable to open SAF file $saf: $!\n");
sysread $saffh, my $magic, 2;
close $saffh;
if ($magic eq "\x1f\x8b") {
	$saffh = new IO::Zlib;
	$saffh->open($saf, "rb");
} else {
	open($saffh, '<', $saf);
}

while(<$saffh>) {
	my @tok = split(/\s+/, $_);
	my $contig = shift @tok;
	my $pos = shift @tok;
	my $maf = 0;
	if ($tok[0] < 0) {
		my $i = 0;
		foreach my $like (@tok) {
			$maf = $i if ($tok[$i] > $tok[$maf]);
			$i++;
		}
	}
	next if ($maf > 0 && $maf < $nminor);
	push @{$contigs{$contig}}, $pos;
}

close $saffh;

foreach my $assembly (@contig_order) {
	if (exists $contigs{$assembly} && scalar @{$contigs{$assembly}} > 0) {
		foreach my $site (@{$contigs{$assembly}}) {
			print "$assembly\t$site\n";
		}
	}
}

exit;
