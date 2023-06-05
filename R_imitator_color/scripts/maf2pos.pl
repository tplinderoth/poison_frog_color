#!/usr/bin/perl

# maf2pos.pl

use warnings;
use strict;
use IO::Zlib;

my $minmaf = 0;

die(qq/
maf2pos.pl <maf file> <minimum maf to retain site [default: $minmaf]>
\n/) if (!@ARGV);

my $maffile = shift @ARGV;
open(my $maffh, '<', $maffile) or die("Unable to open MAF file $maffile: $!\n");
sysread $maffh, my $magic, 2;
close $maffh;
if ($magic eq "\x1f\x8b") {
        $maffh = new IO::Zlib;
        $maffh->open($maffile, "rb");
} else {
        open($maffh, '<', $maffile);
}

$minmaf = shift @ARGV if (@ARGV);
die("ERROR: Minimum MAF < 0") if ($minmaf < 0);

chomp(my @tok = split(/\s+/, <$maffh>));
my $mafidx;
for ($mafidx = 0; $mafidx <= $#tok; $mafidx++) {
	last if ($tok[$mafidx] eq 'knownEM');
}

while(<$maffh>) {
	my @tok = split(/\s+/, $_);
	print ("$tok[0]\t$tok[1]\n") if $tok[$mafidx] >= $minmaf;
}

close $maffh;

exit;
