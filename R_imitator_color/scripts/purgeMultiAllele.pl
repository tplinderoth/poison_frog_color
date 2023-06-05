#!/usr/bin/perl

# purgeMultiAllele.pl
# requires bcftools to be installed and in user's PATH

use warnings;
use strict;

die(qq/purgeMultiAllele.pl <vcf> <ref fasta>\n/) if (!@ARGV && scalar(@ARGV) < 2);
die ("Unable to locate VCF $ARGV[0]\n") if (!-e $ARGV[0]);
die ("Unable to locate FASTA $ARGV[1]\n") if (!-f $ARGV[1]);

my $scriptname = $1 if ($0 =~ /([^\/]+)$/);
my $header_command = "##purgeMultiAlleleCommand=$scriptname @ARGV; Date=" . localtime();
my $bcftools_command = "bcftools norm --atom-overlaps '*' -f $ARGV[1] -m +any $ARGV[0]";

open(my $vcf, '-|', $bcftools_command) or die("Unable to read piped input from '${bcftools_command}': $!\n");

my @deletion; # [chr, start position, end position]
my ($prevchr, $chr, $pos, $ref, $alt) = ('.', '.', 0, '.', '.');
while (<$vcf>) {
	chomp;
	if ($_ =~ /^#/) {
		print STDOUT "$header_command\n" if ($_ =~ /^#CHROM/);
		print STDOUT "$_\n";
		next;
	}

	if ($_ =~ /^(\S+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)/) {
		($chr, $pos, $ref, $alt) = ($1, $2, $3, $4);
	} else {
		die("Incorrect VCF format:\n$_\n");
	}

	undef @deletion if ($chr ne $prevchr);
	die("Multiallelic reference found at ${chr}:${pos}: $ref") if $ref =~ /,/;
	my $reflen = length($ref);
	my $nalt = 0;
	my $indel = ($reflen > 1) ? 1:0;
	foreach my $allele (split(',',$alt)) {
		# find deleted region wrt reference
		my $altlen = length($allele);
		$indel = 1 if ($altlen > 1);
		if ($altlen < $reflen) {
			# deletion
			my $start = $pos + 1; # first deleted position
			my $end = $pos + ($reflen - $altlen); # last deleted base
			if (!@deletion) {
				@deletion = ($chr, $start, $end);
			} elsif ($chr eq $deletion[0] && $start >= $deletion[1] && $end > $deletion[2]) {
				$deletion[1] = $start if ($start > $deletion[2]);
				$deletion[2] = $end;
			}
		}
		$nalt++;
	}

	# check if site is multiallelic or is an indel
	if ($nalt == 1 && !$indel && $alt =~ /[A|C|G|T|\.]{1}/g && $ref =~ /[A|C|G|T]{1}/g) {
		# print site if not spanning a deletion
		print "$_\n" unless (@deletion && $chr eq $deletion[0] && $pos >= $deletion[1] && $pos <= $deletion[2]);
	}

	$prevchr = $chr;
}

close $vcf;

exit;
