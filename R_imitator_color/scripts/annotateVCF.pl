#!/usr/bin/perl

# annotateVCF.pl

use warnings;
use strict;
use IO::Zlib;
use Getopt::Long;

die(qq/
annotateVCF.pl [annotation] [vcf]

Annotations:
--BandMisMap     ngsParalog output file for banded population
--StripeMisMap   ngsParalog output file for striped population
--AdmixMisMap    ngsparalog output file for admixed population
--nsample        TSV file where each row has (1) sample name in VCF, (2) pop ID 'band','stripe','admix','variabilis', or 'summersi' 
\n/) if (!@ARGV || scalar @ARGV < 2 || $ARGV[0] eq '-help' || $ARGV[0] eq '--help');

my $command = "##annotateVCFCommand=" . join(' ',@ARGV) . "\n";

my $vcf = pop @ARGV if (@ARGV);
die("Unable to locate VCF $vcf: $!\n") if (!-f $vcf);

my ($band_mismap, $stripe_mismap, $admix_mismap, $popfile) = (undef, undef, undef, undef);

GetOptions('BandMisMap=s' => \$band_mismap, 'StripeMisMap=s' => \$stripe_mismap, 'AdmixMisMap=s' => \$admix_mismap,
'nsample=s' => \$popfile);

die("Unable to locate banded ngsParalog file $band_mismap: $!\n") if ($band_mismap && !-f $band_mismap);
die("Unable to locate striped ngsParalog file $stripe_mismap: $!\n") if ($stripe_mismap && !-f $stripe_mismap);
die("Unable to locate admixed ngsParalog file $admix_mismap: $!\n") if ($admix_mismap && !-f $admix_mismap);
die("Unable to locate sample population file $popfile: $!\n") if ($popfile && !-f $popfile);

my $countdata = 0;

# parse popID file
my (%popmap, %popidx);
my @pops;

if ($popfile) {
	open(my $popfh, '<', $popfile) or die("Unable to open --ns population file $popfile: $!\n");

	while (<$popfh>) {
		chomp;
		my ($sample, $pop) = ($1, $2) if ($_ =~ /^(\S+)\s+(\S+)/);
		$popmap{$sample} = $pop;
	}
	close $popfh;

	chomp(my @vcfhead = split(/\s+/, my $head = `bcftools view -h $vcf | tail -n -1`));

	for (my $i=9; $i <= $#vcfhead; $i++) {
		if (exists $popmap{$vcfhead[$i]}) {
			my $pop = $popmap{$vcfhead[$i]};
			push @pops, $pop if !exists($popidx{$pop});
			push @{$popidx{$pop}}, $i;
		} else {
			print STDERR "WARNING: $vcfhead[$i] not in $popfile\n";
		}
	}

	$countdata = 1;
}

# parse ngsParalog files

my (%lr, $spcode, @spord);
if ($band_mismap) {
	print STDERR "Parsing banded ngsParalog LRs\n";
	$spcode = 'b';
	open(my $fh, '<', $band_mismap) or die("Unable to open banded ngsParalog file $band_mismap: $!\n");
	parseLR(\%lr, $fh, $spcode);
	close $fh;
	push @spord, $spcode;
}

if ($stripe_mismap) {
        print STDERR "Parsing striped ngsParalog LRs\n";
        $spcode = 's';
        open(my $fh, '<', $stripe_mismap) or die("Unable to open striped ngsParalog file $stripe_mismap: $!\n");
        parseLR(\%lr, $fh, $spcode);
        close $fh;
        push @spord, $spcode;
}

if ($admix_mismap) {
        print STDERR "Parsing admixed ngsParalog LRs\n";
        $spcode = 'a';
        open(my $fh, '<', $admix_mismap) or die("Unable to open admixed ngsParalog file $band_mismap: $!\n");
        parseLR(\%lr, $fh, $spcode);
        close $fh;
        push @spord, $spcode;
}

# open VCF

my $vcffh;
open($vcffh, '<', $vcf) or die("Couldn't open VCF file $vcf: $!\n");
sysread $vcffh, my $magic, 2;
close $vcffh;

if ($magic eq "\x1f\x9d") {
	open($vcffh, '<', $vcf);
} elsif ($magic eq "\x1f\x8b") {
	$vcffh = new IO::Zlib;
	$vcffh->open($vcf, "rb");
} else {
	die("Unrecognized compression type for $vcf\n");
}

# add annotations to VCF header
print STDERR "Processing $vcf\n";

my $vcfline = <$vcffh>;
print STDOUT $vcfline;
print STDOUT $command;
my %seenanno;
my @vcfhead;
while (<$vcffh>) {
	$seenanno{NBand} = 1 if ($_ =~ /ID=NBand,/);
	$seenanno{NStripe} = 1 if ($_ =~ /ID=NStripe,/);
	$seenanno{NAdmix} = 1 if ($_ =~ /ID=NAdmix,/);
	$seenanno{NVariabilis} = 1 if ($_ =~ /ID=NVariabilis,/);
	$seenanno{NSummersi} = 1 if ($_ =~ /ID=NSummersi,/);
	$seenanno{BandMisMap} = 1 if ($_ =~ /ID=BandMisMap,/);
	$seenanno{StripeMisMap} = 1 if ($_ =~ /ID=StripeMisMap,/);
	$seenanno{AdmixMisMap} = 1 if ($_ =~ /ID=AdmixMisMap,/);
	push @vcfhead, $_;
	last if ($_ =~ /^#CHROM/);
}

if ($countdata) {
	print STDOUT "##INFO=<ID=NBand,Number=1,Type=Integer,Description=\"Number of banded Ranitomeya imitator with data\">\n" if (exists $popidx{band} && !exists $seenanno{NBand});
	print STDOUT "##INFO=<ID=NStripe,Number=1,Type=Integer,Description=\"Number of striped Ranitomeya imitator with data\">\n" if (exists $popidx{stripe} && !exists $seenanno{NStripe});
	print STDOUT "##INFO=<ID=NAdmix,Number=1,Type=Integer,Description=\"Number of admixed Ranitomeya imitator with data\">\n" if (exists $popidx{admix} && !exists $seenanno{NAdmix});
	print STDOUT "##INFO=<ID=NVariabilis,Number=1,Type=Integer,Description=\"Number of Ranitomeya variabilis with data\">\n" if (exists $popidx{variabilis} && !exists $seenanno{NVariabilis});
	print STDOUT "##INFO=<ID=NSummersi,Number=1,Type=Integer,Description=\"Number of Ranitomeya summersi with data\">\n" if (exists $popidx{summersi} && !exists $seenanno{NSummersi});
}
print STDOUT "##INFO=<ID=BandMisMap,Number=1,Type=Float,Description=\"Likelihood ratio of mismapping reads for banded imitator\">\n" if ($band_mismap && !exists $seenanno{BandMisMap});
print STDOUT "##INFO=<ID=StripeMisMap,Number=1,Type=Float,Description=\"Likelihood ratio of mismapping reads for striped imitator\">\n" if ($stripe_mismap && !exists $seenanno{StripeMisMap});
print STDOUT "##INFO=<ID=AdmixMisMap,Number=1,Type=Float,Description=\"Likelihood ratio of mismapping reads for admixed imitator\">\n" if ($admix_mismap && !exists $seenanno{AdmixMisMap});

foreach (@vcfhead) {
	print STDOUT $_;
}

# add annotations to sites
$" = "\t";
my @fmtidx = (0, -1); # [GT, DP]

while ($vcfline = <$vcffh>) {
	chomp($vcfline);
	my @tok = split(/\s+/, $vcfline);
	formatFields(\$tok[8], \@fmtidx);

	if ($countdata) {
		if ($fmtidx[1] > 0) {
			if ($popidx{band}) {
				my $nband = countNumData(\@{$popidx{band}}, \@tok, $fmtidx[1]);
				$tok[7] =~ s/NBand=\d+//;
				$tok[7] .= ";NBand=$nband";
			}
			if ($popidx{stripe}) {
				my $nstripe = countNumData(\@{$popidx{stripe}}, \@tok, $fmtidx[1]);
				$tok[7] =~ s/NStripe=\d+//;
				$tok[7] .= ";NStripe=$nstripe";
			}
			if ($popidx{admix}) {
				my $nadmix = countNumData(\@{$popidx{admix}}, \@tok, $fmtidx[1]);
				$tok[7] =~ s/NAdmix=\d+//;
				$tok[7] .= ";NAdmix=$nadmix";
			}
			if ($popidx{variabilis}) {
				my $nvariabilis = countNumData(\@{$popidx{variabilis}}, \@tok, $fmtidx[1]);
				$tok[7] =~ s/NVariabilis=\d+//;
				$tok[7] .= ";NVariabilis=$nvariabilis";
			}
			if ($popidx{summersi}) {
				my $nsummersi = countNumData(\@{$popidx{summersi}}, \@tok, $fmtidx[1]);
				$tok[7] =~ s/NSummersi=\d+//;
				$tok[7] .= ";NSummersi=$nsummersi";
			}
		} else {
			print STDERR "WARNING: No FORMAT/DP at $tok[0] $tok[1]\n";
		}
	}

	if (%lr) {
		$tok[7] =~ s/BandMisMap=[^;]+//;
		$tok[7] =~ s/StripeMisMap=[^;]+//;
		$tok[7] =~ s/AdmixMisMap=[^;]+//;
		if (exists $lr{$tok[0]}{$tok[1]}) {
			$tok[7] .= ";BandMisMap=$lr{$tok[0]}{$tok[1]}{b}" if exists($lr{$tok[0]}{$tok[1]}{b});
			$tok[7] .= ";StripeMisMap=$lr{$tok[0]}{$tok[1]}{s}" if exists($lr{$tok[0]}{$tok[1]}{s});
			$tok[7] .= ";AdmixMisMap=$lr{$tok[0]}{$tok[1]}{a}" if exists($lr{$tok[0]}{$tok[1]}{a});
		}
	}

	$tok[7] =~ s/;{2,}/;/g;
	print "@tok\n";
}

close $vcffh;

exit;

sub countNumData {
	my ($idx, $vcf, $dpfield) = @_;
	my $c = 0;
	foreach my $geno (@$vcf[@$idx]) {
		$c++ if (split(':', $geno))[$dpfield] > 0;
	}
	return $c;
}

sub formatFields {
	my ($format, $fmtidx) =@_;
	for (my $i=1; $i<=$#$fmtidx; $i++) {
		$fmtidx[$i] = -1;
	}
	my $j = 0;
	foreach(split(':',$$format)) {
		$$fmtidx[1] = $j if ($_ eq 'DP');
		$j++;
	}
}

sub parseLR {
	my ($lr, $fh, $spcode) = @_;

	while (<$fh>) {
		my @tok = split(/\s+/, $_);
		$$lr{$tok[0]}{$tok[1]}{$spcode} = $tok[4];
	}
}
