#!/usr/bin/perl
use strict;
use Getopt::Long;
my %opt;
GetOptions (\%opt, "sv:s", "gff:s", "help");
my $help=<<USAGE;
This script is edited to annotate 'sv of longmi genome' against gene 'structure'
Edit by liuyang, 2021-11-02
Note: very slow

--sv : ALL_del_ins_inv.bed.merge.txt
#svID   CHROM   POS     END     SVTYPE  SVLEN   MERGE_SAMPLES_NUM       MERGE_SAMPLES   INS_SEQ REF_SITE
1000001 LGMI4.Chr01     95083   95139   DEL     56      1       BC310.Chr01_108692_108692_56
1000002 LGMI4.Chr01     128681  131953  DEL     3272    2       BC382.Chr01_135253_135253_3272,BC434.Chr01_132856_132856_3272

--gff : Pmlongmi4_features.gff, contain "CDS", "intron", "intergenic", "upATG" and "downTAA" in column 3
Chr01   CoGe    intergenic      1       9831    .       .       .       .
Chr01   CoGe    downTAA 9832    11831   .       -       .       ID=longmi010908_T1;Parent=longmi010908_T1
Chr01   CoGe    CDS     11832   12101   .       -       .       Parent=gene:longmi010908.mRNA1;ID=CDS:longmi010908_T1
Chr01   CoGe    intron  12102   12197   .       -       .       ID=gene:longmi010908;Parent=gene:longmi010908.mRNA1

Run: perl $0 -sv ALL_del_ins_inv.bed.merge.txt -gff Pmlongmi4_features.gff > 

USAGE
if ($opt{help} or keys %opt < 1) {
	print "$help";
	exit();
}

my $sv = $opt{sv};
my $gff = $opt{gff};

open IN_GFF, $gff or die $!;
my @lines = <IN_GFF>;

my %cdsFlag;
my %intronFlag;
open IN_SV, $sv or die $!;
while (<IN_SV>) {
	chomp;
	next if /^#/;
	my ($svID,$refChr,$refStart,$refEnd,$svType,$svLen,$samples) = (split)[0,1,2,3,4,5,7];  # Note: $refChr = LGMI4.Chr01
	my $refStart = $refStart + 1;  # bed -> gff
	
	# 
	my ($svMid, $svHalfLen);
	if ($refStart == $refEnd + 1) {
		$svMid = $refStart;
		$svHalfLen = 0;
	} else {
		$svMid = abs($refStart - $refEnd) / 2 + $refStart;  # midpoint of sv
		$svHalfLen = (abs($refStart - $refEnd) + 1) / 2;  # half length of sv
	}
	
	# check if overlap
	my ($type,$start,$end,$attri);  # for features lines
	foreach (@lines) {  # each line in features gff
		chomp;
		my ($chr,$type,$start,$end,$attri) = (split)[0,2,3,4,8];  # bed format
		next unless $refChr =~ /$chr/;  # Note: $refChr = LGMI4.Chr01
		#my $id = $1 if ($attri =~ /(longmi\d+)/);
		my $id;
		$id = $1 if ($attri =~ /(longmi\d+.mRNA\d+)/ or $attri =~ /(longmi\d+_T\d+)/);
		
		if ( ($svMid >= ($start+1) and $svMid <= $end) or (abs($svMid - ($start+1)) <= $svHalfLen) or (abs($svMid - $end) <= $svHalfLen) ) {
			if ($type =~ /CDS/) {  # CDS
				unless ($cdsFlag{$id} == 1) {
					print "$svID\t$refChr\t$refStart\t$refEnd\t$svLen\t$svType\t$type\t$id\t$samples\n";
					$cdsFlag{$id} = 1;
				}
			} elsif ($type =~ /intron/) {  # intron
				unless ($intronFlag{$id} == 1) {
					print "$svID\t$refChr\t$refStart\t$refEnd\t$svLen\t$svType\t$type\t$id\t$samples\n";
					$intronFlag{$id} = 1;
				}
			} elsif ($type =~ /up|down/) {  # upstream and downstream
				print "$svID\t$refChr\t$refStart\t$refEnd\t$svLen\t$svType\t$type\t$id\t$samples\n";
			} elsif ($type =~ /intergenic/) {  # intergenic
				print "$svID\t$refChr\t$refStart\t$refEnd\t$svLen\t$svType\t$type\t.\t$samples\n";
			}
		}
		
		#my $first = shift @lines if ($refStart > $end);
		
		#next if ($refEnd < $start);
	}
}

close IN_SV;
close IN_GFF;

