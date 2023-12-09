#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
my %opt;
GetOptions (\%opt, "gff:s", "sv:s", "svtype:s", "help");
my $help=<<USAGE;
This script is edited to filer SVs which overlapping CEN region. 
Edited by liuyang, 2022-01-20

--gff : ALL.cent.gff
Chr01   BC027   CEN     7364496 7365259 .       .       .       .
Chr01   BC027   CEN     11837819        11837974        .       .       .       .

--sv : All.del.sort.51bp.nogap.bed
LGMI4.Chr01  95083  95139  DEL_PAV  BC310.Chr01  108692  108692  56  +  54396  40873  AGGGACACCATGCCCTTTAG  ATATATATATATAACTCGTT  AGGGACACCATGCCCTTTAGATATATATATATAACTCGTT  ATATATATATATATATATATATATATATATATATATATATATATATATATATATATA  A

--sv : All.ins.sort.51bp.nogap.bed
BC382.Chr01  16641  16692  DEL_PAV  LGMI4.Chr01  10083  10083  51  +  3  6137  AGTAGATCGTTGGCCTGCAG  TTGGATTCCCCTCTCACCCC  GCAGTAGATCGTTGGCCCAGTTGGATTCCCCTCTCACCCC  GTATATTATTTGATAGAGAATGCACCCGGTGTTGATATCATATATGACCCTT  T

--svtype : "DEL"/"INS"

Note: Time-consuming: ~10 miniutes.

Run: perl $0 -gff ALL.cen.gff -sv All.del.sort.51bp.nogap.bed -svtype "DEL" > out

USAGE
if ($opt{help} or keys %opt < 1) {
	print "$help";
	exit();
}

my $gff = $opt{gff};
my $sv = $opt{sv};
my $type = $opt{svtype};

######
open IN_GFF, $gff or die $!;
my @lines = <IN_GFF>;

######
open IN_SV, $sv or die $!;
while (<IN_SV>) {
	next if /^#/;
	chomp;
	
	my $line = $_;
	my $flag = 0;  # for determining to filter or not
	
	### for DEL
	if ($type =~ /DEL/) {
		my ($refChr,$refStart,$refEnd) = (split)[0,1,2];
		$refStart = $refStart + 1;  # bed -> gff
		my $ref = "Pmlongmi4" if $refChr =~ /LGMI4/;
		my $chrRef = $1 if $refChr =~ /(Chr\d+)/;
		# check if overlap
		my ($refMid, $refHalfLen);
		$refMid = abs($refEnd - $refStart) / 2 + $refStart;  # 'gff' format positon
		$refHalfLen = abs($refEnd - $refStart + 1) / 2;
		#
		foreach (@lines) {
			chomp;
			my ($start,$end) = (split)[3,4];  # 'gff' format positon
			if ($_ =~ /$ref/ && $_ =~ /$chrRef/) {
				if ( ($refMid >= $start and $refMid <= $end) or (abs($refMid - $start) <= $refHalfLen) or (abs($refMid - $end) <= $refHalfLen) ) {
					$flag = 1;  # if overlap
					next;
				}
			}
		}
	}
	
	### for INS
	if ($type =~ /INS/) {
		my ($accChr,$accStart,$accEnd) = (split)[0,1,2];
		$accStart = $accStart + 1;  # bed -> gff
		my $acc = $1 if $accChr =~ /(BC\d+)/;
		my $chrAcc = $1 if $accChr =~ /(Chr\d+)/;
		# check if overlap
		my ($accMid, $accHalfLen);
		$accMid = abs($accEnd - $accStart) / 2 + $accStart;  # 'gff' format positon
		$accHalfLen = abs($accEnd - $accStart + 1) / 2;
		#
		foreach (@lines) {
			chomp;
			my ($start,$end) = (split)[3,4];  # 'gff' format positon
			if ($_ =~ /$acc/ && $_ =~ /$chrAcc/) {
				if ( ($accMid >= $start and $accMid <= $end) or (abs($accMid - $start) <= $accHalfLen) or (abs($accMid - $end) <= $accHalfLen) ) {
					$flag = 1;  # if overlap
					next;
				}
			}
		}
	}
	
	### print
	print "$line\n" if $flag == 0;	
}

close IN_SV;
close IN_GFF;
