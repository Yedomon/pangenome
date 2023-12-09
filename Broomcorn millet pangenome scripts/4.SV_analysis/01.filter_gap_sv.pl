#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
my %opt;
GetOptions (\%opt, "bed:s", "sv:s", "svtype:s", "help");
my $help=<<USAGE;
This script is edited to filer SVs which overlapping gap region. 
Edit by liuyang, 2021-10-03

--bed : ALL.gap.bed
BC027  Chr01_RagTag  12096597  12096697
Pmlongmi4  Chr18  30754658  30755157

--sv : All.del.sort.bed
#CHROM  POS  END  ID  SVTYPE  SVLEN  MERGE_SAMPLES_NUM  MERGE_SAMPLES  QUERY  Q_POS  Q_END  FIVE  THREE
LGMI4.Chr01  20124  20126  DEL_PAV  BC264.Chr01  21453  21453  2  +  10329  20559  CGATGATTCAGTTCATTCCT  ACACACACACACACACACAC  CGATGATTCAGTTCATTCCT  ACACACACACACACACACAC  ACA  A

--sv : All.ins.sort.bed
BC382.Chr01  16636  16638  DEL_PAV  LGMI4.Chr01  10080  10080  2  +  279  3  ATAGCAGTAGATCGTTGGCC  CAGGTATATTATTTGATAGA  ATAGCAGTAGATCGTTGGCC  CAGTTGGATTCCCCTCTCAC  TGC  C

--sv : All.inv.reference.sort.bed
LGMI4.Chr01  3768955  3793941  INV  BC498.Chr01  3787214  3812176  24986  +  10126  13714  GCCTCTTGAGAATCGGCTAA  ATGGCATGCATTTTTTGGTG  GCCCCTCGAGAATCGGCTAA  ATGACATGCATTTTTTGGTG  CCTTCCAATC...ACCTTTGGAT  GGCCCTC...GAGTAAAACA

--sv : All.trl.reference.sort.bed
LGMI4.Chr01     27120951        27337465        BC136.Chr01     25622357        25919641        +       216514  297284  IntraChr        syn

--svtype : "DEL"/"INS"/"INV"/"TRL"

Note: the start is start from 0.
Note: Time-consuming: ~10 miniutes.

Run: perl $0 -bed ALL.gap.bed -sv All.del.sort.bed -svtype "DEL" > out

USAGE
if ($opt{help} or keys %opt < 1) {
	print "$help";
	exit();
}

my $bed = $opt{bed};
my $sv = $opt{sv};
my $type = $opt{svtype};

######
open IN_BED, $bed or die $!;
my @lines = <IN_BED>;

######
open IN_SV, $sv or die $!;
while (<IN_SV>) {
	next if /^#/;
	chomp;
	my ($refChr,$refStart,$refEnd,$accChr,$accStart,$accEnd,$svLen);
	if ($type =~ /DEL|INV/) {
		($refChr,$refStart,$refEnd,$accChr,$accStart,$accEnd,$svLen) = (split)[0,1,2,4,5,6,7];
	} elsif ($type =~ /INS/) {
		($accChr,$accStart,$accEnd,$refChr,$refStart,$refEnd,$svLen) = (split)[0,1,2,4,5,6,7];
	} elsif ($type =~ /TRL/) {
		($refChr,$refStart,$refEnd,$accChr,$accStart,$accEnd,$svLen) = (split)[0,1,2,3,4,5,7];
	}
	
	next unless $svLen > 50;  # filter sv length
	
	my $line = $_;
	my $flag = 0;  # for determining to filter or not
	
	## get IDs
	my $ref = "Pmlongmi4" if $refChr =~ /LGMI4/;
	my $chrRef = $1 if $refChr =~ /(Chr\d+)/;
	my $acc = $1 if $accChr =~ /(BC\d+)/;
	my $chrAcc = $1 if $accChr =~ /(Chr\d+)/;
	
	### for reference --------------------
	my ($refMid, $refHalfLen);
	if ($refStart == $refEnd) {  # it is 'bed' format position
		$refMid = $refStart + 1;  # change to 'gff' format positon
		$refHalfLen = 0;
	} else {
		$refMid = abs($refEnd - $refStart + 1) / 2 + $refStart;  # 'gff' format positon
		$refHalfLen = abs($refEnd - $refStart) / 2;
	}
	# searching sv overlap with gap
	foreach (@lines) {
		chomp;
		my ($start, $end) = (split)[2,3];  # 'bed' format positon
		if ($_ =~ /$ref/ && $_ =~ /$chrRef/) {
			if ( ($refMid >= ($start+1) and $refMid <= $end) or (abs($refMid - ($start+1)) <= $refHalfLen) or (abs($refMid - $end) <= $refHalfLen) ) {
				$flag = 1;  # if overlap
				next;
			}
			if ($type =~ /INV|TRL/) {  # just for INV, extending 10 bp each sides
				if ( (abs($refMid - ($start+1)) <= $refHalfLen + 10) or (abs($refMid - $end) <= $refHalfLen + 10) ) {
					$flag = 1;
					next;
				}
			}
		}
	}
	
	### for samples --------------------
	my ($accMid, $accHalfLen);
	if ($accStart == $accEnd) {  # it is 'bed' format position
		$accMid = $accStart + 1;  # 'gff' format positon
		$accHalfLen = 0;
	} else {
		$accMid = abs($accEnd - $accStart + 1) / 2 + $accStart;  # 'gff' format positon
		$accHalfLen = abs($accEnd - $accStart) / 2;
	}
	# searching sv overlap with gap
	foreach (@lines) {
		chomp;
		my ($start, $end) = (split)[2,3];  # bed format
		if ($_ =~ /$acc/ && $_ =~ /$chrAcc/) {
			if ( ($accMid >= ($start+1) and $accMid <= $end) or (abs($accMid - ($start+1)) <= $accHalfLen) or (abs($accMid - $end) <= $accHalfLen) ) {
				$flag = 1;  # if overlap
				next;
			}
			if ($type =~ /INV|TRL/) {  # just for INV, extending 10 bp each sides
				if ( (abs($accMid - ($start+1)) <= $accHalfLen + 10) or (abs($accMid - $end) <= $accHalfLen + 10) ) {
					$flag = 1;
					next;
				}
			}
		}
	}
	
	### print
	print "$line\n" if $flag == 0;
}

close IN_SV;
close IN_BED;
