#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum0);
use Statistics::Test::WilcoxonRankSum;  #https://metacpan.org/pod/Statistics::Test::WilcoxonRankSum
#use Data::Dumper;

my %opt;
GetOptions (\%opt, "func:s", "exp:s", "anno:s", "help");
my $help=<<USAGE;
This script is edited to get a table of "sv overlapping genes (with rice functional genes if exists)". 
Modified for the new RNAseq data.
Edited by liuyang, 2022-12-18.

--func : Longmi_gene2rice-.txt
longmi000014    LOC_Os01g13770:OsTPT1 LOC_Os05g15160:OsTPT2
longmi000073    LOC_Os01g12890:CCP1|DFO1|EMF1|OsEMF1

--exp : ALL.RPKM.NEW.txt
Leaf_BC027-1    longmi028140    369     0       20812292         0
Leaf_BC027-1    longmi028141    369     0       20812292         0

--anno : ALL.PAV.anno.txt.edit.newid
LGMI4.Chr01_100001      LGMI4.Chr01     95084   95139   56      DEL     intron  longmi010923    BC310.Chr01_108692_108692_56
LGMI4.Chr01_100001      LGMI4.Chr01     95084   95139   56      DEL     downTAA longmi010923    BC310.Chr01_108692_108692_56

Run: perl $0 -func Longmi_gene2rice-.txt -exp ALL.FPKM.NEW.txt -anno ALL.PAV.anno.txt.edit.newid > out

USAGE
if ($opt{help} or keys %opt < 1) {
	print "$help";
	exit();
}

my $func = $opt{func};
my $exp = $opt{exp};
my $anno = $opt{anno};

# header
print "#longmiGene\tsvPos\tsvClass\toverlapFeature\taccNum\talt/ref(G)\tfc(G)\tpvalue(G)\talt/ref(R)\tfc(R)\tpvalue(R)\taccSV\triceGene\n";

######
open IN_FUNC, $func or die $!;
my %FUNC;
while (<IN_FUNC>) {
	chomp;
	# longmi000014    LOC_Os01g13770:OsTPT1 LOC_Os05g15160:OsTPT2
	my ($lm, $rice) = split/\t/;
	$FUNC{$lm} = $rice;  # $FUNC{longmi000026} = LOC_Os01g13540:OsNIN-like3
}

######
open IN_EXP, $exp or die $!;
my %FPKM;
while (<IN_EXP>) {
	chomp;
	# Leaf_BC027-1    longmi028140    369     0       20812292         0
	my ($acc,$lmgene,$fpkm) = (split)[0,1,5];
	$FPKM{$acc}{$lmgene} = $fpkm;  # $FPKM{Leaf_BC027-1}{longmi028140} = 0;
}

# foreach my $key1 (keys %FPKM) {
	# my $hash2 = $FPKM{$key1};
	# foreach my $key2 (keys %$hash2) {
		# print "$key1\t$key2\t$hash2->{$key2}\n";
	# }
# }
#print "test\t$FPKM{G434A}{longmi025893}\n";

######
open IN_ANNO, $anno or die $!;
my %store;  # store printing lines
my $final_feature = "";  # for printing
while (<IN_ANNO>) {
	chomp;
	my ($svID,$refChr,$refStart,$refEnd,$svLen,$svClass,$feature,$overlapGene,$samples) = split;  # $samples= BC310.Chr01_108692_108692_56,...

	if ($feature =~ /CDS|intron|upATG|downTAA/) {
		my $sv = $refChr . "_" . $refStart . "_" . $refEnd . "_" . $svLen;  # LGMI.Chr05_78169_92491_14323
		my $gene = $1 if $overlapGene =~ /(longmi\d+)/;
		
		# counting SV number in each group ----------------
		my (@C1, @C2, @C3, @W);
		push (@C1, $1) while $samples =~ /BC040|BC100/mg;
		  my $C1 = $#C1 + 1;
		push (@C2, $1) while $samples =~ /BC311|BC494|BC498/mg;
		  my $C2 = $#C2 + 1;
		push (@C3, $1) while $samples =~ /BC475|BC477|BC244|BC027|BC048|BC136|BC170|BC188|BC204|BC217|BC235|BC264|BC292|BC310|BC315|BC332|BC350|BC360|BC362/mg;
		  my $C3 = $#C3 + 1;
		push (@W, $1) while $samples =~ /BC328|BC382|BC398|BC404|BC407|BC418|BC426|BC434/mg;
		  my $W = $#W + 1;
		
	### calculating mean FPKM and foldchange and pvalue ----------------
	### two tissues, three replicates.
	# all accessions
		my @allAcc = ("BC027","BC040","BC048","BC100","BC136","BC170","BC188","BC204","BC217","BC235","BC244","BC264","BC292","BC310","BC311","BC315","BC328","BC332","BC350","BC360","BC362","BC382","BC398","BC404","BC407","BC418","BC426","BC434","BC477","BC475","BC494","BC498");
	# accessions with alternative genotype
		my @altAcc = $samples =~ /BC\d+/mg;
		my $altNum = $#altAcc + 1;  # alt accession number
	# get accessions with reference genotype
		my %hash_altAcc = map {$_ => 1} @altAcc;
		my @refAcc = grep {! $hash_altAcc{$_}} @allAcc;
		my $refNum = $#refAcc + 1;  # ref accession number

	# calculating altAcc mean FPKM (between replicates) ----------------
		my (@array_altFPKMmeanG, @array_altFPKMmeanR);
		foreach (@altAcc) {
		# make RNAseq sample IDs
			my $No = $1 if $_ =~ /(BC\d+)/;
			my $G1 = "Leaf_".$No."-1";  # Leaf_BC027-1
			my $G2 = "Leaf_".$No."-2";
			my $G3 = "Leaf_".$No."-3";
			my $R1 = "Root_".$No."-1";
			my $R2 = "Root_".$No."-2";
			my $R3 = "Root_".$No."-3";
		# array, store the mean FPKM of each alt-accession
			my $altFPKMmeanGG1 = ( $FPKM{$G1}{$gene} + $FPKM{$G2}{$gene} + $FPKM{$G3}{$gene}) / 3;
			my $altFPKMmeanRR1 = ( $FPKM{$R1}{$gene} + $FPKM{$R2}{$gene} + $FPKM{$R3}{$gene}) / 3;
			push @array_altFPKMmeanG, $altFPKMmeanGG1;
			push @array_altFPKMmeanR, $altFPKMmeanRR1;
		}
		
		# calculating refAcc mean FPKM (between replicates) ----------------
		my (@array_refFPKMmeanG, @array_refFPKMmeanR);
		foreach (@refAcc) {
		# make RNAseq id
			my $No = $1 if $_ =~ /(BC\d+)/;
			my $G1 = "Leaf_".$No."-1";  # Leaf_BC027-1
			my $G2 = "Leaf_".$No."-2";
			my $G3 = "Leaf_".$No."-3";
			my $R1 = "Root_".$No."-1";
			my $R2 = "Root_".$No."-2";
			my $R3 = "Root_".$No."-3";
		# array, store the FPKM of each ref-accession
			my $refFPKMmeanGG1 = ( $FPKM{$G1}{$gene} + $FPKM{$G2}{$gene} + $FPKM{$G3}{$gene}) / 3;
			my $refFPKMmeanRR1 = ( $FPKM{$R1}{$gene} + $FPKM{$R2}{$gene} + $FPKM{$R3}{$gene}) / 3;
			push @array_refFPKMmeanG, $refFPKMmeanGG1;
			push @array_refFPKMmeanR, $refFPKMmeanRR1;
		}
		
		# calculating altAcc mean FPKM (between accessions) ----------------
		my $altFPKMmeanG = 0.00001;  # can't be 0
		my $altFPKMmeanR = 0.00001;
		my $altFPKMsumG = sum0 @array_altFPKMmeanG;
		my $altFPKMsumR = sum0 @array_altFPKMmeanR;
		if ($altNum > 0) {
			$altFPKMmeanG = $altFPKMsumG / $altNum;
			$altFPKMmeanR = $altFPKMsumR / $altNum;
		}
		
		# calculating refAcc mean FPKM (between accession) ----------------
		my $refFPKMmeanG = 0.00001;
		my $refFPKMmeanR = 0.00001;
		my $refFPKMsumG = sum0 @array_refFPKMmeanG;
		my $refFPKMsumR = sum0 @array_refFPKMmeanR;
		if ($refNum > 0) {
			$refFPKMmeanG = $refFPKMsumG / $refNum;
			$refFPKMmeanR = $refFPKMsumR / $refNum;
		}
		
		# # log2 fold change ----------------
		# my $log2fcG = log2($altFPKMmeanG + 1) - log2($refFPKMmeanG + 1);
		# my $log2fcR = log2($altFPKMmeanR + 1) - log2($refFPKMmeanR + 1);
		
		# fold change
		my $fcG = ($altFPKMmeanG + 1) / ($refFPKMmeanG + 1);
		my $fcR = ($altFPKMmeanR + 1) / ($refFPKMmeanR + 1);
		#$fcG = 1 / $fcG if $fcG < 1;
		#$fcR = 1 / $fcR if $fcR < 1;
		
		# pvalue (Mann-Whitney U test) ----------------
		my $probG = 1;
		my $probR = 1;
		if ($altFPKMsumG > 0 && $altFPKMsumR > 0 && $refFPKMsumG > 0 && $refFPKMsumR > 0) {
			my $wilcox_test_G = Statistics::Test::WilcoxonRankSum->new();
			my $wilcox_test_R = Statistics::Test::WilcoxonRankSum->new();
			$wilcox_test_G->load_data(\@array_altFPKMmeanG,\@array_refFPKMmeanG);
			$wilcox_test_R->load_data(\@array_altFPKMmeanR,\@array_refFPKMmeanR);
			$probG = $wilcox_test_G->probability();
			$probR = $wilcox_test_R->probability();
		}
		
		# format print
		$altFPKMmeanG = sprintf "%.4f", $altFPKMmeanG;
		$altFPKMmeanR = sprintf "%.4f", $altFPKMmeanR;
		$refFPKMmeanG = sprintf "%.4f", $refFPKMmeanG;
		$refFPKMmeanR = sprintf "%.4f", $refFPKMmeanR;
		$fcG = sprintf "%.4f", $fcG;
		$fcR = sprintf "%.4f", $fcR;
		$probG = sprintf "%f", $probG;
		$probR = sprintf "%f", $probR;
		
		# printing line to the hash ----------------
		if (exists $store{$gene}) {
			if ($final_feature !~ /$feature/) {
				$final_feature = $final_feature . " " . $feature;
			}
		} else {
			$final_feature = $feature;
		}
		
		$store{$gene} = "$gene\t$sv\t$svClass\t$final_feature\tC1:$C1 C2:$C2 C3:$C3 W:$W\t$altFPKMmeanG:$refFPKMmeanG\t$fcG\t$probG\t$altFPKMmeanR:$refFPKMmeanR\t$fcR\t$probR\t$samples";
		
		# # feature coding
		# $store{$gene} =~ s/CDS/1/;
		# $store{$gene} =~ s/intron/0/;
		# $store{$gene} =~ s/upATG/-2000/;
		# $store{$gene} =~ s/downTAA/2000/;

	}
}

######
foreach my $k (sort {$a cmp $b} keys %store) {
	if ($FUNC{$k}) {
		print "$store{$k}\t$FUNC{$k}\n";
	} else {
		print "$store{$k}\tNA\n";
	}
}

######
close IN_ANNO;

#------------------
sub log2 {
   my $n = shift;
   return log($n)/log(2);
}
