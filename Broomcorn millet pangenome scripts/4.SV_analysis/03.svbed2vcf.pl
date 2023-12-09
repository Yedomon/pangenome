#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
my %opt;
GetOptions (\%opt, "sv:s", "svtype:s", "help");
my $help=<<USAGE;
This script is edited to combine sv
Edit by liuyang, 2021-10-06

--sv : All.del.sort.51bp.nogap.bed.combine.txt
#svID   CHROM   POS     END     SVTYPE  SVLEN   MERGE_SAMPLES_NUM       MERGE_SAMPLES
1000001 LGMI4.Chr01     95083   95139   DEL     56      1       BC310.Chr01_108692_108692_56
1000002 LGMI4.Chr01     128681  131953  DEL     3272    2       BC382.Chr01_135253_135253_3272,BC434.Chr01_132856_132856_3272

--sv : 
#svID   CHROM   POS     END     SVTYPE  SVLEN   MERGE_SAMPLES_NUM       MERGE_SAMPLES   INS_SEQ REF_SEQ
1000001 LGMI4.Chr01     10083   10083   INS     51      1       BC382.Chr01_16641_16692_51      GTATATTATTTGATAGAGAATGCACCCGGTGTTGATATCATATATGACCCTT    T

--svtype : "DEL"/"INS"/"INV"
Run: perl $0 -sv All.del.sort.51bp.nogap.bed.merge.txt -svtype "DEL" > out

USAGE
if ($opt{help} or keys %opt < 1) {
	print "$help";
	exit();
}

my $sv = $opt{sv};
my $svtype = $opt{svtype};

my $anno = 
"##fileformat=VCFv4.2\n" . 
"##fileDate=2021-10-26\n" . 
"##source=Liuyang_20211026\n" . 
"##reference=LGMI4\n" . 
"##contig=<ID=LGMI4.Chr01,length=69183459,assembly=LGMI4,md5=XXXXXXXX,species=\"Guzi\",taxonomy=x>\n" . 
"##contig=<ID=LGMI4.Chr02,length=61153219,assembly=LGMI4,md5=XXXXXXXX,species=\"Guzi\",taxonomy=x>\n" . 
"##contig=<ID=LGMI4.Chr03,length=57970102,assembly=LGMI4,md5=XXXXXXXX,species=\"Guzi\",taxonomy=x>\n" . 
"##contig=<ID=LGMI4.Chr04,length=56286655,assembly=LGMI4,md5=XXXXXXXX,species=\"Guzi\",taxonomy=x>\n" . 
"##contig=<ID=LGMI4.Chr05,length=54126031,assembly=LGMI4,md5=XXXXXXXX,species=\"Guzi\",taxonomy=x>\n" . 
"##contig=<ID=LGMI4.Chr06,length=52839179,assembly=LGMI4,md5=XXXXXXXX,species=\"Guzi\",taxonomy=x>\n" . 
"##contig=<ID=LGMI4.Chr07,length=51234605,assembly=LGMI4,md5=XXXXXXXX,species=\"Guzi\",taxonomy=x>\n" . 
"##contig=<ID=LGMI4.Chr08,length=48259421,assembly=LGMI4,md5=XXXXXXXX,species=\"Guzi\",taxonomy=x>\n" . 
"##contig=<ID=LGMI4.Chr09,length=45112342,assembly=LGMI4,md5=XXXXXXXX,species=\"Guzi\",taxonomy=x>\n" . 
"##contig=<ID=LGMI4.Chr10,length=44648547,assembly=LGMI4,md5=XXXXXXXX,species=\"Guzi\",taxonomy=x>\n" . 
"##contig=<ID=LGMI4.Chr11,length=43177482,assembly=LGMI4,md5=XXXXXXXX,species=\"Guzi\",taxonomy=x>\n" . 
"##contig=<ID=LGMI4.Chr12,length=42466157,assembly=LGMI4,md5=XXXXXXXX,species=\"Guzi\",taxonomy=x>\n" . 
"##contig=<ID=LGMI4.Chr13,length=40720392,assembly=LGMI4,md5=XXXXXXXX,species=\"Guzi\",taxonomy=x>\n" . 
"##contig=<ID=LGMI4.Chr14,length=38490750,assembly=LGMI4,md5=XXXXXXXX,species=\"Guzi\",taxonomy=x>\n" . 
"##contig=<ID=LGMI4.Chr15,length=34360906,assembly=LGMI4,md5=XXXXXXXX,species=\"Guzi\",taxonomy=x>\n" . 
"##contig=<ID=LGMI4.Chr16,length=33613985,assembly=LGMI4,md5=XXXXXXXX,species=\"Guzi\",taxonomy=x>\n" . 
"##contig=<ID=LGMI4.Chr17,length=32993148,assembly=LGMI4,md5=XXXXXXXX,species=\"Guzi\",taxonomy=x>\n" . 
"##contig=<ID=LGMI4.Chr18,length=32237550,assembly=LGMI4,md5=XXXXXXXX,species=\"Guzi\",taxonomy=x>\n" . 
"##INFO=<ID=BKPTID,Number=.,Type=String,Description=\"ID of the assembled alternate allele in the assembly file\">\n" . 
"##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">\n" . 
"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n" . 
"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n" . 
"##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description=\"Length of base pair identical micro-homology at event breakpoints\">\n" . 
"##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical micro-homology at event breakpoints\">\n" . 
"##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n" . 
"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n" . 
"##ALT=<ID=DEL,Description=\"Deletion\">\n" . 
"##ALT=<ID=DEL:ME:ALU,Description=\"Deletion of ALU element\">\n" . 
"##ALT=<ID=DEL:ME:L1,Description=\"Deletion of L1 element\">\n" . 
"##ALT=<ID=DUP,Description=\"Duplication\">\n" . 
"##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">\n" . 
"##ALT=<ID=INS,Description=\"Insertion of novel sequence\">\n" . 
"##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element\">\n" . 
"##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">\n" . 
"##ALT=<ID=INV,Description=\"Inversion\">\n" . 
"##ALT=<ID=CNV,Description=\"Copy number variable region\">\n" . 
"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" . 
"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n" . 
"##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">\n" . 
"##FORMAT=<ID=CNQ,Number=1,Type=Float,Description=\"Copy number genotype quality for imprecise events\">";

my $head = 
"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tBC027\tBC040\tBC048\tBC100\tBC136\tBC170\tBC188\tBC204\tBC217\tBC235\tBC246\tBC264\tBC292\tBC310\tBC311\tBC315\tBC328\tBC332\tBC350\tBC357\tBC362\tBC382\tBC396\tBC404\tBC407\tBC418\tBC426\tBC434\tBC476\tBC479\tBC494\tBC498";

print "$anno\n$head\n";

my %hash;
my @allAcc = ("BC027","BC040","BC048","BC100","BC136","BC170","BC188","BC204","BC217","BC235","BC246","BC264","BC292","BC310","BC311","BC315","BC328","BC332","BC350","BC357","BC362","BC382","BC396","BC404","BC407","BC418","BC426","BC434","BC476","BC479","BC494","BC498");

######
open IN_SV, $sv or die $!;
while (<IN_SV>) {
	chomp;
	next if $_ =~ /^#/;
	
	if ($svtype =~ /DEL|INV/) {  # for DEL/INV
		my ($svID,$refChr,$refStart,$refEnd,$svType,$svLen,$mergeNum,$samples) = split;
		## get accessions with alternative genotype
		my @altAcc;
		push (@altAcc, $1) while $samples =~ /(BC\d+)/mg;
		## get accessions with reference genotype
		my %hash_altAcc = map {$_ => 1} @altAcc;
		my @refAcc = grep {! $hash_altAcc{$_}} @allAcc;
		
		$hash{$_} = "1|1:22" foreach(@altAcc);
		$hash{$_} = "0|0:22" foreach(@refAcc);
		
		my $vcf = "";
		foreach my $acc (sort { $a cmp $b } @allAcc) {
			if ($vcf eq "") {  # the first element
				$vcf = $hash{$acc};
			} else {
				$vcf = $vcf . "\t" . $hash{$acc};
			}
		}
		print "$refChr\t$refStart\tSV$refStart\tN\t<DEL>\t22\tPASS\tSVTYPE=DEL;END=$refEnd;SVLEN=$svLen\tGT:GQ\t$vcf\n" if $svtype =~ /DEL/;
		print "$refChr\t$refStart\tSV$refStart\tN\t<INV>\t22\tPASS\tSVTYPE=INV;END=$refEnd;SVLEN=$svLen\tGT:GQ\t$vcf\n" if $svtype =~ /INV/;
	}
	
	elsif ($svtype =~ /INS/) {  # for INS
		my ($svID,$refChr,$refStart,$refEnd,$svType,$svLen,$mergeNum,$samples,$insSeq,$refSeq) = split;
		## get accessions with alternative genotype
		my @altAcc;
		push (@altAcc, $1) while $samples =~ /(BC\d+)/mg;
		## get accessions with reference genotype
		my %hash_altAcc = map {$_ => 1} @altAcc;
		my @refAcc = grep {! $hash_altAcc{$_}} @allAcc;
		
		$hash{$_} = "1|1:22" foreach(@altAcc);
		$hash{$_} = "0|0:22" foreach(@refAcc);
		
		my $vcf = "";
		foreach my $acc (sort { $a cmp $b } @allAcc) {
			if ($vcf eq "") {  # the first element
				$vcf = $hash{$acc};
			} else {
				$vcf = $vcf . "\t" . $hash{$acc};
			}
		}
		print "$refChr\t$refStart\tSV$refStart\t$refSeq\t$insSeq\t22\tPASS\tSVTYPE=INS;SVLEN=$svLen\tGT:GQ\t$vcf\n";
	}

}

close IN_SV;
