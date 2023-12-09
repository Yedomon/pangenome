#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(max);
my %opt;
GetOptions (\%opt, "sv:s", "svtype:s", "help");
my $help=<<USAGE;
This script is edited to merge sv
Edit by liuyang, 2021-10-28

!!!Note: This input file must be sorted in advance, by reference chromosome and position (and sample).
--sv : All.del.sort.51bp.nogap.nocen.bed
LGMI4.Chr01  95083  95139  DEL_PAV  BC310.Chr01  108692  108692  56  +  54396  40873  AGGGACACCATGCCCTTTAG  ATATATATATATAACTCGTT        AGGGACACCATGCCCTTTAG    ATATATATATATAACTCGTT    ATATATATATATATATATATATATATATATATATATATATATATATATATATATATA       A

--sv : All.ins.sort.51bp.nogap.nocen.bed
BC382.Chr01  16641  16692  DEL_PAV  LGMI4.Chr01  10083  10083  51  +  3  6137  AGTAGATCGTTGGCCTGCAG  TTGGATTCCCCTCTCACCCC  GCAGTAGATCGTTGGCCCAGTTGGATTCCCCTCTCACCCC  GTATATTATTTGATAGAGAATGCACCCGGTGTTGATATCATATATGACCCTT  T
Note: very slow, should run by qsub

--sv : All.inv.reference.sort.51bp.nogap.bed

--sv : All.trl.reference.sort.51bp.nogap.bed
LGMI4.Chr01  34006127  34056659  BC404.Chr09  28707008  28756572  -  50532  49564  InterChr  nonSyn

--svtype : "DEL"/"INS"/"INV"/"TRL"

Run: perl $0 -sv All.del.sort.51bp.nogap.bed -svtype "DEL" > out

USAGE
if ($opt{help} or keys %opt < 1) {
	print "$help";
	exit();
}

my $sv = $opt{sv};
open IN_SV, $sv or die $!;
my $svtype = $opt{svtype};  # different rules for DEL / INS

my %hash;  # store the line
my $svID;
my $svNo;
$svNo = 100000 if $svtype =~ /DEL/;
$svNo = 200000 if $svtype =~ /INS/;
$svNo = 300000 if $svtype =~ /INV/;
$svNo = 400000 if $svtype =~ /TRL/;
my $refChrFlag = "LGMI4.Chr01";  # for initializing
my $refStartFlag = 0;  # store previous position
my $refEndFlag = 0;  # store previous position
my ($samID, $mergeSam, $mergeNum);
# The following is for PAV_INS
my $insLenFlag = 1;
my $seqAltFlag = "X";
my $seqRefFlag = "X";
# for reading line
my ($refChr,$refStart,$refEnd,$samChr,$samStart,$samEnd);  # for all
my ($len,$seq);  # for DEL
my ($insLen,$strand,$seqAlt,$seqRef);  # for INS

### 
while (<IN_SV>) {
	chomp;
	if ($svtype =~ /DEL|INV|TRL/) {  # different rules between DEL/INV and INS
	#LGMI4.Chr01	95083	95139	DEL_PAV	BC310.Chr01	108692	108692	56	+	54396	40873	AGGGACACCATGCCCTTTAG	ATATATATATATAACTCGTT	AGGGACACCATGCCCTTTAG	ATATATATATATAACTCGTT	ATATATATATATATATATATATATATATATATATATATATATATATATATATATATA	A
		if ($svtype =~ /DEL|INV/) {
			($refChr,$refStart,$refEnd,$samChr,$samStart,$samEnd,$len,$seq) = (split)[0,1,2,4,5,6,7,15];
		} else {
			($refChr,$refStart,$refEnd,$samChr,$samStart,$samEnd,$len) = (split)[0,1,2,3,4,5,7];
		}
		
		if ($refChr =~ /$refChrFlag/) {  # make sure the chr number
			if ($refEndFlag-$refStart > 0) {  # if overlap
				my $refEndmax = ($refEndFlag > $refEnd) ? $refEndFlag : $refEnd;
				my $refEndmin = ($refEndFlag > $refEnd) ? $refEnd : $refEndFlag;
				## 
				if (($refEndmin-$refStart)/($refEndmax-$refStartFlag) > 0.9) {  # overlap ratio > 0.9
					$samID = "$samChr\_$samStart\_$samEnd\_$len";
					$mergeSam = "$mergeSam,$samID";
					$mergeNum ++;
					$svID = $refChr . "_" . $svNo;
					$hash{$svID} = "$refChr\t$refStart\t$refEnd\t$svtype\t$len\t$mergeNum\t$mergeSam\t$seq";
				} else {  # overlap ratio <= 0.9
					subPrintDEL();
				}
			} else {  # if no overlap
				subPrintDEL();
			}
		} else {  # initialize the chr number
			$svNo ++;
			$samID = "$samChr\_$samStart\_$samEnd\_$len";
			$mergeSam = $samID;
			$mergeNum = 1;
			$svID = $refChr . "_" . $svNo;
			$hash{$svID} = "$refChr\t$refStart\t$refEnd\t$svtype\t$len\t$mergeNum\t$mergeSam\t$seq";
			$refChrFlag = $refChr;
			$refStartFlag = $refStart;
			$refEndFlag = $refEnd;
		}
	}
	
	elsif ($svtype =~ /INS/) {  # rules for insertion
	#BC382.Chr01	16641	16692	DEL_PAV	LGMI4.Chr01	10083	10083	51	+	3	6137	AGTAGATCGTTGGCCTGCAG	TTGGATTCCCCTCTCACCCC	GCAGTAGATCGTTGGCCCAGTTGGATTCCCCTCTCACCCC	GTATATTATTTGATAGAGAATGCACCCGGTGTTGATATCATATATGACCCTT	T
		($samChr,$samStart,$samEnd,$refChr,$refStart,$refEnd,$insLen,$strand,$seqAlt,$seqRef) = (split)[0,1,2,4,5,6,7,8,15,16];
		
		# strand "-" to "+"
		if ($strand eq "-") {
			$seqAlt = rc($seqAlt);
			$seqRef = rc($seqRef);
		}
		
		if ($refChr =~ /$refChrFlag/) {  # make sure the chr, LGMI4.Chr01
			if (abs($refStart-$refStartFlag) < 10) {  # distance between insertion site
				
				# cal coverage
				my $shorterLen = ( length($seqAlt) < length($seqAltFlag) ) ? length($seqAlt) : length($seqAltFlag);
				#my $longerLen  = ( length($seqAlt) > length($seqAltFlag) ) ? length($seqAlt) : length($seqAltFlag);
				#my $coverage = $shorterLen / $longerLen;
				
				# align
				my $seqs = align($seqAlt, $seqAltFlag);  # $seqs = "AATTTTGGGGCCCC	AA-TTT--GGCCCC"
				
				# cal identity
				my $identity;
				my ($seq_1,$seq_2) = split(/\t/, $seqs);
				if (length($seq_1) > 0) {
					$identity = calIdentity($seq_1, $seq_2, $shorterLen);
				} else {
					$identity = 0;
				}
				
				###
				if ($identity > 0.8) {  # 
				#if ($identity > 0.5 && $coverage > 0.5) {  # 
					if ($insLen < $insLenFlag) {  # keep the shorter INS
						$samID = "$samChr\_$samStart\_$samEnd\_$insLen";
						$mergeSam = "$samID,$mergeSam";
						$mergeNum ++;
						# refresh the flags
						$refStartFlag = $refStart;
						$refEndFlag = $refEnd;
						$insLenFlag = $insLen;
						$seqAltFlag = $seqAlt;
						$seqRefFlag = $seqRef;
					} else {
						$samID = "$samChr\_$samStart\_$samEnd\_$insLen";
						$mergeSam = "$mergeSam,$samID";
						$mergeNum ++;
						# keep the last ones
						$refStart = $refStartFlag;
						$refEnd = $refEndFlag;
						$insLen = $insLenFlag;
						$seqAlt = $seqAltFlag;
						$seqRef = $seqRefFlag;
					}
					$svID = $refChr . "_" . $svNo;
					$hash{$svID} = "$refChr\t$refStart\t$refEnd\t$svtype\t$insLen\t$mergeNum\t$mergeSam\t$seqAlt\t$seqRef";
					
				} else {  # if identity <= 0.5
					subPrintINS();
				}
				
			} else {  # if distance of INS site > 10
				subPrintINS();
			}
			
		} else {  # initialize the chr
			$svNo ++;
			$samID = "$samChr\_$samStart\_$samEnd\_$insLen";
			$mergeSam = $samID;
			$mergeNum = 1;
			$svID = $refChr . "_" . $svNo;
			$hash{$svID} = "$refChr\t$refStart\t$refEnd\t$svtype\t$insLen\t$mergeNum\t$mergeSam\t$seqAlt\t$seqRef";
			# refresh the flags
			$refChrFlag = $refChr;
			$refStartFlag = $refStart;
			$refEndFlag = $refEnd;
			$insLenFlag = $insLen;
			$seqAltFlag = $seqAlt;
			$seqRefFlag = $seqRef;
		}
	}
}

print "#svNo\tCHROM\tPOS\tEND\tSVTYPE\tSVLEN\tMERGE_SAMPLES_NUM\tMERGE_SAMPLES\tINS_SEQ\tREF_SEQ\n";

foreach my $id (sort {$a cmp $b} keys %hash) {  # print in order
	print "$id\t$hash{$id}\n";
}

close IN_SV;

### The followings are subroutines
#------------------
sub subPrintDEL {
	$svNo ++;
	$samID = "$samChr\_$samStart\_$samEnd\_$len";
	$mergeSam = $samID;
	$mergeNum = 1;
	$svID = $refChr . "_" . $svNo;
	$hash{$svID} = "$refChr\t$refStart\t$refEnd\t$svtype\t$len\t$mergeNum\t$mergeSam\t$seq";
	$refStartFlag = $refStart;
	$refEndFlag = $refEnd;
}

#------------------
sub subPrintINS {
	$svNo ++;
	$samID = "$samChr\_$samStart\_$samEnd\_$insLen";
	$mergeSam = $samID;
	$mergeNum = 1;
	$svID = $refChr . "_" . $svNo;
	$hash{$svID} = "$refChr\t$refStart\t$refEnd\t$svtype\t$insLen\t$mergeNum\t$mergeSam\t$seqAlt\t$seqRef";
	# refresh the flags
	$refStartFlag = $refStart;
	$refEndFlag = $refEnd;
	$insLenFlag = $insLen;
	$seqAltFlag = $seqAlt;
	$seqRefFlag = $seqRef;
}

#------------------
sub align {
	# Smith-Waterman Algorithm
	# origin: https://github.com/RodenLuo/Smith-Waterman-in-Perl
	
=example
	seq1 : AAAATTTTGGGGCCCC
	seq2 : ATAATTTGGCCCC	
	out  : 
		AATTTTGGGGCCCC
		AA-TTT--GGCCCC
=cut
	
	# get sequences
	my ($seq1, $seq2) = @_;

	# scoring scheme
	my $MATCH     =  1;  # +1 for letters that match
	my $MISMATCH  = -1;  # -1 for letters that mismatch
	my $GAP       = -1;  # -1 for any gap

	# initialization
	my @matrix;
	$matrix[0][0]{score}   = 0;
	$matrix[0][0]{pointer} = "none";
	for (my $j = 1; $j <= length($seq1); $j++) {
		$matrix[0][$j]{score}   = 0;
		$matrix[0][$j]{pointer} = "none";
	}
	for (my $i = 1; $i <= length($seq2); $i++) {
		$matrix[$i][0]{score}   = 0;
		$matrix[$i][0]{pointer} = "none";
	}

	# fill
	my $max_i     = 0;
	my $max_j     = 0;
	my $max_score = 0;


	for (my $i = 1; $i <= length($seq2); $i++) {
		for (my $j = 1; $j <= length($seq1); $j++) {
			my ($diagonal_score, $left_score, $up_score);
			
			# calculate match score
			my $letter1 = substr($seq1, $j-1, 1);
			my $letter2 = substr($seq2, $i-1, 1);      
			if ($letter1 eq $letter2) {
				$diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
			} else {
				$diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
			}
			
			# calculate gap scores
			$up_score   = $matrix[$i-1][$j]{score} + $GAP;
			$left_score = $matrix[$i][$j-1]{score} + $GAP;
			
			if ($diagonal_score <= 0 and $up_score <= 0 and $left_score <= 0) {
				$matrix[$i][$j]{score}   = 0;
				$matrix[$i][$j]{pointer} = "none";
				next;  # terminate this iteration of the loop
			}
			
			# choose best score
			if ($diagonal_score >= $up_score) {
				if ($diagonal_score >= $left_score) {
					$matrix[$i][$j]{score}   = $diagonal_score;
					$matrix[$i][$j]{pointer} = "diagonal";
				} else {
					$matrix[$i][$j]{score}   = $left_score;
					$matrix[$i][$j]{pointer} = "left";
				}
			} else {
				if ($up_score >= $left_score) {
					$matrix[$i][$j]{score}   = $up_score;
					$matrix[$i][$j]{pointer} = "up";
				} else {
					$matrix[$i][$j]{score}   = $left_score;
					$matrix[$i][$j]{pointer} = "left";
				}
			}
			
			# set maximum score
			if ($matrix[$i][$j]{score} > $max_score) {
				$max_i     = $i;
				$max_j     = $j;
				$max_score = $matrix[$i][$j]{score};
			}
		}
	}

	# trace-back
	my $align1 = "";
	my $align2 = "";

	my $j = $max_j;
	my $i = $max_i;

	while (1) {
		last if $matrix[$i][$j]{pointer} eq "none";
		
		if ($matrix[$i][$j]{pointer} eq "diagonal") {
			$align1 .= substr($seq1, $j-1, 1);
			$align2 .= substr($seq2, $i-1, 1);
			$i--; $j--;
		}
		elsif ($matrix[$i][$j]{pointer} eq "left") {
			$align1 .= substr($seq1, $j-1, 1);
			$align2 .= "-";
			$j--;
		}
		elsif ($matrix[$i][$j]{pointer} eq "up") {
			$align1 .= "-";
			$align2 .= substr($seq2, $i-1, 1);
			$i--;
		}
	}

	$align1 = reverse $align1;
	$align2 = reverse $align2;

	#print "$align1\n";
	#print "$align2\n";
	
	return "$align1\t$align2";
}

#------------------
sub calIdentity {
	my ($seq1, $seq2, $shortlen) = @_;
	
	# calulate the identity (the basic algorithm is from http://bbs.chinaunix.net/forum.php?mod=viewthread&tid=3621368)
	my $str = $seq1 ^ $seq2;
	my $num = $str =~ s/\0/*/g;
	my $ident = $num / $shortlen;  #
	
	return $ident;
}

#------------------
sub rc {
	my ($seq) = @_;
	my @seqs = split(//, $seq);
	@seqs = reverse @seqs;  # reverse seq
	my $seq_r = join "", @seqs;
	$seq_r =~ tr/agctAGCT/tcgaTCGA/;  # complement seq
	my $seq_rc .= $seq_r;
	
	return $seq_rc;
}
