#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my %opt;
GetOptions (\%opt, "fa:s", "group:s", "help");
my $help=<<USAGE;
The script is edited to 
By liuyang.
Run: perl $0 -fa ALL.BC.cds.fa -group BCnr100_Orthogroups.num.txt > out

--fa : 
>BC027_01G000100.1
ATGGCTCTCAACTACCTGGCCGTCGCCGCCATCAACTTCGTCGCCGCGCTGCTCTCCATCCCCATCATCG
CCGCCGGCATCTGGCTGTCTACGCAGCCCGACAACGCCTGCGTCCAGATCCTCCAGTGGCCTGTCATTGC

--group : 
32	3	5	16	8	OG0000000	BC027_01G032000.1 BC027_01G032400.1 BC027_01G303600.1
32	3	5	16	8	OG0000001	BC027_14G267300.1 BC027_UnG112200.1 BC027_UnG122300.1

USAGE

if ($opt{help} or keys %opt < 1) {
	print "$help";
	exit();
}

######
my %hash_seq;
getfastaseq($opt{fa});

######
my $group = $opt{group};
open GROUP, $group or die $!;
while (<GROUP>) {
	chomp;
	my ($num,$familyID,$genes) = (split/\t/)[0,5,6];
	
	# 
	my (@geneList, @accList);
	push (@geneList, $1) while ($genes =~ /(BC\S+)/mg);  # BC027_01G032000.1
	push (@accList, $1) while ($genes =~ /(BC\d+)/mg);  # BC027

	# delete the repeat elements
	my %count;
	my @nrAcc = grep { ++$count{$_} == 1 } @accList;
	my $nrNum = $#nrAcc + 1;  #
	
	# 
	my %hash_gene;
	foreach (@geneList) {
		my $accID = $1 if ($_ =~ /(BC\d+)/);
		$hash_gene{$accID} = $_ unless exists $hash_gene{$accID};
	}
	
	my $class;
	if ($num == 32) {
		$class = "pan_core.";
	} elsif ($num == 30 or $num == 31) {
		$class = "pan_softcore.";
	} elsif ($num >= 2 && $num <= 29) {
		$class = "pan_dispensable.";
	} elsif ($num == 1) {
		$class = "pan_private.";
	}
	
	# print seq file
	open OUT, ">", $class . $familyID . ".fa" or die $!;
	foreach my $k (sort {$a cmp $b} keys %hash_gene) {  # print in order
		my $id = $hash_gene{$k};
		print OUT ">$id\n$hash_seq{$id}\n";
	}
	close OUT;
}

close GROUP;


####################
sub getfastaseq
{
$/=">";
#my %hash_seq;
my ($file) = @_;
open IN, "$file" or die "$!";
while (<IN>) {
	#next if (length $_ < 2);
	my @unit = split("\n", $_);
	my $temp = shift @unit;
	my @temp1 = split(" ", $temp);
	my $head = $temp1[0];
	my $seq = join("\n", @unit);
	$seq =~ s/\>//g;
	$seq =~ s/\n//g;
	$seq =~ s/\s//g;
	$hash_seq{$head} = $seq;
}
close IN;
$/="\n";
return \%hash_seq;
}
