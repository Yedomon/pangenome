#use warnings;
#use strict;

open(FILE1,"$ARGV[0]");


$tg=join('',<FILE1>);
@yh=split('\n',$tg);

foreach $lk(@yh){
	#	@lk1=split('\t',$lk);
system("awk \'\$2~/^$lk$\/\' uniq_list_of_genes_in_pair_Os_Oc.txt > temp10.txt");




#system("awk \'$lk1[1]==\$1\{print \$2\}\' rcb_os_oc3.txt > temp2.txt");

#system("awk \'\!x\[\$0\]\+\+\' temp1.txt > temp3.txt");
#system("awk \'\!x\[\$0\]\+\+\' temp2.txt > temp4.txt");

#system("awk \'NR\=\=FNR\{a\[\$0\]\+\+\;next\} a\[\$0\]\' temp3.txt temp4.txt > temp5.txt");
$gene=`paste -s -d '\t' temp10.txt`;

chomp($gene);

print "$gene\n";

	     
}




  


