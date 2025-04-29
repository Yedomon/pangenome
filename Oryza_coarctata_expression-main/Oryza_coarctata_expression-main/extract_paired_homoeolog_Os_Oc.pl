#use warnings;
#use strict;

open(FILE1,"$ARGV[0]");


$tg=join('',<FILE1>);
@yh=split('\n',$tg);

foreach $lk(@yh){
        @lk1=split('\t',$lk);
system("awk \'$lk1[0]==\$1\{print \$2\}\' uniq_list_of_genes_in_pair_Os_Oc.txt > temp1.txt");
system("awk \'$lk1[1]==\$1\{print \$2\}\' uniq_list_of_genes_in_pair_Os_Oc.txt > temp2.txt");

system("awk \'\!x\[\$0\]\+\+\' temp1.txt > temp3.txt");
system("awk \'\!x\[\$0\]\+\+\' temp2.txt > temp4.txt");

system("awk \'NR\=\=FNR\{a\[\$0\]\+\+\;next\} a\[\$0\]\' temp3.txt temp4.txt > temp5.txt");
$gene=`paste -s -d ',' temp5.txt`;

#print "$gene\n";

chomp($gene);
print "g$lk1[0]\tg$lk1[1]\t$gene\n";

             
}

