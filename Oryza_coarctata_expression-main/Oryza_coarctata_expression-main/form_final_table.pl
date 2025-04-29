#use warnings;
#use strict;

open(FILE1,"$ARGV[0]");
open(FILE2,"$ARGV[1]");
$subgenome=$ARGV[2];
   while ($line1 = <FILE1>)
   {
       chomp($line1);
       @rf=split('\t',$line1);
     
            $hash{"$id"} = $fasta;
            $id = $rf[0];
            $fasta = $rf[1];
             }
$hash{"$id"} = $fasta;  

$tg=join('',<FILE2>);
@yh=split('\n',$tg);


print "Gene_Id\tClass\tlog2(TPM+1)\n";

foreach $lk(@yh){
	@lk1=split('\t',$lk);

if(exists($hash{$lk1[0]})){
       $pheno1=$hash{$lk1[0]};
          
                
   $pheno2 = log($pheno1)/log(2)
                
                 }


print "$lk1[0]\t$lk1[1]\t$pheno2\n";
	     
}

