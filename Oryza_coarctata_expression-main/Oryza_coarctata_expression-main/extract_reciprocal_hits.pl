#use warnings;
#use strict;

open(FILE1,"$ARGV[0]");
open(FILE2,"$ARGV[1]");

   while ($line1 = <FILE1>)
   {
       chomp($line1);
       @rf=split('\t',$line1);
     
            $hash{"$id"} = $fasta;
            $id = $rf[0];
            $fasta = $line1;
             }
$hash{"$id"} = $fasta;  


$tg=join('',<FILE2>);
@yh=split('\n',$tg);

foreach $lk(@yh){
	@lk1=split('\t',$lk);

if(exists($hash{$lk1[0]})){
       $pheno=$hash{$lk1[0]};
          
   print "$lk\t$pheno\n";

           }

	     
}




  


