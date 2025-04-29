#! usr/bin/perl -w
open (FILE1, "$ARGV[0]") or die $!;
     
$file3=join('',<FILE1>);
@file4=split('\n',$file3);
foreach $line(@file4){

@line1=split('\t',$line);

if($line1[0]=~/(S0\_[0-9]*[a-z]*)\_[0-9]*.fq.gz/){

	$id=$1;
}
print "mkdir $id\n";

print "fastqc \-o $id $line1[0] $line1[1]\n";



}

               


