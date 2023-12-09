#!/bin/bash
#PBS -l nodes=1:ppn=48
#PBS -l walltime=50:00:00
#PBS -d ./
#PBS -j oe

CPU=$PBS_NP
if [ ! $CPU ]; then CPU=2; fi

N=$PBS_ARRAYID
if [ ! $N ]; then N=1; fi

#----------------
hifiasm=/public1/home/liuyang/program/hifiasm-0.14.2/hifiasm
reads=`ls *fastq.gz | head -n $N | tail -n 1`
prefix=${reads%%.*}
#prefix=${reads%.fastq.gz}

start=`date +%s`
start_date=`date -d @"$start" "+%Y-%m-%d %H:%M:%S"`
echo "<<<<< Start= $start_date"
echo "<<<<< CPU= $CPU"
echo "<<<<< FILE= $reads"
echo "<<<<< FILE_PREFIX= $prefix"

#----------------
if [ ! -e ${prefix}.asm ]; then
	echo "<<<<< Running hifiasm ..."
	echo "Command: \
	$hifiasm -o ${prefix}.asm -t $CPU -l 2 -u ${prefix}.fastq.gz 2> ${prefix}.asm.log "
	#
	$hifiasm -o ${prefix}.asm -t $CPU -l 2 -u ${prefix}.fastq.gz 2> ${prefix}.asm.log
	#
	awk '/^S/{print ">"$2;print $3}' ${prefix}.asm.p_ctg.gfa > ${prefix}.asm.p_ctg.fa
	#
	perl -e 'print "-- Nxx Summary --\n";$min=1e10;$max=0;my ($len,$total)=(0,0);my @x;while(<>){if(/>/){$contig++;}if(/^[\>\@]/){if($len>0){$total+=$len;push@x,$len;};$len=0;}else{s/\s//g;$len+=length($_);}}if ($len>0){$total+=$len;$mean=$total/$contig;push @x,$len;}@x=sort{$b<=>$a}@x;my ($count,$flag)=(0,0);for (my $j=0;$j<@x;$j++){$count+=$x[$j];if (($count>=$total*0.5)&&($flag==0)){print "N50\t$x[$j] bp\n";$flag=1;}elsif (($count>=$total*0.9)&&($flag==1)){print "N90\t$x[$j] bp\n";$flag=2;}elsif (($count>=$total*0.95)&&($flag==2)){print "N95\t$x[$j] bp\n";$flag=3;}elsif ($count==$total){print "Contig number\t$contig\nTotal length\t$total bp\nMean length\t$mean bp\nMax length\t$x[0] bp\nMin length\t$x[$j] bp\n";exit;}}' ${prefix}.asm.p_ctg.fa > ${prefix}.asm.p_ctg.fa.N50
fi

#----------------
end=`date +%s`
end_date=`date -d @"$end" "+%Y-%m-%d %H:%M:%S"`
runtime=$((end-start))
h=$(($runtime/3600))
hh=$(($runtime%3600))
m=$(($hh/60))
s=$(($hh%60))
echo "<<<<< CPU= $CPU"
echo "<<<<< Start= $start_date"
echo "<<<<< End= $end_date"
echo "<<<<< Run time= $h:$m:$s"
echo "<<<<< Done!!!"

