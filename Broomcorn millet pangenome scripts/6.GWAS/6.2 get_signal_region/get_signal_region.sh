if [ -d $1 ];then            #$1 is a folder containing significant gwas signal results for the trait across all environments
python merge_sig.py $1 
fi
if [ -s "$1".txt ];then 
   python get_signal_block.py $1.txt >./$1.sig.block
else
   echo "No significant signal was found for this trait....."
fi

