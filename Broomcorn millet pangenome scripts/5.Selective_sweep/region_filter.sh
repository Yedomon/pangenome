#!/usr/bin/bash

top_size=(5)
pops=("C1_vs_C2" "C1_vs_C3" "C1_vs_W" "C2_vs_W" "C3_vs_W" "C_vs_W")

for i in ${top_size[*]};do
    for j in ${pops[*]};do
        # filter by xpclr, fst, pi
        bedtools intersect -a xpclr_top/xpclr_top${i}_merged/${j}.win_0.002.regions.txt \
        -b fst_top/fst_top${i}_merged/${j}.fst.top${i}.win_20000.regions.txt > \
        regions_tmp/${j}.xpclr.fst.top${i}.win_20000.regions.txt
        bedtools intersect -a regions_tmp/${j}.xpclr.fst.top${i}.win_20000.regions.txt \
        -b pi_top/pi_top50_merged/${j}.top50.pi.win_20000.regions.txt > \
        regions/${j}.top${i}.20kb.regions.txt

    done
done

