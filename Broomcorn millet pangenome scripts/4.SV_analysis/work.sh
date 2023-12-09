### filter gap
perl 01.filter_gap_sv.pl -bed ALL.gap.bed -sv All.del.sort.bed -svtype "DEL" > All.del.sort.51bp.nogap.bed
perl 01.filter_gap_sv.pl -bed ALL.gap.bed -sv All.ins.sort.bed -svtype "INS" > All.ins.sort.51bp.nogap.bed
perl 01.filter_gap_sv.pl -bed ALL.gap.bed -sv All.inv.reference.sort.bed -svtype "INV" > All.inv.reference.sort.51bp.nogap.bed
perl 01.filter_gap_sv.pl -bed ALL.gap.bed -sv All.trl.reference.sort.bed -svtype "TRL" > All.trl.reference.sort.51bp.nogap.bed

### filter cen
perl 01.filter_cen_sv.pl -gff ALL.cent.gff -sv All.del.sort.51bp.nogap.bed -svtype "DEL" > All.del.sort.51bp.nogap.nocen.bed 
perl 01.filter_cen_sv.pl -gff ALL.cent.gff -sv All.ins.sort.51bp.nogap.bed -svtype "INS" > All.ins.sort.51bp.nogap.nocen.bed

### merge (sort)
# DEL
cat All.del.sort.51bp.nogap.nocen.bed |sort -k1,1 -k2,2n -k3,3n -k8,8n -k16,16 > All.del.sort.51bp.nogap.nocen.bed.sort
perl 02.merge_sv_smith-identity_seq.pl -sv All.del.sort.51bp.nogap.nocen.bed.sort -svtype "DEL" > All.del.sort.51bp.nogap.nocen.sort.merge.txt
# INS
mkdir All.ins.sort.51bp.nogap.nocen.bed.sort.chr
perl 02.merge_sv_smith-identity_seq.pl -sv All.ins.sort.51bp.nogap.nocen.bed.sort -svtype "INS" > All.ins.sort.51bp.nogap.nocen.sort.merge.txt
# INV
perl 02.merge_sv.pl -sv All.inv.reference.sort.51bp.nogap.bed -svtype "INV" > All.inv.reference.sort.51bp.nogap.merge.txt
# TRL
perl 02.merge_sv_smith-identity_seq.pl -sv All.trl.reference.sort.51bp.nogap.bed -svtype "TRL" > All.trl.reference.sort.51bp.nogap.bed.merge.txt

### bed2vcf
perl 03.svbed2vcf.pl -sv All.del.sort.51bp.nogap.nocen.merge.txt -svtype "DEL" > All.del.sort.51bp.nogap.merge.vcf
perl 03.svbed2vcf.pl -sv All.ins.sort.51bp.nogap.nocen.merge.txt -svtype "INS" > All.ins.sort.51bp.nogap.merge.vcf

### get anno
perl 04.sv_feature_anno_LM_new.pl -sv All.del.sort.51bp.nogap.nocen.merge.txt -gff Pmlongmi4_features.gff > All.del.anno.txt
perl 04.sv_feature_anno_LM_new.pl -sv All.ins.sort.51bp.nogap.nocen.merge.txt -gff Pmlongmi4_features.gff > All.ins.anno.txt

### count fold change between sv-gene and nonsv-gene
perl 05.sv_gene_anno_FCnoLog_p_4newid.pl -func Longmi_gene2rice-.txt -exp ALL.FPKM.NEW.txt -anno ALL.PAV.anno.txt.edit.newid > ALL.PAV.anno.txt.edit.newid.FC
