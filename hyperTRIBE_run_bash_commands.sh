#ECT2 project run
for fn in `cat /binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/pipeline/files/prefix_names_missing_shoots.txt`
do
nice bash /binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/pipeline/CODE_HT/trim_and_align_edited.sh $fn
done

samtools index E2T_St1.sort.bam
samtools index E2T_St2.sort.bam

#FROM SORTED BAM FILES START HERE: start with samtools mpileup
#mpileup + compile counts for roots and shoots separately:
samtools mpileup --max-depth 50000 -Q 30 --skip-indels -f /binf-isilon/alab/projects/ECT2_TC/annotations/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta/genome.fa E2T_Rc11.sort.bam E2T_Rc12.sort.bam E2T_Rc13.sort.bam E2T_Rc14.sort.bam E2T_Rc15.sort.bam E2T_Re6.sort.bam  E2T_Re7.sort.bam  E2T_Re8.sort.bam  E2T_Re9.sort.bam E2T_Re10.sort.bam E2T_Rt1.sort.bam E2T_Rt2.sort.bam E2T_Rt3.sort.bam E2T_Rt4.sort.bam E2T_Rt5.sort.bam | perl /home/sarah/R_code/ECT2/hyperTRIBE_mpileup2bases.pl> baseCounts_roots_hyperTRIBE.txt &
samtools mpileup --max-depth 50000 -Q 30 --skip-indels -f /binf-isilon/alab/projects/ECT2_TC/annotations/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta/genome.fa E2T_Sc11.sort.bam E2T_Sc12.sort.bam  E2T_Sc13.sort.bam  E2T_Sc14.sort.bam  E2T_Sc15.sort.bam E2T_Se6.sort.bam  E2T_Se7.sort.bam  E2T_Se8.sort.bam  E2T_Se9.sort.bam  E2T_Se10.sort.bam E2T_St1.sort.bam  E2T_St2.sort.bam  E2T_St3.sort.bam  E2T_St4.sort.bam  E2T_Se5.sort.bam | perl /home/sarah/R_code/ECT2/hyperTRIBE_mpileup2bases.pl> baseCounts_shoots_hyperTRIBE.txt &
