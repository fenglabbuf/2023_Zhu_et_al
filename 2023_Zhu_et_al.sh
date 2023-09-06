### align fastq by STAR ###
# create refGene index on hg38
STAR --runThreadN 5 --runMode genomeGenerate \
     --genomeDir STAR/index \
     --genomeFastaFiles STAR/hg38.fasta \
     --sjdbGTFfile STAR/hg38.refGene.gtf \
     --sjdbOverhang 49 \
     --genomeSAsparseD 2

# align fastq and generate splice junction count files
STAR --runThreadN 5 \
     --genomeDir STAR/index \
     --readFilesCommand zcat \
     --readFilesIn Anp_D3_1F.fastq.gz Anp_D3_1R.fastq.gz \
     --outFileNamePrefix "STAR/BAM" \
     --genomeSAsparseD 2 \
     --outFilterType BySJout \
     --outFilterMultimapNmax 1 \
     --alignSJoverhangMin 8 \
     --alignSJDBoverhangMin 1 \
     --outFilterMismatchNmax 999 \
     --alignIntronMin 20 \
     --alignIntronMax 1000000 \
     --alignMatesGapMax 1000000 \
     --outSAMstrandField intronMotif \
     --outSAMtype BAM SortedByCoordinate \
     --limitBAMsortRAM 31000000000 \
     --quantMode TranscriptomeSAM


### use leafcutter to calculate statistics on splice junctions
# reformat SJ.out.tab to .junc
awk 'BEGIN{FS="\t";OFS="\t"}{print $1, $2-1, $3, "JUNCBJ"NR, $7, ($4 == 1)? "+":"-", $2-1, $3, "255,0,0", "2", "0,0", "0,0"}' $i > "${i}"".junc"

# generate junction count matrix for all samples
regtools junctions extract -a 8 -m 50 -M 500000 -s 1 $i -o "${i}"".junc"

python ~/Desktop/leafcutter/scripts/leafcutter_cluster_regtools_py3.py -j juncfiles.txt \
	   -o leafcutterjcount -l 500000 -m 20

# compare between samples and calculate statistics
Rscript ~/Desktop/leafcutter/scripts/leafcutter_ds.R \
        leafcutterjcount_perind_numers.counts.gz \
        "group""$i"".txt" -o $i \
        -p 4 -i 3

Rscript ~/Desktop/leafcutter/leafviz/prepare_results.R \
        -o "$i"".RData" \
        -m "$i"".txt" \
        -f 0.1 \
        leafcutterjcount_perind_numers.counts.gz \
        "$i""_cluster_significance.txt" \
        "$i""_effect_sizes.txt" \
        ~/Desktop/leafcutter/leafviz/leafanno