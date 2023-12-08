curr_dir=$(pwd)
bamtofastq=/opt/lib/biobambam2-2.0.180/bin/bamtofastq
star_index_dir=/opt/ref-genomes/human_genome_star
#assembly=$star_index_dir/Homo_sapiens.GRCh38.dna.primary_assembly.fa
#annotation=$star_index_dir/Homo_sapiens.GRCh38.103.gtf
arriba=/opt/lib/arriba_v2.1.0/arriba 
arriba_db=/opt/lib/arriba_v2.1.0/database
featureCounts=/opt/lib/subread-2.0.2-source/bin/featureCounts
threads=12
out_dir=$curr_dir


for x in $(awk 'BEGIN {FS="\t"}; NR>0 && NR<2 {print $1}' list.txt);
do
        echo $x;
#       $bamtofastq collate=1 exclude=QCFAIL,SECONDARY,SUPPLEMENTARY filename=$x gz=1 level=1  outputperreadgroup=1 outputperreadgroupprefix=${x%.*} inputformat=bam outputdir=$curr_dir

echo "## Running STAR for PE reads ##"

STAR \
    --runThreadN $threads \
    --genomeDir $star_index_dir \
    --genomeLoad NoSharedMemory \
    --readFilesIn ${x%.*}'_default_1.fq' ${x%.*}'_default_2.fq' \
    --readFilesCommand zcat \
    --outStd BAM_Unsorted \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped Within \
    --outBAMcompression 7 \
    --outFilterMultimapNmax 50 \
    --peOverlapNbasesMin 10 \
    --alignSplicedMateMapLminOverLmate 0.5 \
    --alignSJstitchMismatchNmax 5 -1 5 5 \
    --chimSegmentMin 10 \
    --chimOutType WithinBAM HardClip \
    --chimJunctionOverhangMin 10 \
    --chimScoreDropMax 30 \
    --chimScoreJunctionNonGTAG 0 \
    --chimScoreSeparation 1 \
    --chimSegmentReadGapMax 3 \
    --chimMultimapNmax 50 \
    > ${x%.*}'_aligned.out.PE.bam' 

echo "## Running STAR for o1 reads ##"

STAR \
    --runThreadN $threads \
    --genomeDir $star_index_dir \
    --genomeLoad NoSharedMemory \
    --readFilesIn ${x%.*}'_default_o1.fq' \
    --readFilesCommand zcat \
    --outStd BAM_Unsorted \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped Within \
    --outBAMcompression 7 \
    --outFilterMultimapNmax 50 \
    --peOverlapNbasesMin 10 \
    --alignSplicedMateMapLminOverLmate 0.5 \
    --alignSJstitchMismatchNmax 5 -1 5 5 \
    --chimSegmentMin 10 \
    --chimOutType WithinBAM HardClip \
    --chimJunctionOverhangMin 10 \
    --chimScoreDropMax 30 \
    --chimScoreJunctionNonGTAG 0 \
    --chimScoreSeparation 1 \
    --chimSegmentReadGapMax 3 \
    --chimMultimapNmax 50 \
    > ${x%.*}'_aligned.out.SE1.bam'

echo "## Running STAR for o2 reads ##"

STAR \
    --runThreadN $threads \
    --genomeDir $star_index_dir \
    --genomeLoad NoSharedMemory \
    --readFilesIn ${x%.*}'_default_o2.fq' \
    --readFilesCommand zcat \
    --outStd BAM_Unsorted \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped Within \
    --outBAMcompression 7 \
    --outFilterMultimapNmax 50 \
    --peOverlapNbasesMin 10 \
    --alignSplicedMateMapLminOverLmate 0.5 \
    --alignSJstitchMismatchNmax 5 -1 5 5 \
    --chimSegmentMin 10 \
    --chimOutType WithinBAM HardClip \
    --chimJunctionOverhangMin 10 \
    --chimScoreDropMax 30 \
    --chimScoreJunctionNonGTAG 0 \
    --chimScoreSeparation 1 \
    --chimSegmentReadGapMax 3 \
    --chimMultimapNmax 50 \
    > ${x%.*}'_aligned.out.SE2.bam'

echo "## merging BAM files ##"

# merge bam files
samtools merge -1 --threads $threads ${x%.*}'_aligned.out.bam' ${x%.*}'_aligned.out.PE.bam' ${x%.*}'_aligned.out.SE1.bam' ${x%.*}'_aligned.out.SE2.bam'


# tee ${x%.*}'_aligned.out.bam' |
echo "## Running arriba for fusion detection ##"

$arriba \
    -x ${x%.*}'_aligned.out.bam' \
    -o ${x%.*}'_fusions.tsv' \
    -O ${x%.*}'_fusions.discarded.tsv' \
    -a $star_index_dir/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -g $star_index_dir/Homo_sapiens.GRCh38.103.gtf \
    -b $arriba_db/blacklist_hg38_GRCh38_v2.1.0.tsv.gz \
    -k $arriba_db/known_fusions_hg38_GRCh38_v2.1.0.tsv.gz \
    -t $arriba_db/known_fusions_hg38_GRCh38_v2.1.0.tsv.gz \
    -p $arriba_db/protein_domains_hg38_GRCh38_v2.1.0.gff3 \
    -X 


# moving files to the proper location
echo "## moving files to the Final location"

mkdir $out_dir/STAR_out
mv ${x%.*}'_aligned.out.bam' $out_dir/STAR_out

mkdir $out_dir/arriba_out
mv ${x%.*}'_fusions.tsv' ${x%.*}'_fusions.discarded.tsv' $out_dir/arriba_out

# cleaning up
# rm ${x%.*}'_default_1.fq' ${x%.*}'_default_2.fq' ${x%.*}'_default_o1.fq' ${x%.*}'_default_o2.fq'

# running featureCounts 
# $featureCounts -p -t exon -g gene_id -a $star_index_dir/Homo_sapiens.GRCh38.103.gtf -o counts.txt mapping_results_PE.bam

done



