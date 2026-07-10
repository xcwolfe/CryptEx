#!/bin/bash
# CRYPTEX - Optimized for Yale Bouchet HPC Cluster (Slurm)
# Jack Humphrey / UCL / Adapted for Dong Lab Pipeline Standard
set -euo pipefail

## Default base variables
oFolder="/home/zw529/donglab/pipelines/modules/rnaseq/bin/CryptEx"
Step2_master="${oFolder}/Step2_master.sh"

# R scripts for CryptEx
Rbin=$(which R || echo "R")
dexseqFinalProcessR="${oFolder}/dexseq/forked_dexseq_pipeline_v2.R"
countPrepareR="${oFolder}/dexseq/forked_counts_prepare_pipeline.R"
deseqFinalProcessR="${oFolder}/dexseq/forked_deseq2_pipeline.R"
R_support_chopper="${oFolder}/support_frame_chopper.R"

# R scripts for downstream analyses
R_splice_junction_analyzer="${oFolder}/downstream_analyses/splice_junction_analyzer.R"
R_functional_enrichment="${oFolder}/downstream_analyses/Functional_Enrichment.R"

# Explicit path to verified DEXSeq cluster script location
pycount="/home/zw529/R/x86_64-pc-linux-gnu-library/4.4/DEXSeq/python_scripts/dexseq_count.py"

# ==============================================================================
# DEFAULT BASE COHORT ASSIGNMENTS (Optimized for Zach's TargetALS Human Pipelines)
# ==============================================================================
species="human"
protein="TDP43"
support="/home/zw529/donglab/data/target_ALS/CryptEx/Motor_Cortex/Motor_Cortex_support.tab"
submit="no"
splice_extractor="no"
gff_creator="yes"
read_counter="yes"
DEXSeq="no"
DESeq="no"
cohort_merger="no"
strict="no"
strict_num=500
paired="no"
stranded="no"
splice_junction_analyzer="no"
functional_enrichment="no"
intron_retainer="no"
intron_DEXSeq="no"
hold_Step1="no"
IGV="no"
intron_BED_passed="" # Tracking holder for dynamic override

# Argument parsing matrix
until [ -z "${1:-}" ]; do
    case "$1" in
        --species) shift; species=$1 ;;
        --annotation_file) shift; annotation_file=$1 ;;
        --protein) shift; protein=$1 ;;
        --submit) shift; submit=$1 ;;
        --support) shift; support=$1 ;;
        --gff) shift; gff=$1 ;;
        --intron_BED) shift; intron_BED_passed=$1 ;; # Added parser check
        --splice_extractor) shift; splice_extractor=$1 ;;
        --gff_creator) shift; gff_creator=$1 ;;
        --read_counter) shift; read_counter=$1 ;;
        --DEXSeq) shift; DEXSeq=$1 ;;
        --DESeq) shift; DESeq=$1 ;;
        --intron_retainer) shift; intron_retainer=$1 ;;
        --intron_DEXSeq) shift; intron_DEXSeq=$1 ;;
        --cohort_merger) shift; cohort_merger=$1 ;;
        --strict) shift; strict=$1 ;;
        --strict_num) shift; strict_num=$1 ;; # Kept synchronized
        --IGV) shift; IGV=$1 ;;
        --splice_junction_analyzer) shift; splice_junction_analyzer=$1 ;;
        --functional_enrichment) shift; functional_enrichment=$1 ;;
        --hold_Step1) shift; hold_Step1=$1 ;;
        --paired) shift; paired=$1 ;;
        --stranded) shift; stranded=$1 ;;
        --outdir) shift; oFolder=$1 ;;
        *) echo "Unrecognised argument: $1"; exit 1 ;;
    esac
    shift
done

# Resolve execution pathing layout
if [ -z "${oFolder:-}" ]; then
    oFolder="/home/zw529/donglab/pipelines/modules/rnaseq/bin/CryptEx"
fi

results="${oFolder}/${protein}_${species}"
reference="/home/zw529/donglab/pipelines/modules/rnaseq/bin/CryptEx/reference"
clusterFolder="${results}/cluster"

for folder in "$oFolder" "$results" "$reference" "$clusterFolder" "${clusterFolder}/out" "${clusterFolder}/error" "${clusterFolder}/R" "${clusterFolder}/submission"; do
    mkdir -p "$folder"
done

# Absolute Gencode v49 Annotation Anchors
ref_dir="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode"
gff_base="${ref_dir}/gencode.v49.annotation"
annotation_file="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/STAR/geneInfo.tab"

exon_GFF="${gff_base}_exons_only.gff"

# Use passed --intron_BED parameter if present, otherwise fall back to default
if [ -n "$intron_BED_passed" ]; then
    intron_BED="$intron_BED_passed"
else
    intron_BED="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/STAR/introns.sorted.bed"
fi

intron_tweaked_GFF="${oFolder}/reference/human_introns_for_HTseq.gff"

if [ ! -e "$intron_tweaked_GFF" ]; then
    mkdir -p "${oFolder}/reference"
    if [ -f "${gff_base}.gtf" ]; then
        sed 's/intron/exonic/g' "${gff_base}.gtf" > "$intron_tweaked_GFF"
    else
        sed 's/intron/exonic/g' "${ref_dir}/gencode.v49.annotation.gtf" > "$intron_tweaked_GFF"
    fi
fi
intron_GFF="$intron_tweaked_GFF"

## Central Log File Allocation
report_file="${results}/report.txt"
if [[ "$strict" == "yes" ]]; then
    report_file="${results}/report_strict${strict_num}.txt"
fi

echo "CryptEx
Started at: $(date)
--Protein:  $protein
--Species:  $species
--Step1:    $splice_extractor
--hold_Step1:   $hold_Step1
--Step2:    $gff_creator
--Step3:    $read_counter
--Step4:    $DEXSeq
--Strict:   $strict $strict_num
Support file:
" >> "$report_file"
cat "$support" >> "$report_file"

#################################
### STEP 1: SPLICE EXTRACTION ###
#################################

Step1_jobID="Step1_${protein}_${species}"
if [[ "$strict" == "yes" ]]; then 
    Step1_jobID="Step1_${protein}_${species}_strict_${strict_num}"
fi

Step1_jobscript="${clusterFolder}/submission/Step1_${protein}_${species}.sh"
Step1_taskfile="${clusterFolder}/submission/Step1_${protein}_${species}_tasks.txt"
rm -f "$Step1_taskfile"

sample_num_1=$(awk 'NR > 1' "$support" | wc -l)

if [[ "$splice_extractor" == "yes" ]]; then
    echo "creating job scripts and task file for spliced read extraction" 

    echo "#!/bin/bash
#SBATCH --mem=6G
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=1
#SBATCH --job-name=${Step1_jobID}
#SBATCH --array=1-${sample_num_1}%30

if [[ \"\$SLURM_ARRAY_TASK_ID\" == \"1\" ]]; then
    echo \"Step1 started at \$(date +%H:%M:%S)\" >> $report_file
fi    

TARGET_SCRIPT=\$(sed -n \"\${SLURM_ARRAY_TASK_ID}p\" \"$Step1_taskfile\")
bash \"\$TARGET_SCRIPT\"
" > "$Step1_jobscript"

    awk 'NR > 1 {print $1,$2,$3}' "$support" | while read -r sample bam dataset; do
        splicefolder="${results}/${dataset}/splice_extraction"
        mkdir -p "$splicefolder"
        output="${splicefolder}/${dataset}_${sample}"
        sample_jobscript="${clusterFolder}/submission/splice_extract_${protein}_${species}_${dataset}_${sample}.sh"
        
        echo "#!/bin/bash
samtools view -h -F 256 $bam | awk '\$1 ~ /@/ || \$6 ~ /N/' | samtools view -bh - > ${output}.spliced.bam 
bedtools intersect -a ${output}.spliced.bam -b ${exon_GFF} > ${output}.spliced.exons.bam
bedtools bamtobed -i ${output}.spliced.exons.bam -split | sort -k1,1 -k2,2n > ${output}.spliced.bed
bedtools intersect -a ${output}.spliced.bed -b ${exon_GFF} -v > ${output}.spliced.introns.bed
echo \"Step 1 finished for $sample at \$(date +%H:%M:%S)\" >> $report_file 
" > "$sample_jobscript"

        echo "$sample_jobscript" >> "$Step1_taskfile"
    done
fi

############################################
### STEP 2: BED MERGING AND GFF CREATION ###
############################################

Step2_jobscript="${clusterFolder}/submission/Step2_GFF_creator_${protein}_${species}.sh"
if [[ "$strict" == "yes" ]]; then 
    Step2_jobscript="${clusterFolder}/submission/Step2_GFF_creator_${protein}_${species}_strict_${strict_num}.sh"
fi
rm -f "$Step2_jobscript"

Step2_jobID="Step2_${protein}_${species}"
if [[ "$strict" == "yes" ]]; then 
    Step2_jobID="Step2_${protein}_${species}_strict_${strict_num}"
fi

Step2_taskfile="${clusterFolder}/submission/Step2_${protein}_${species}_tasks.txt"
rm -f "$Step2_taskfile"

sample_num_2=$(awk 'NR > 1 {print $3}' "$support" | uniq | wc -l)

if [ "$gff_creator" = "yes" ]; then
    echo "creating job script and task file for GFF creation"

    echo "#!/bin/bash
#SBATCH --mem=24G
#SBATCH --time=9:00:00
#SBATCH --cpus-per-task=4
#SBATCH --job-name=${Step2_jobID}

if [[ \"\$SLURM_ARRAY_TASK_ID\" == \"1\" ]]; then
    echo \"Step2 started at \$(date +%H:%M:%S)\" >> $report_file
fi      

TARGET_SCRIPT=\$(sed -n \"\${SLURM_ARRAY_TASK_ID}p\" \"$Step2_taskfile\")
bash \"\$TARGET_SCRIPT\"
" > "$Step2_jobscript"

    for dataset in $(awk 'NR > 1 {print $3}' "$support" | uniq); do
        mkdir -p "${results}/${dataset}/GFF"
        spliced_beds="${results}/${dataset}/splice_extraction/"
        output="${results}/${dataset}/GFF/${protein}_${species}_${dataset}"
        step2_dataset_script="${clusterFolder}/submission/GFF_creator_${protein}_${species}_${dataset}.sh"

        echo "#!/bin/bash
bash $Step2_master --dataset ${dataset} --output ${output} --intron_BED ${intron_BED} --exon_GFF ${exon_GFF} --strict ${strict} --strict_num ${strict_num} --spliced_beds ${spliced_beds}
" > "$step2_dataset_script"

        echo "$step2_dataset_script" >> "$Step2_taskfile"
    done
fi

##############################################################
## STEP 3: READ COUNTING FOR EACH BAM WITH THE NEW GFF FILE ###
##############################################################

if [[ "$read_counter" = "yes" ]]; then 

    Step3_master_jobscript="${clusterFolder}/submission/Step3_count_${protein}_${species}.sh"
    Step3_taskfile="${clusterFolder}/submission/Step3_${protein}_${species}_tasks.txt"
    rm -f "$Step3_taskfile"

    step3_num=$(awk 'NR > 1' "$support" | wc -l)
    Step3_jobID="Step3_${protein}_${species}"
    if [[ "$strict" == "yes" ]]; then 
        Step3_jobID="Step3_${protein}_${species}_strict_${strict_num}"
    fi

    echo "#!/bin/bash
#SBATCH --mem=24G
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=2
#SBATCH --job-name=$Step3_jobID
#SBATCH --array=1-${step3_num}%30

if [[ \"\$SLURM_ARRAY_TASK_ID\" == \"1\" ]]; then
    echo \"Step3 started at \$(date +%H:%M:%S)\" >> $report_file
fi

TARGET_SCRIPT=\$(sed -n \"\${SLURM_ARRAY_TASK_ID}p\" \"$Step3_taskfile\")
bash \"\$TARGET_SCRIPT\"
" > "$Step3_master_jobscript"

    awk 'NR > 1 {print $1,$2,$3,$4}' "$support" | while read -r sample bam dataset condition; do
        if [[ "$strict" == "no" ]]; then
            GFF="${results}/${dataset}/GFF/${protein}_${species}_${dataset}.total.cryptics.gff"
        else
            GFF="${results}/${dataset}/GFF/${protein}_${species}_${dataset}.strict.${strict_num}.total.cryptics.gff"
        fi

        info_table="${oFolder}/support/${protein}_${species}_info.tab"
        if [ -e "$info_table" ]; then
            paired_val=$(awk -v d="$dataset" '\$1==d {print \$2}' "$info_table")
            stranded_val=$(awk -v d="$dataset" '\$1==d {print \$3}' "$info_table")
        else
            paired_val=$paired
            stranded_val=$stranded
        fi

        countFolder="${results}/${dataset}/counts"
        Step3_sample_jobscript="${clusterFolder}/submission/count_${dataset}_${sample}.sh"

        if [[ "$strict" = "yes" ]]; then
            countFolder="${results}/${dataset}/strict_${strict_num}/counts"
            Step3_sample_jobscript="${clusterFolder}/submission/count_${dataset}_${sample}_strict_${strict_num}.sh"
        fi        
        output="${countFolder}/${sample}_dexseq_counts.txt"
        mkdir -p "$countFolder"

        echo "#!/bin/bash
# Execute counting and save the raw file natively
python $pycount --stranded no -p ${paired_val} -f bam -r pos $GFF $bam ${output}

# Immediately strip out the 2-column HTSeq summary rows in-place
sed -i '/^_[a-z]/d' ${output}

echo "Step 3 finished for $sample at \$(date +%H:%M:%S)" >> $report_file 
" > "$Step3_sample_jobscript"
        
        echo "$Step3_sample_jobscript" >> "$Step3_taskfile"
    done
    echo "creating job scripts for read counting"
fi

####################
## STEP 4: DEXSeq ##
####################

sanity_check=$(awk 'NR == 1 {print $3}' "$support")

if [ "$DEXSeq" = "yes" ] && [ "$sanity_check" = "dataset" ]; then

    Step4_master_jobscript="${clusterFolder}/submission/Step4_DEXSeq_${protein}_${species}.sh"
    if [[ "$strict" == "yes" ]]; then 
        Step4_master_jobscript="${clusterFolder}/submission/Step4_DEXSeq_${protein}_${species}_strict_${strict_num}.sh"
    fi
    dataset_num=$(awk 'NR > 1 {print $3}' "$support" | uniq | wc -l)
    Step4_jobID="Step4_${protein}_${species}"
    if [[ "$strict" == "yes" ]]; then 
        Step4_jobID="Step4_${protein}_${species}_strict_${strict_num}"
    fi

    sample_num=$(awk 'NR > 1' "$support" | wc -l)

    echo "#!/bin/bash
#SBATCH --mem=24G
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=${clusterFolder}/out/Step4_%A_%a.out
#SBATCH --error=${clusterFolder}/error/Step4_%A_%a.err
#SBATCH --job-name=${Step4_jobID}
#SBATCH --chdir=${oFolder}
#SBATCH --array=1-${dataset_num}%30

module load R

if [[ \"\$SLURM_ARRAY_TASK_ID\" == \"1\" ]]; then
    echo \"Step4 started at \$(date +%H:%M:%S)\" >> $report_file
fi
jobs=\"" > "$Step4_master_jobscript"

    for dataset in $(awk 'NR > 1 {print $3}' "$support" | uniq); do
        DEXSeqFolder="${results}/${dataset}/dexseq"
        mkdir -p "$DEXSeqFolder"
        support_frame="${DEXSeqFolder}/${dataset}_dexseq_frame.tab"
        jobscript="${clusterFolder}/submission/DEXSeq_${dataset}.sh"
        GFF="${results}/${dataset}/GFF/${protein}_${species}_${dataset}.total.cryptics.gff"
        keepSex=TRUE
        keepDups=FALSE
        cryptic=TRUE
        iFolder="${results}/${dataset}"
        error_file="${clusterFolder}/R/dexseq_${dataset}.out"

        if [ "$strict" = "yes" ]; then
            iFolder="${results}/${dataset}/strict_${strict_num}"
            DEXSeqFolder="${iFolder}/dexseq"
            mkdir -p "$DEXSeqFolder"
            GFF="${results}/${dataset}/GFF/${protein}_${species}_${dataset}.strict.${strict_num}.total.cryptics.gff"
            jobscript="${clusterFolder}/submission/DEXSeq_${dataset}_strict_${strict_num}.sh"
            support_frame="${DEXSeqFolder}/${dataset}_dexseq_frame.tab"
            error_file="${clusterFolder}/R/dexseq_${dataset}_strict_${strict_num}.out"
        fi

        awk -v dataset="$dataset" 'NR == 1 {print $0} $3 == dataset {print $0}' "$support" > "$support_frame"
        bam_list=$(awk -v results="$results" -v dataset="$dataset" 'NR>1 {print results"/"dataset"/splice_extraction/"dataset"_"$1".spliced.exons.bam"}' "$support_frame" | tr '\n' '\t')
        sample_list=$(awk 'NR>1{print $1}' "$support_frame" | tr '\n' '\t')

        echo "#!/bin/bash
${Rbin}script ${dexseqFinalProcessR} --cryptic ${cryptic} --gff ${GFF} --keep.sex ${keepSex} --keep.dups ${keepDups} --support.frame ${support_frame} --code ${dataset} --annotation.file ${annotation_file} --iFolder ${iFolder} > ${error_file} 2>&1
cat ${error_file} >> ${report_file}

for bam in $bam_list; do
    samtools index \$bam
done

for i in \$(ls ${DEXSeqFolder}/); do
    comparison=\$(basename \$i)
    if [ -e ${DEXSeqFolder}/\$comparison/${dataset}_\${comparison}_CrypticExons.bed ]; then
        bedtools multicov -bed ${DEXSeqFolder}/\$comparison/${dataset}_\${comparison}_CrypticExons.bed -bams $bam_list -split > ${DEXSeqFolder}/\$comparison/${dataset}_\${comparison}_SJ_analysis.bed 
        echo -e \"chr\tstart\tend\tgene.id\tstrand\tintron.id\tlog2FC\tFDR\t$sample_list\" > ${DEXSeqFolder}/\$comparison/header       
        cat ${DEXSeqFolder}/\$comparison/header ${DEXSeqFolder}/\$comparison/${dataset}_\${comparison}_SJ_analysis.bed > tmp 
        mv tmp ${DEXSeqFolder}/\$comparison/${dataset}_\${comparison}_SJ_analysis.bed  
    fi
done
" > "$jobscript"

        echo "$jobscript" >> "$Step4_master_jobscript"
    done

    echo "\"
script=\$(echo \$jobs | cut -f\$SLURM_ARRAY_TASK_ID -d \" \")
bash \$script
" >> "$Step4_master_jobscript"
    echo "creating job scripts for DEXSeq testing"
fi

############################
### STEP 4b: SJ ANALYZER ###
############################

if [[ "$splice_junction_analyzer" == "yes" ]]; then
    Step4b_master_jobscript="${clusterFolder}/submission/Step4b_SJ_analyzer_${protein}_${species}.sh"
    dataset_num=$(awk 'NR > 1 {print $3}' "$support" | uniq | wc -l)
    Step4b_jobID="Step4b_${protein}_${species}"

    echo "#!/bin/bash
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=${clusterFolder}/out/Step4b_%A_%a.out
#SBATCH --error=${clusterFolder}/error/Step4b_%A_%a.err
#SBATCH --job-name=${Step4b_jobID}
#SBATCH --chdir=${oFolder}
#SBATCH --array=1-${dataset_num}%20

module load R

if [[ \"\$SLURM_ARRAY_TASK_ID\" == \"1\" ]]; then
    echo \"Step4b started at \$(date +%H:%M:%S)\" >> $report_file
fi
jobs=\"" > "$Step4b_master_jobscript"

    for dataset in $(awk 'NR > 1 {print $3}' "$support" | uniq); do
        DEXSeqFolder="${results}/${dataset}/strict_${strict_num}/dexseq"
        if [ ! -e "$DEXSeqFolder" ]; then
            DEXSeqFolder="${results}/${dataset}/dexseq"
        fi
        outFolder="${results}/${dataset}/splice_junction_analysis"
        mkdir -p "$outFolder"
        jobscript="${clusterFolder}/submission/SJ_analyzer_${dataset}.sh"
        
        echo "#!/bin/bash
# SPLICE JUNCTION ANALYSIS FOR ${dataset}" > "$jobscript"

        if [ -d "${results}/${dataset}/strict_500/dexseq/" ]; then
            for i in $(ls "${results}/${dataset}/strict_500/dexseq/"); do
                if [[ "$i" =~ .*[\.][a-z]+ ]]; then continue; fi
                condition_names=$(basename "$i")
                
                if [ -e "${DEXSeqFolder}/${condition_names}/${dataset}_${condition_names}_SignificantExons.csv" ]; then
                    dexseq_res="${DEXSeqFolder}/${condition_names}/${dataset}_${condition_names}_SignificantExons.csv"
                else 
                    continue
                fi
                
                support_frame="${outFolder}/${dataset}_support_frame.tab"
                error_file="${clusterFolder}/R/SJ_analyzer_${dataset}.out"
                awk -v dataset="$dataset" 'NR == 1 {print $0} $3 == dataset {print $0}' "$support" > "$support_frame"

                echo "${Rbin}script ${R_splice_junction_analyzer} --support.frame ${support_frame} --code ${dataset} --species $species --condition.names ${condition_names} --dexseq.res ${dexseq_res} --outFolder ${outFolder} > ${error_file} 2>&1" >> "$jobscript"
            done
        fi
        echo "$jobscript" >> "$Step4b_master_jobscript"
    done

    echo "\"
script=\$(echo \$jobs | cut -f\$SLURM_ARRAY_TASK_ID -d \" \")
bash \$script
" >> "$Step4b_master_jobscript"
fi

##################################
## STEP 4c: FUNCTIONAL ENRICHMENT
##################################

if [[ "$functional_enrichment" == "yes" ]]; then
    Step4c_master_jobscript="${clusterFolder}/submission/Step4c_Functional_Enrichment_${protein}_${species}.sh"
    dataset_num=$(awk 'NR > 1 {print $3}' "$support" | uniq | wc -l)
    Step4c_jobID="Step4c_${protein}_${species}"

    echo "#!/bin/bash
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=${clusterFolder}/out/Step4c_%A_%a.out
#SBATCH --error=${clusterFolder}/error/Step4c_%A_%a.err
#SBATCH --job-name=${Step4c_jobID}
#SBATCH --chdir=${oFolder}
#SBATCH --array=1-${dataset_num}%30

module load R

if [[ \"\$SLURM_ARRAY_TASK_ID\" == \"1\" ]]; then
    echo \"Step4c started at \$(date +%H:%M:%S)\" >> $report_file
fi
jobs=\"" > "$Step4c_master_jobscript"

    for dataset in $(awk 'NR > 1 {print $3}' "$support" | uniq); do
        SJAnalysisFolder="${results}/${dataset}/splice_junction_analysis"
        outFolder="${results}/${dataset}"
        jobscript="${clusterFolder}/submission/functional_enrichment_${dataset}.sh"
        
        echo "#!/bin/bash
# ENRICHMENT ANALYSIS FOR ${dataset}" > "$jobscript"

        if [ -d "${results}/${dataset}/strict_500/dexseq/" ]; then
            for i in $(ls "${results}/${dataset}/strict_500/dexseq/"); do
                if [[ "$i" =~ .*[\.][a-z]+ ]]; then continue; fi
                condition_names=$(basename "$i")
                SJAnalysis_res="${SJAnalysisFolder}/${dataset}_${condition_names}_splicing_analysis.tab"
                support_frame="${outFolder}/${dataset}_support_frame.tab"
                error_file="${clusterFolder}/R/Functional_Enrichment_${dataset}.out"
                awk -v dataset="$dataset" 'NR == 1 {print $0} $3 == dataset {print $0}' "$support" > "$support_frame"

                echo "${Rbin}script ${R_functional_enrichment} --support.frame ${support_frame} --code ${dataset} --species $species --condition.names ${condition_names} --outFolder ${outFolder} > ${error_file} 2>&1" >> "$jobscript"
            done
        fi
        echo "$jobscript" >> "$Step4c_master_jobscript"
    done

    echo "\"
script=\$(echo \$jobs | cut -f\$SLURM_ARRAY_TASK_ID -d \" \")
bash \$script
" >> "$Step4c_master_jobscript"
fi

#######################
### STEP 5: DESeq #####
#######################

sanity_check=$(awk 'NR == 1 {print $3}' "$support")

if [ "$DESeq" = "yes" ] && [ "$sanity_check" = "dataset" ]; then
    Step5_master_jobscript="${clusterFolder}/submission/Step5_DESeq_${protein}_${species}.sh"
    dataset_num=$(awk 'NR > 1 {print $3}' "$support" | uniq | wc -l)
    Step5_jobID="Step5_${protein}_${species}"

    echo "#!/bin/bash
#SBATCH --mem=16G
#SBATCH --time=9:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=${clusterFolder}/out/Step5_%A_%a.out
#SBATCH --error=${clusterFolder}/error/Step5_%A_%a.err
#SBATCH --job-name=${Step5_jobID}
#SBATCH --chdir=${oFolder}
#SBATCH --array=1-${dataset_num}%30

module load R

if [[ \"\$SLURM_ARRAY_TASK_ID\" == \"1\" ]]; then
    echo \"Step5 started at \$(date +%H:%M:%S)\" >> $report_file
fi
jobs=\"" > "$Step5_master_jobscript"

    for dataset in $(awk 'NR > 1 {print $3}' "$support" | uniq); do
        outFolder="${results}/${dataset}/expression"
        mkdir -p "$outFolder"
        support_frame="${outFolder}/${dataset}_deseq_frame.tab"
        jobscript="${clusterFolder}/submission/DESeq_${dataset}.sh"
        
        if [[ "$species" == "mouse" ]]; then
            GFF="/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Mouse/Mus_musculus.GRCm38.82_fixed.gff"
        elif [[ "$species" == "human" ]]; then
            GFF="/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Human_hg38/Homo_sapiens.GRCh38.82_fixed.gff"
        fi
        keepSex=TRUE; keepDups=FALSE; cryptic=TRUE

        dexseq_counts="${results}/${dataset}/counts"
        if [ ! -e "$dexseq_counts" ]; then
            dexseq_counts="${results}/${dataset}/strict_${strict_num}/counts"
        fi
        new_countFolder="${outFolder}/dexseq"
        mkdir -p "$new_countFolder"

        awk -v dataset="$dataset" 'BEGIN{OFS="\t"} NR == 1 {print $0} $3 == dataset {print $0}' "$support" > "$support_frame"
        Rscript "${R_support_chopper}" "${support_frame}"

        echo "#!/bin/bash
for countfile in \$(ls $dexseq_counts); do 
    awk '\$1 !~ /i/' ${dexseq_counts}/\$countfile > $new_countFolder/\$countfile
done
${Rbin}script ${countPrepareR} --gff ${GFF} --keep.dups ${keepDups} --support.frame ${support_frame} --code ${dataset} --annotation.file ${annotation_file} --iFolder ${outFolder} > ${clusterFolder}/R/prepare_counts_${dataset}.out 2>&1
cat ${clusterFolder}/R/prepare_counts_${dataset}.out >> $report_file

${Rbin}script ${deseqFinalProcessR} --keep.sex ${keepSex} --support.frame ${support_frame} --keep.dups ${keepDups} --code ${dataset} --annotation.file ${annotation_file} --iFolder ${outFolder} > ${clusterFolder}/R/deseq_${dataset}.out 2>&1
cat ${clusterFolder}/R/deseq_${dataset}.out >> $report_file
" > "$jobscript"

        echo "$jobscript" >> "$Step5_master_jobscript"
    done

    echo "\"
script=\$(echo \$jobs | cut -f\$SLURM_ARRAY_TASK_ID -d \" \")
bash \$script
" >> "$Step5_master_jobscript"
fi

##############################
### FINAL STEP: SUBMISSION ###
##############################

if [[ "$submit" == "yes" ]]; then
    hold=""
    
    # Standard log and directory parameters passed explicitly
    LOG_ARGS="--chdir=${oFolder} --output=${clusterFolder}/out/%x_%A_%a.out --error=${clusterFolder}/error/%x_%A_%a.err"
    
    if [[ "$splice_extractor" == "yes" ]]; then
        JOB_OUT=$(sbatch $LOG_ARGS --error="${report_file}" "$Step1_jobscript")
        Step1_jobID=$(echo "$JOB_OUT" | awk '{print $NF}')
        hold="--dependency=afterok:$Step1_jobID"
        echo "Submitted Step 1: $Step1_jobID"
    fi
    
    if [[ "$gff_creator" == "yes" ]]; then
        if [[ "$hold_Step1" == "yes" ]]; then
            Step1_jobID="Step1_${protein}_${species}"          
            hold="--dependency=afterok:$Step1_jobID"
        fi
        JOB_OUT=$(sbatch $hold $LOG_ARGS --error="${report_file}" "$Step2_jobscript")
        Step2_jobID=$(echo "$JOB_OUT" | awk '{print $NF}')
        hold="--dependency=afterok:$Step2_jobID"
        echo "Submitted Step 2: $Step2_jobID"
    fi
    
    if [[ "$read_counter" == "yes" ]]; then
        JOB_OUT=$(sbatch $hold $LOG_ARGS "$Step3_master_jobscript")
        Step3_jobID=$(echo "$JOB_OUT" | awk '{print $NF}')
        hold="--dependency=afterok:$Step3_jobID"
        echo "Submitted Step 3: $Step3_jobID"
    fi
    
    if [[ "$DESeq" == "yes" ]]; then
        sbatch $hold "$Step5_master_jobscript"
        echo "Submitted Step 5 (DESeq)"
    fi
    
    if [[ "$DEXSeq" == "yes" ]]; then
        JOB_OUT=$(sbatch $hold "$Step4_master_jobscript")
        Step4_jobID=$(echo "$JOB_OUT" | awk '{print $NF}')
        echo "Submitted Step 4 (DEXSeq): $Step4_jobID"
        
        if [[ "$splice_junction_analyzer" == "yes" ]]; then
            JOB_OUT=$(sbatch --dependency=afterok:$Step4_jobID "$Step4b_master_jobscript")
            Step4b_jobID=$(echo "$JOB_OUT" | awk '{print $NF}')
            echo "Submitted Step 4b (SJ Analyzer): $Step4b_jobID"
            
            if [[ "$functional_enrichment" == "yes" ]]; then
                sbatch --dependency=afterok:$Step4b_jobID "$Step4c_master_jobscript"
                echo "Submitted Step 4c (Functional Enrichment)"
            fi
        fi
    fi
fi

exit 0
