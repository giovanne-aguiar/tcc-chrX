#!/bin/bash
set -euo pipefail

# ======== CONFIG ========
REF="/home/ec2-user/SSD/pipe/REF/H_sapiens38_noalt.fasta"
DBSNP="/home/ec2-user/SSD/pipe/REF/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
HAPMAP="/home/ec2-user/SSD/pipe/REF/hapmap_3.3.hg38.vcf.gz"
OMNI="/home/ec2-user/SSD/pipe/REF/1000G_omni2.5.hg38.vcf.gz"
G1000="/home/ec2-user/SSD/pipe/REF/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
MILLS="/home/ec2-user/SSD/pipe/REF/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
PLOIDY="/home/ec2-user/SSD/pipe/REF/ploidy_regions_XY.bed"
INTERVALS="/home/ec2-user/SSD4T/newArray_EPIGEN/intervals.bed"


THREADS=8

# Lista de amostras: ID sexo (F/M)
SAMPLES="samples.txt"

# Pastas
OUTRAW="/home/ec2-user/SSD4T/teste/rawvcf"
OUTREC="/home/ec2-user/SSD4T/teste/recalibrated"
OUTRES="/home/ec2-user/SSD4T/teste/results"
LOGS="/home/ec2-user/SSD4T/teste/logs"
JOINT_DIR="/home/ec2-user/SSD4T/teste/joint"
INDIVIDUAL_DIR="/home/ec2-user/SSD4T/teste/individual"
COMPARE_SCRIPT="/home/ec2-user/SSD4T/comparev2.0.py"
LOGGERAL="$LOGS/compare_all.log"

mkdir -p "$OUTRAW" "$OUTREC" "$OUTRES" "$LOGS" "$JOINT_DIR" "$INDIVIDUAL_DIR"

# zera log geral
> "$LOGGERAL"

# ======== FUNÇÕES ========

call_variants_array() {
    local sid="$1"
    local bam="$2"
    local sex="$3"
    local gvcf_out="$OUTRAW/${sid}.g.vcf.gz"

    echo ">> [${sid}] HaplotypeCaller (gVCF) - sex=${sex}"

    
    if [[ "$sex" == "M" ]]; then
        gatk --java-options "-Xmx4g" HaplotypeCaller \
          -R "$REF" -I "$bam" -O "$gvcf_out" \
          -ERC GVCF \
          --native-pair-hmm-threads "$THREADS" \
          --ploidy-regions "$PLOIDY" \
          --output-mode EMIT_ALL_CONFIDENT_SITES \
          -L "$INTERVALS"
    else
        gatk --java-options "-Xmx4g" HaplotypeCaller \
          -R "$REF" -I "$bam" -O "$gvcf_out" \
          -ERC GVCF \
          --native-pair-hmm-threads "$THREADS" \
          --output-mode EMIT_ALL_CONFIDENT_SITES \
          -L "$INTERVALS"
    fi

}

joint_genotyping() {
    local joint_vcf="$JOINT_DIR/joint.vcf.gz"
    local combined_gvcf="$JOINT_DIR/combined.g.vcf.gz"

    echo ">>> Combinando todos os gVCFs individuais." >&2

    gatk --java-options "-Xmx20g" CombineGVCFs \
      -R "$REF" \
      $(find "$OUTRAW" -name "*.g.vcf.gz" -printf "-V %p ") \
      -O "$combined_gvcf"

    echo ">>> Executando joint genotyping para todas as amostras do cromossomo X." >&2

    # roda GenotypeGVCFs com o gVCF combinado
    gatk --java-options "-Xmx20g" GenotypeGVCFs \
      -R "$REF" \
      -V "$combined_gvcf" \
      -O "$joint_vcf"

    echo "$joint_vcf"

}

recalibrate_joint() {
    local joint_vcf="$1"
    local base_name="joint_chrX"
    local rec_snp="$OUTREC/${base_name}_snp.vcf.gz"
    local rec_indel="$OUTREC/${base_name}_indel.vcf.gz"
    FILTERED_JOINT_VCF="$OUTREC/${base_name}_filtered.vcf.gz"

    echo ">>> Executando VQSR no VCF joint" >&2

    # Quick check: número de variantes total (para avaliar se VQSR faz sentido)
    local nv
    nv=$(bcftools view -H "$joint_vcf" | wc -l || true)
    echo "Número total de linhas (variants) no joint VCF: $nv" >&2
    if (( nv < 50000 )); then
        echo "AVISO: número de variantes ($nv) abaixo do recomendado para VQSR robusto." >&2
        echo "Considere usar hard-filters se o número for muito pequeno." >&2
    fi

    # SNP recalibration
    gatk VariantRecalibrator \
      -R "$REF" -V "$joint_vcf" \
      --resource:hapmap,known=false,training=true,truth=true,prior=15.0 "$HAPMAP" \
      --resource:omni,known=false,training=true,truth=false,prior=12.0 "$OMNI" \
      --resource:1000G,known=false,training=true,truth=false,prior=10.0 "$G1000" \
      --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "$DBSNP" \
      -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
      -mode SNP \
      -O "$OUTREC/${base_name}_snp.recal" \
      --tranches-file "$OUTREC/${base_name}_snp.tranches"

    gatk ApplyVQSR -R "$REF" -V "$joint_vcf" \
      -O "$rec_snp" \
      --truth-sensitivity-filter-level 99.0 \
      --tranches-file "$OUTREC/${base_name}_snp.tranches" \
      --recal-file "$OUTREC/${base_name}_snp.recal" \
      -mode SNP

    # INDEL recalibration (aplica sobre rec_snp)
    gatk VariantRecalibrator \
      -R "$REF" -V "$rec_snp" \
      -resource:mills,known=false,training=true,truth=true,prior=12.0 "$MILLS" \
      -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "$DBSNP" \
      -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
      -mode INDEL \
      -O "$OUTREC/${base_name}_indel.recal" \
      --tranches-file "$OUTREC/${base_name}_indel.tranches"

    gatk ApplyVQSR -R "$REF" -V "$rec_snp" \
      -O "$rec_indel" \
      --truth-sensitivity-filter-level 99.0 \
      --tranches-file "$OUTREC/${base_name}_indel.tranches" \
      --recal-file "$OUTREC/${base_name}_indel.recal" \
      -mode INDEL

    # Filtrar apenas PASS e indexar
    bcftools view -f PASS -Oz -o "$FILTERED_JOINT_VCF" "$rec_indel"
    tabix -p vcf "$FILTERED_JOINT_VCF"

}

extract_individual_samples() {
    local joint_filtered_vcf="$1"
    echo ">>> Separando amostras individuais do VCF joint filtrado" >&2

    while read -r SID SEX BAM ARRAYVCF; do
        local individual_vcf="$INDIVIDUAL_DIR/${SID}_chrX_filtered.vcf.gz"
        echo " - extraindo amostra $SID" >&2
        bcftools view -s "$SID" -Oz -o "$individual_vcf" "$joint_filtered_vcf"
        tabix -p vcf "$individual_vcf"
    done < "$SAMPLES"
}

run_comparisons() {
    echo ">>> Executando comparações individuais"

    while read -r SID SEX BAM ARRAYVCF; do
        local individual_vcf="$INDIVIDUAL_DIR/${SID}_chrX_filtered.vcf.gz"
        if [[ ! -f "$individual_vcf" ]]; then
            echo "Aviso: vcf individual não encontrado para $SID, pulando: $individual_vcf"
            continue
        fi

        echo ">>> Comparando array vs joint-filtered para $SID"
        mkdir -p "$OUTRES/${SID}_joint_analysis"
        python3 "$COMPARE_SCRIPT" "$ARRAYVCF" "$individual_vcf" "$OUTRES/${SID}_joint_analysis"

        # agregar logs
        echo "=== ANÁLISE JOINT-FILTERED - $SID ($SEX) ===" >> "$LOGGERAL"
        if [[ -f "$OUTRES/${SID}_joint_analysis/compare.log" ]]; then
            cat "$OUTRES/${SID}_joint_analysis/compare.log" >> "$LOGGERAL"
        fi
        echo -e "\n==============================\n" >> "$LOGGERAL"
    done < "$SAMPLES"
}

# ======== MAIN LOOP ========


echo "=== INICIANDO PIPELINE COM JOINT CALLING - CROMOSSOMO X ==="
echo "Data/Hora: $(date)"
echo "========================================================="

# ETAPA 1: Chamada individual de variantes no cromossomo X (gera gVCFs)
echo ""
echo "ETAPA 1: Chamada de variantes individual (gVCFs) - cromossomo X"
echo "---------------------------------------------------------------"

while read -r SID SEX BAM ARRAYVCF; do
    echo ">>> Chamando variantes para $SID ($SEX)"
    call_variants_array "$SID" "$BAM" "$SEX"
done < "$SAMPLES"

# ETAPA 2: Montar sample_map e executar joint genotyping
echo ""
echo "ETAPA 2: Joint Genotyping"
echo "-------------------------"

JOINT_VCF=$(joint_genotyping)
echo "Joint VCF criado: $JOINT_VCF"

# ETAPA 3: Aplicar VQSR no VCF joint
echo ""
echo "ETAPA 3: Aplicando VQSR no VCF joint"
echo "------------------------------------"

recalibrate_joint "$JOINT_VCF"
echo "VCF joint filtrado: $FILTERED_JOINT_VCF"

# ETAPA 4: Separar amostras individuais
echo ""
echo "ETAPA 4: Separando amostras individuais"
echo "-------------------------------------"

extract_individual_samples "$FILTERED_JOINT_VCF"

# ETAPA 5: Executar comparações
echo ""
echo "ETAPA 5: Executando comparações"
echo "-------------------------------"

run_comparisons



echo "=========================================================="
echo "=== PIPELINE FINALIZADO ==="
echo "Data/Hora: $(date)"
echo "Log geral disponível em: $LOGGERAL"
echo ""
echo "=== ESTRUTURA DO PIPELINE ==="
echo "Arquivos principais:"
echo "  - gVCFs (raw): $OUTRAW/"
echo "  - Joint VCF: $JOINT_DIR/joint_chrX.vcf.gz"
echo "  - Joint VCF filtrado: $FILTERED_JOINT_VCF"
echo "  - Amostras individuais extraídas: $INDIVIDUAL_DIR/"
echo "  - Comparações: $OUTRES/"
echo ""
echo "Cada pasta de resultado contém:"
echo "  - compare.log (log detalhado)"
echo "  - summary.txt (resumo)"
echo "  - erros.txt (diferenças de genótipo)"
echo "  - *.txt (outros arquivos de análise)"
echo "=========================================================="
