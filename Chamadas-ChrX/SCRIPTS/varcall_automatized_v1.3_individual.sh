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


THREADS=6

# Lista de amostras: ID sexo (F/M)
SAMPLES="samples.txt"

# Pastas
OUTRAW="/home/ec2-user/SSD4T/varcall_sjoint/rawvcf"
OUTREC="/home/ec2-user/SSD4T/varcall_sjoint/recalibrated"
OUTRES="/home/ec2-user/SSD4T/varcall_sjoint/results"
LOGS="/home/ec2-user/SSD4T/varcall_sjoint/logs"
COMPARE_SCRIPT="/home/ec2-user/SSD4T/comparev2.0.py"
LOGGERAL="$LOGS/G04.log"


mkdir -p "$OUTRAW" "$OUTREC" "$OUTRES" "$LOGS"

# zera log geral
> "$LOGGERAL"

# ======== FUNÇÕES ========

call_variants_array() {
    local sid="$1"
    local bam="$2"
    local array_vcf="$3"
    local vcf_array="$OUTRAW/${sid}_array.vcf.gz"

    gatk --java-options "-Xmx4g -Xms4g" HaplotypeCaller \
      -R "$REF" -I "$bam" -O "$vcf_array" \
      --native-pair-hmm-threads "$THREADS" \
      --output-mode EMIT_ALL_CONFIDENT_SITES \
      -L "$INTERVALS"
}

call_variants_ploidy() {
    local sid="$1"
    local bam="$2"
    local array_vcf="$3"
    local vcf_ploidy="$OUTRAW/${sid}_ploidy.vcf.gz"

    gatk --java-options "-Xmx4g -Xms4g" HaplotypeCaller \
      -R "$REF" -I "$bam" -O "$vcf_ploidy" \
      --native-pair-hmm-threads "$THREADS" \
      --ploidy-regions "$PLOIDY" \
      --output-mode EMIT_ALL_CONFIDENT_SITES \
      -L "$INTERVALS"
}

recalibrate() {
    local base="$1"     # ex.: D03_array  ou  D03_ploidy
    local vcf_in="$2"   # caminho do VCF a recalibrar

    local rec_snp="$OUTREC/${base}_snp.vcf.gz"
    local rec_indel="$OUTREC/${base}_indel.vcf.gz"

    # SNP
    gatk VariantRecalibrator \
      -R "$REF" -V "$vcf_in" \
      --resource:hapmap,known=false,training=true,truth=true,prior=15.0 "$HAPMAP" \
      --resource:omni,known=false,training=true,truth=false,prior=12.0 "$OMNI" \
      --resource:1000G,known=false,training=true,truth=false,prior=10.0 "$G1000" \
      --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "$DBSNP" \
      -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
      -mode SNP \
      -O "$OUTREC/${base}_snp.recal" \
      --tranches-file "$OUTREC/${base}_snp.tranches"

    gatk ApplyVQSR -R "$REF" -V "$vcf_in" \
      -O "$rec_snp" \
      --truth-sensitivity-filter-level 99.0 \
      --tranches-file "$OUTREC/${base}_snp.tranches" \
      --recal-file "$OUTREC/${base}_snp.recal" \
      -mode SNP

    # INDEL
    gatk VariantRecalibrator \
      -R "$REF" -V "$vcf_in" \
      -resource:mills,known=false,training=true,truth=true,prior=12.0 "$MILLS" \
      -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "$DBSNP" \
      -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
      -mode INDEL \
      -O "$OUTREC/${base}_indel.recal" \
      --tranches-file "$OUTREC/${base}_indel.tranches"

    gatk ApplyVQSR -R "$REF" -V "$vcf_in" \
      -O "$rec_indel" \
      --truth-sensitivity-filter-level 99.0 \
      --tranches-file "$OUTREC/${base}_indel.tranches" \
      --recal-file "$OUTREC/${base}_indel.recal" \
      -mode INDEL
}

merge_and_filter() {
    local base="$1"   # ex.: D03_array  ou  D03_ploidy
    local rec_snp="$OUTREC/${base}_snp.vcf.gz"
    local rec_indel="$OUTREC/${base}_indel.vcf.gz"
    local merged="$OUTRES/${base}.merged.vcf.gz"
    local pass_vcf="$OUTRES/${base}.PASS.vcf.gz"

    gatk MergeVcfs \
      -I "$rec_snp" \
      -I "$rec_indel" \
      -O "$merged"

    bcftools view -f PASS -Oz -o "$pass_vcf" "$merged"
    tabix -p vcf "$pass_vcf"
}

# ======== MAIN LOOP ========

while read -r SID SEX BAM ARRAYVCF; do
    echo ">>> Rodando $SID ($SEX)"
    
    # Criar flag de sucesso para esta amostra
    SAMPLE_SUCCESS=true
    
    # Desabilitar exit-on-error temporariamente para esta amostra
    set +e

    # nomes fixos para esta amostra
    VCF_ARRAY="$OUTRAW/${SID}_array.vcf.gz"
    VCF_PLOIDY="$OUTRAW/${SID}_ploidy.vcf.gz"

    if [[ "$SEX" == "M" ]]; then
        echo ">>> Processando homem: baseline (array) e ploidia-aware"

        # homens: baseline (array) e ploidia-aware
        if ! call_variants_array "$SID" "$BAM" "$ARRAYVCF"; then
            echo "ERRO: Falha na chamada de variantes array para $SID" >&2
            SAMPLE_SUCCESS=false
        fi

        if [[ "$SAMPLE_SUCCESS" == "true" ]] && ! call_variants_ploidy "$SID" "$BAM" "$ARRAYVCF"; then
            echo "ERRO: Falha na chamada de variantes ploidy para $SID" >&2
            SAMPLE_SUCCESS=false
        fi

        # recalibração para os dois
        if [[ "$SAMPLE_SUCCESS" == "true" ]] && ! recalibrate "${SID}_array" "$VCF_ARRAY"; then
            echo "ERRO: Falha na recalibração array para $SID" >&2
            SAMPLE_SUCCESS=false
        fi

        if [[ "$SAMPLE_SUCCESS" == "true" ]] && ! recalibrate "${SID}_ploidy" "$VCF_PLOIDY"; then
            echo "ERRO: Falha na recalibração ploidy para $SID" >&2
            SAMPLE_SUCCESS=false
        fi

        if [[ "$SAMPLE_SUCCESS" == "true" ]] && ! merge_and_filter "${SID}_array"; then
            echo "ERRO: Falha no merge array para $SID" >&2
            SAMPLE_SUCCESS=false
        fi

        if [[ "$SAMPLE_SUCCESS" == "true" ]] && ! merge_and_filter "${SID}_ploidy"; then
            echo "ERRO: Falha no merge ploidy para $SID" >&2
            SAMPLE_SUCCESS=false
        fi

        # rodar comparador apenas se tudo deu certo
        if [[ "$SAMPLE_SUCCESS" == "true" ]]; then
            echo ">>> Comparando array vs baseline para $SID"
            if python3 "$COMPARE_SCRIPT" "$ARRAYVCF" "$OUTRES/${SID}_array.PASS.vcf.gz" "$OUTRES/${SID}_array_analysis"; then
                echo ">>> Comparando array vs ploidy para $SID"
                python3 "$COMPARE_SCRIPT" "$ARRAYVCF" "$OUTRES/${SID}_ploidy.PASS.vcf.gz" "$OUTRES/${SID}_ploidy_analysis"
                
                # concatenar AMBOS os logs no log geral
                echo "=== ANÁLISE BASELINE (ARRAY) - $SID ===" >> "$LOGGERAL"
                cat "$OUTRES/${SID}_array_analysis/compare.log" >> "$LOGGERAL"
                echo -e "\n" >> "$LOGGERAL"

                echo "=== ANÁLISE PLOIDY - $SID ===" >> "$LOGGERAL"
                cat "$OUTRES/${SID}_ploidy_analysis/compare.log" >> "$LOGGERAL"
                echo -e "\n==============================\n" >> "$LOGGERAL"
                
                echo ">>> Amostra $SID processada completamente"
            else
                echo "ERRO: Falha nas comparações para $SID" >&2
                SAMPLE_SUCCESS=false
            fi
        fi

    elif [[ "$SEX" == "F" ]]; then
        echo ">>> Processando mulher: apenas baseline (diploide)"

        # mulheres: só baseline (diploide)
        if ! call_variants_array "$SID" "$BAM" "$ARRAYVCF"; then
            echo "ERRO: Falha na chamada de variantes para $SID" >&2
            SAMPLE_SUCCESS=false
        fi

        if [[ "$SAMPLE_SUCCESS" == "true" ]] && ! recalibrate "${SID}_array" "$VCF_ARRAY"; then
            echo "ERRO: Falha na recalibração para $SID" >&2
            SAMPLE_SUCCESS=false
        fi

        if [[ "$SAMPLE_SUCCESS" == "true" ]] && ! merge_and_filter "${SID}_array"; then
            echo "ERRO: Falha no merge para $SID" >&2
            SAMPLE_SUCCESS=false
        fi

        # rodar comparador só contra baseline
        if [[ "$SAMPLE_SUCCESS" == "true" ]]; then
            echo ">>> Comparando array vs baseline para $SID"
            if python3 "$COMPARE_SCRIPT" "$ARRAYVCF" "$OUTRES/${SID}_array.PASS.vcf.gz" "$OUTRES/${SID}_analysis"; then
                # concatenar no log geral
                echo "=== ANÁLISE BASELINE (ARRAY) - $SID ===" >> "$LOGGERAL"
                cat "$OUTRES/${SID}_analysis/compare.log" >> "$LOGGERAL"
                echo -e "\n==============================\n" >> "$LOGGERAL"
                
                echo ">>> Amostra $SID processada completamente"
            else
                echo "ERRO: Falha na comparação para $SID" >&2
                SAMPLE_SUCCESS=false
            fi
        fi
    fi

    # Log do resultado da amostra
    if [[ "$SAMPLE_SUCCESS" == "false" ]]; then
        echo "!!! AMOSTRA $SID FALHOU - continuando para próxima..." >&2
        echo "=== AMOSTRA FALHADA: $SID ($SEX) ===" >> "$LOGGERAL"
        echo "Erro durante processamento - verifique logs acima" >> "$LOGGERAL"
        echo -e "\n==============================\n" >> "$LOGGERAL"
    fi

    # Reabilitar exit-on-error para o próximo loop
    set -e

done < "$SAMPLES"



echo ""
echo "=== PIPELINE FINALIZADO ==="
echo "Data/Hora: $(date)"
echo "Log geral disponível em: $LOGGERAL"
echo ""
echo "=== ESTRUTURA DE RESULTADOS ==="
echo "Para cada amostra masculina você terá:"
echo "  - ${SID}_array_analysis/    (comparação baseline)"
echo "  - ${SID}_ploidy_analysis/   (comparação ploidy-aware)"
echo ""
echo "Para cada amostra feminina você terá:"
echo "  - ${SID}_analysis/    (comparação baseline)"
echo ""
echo "Cada pasta contém:"
echo "  - compare.log (log detalhado)"
echo "  - summary.txt (resumo)"
echo "  - erros.txt (diferenças de genótipo)"
echo "  - *.txt (outros arquivos de análise)"
echo "================================"
echo ""
echo "=== RESUMO FINAL ==="
echo "Verificando amostras processadas vs falhadas:"

SUCCESS_COUNT=0
FAILED_COUNT=0

while read -r SID SEX BAM ARRAYVCF; do
    if [[ -f "$OUTRES/${SID}_array_analysis/summary.txt" ]] || [[ -f "$OUTRES/${SID}_analysis/summary.txt" ]]; then
        echo "✓ $SID - SUCESSO"
        ((SUCCESS_COUNT++))
    else
        echo "✗ $SID - FALHOU"
        ((FAILED_COUNT++))
    fi
done < "$SAMPLES"

echo ""
echo "TOTAL: $((SUCCESS_COUNT + FAILED_COUNT)) amostras"
echo "SUCESSO: $SUCCESS_COUNT"
echo "FALHAS: $FAILED_COUNT"
echo "========================"
