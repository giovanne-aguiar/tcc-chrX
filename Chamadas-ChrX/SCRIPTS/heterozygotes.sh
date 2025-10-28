#!/bin/bash
set -euo pipefail

# Diretório base dos resultados
BASE_DIR="/home/ec2-user/SSD4T/teste_default/results"
OUTPUT_FILE="/home/ec2-user/SSD4T/teste_default/heterozygotes/heterozygotes_joint_array.txt"

# Lista de amostras
SAMPLES_FILE="samples.txt"

# Criar cabeçalho do arquivo de saída
echo -e "sample_id\tsex\tchr\tpos\tref\talt\tgenotype" > "$OUTPUT_FILE"

echo "Iniciando extração de heterozigotos..."

# Loop através das amostras
while read -r SID SEX BAM ARRAYVCF; do
    echo "Processando amostra: $SID"
    
    # Caminho para o arquivo 0002.vcf.gz
    VCF_FILE="${BASE_DIR}/${SID}_joint_analysis/0002.vcf.gz"
    
    # Verificar se o arquivo existe
    if [[ ! -f "$VCF_FILE" ]]; then
        echo "AVISO: Arquivo não encontrado para $SID: $VCF_FILE"
        continue
    fi
    
    # Extrair posições com genótipos heterozigotos
    # Filtrar apenas genótipos 0/1, 1/0, 0|1, 1|0 (heterozigotos)
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' "$VCF_FILE" | \
    awk -v sid="$SID" -v sex="$SEX" '
    BEGIN { OFS="\t" }
    {
        # Verificar se é heterozigoto
        gt = $5
        if (gt == "0/1" || gt == "1/0" || gt == "0|1" || gt == "1|0") {
            print sid, sex, $1, $2, $3, $4, gt
        }
    }' >> "$OUTPUT_FILE"
    
    echo "  - $SID processado"
    
done < "$SAMPLES_FILE"

echo "Extração concluída!"
echo "Arquivo gerado: $OUTPUT_FILE"

# Mostrar estatísticas básicas
echo ""
echo "=== ESTATÍSTICAS ==="
total_variants=$(tail -n +2 "$OUTPUT_FILE" | wc -l)
echo "Total de variantes heterozigotas: $total_variants"

echo ""
echo "Heterozigotos por amostra:"
tail -n +2 "$OUTPUT_FILE" | cut -f1 | sort | uniq -c | sort -nr
