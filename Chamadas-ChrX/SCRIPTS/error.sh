#!/bin/bash
set -euo pipefail

# Arquivo de saída consolidado
OUTPUT_FILE="/home/ec2-user/SSD4T/all_errors_consolidated.txt"

# Lista de amostras
SAMPLES_FILE="samples.txt"

# Diretórios base para cada metodologia
JOINT_PLOIDY_DIR="/home/ec2-user/SSD4T/teste/results"
JOINT_DEFAULT_DIR="/home/ec2-user/SSD4T/teste_default/results"

# Criar cabeçalho do arquivo consolidado
echo -e "sample_id\tsex\tmethod\tpos_array\trsid_array\tref_array\talt_array\tgt_array\tpos_wgs\trsid_wgs\tref_wgs\talt_wgs\tgt_wgs" > "$OUTPUT_FILE"

echo "Iniciando consolidação dos arquivos de erros..."

# Função para processar um diretório de metodologia
process_methodology() {
    local base_dir="$1"
    local method_name="$2"
    local suffix="$3"  # sufixo específico para cada metodologia (ex: "_joint_analysis")
    
    echo "Processando metodologia: $method_name"
    
    while read -r SID SEX BAM ARRAYVCF; do
        local errors_file="${base_dir}/${SID}_joint_analysis/erros.txt"
        
        if [[ -f "$errors_file" ]]; then
            # Verificar se arquivo não está vazio (além do cabeçalho)
            local line_count=$(wc -l < "$errors_file")
            if [[ $line_count -gt 1 ]]; then
                echo "  Processando $SID ($method_name): $(($line_count - 1)) erros"
                
                # Pular cabeçalho (primeira linha) e adicionar informações da amostra
                tail -n +2 "$errors_file" | \
                awk -v sid="$SID" -v sex="$SEX" -v method="$method_name" '
                BEGIN { OFS="\t" }
                {
                    print sid, sex, method, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10
                }' >> "$OUTPUT_FILE"
            else
                echo "  $SID ($method_name): sem erros"
            fi
        else
            echo "  AVISO: Arquivo de erros não encontrado para $SID em $method_name: $errors_file"
        fi
    done < "$SAMPLES_FILE"
}

# Processar cada metodologia
echo ""
process_methodology "$JOINT_PLOIDY_DIR" "joint_ploidy" "_joint_analysis"
echo ""
process_methodology "$JOINT_DEFAULT_DIR" "joint_default" "_joint_analysis"

echo ""
echo "Consolidação concluída!"
echo "Arquivo gerado: $OUTPUT_FILE"

# Estatísticas gerais
echo ""
echo "=== ESTATÍSTICAS GERAIS ==="
total_errors=$(tail -n +2 "$OUTPUT_FILE" | wc -l)
echo "Total de erros consolidados: $total_errors"

echo ""
echo "Script finalizado com sucesso!"
echo "Use o arquivo '$OUTPUT_FILE' para análises posteriores em R/Python."
