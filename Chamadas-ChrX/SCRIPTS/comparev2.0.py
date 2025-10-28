#!/usr/bin/env python3
import sys
import gzip
import os

# =============================
# ARGUMENTOS
# =============================
if len(sys.argv) != 4:
    print(f"Uso: {sys.argv[0]} <vcf_array> <vcf_wgs> <output_dir>")
    sys.exit(1)

vcf_array = sys.argv[1]   # VCF do array (chip)
vcf_wgs   = sys.argv[2]   # VCF do WGS (ploidy ou baseline)
output    = sys.argv[3]   # pasta de saída (ex.: results/D03_array_analysis)

# Criar pasta de saída
os.makedirs(output, exist_ok=True)

# Configurar log
logfile = os.path.join(output, "compare.log")
sys.stdout = open(logfile, "w")   # redireciona prints para log

print(f"Rodando comparação completa:")
print(f"Array = {vcf_array}")
print(f"WGS   = {vcf_wgs}")
print(f"Saída = {output}")
print("=" * 60)

# =============================
# FUNÇÕES AUXILIARES
# =============================
def open_vcf(path):
    """Abre VCF comprimido ou não"""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def load_variants(vcf_file):
    """Carrega variantes básicas (chr, pos, ref, alt)"""
    variants = set()
    with open_vcf(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            chrom, pos, ref, alt = parts[0], parts[1], parts[3], parts[4]
            variants.add((chrom, pos, ref, alt))
    return variants

def corrigir_genotipo(lista):
    """Função para corrigir genótipos haploides → diploides (ex: '0' → '0/0', '1' → '1/1')"""
    for linha in lista:
        gt_info = linha[-1].split(':')  # divide no ':'
        gt = gt_info[0]
        
        # Corrigir genótipos haploides para diploides
        if gt == '0':
            gt_info[0] = '0/0'  # homozigoto referência
        elif gt == '1':
            gt_info[0] = '1/1'  # homozigoto alternativo
        
        linha[-1] = ':'.join(gt_info)  # junta de volta
    return lista

# =============================
# ANÁLISE BÁSICA DE VARIANTES
# =============================
print("ETAPA 1: Análise básica de variantes")
print("-" * 40)

try:
    array_vars = load_variants(vcf_array)
    wgs_vars = load_variants(vcf_wgs)
    
    only_array = array_vars - wgs_vars
    only_wgs = wgs_vars - array_vars
    shared = array_vars & wgs_vars
    
    print(f"Total variantes Array = {len(array_vars)}")
    print(f"Total variantes WGS   = {len(wgs_vars)}")
    print(f"Compartilhadas        = {len(shared)}")
    print(f"Exclusivas Array      = {len(only_array)}")
    print(f"Exclusivas WGS        = {len(only_wgs)}")
    
    # Salvar arquivos básicos
    with open(os.path.join(output, "only_array.txt"), "w") as f:
        for v in sorted(only_array):
            f.write("\t".join(v) + "\n")
    
    with open(os.path.join(output, "only_wgs.txt"), "w") as f:
        for v in sorted(only_wgs):
            f.write("\t".join(v) + "\n")
    
    with open(os.path.join(output, "shared.txt"), "w") as f:
        for v in sorted(shared):
            f.write("\t".join(v) + "\n")
            
    print("Arquivos básicos salvos com sucesso!")
    
except Exception as e:
    print(f"ERRO na análise básica: {e}")
    sys.exit(1)

# =============================
# ANÁLISE DETALHADA COM BCFTOOLS
# =============================
print("\n" + "=" * 60)
print("ETAPA 2: Análise detalhada com bcftools")
print("-" * 40)

try:
    # Passo 1: Remover missing do array
    array_no_missing = os.path.join(output, "array_noMissing.vcf.gz")
    cmd = f"bcftools view -e 'GT=\"./.\"' -Oz -o {array_no_missing} {vcf_array} && tabix {array_no_missing}"
    print("Excluindo missing do Array:")
    print(cmd)
    result = os.system(cmd)
    if result != 0:
        print(f"ERRO ao remover missing: código {result}")
        sys.exit(1)

    # Passo 1.5: Remover missing do WGS
    wgs_no_missing = os.path.join(output, "wgs_noMissing.vcf.gz")
    cmd_wgs = f"bcftools view -e 'GT=\"./.\" || GT=\".\"' -Oz -o {wgs_no_missing} {vcf_wgs} && tabix {wgs_no_missing}"
    print("Excluindo missing do WGS:")
    print(cmd_wgs)
    result = os.system(cmd_wgs)
    if result != 0:
        print(f"ERRO ao remover missing do WGS: código {result}")
        sys.exit(1)
    
    # Passo 2: Intersecção com bcftools
    print("\nFazendo intersecção:")
    isec_cmd = f"bcftools isec {array_no_missing} {wgs_no_missing} -p {output}"
    print(isec_cmd)
    result = os.system(isec_cmd)
    if result != 0:
        print(f"ERRO na intersecção: código {result}")
        sys.exit(1)
    
    # Passo 3: Comprimir arquivos de intersecção
    print("\nComprimindo arquivos de intersecção:")
    bgzip_cmd = f"bgzip -f {output}/0002.vcf && bgzip -f {output}/0003.vcf"
    print(bgzip_cmd)
    result = os.system(bgzip_cmd)
    if result != 0:
        print(f"ERRO na compressão: código {result}")
        sys.exit(1)
        
except Exception as e:
    print(f"ERRO no processamento bcftools: {e}")
    sys.exit(1)

# =============================
# COMPARAÇÃO DE GENÓTIPOS
# =============================
print("\n" + "=" * 60)
print("ETAPA 3: Comparação detalhada de genótipos")
print("-" * 40)

try:
    # Carregar VCFs processados
    vcf1_path = os.path.join(output, "0002.vcf.gz")
    vcf2_path = os.path.join(output, "0003.vcf.gz")
    
    # Verificar se arquivos existem
    if not os.path.exists(vcf1_path) or not os.path.exists(vcf2_path):
        print("AVISO: Arquivos de intersecção não encontrados. Análise de genótipos não será realizada.")
        print("Isso pode acontecer se não houver variantes compartilhadas.")
    else:
        # Carregar e processar VCFs
        with gzip.open(vcf1_path, 'rt') as f:
            vcf1 = [line.rstrip().split() for line in f if not line.startswith('##')]
        
        with gzip.open(vcf2_path, 'rt') as f:
            vcf2 = [line.rstrip().split() for line in f if not line.startswith('##')]
        
        print(f"Carregados {len(vcf1)} linhas do array e {len(vcf2)} linhas do WGS")
        
        if len(vcf1) <= 1 or len(vcf2) <= 1:
            print("AVISO: Um dos VCFs está vazio ou contém apenas header")
            print("Total de variantes compartilhadas: 0")
            print("Taxa de concordância: N/A")
        else:
            # Corrigir genótipos
            vcf1 = corrigir_genotipo(vcf1)
            vcf2 = corrigir_genotipo(vcf2)
            
            # Comparar variantes
            ref_alt_diff = []
            geno_diff = []
            
            print(f"Comparando {len(vcf1)-1} variantes compartilhadas...")
            
            for i, j in zip(vcf1[1:], vcf2[1:]):
                info1 = i[1:5] + [i[-1].split(':')[0]]
                info2 = j[1:5] + [j[-1].split(':')[0]]  # fiz isso pq tem outras infos além do genótipo
                
                if info1[2:4] != info2[2:4]:
                    ref_alt_diff.append([info1[2:4], info2[2:4]])  # caso a ref/alt não bata
                elif info1[2:4] == info2[2:4]:
                    if info1[-1] != info2[-1]:
                        geno_diff.append([info1, info2])
            
            # Normalizar genótipos para comparação
            for par in geno_diff:
                for entrada in par:
                    entrada[-1] = entrada[-1].replace('|', '/')
                    entrada[-1] = entrada[-1].replace('0/1', '1/0')
            
            # Identificar erros reais
            erro = []
            for i in geno_diff:
                if i[0][-1] != i[1][-1]:
                    erro.append(i)
            
            # Salvar resultados detalhados
            header = ['POS_array', 'RSID_array', 'REF_array', 'ALT_array', 'GT_array',
                      'POS_wgs', 'RSID_wgs', 'REF_wgs', 'ALT_wgs', 'GT_wgs']
            
            # Salvar erros de genótipo
            with open(os.path.join(output, 'erros.txt'), 'w') as escreve:
                escreve.write('\t'.join(header) + '\n')  # cabeçalho
                for par in erro:
                    linha = '\t'.join(par[0] + par[1])   # concatena as duas sublistas
                    escreve.write(linha + '\n')
            
            # Salvar diferenças REF/ALT se existirem
            if len(ref_alt_diff) != 0:
                print('Posições de REF e ALT não batem para alguns genótipos. Elas foram escritas no arquivo ref_alt_diff.txt')
                with open(os.path.join(output, 'ref_alt_diff.txt'), 'w') as escreve:
                    escreve.write('\t'.join(header) + '\n')  # cabeçalho
                    for par in ref_alt_diff:
                        linha = '\t'.join(par[0] + par[1])   # concatena as duas sublistas
                        escreve.write(linha + '\n')
            
            # Calcular e exibir estatísticas finais
            total_shared = len(vcf1) - 1  # -1 por causa do header
            total_errors = len(erro)
            total_ref_alt_diff = len(ref_alt_diff)
            concordance_rate = ((total_shared - total_errors) / total_shared * 100) if total_shared > 0 else 0
            error_rate = (total_errors / total_shared * 100) if total_shared > 0 else 0
            
            print("\n" + "=" * 60)
            print("RESULTADOS FINAIS:")
            print("-" * 40)
            print(f"Total de variantes compartilhadas: {total_shared}")
            print(f"Diferenças de genótipo: {total_errors}")
            print(f"Diferenças REF/ALT: {total_ref_alt_diff}")
            print(f"Taxa de concordância: {concordance_rate:.2f}%")
            print(f"Taxa de erro: {error_rate:.2f}%")
            print("=" * 60)
            
            # Salvar resumo em arquivo separado
            with open(os.path.join(output, 'summary.txt'), 'w') as f:
                f.write("RESUMO DA COMPARAÇÃO\n")
                f.write("=" * 40 + "\n")
                f.write(f"Array VCF: {vcf_array}\n")
                f.write(f"WGS VCF: {vcf_wgs}\n")
                f.write(f"Data/Hora: {os.popen('date').read().strip()}\n")
                f.write("\nESTATÍSTICAS BÁSICAS:\n")
                f.write(f"Total variantes Array: {len(array_vars)}\n")
                f.write(f"Total variantes WGS: {len(wgs_vars)}\n")
                f.write(f"Variantes compartilhadas: {len(shared)}\n")
                f.write(f"Exclusivas Array: {len(only_array)}\n")
                f.write(f"Exclusivas WGS: {len(only_wgs)}\n")
                f.write("\nANÁLISE DE GENÓTIPOS:\n")
                f.write(f"Total compartilhadas analisadas: {total_shared}\n")
                f.write(f"Diferenças de genótipo: {total_errors}\n")
                f.write(f"Diferenças REF/ALT: {total_ref_alt_diff}\n")
                f.write(f"Taxa de concordância: {concordance_rate:.2f}%\n")
                f.write(f"Taxa de erro: {error_rate:.2f}%\n")

except Exception as e:
    print(f"ERRO na comparação de genótipos: {e}")
    import traceback
    print("Detalhes do erro:")
    traceback.print_exc()

print("\n" + "=" * 60)
print("ANÁLISE CONCLUÍDA!")
print(f"Todos os resultados foram salvos em: {output}")
print("Arquivos gerados:")
print("- compare.log (este log)")
print("- summary.txt (resumo completo)")
print("- only_array.txt (variantes exclusivas do array)")
print("- only_wgs.txt (variantes exclusivas do WGS)")
print("- shared.txt (variantes compartilhadas)")
print("- erros.txt (diferenças de genótipo)")
print("- ref_alt_diff.txt (diferenças REF/ALT, se houver)")
print("=" * 60)
