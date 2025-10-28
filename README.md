# Avaliação de metodologias para chamada de variantes a partir de dados de Sequenciamento de Nova Geração (NGS) do cromossomo X humano

**Autor:** Giovanne de Aguiar Barra

**Orientadores:** Rafael Tou
			  Giovanna Giudicelli
			  Eduardo Tarazona
			  
**Projeto:** Trabalho de Conclusão de Curso (TCC)

## Descrição

Este repositório contém os scripts, os principais arquivos de referência utilizados, o pipeline utilizado e os resultados para a avaliação de metodologias de chamada de variantes genéticas no cromossomo X humano. A análise foi realizada a partir de dados de Sequenciamento de Nova Geração (NGS) de 30 amostras do projeto EPIGEN-Brasil.

O principal objetivo foi testar e comparar abordagens de bioinformática para lidar com as particularidades da genotipagem do cromossomo X.

## As seguintes ferramentas foram utilizadas para conduzir as análises:
* bcftools (versão 1.21)
* samtools (versão 1.17)
* fastqc (versão 0.12.1)
* multiqc (versão 1.29)
* Trimmomatic (versão 0.39)
* BWA (versão 0.7.15-r1140)
* GATK (versão 4.6.1.0)

* Python (versão 3.10.19)


## Estrutura do Repositório

-   `/REFERENCE`: Contém todas os arquivos referência utilizados.
-   `/RESULTS`: Contém exemplos das tabelas geradas e do log das chamadas.
-   `/SCRIPTS`: Contém todos os scripts e comandos utilizados.