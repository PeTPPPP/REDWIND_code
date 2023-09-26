# General

1. Clone repo
2. Unzip '.7z' files in wdir and resource/Genes2Drug/9657Genes2Drugs

# Gene2Drug
This section describes prediction of protein binding prediction of drugs using differentially expressed genes in the transcriptomic experiment

## Input-Output
- Input: A list of gens at ./data/gene_list.txt
  - Format: Each line is a GeneName

- Output: Potentially associated drugs with each gen at ./wdir/Gene2Drugs.txt
  - Format: Each line contains predicted associated scores of drug with a gene: 
    - GeneName||drug1Id|score drug2Id|score...
## Running
### Training
```shell
    cd HyperAttentionDTI
    python HyperAttentionDTI_main.py
    cd ..
```
### Creating a dataset for prediction from a given gene list
```shell
   cd wdir
   unzip UniprotReviewedMapping2.zip
   cd ..
   python main.py -c
```

### Prediction
```shell
    python main.py -p

```
### Exporting predicted gene-drug scores
```shell
    python main.py -x

```

# Pathway2Drug

This section describes prediction of protein binding prediction of drugs to pathways using differentially expressed genes of transcriptomic experiment

## Input-Output
- Input: A Genes2Drugs.txt containing gene-drug predicted binding scores
  - Format: Potentially associated drugs with each gen at ./wdir/Gene2Drugs.txt
- Input2: ALL_pathway_gene_list_total.txt containing pathway-gene associations
  - Format: One pathway per line PathwayName&nbsp;   gene1&nbsp;   gene2 ...
- Input3: DEG_results_S1_all_SAAS.csv and DEG_results_S2_all_SAAS.csv containing expression values for all genes as output from edgeR analysis at the two time points
  - Format: ID&nbsp;   logFC&nbsp;   logCPM&nbsp;   F&nbsp;   PValue&nbsp;   names&nbsp;   symbol&nbsp;   entrezid&nbsp;   FDR
- Input4: Top10Stats.txt containing the number of bindings in the top genes of drugs (for filtering)
  - Format: Drug name&nbsp;   Number of top10 occurences

- Output: significant_results_0.0005.txt containing significant L2 scores for given drug and pathway combinations
  - Format: Each line contains predicted associated scores of drug with a pathway: 
    - DrugName&nbsp;   Pathway&nbsp;   score type (L2)&nbsp;   empirical p-value
    
## Running
### Calculating scores and significance
```shell
    python main.py -w

```


