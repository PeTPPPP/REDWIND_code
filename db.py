import numpy as np
import data_factory.extraction as extraction
import params
import utils
import main_internal.matrix_calculation_class as internal
from main_internal.stats_int import fastCalPValue

def ttest(dc):
    drug_id_dict = dc
    drug_name_list = [k for (k, _) in utils.sort_dict(drug_id_dict)[::-1]]
    print(drug_name_list)


def extract_gene_pathway_associations(file_in, gene_id_dict):
    fin = open(file_in)
    nrow = len(gene_id_dict)
    arr = np.zeros((nrow, 0))
    score = 1
    pwys = []
    missing_genes = dict()
    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip().split("\t")
        pwy_name = line[0]
        pwys.append(pwy_name)
        genes = line[1:]
        _vec = np.zeros(nrow)
        for gene in genes:
            if gene in gene_id_dict.keys():
                gene_id = gene_id_dict[gene]
                _vec[gene_id] = score
            else:
                missing_genes.update({gene: "NA"})
        _vec = _vec.reshape(len(_vec), 1)
        arr = np.append(arr, _vec, axis=1)
    print("Warning: gene not found in binding prediction:", list(missing_genes.keys()))
    print(arr)
    return arr, pwys


if __name__ == "__main__":
    pass
