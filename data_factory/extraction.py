import numpy as np
import statistics

import pandas as pd

import params
import utils
import numpy
BINDING_THRESHOLD = params.BINDING_THRESHOLD
## Extracting gene expressions from files

def extract_DEGgenes(S1_file_in):
    fin = open(S1_file_in)
    gene_expression_dict = dict()
    cc = 0
    while True:

        line = fin.readline()
        if line == "":
            break
        cc += 1
        if cc == 1:
            continue
        else:
            line = line.split(sep="\t")

            gene_key = line[6]
            if gene_key == "":
                continue
            else:
                gene_expr_values = float(line[1])
                # print(gene_expr_values, type(gene_expr_values))
                if not params.GENEXP_LOG2:
                   gene_expr_values = np.power(2, gene_expr_values)

            gene_expression_dict.update({gene_key: gene_expr_values})
        print("Processing line: %s" % cc, end='\r')
    fin.close()
    return gene_expression_dict

def complement(S1_file_in, S2_file_in, mean=params.use_mean_of_expression_values):
    dcS1 = extract_DEGgenes(S1_file_in)
    dcS2 = extract_DEGgenes(S2_file_in)
    dcS1 = utils.expand_dict_by_complement_values(dcS1, dcS2, "NA")
    dcS3 = dict((k, [v] + ([dcS2[k]] if k in dcS2 else ["NA"])) for (k, v) in dcS1.items())
    if mean == True:
        for key in dcS3.keys():
            gene_expression_values = dcS3[key]
            if "NA" in gene_expression_values:
                output = 0
            else:
                gene_expression_values = [float(item) for item in dcS3[key]]
                if not params.USE_ABS_EXPR:
                    output = statistics.mean(gene_expression_values)
                else:
                    output = sum(gene_expression_values)

            # gene_expression_values.append(output)
            dcS3.update({key: output})
    else:
        for key in dcS3.keys():
            gene_expression_values = dcS3[key]
            if "NA" in gene_expression_values:
                output = 0
            else:
                gene_expression_values = [float(item) for item in dcS3[key]]
                output = sum(gene_expression_values)

            # gene_expression_values.append(output)
            dcS3.update({key: output})
    return dcS3


## Extracting gene-drug binding predictions
def restrict_drugs_to_selective_ones(file_in, file_in2, file_out):
    bfin = open(file_in, "r")
    rfin = open(file_in2, "r")
    restr_ls = []
    while True:
        line = rfin.readline()
        if line == "":
            break
        line = line.strip().split("\t")
        # print(line[1])
        if int(line[1]) > params.number_of_genes_to_bind_in_tops:
            restr_ls.append(line[0])
    # print(restr_ls)
    rfin.close()
    fout = open(file_out, "w")
    while True:
        line = bfin.readline()
        if line == "":
            break
        else:
            line = line.strip().split(sep="||")
            gene_name = line[0]
            # print(gene_name)
            binding_inf = line[1]
            binding_pairs = binding_inf.split("\t")
            binding_pairs_correct = []
            for bdp in binding_pairs:
                bdp_splitted = bdp.split("|")
                # print(bdp_splitted[0])
                if bdp_splitted[0] in restr_ls:
                    pass
                else:
                    binding_pairs_correct.append(bdp)

            binding_inf = "\t".join(binding_pairs_correct)
            line = gene_name + "||" + binding_inf + "\n"

            fout.write(line)
            # print(binding_pairs)
    fout.close()

def extract_genelist_druglist_frombinding(file_in):
    fin = open(file_in)
    gene_dict = dict()
    drug_dict = dict()
    gene_id_binding = dict()

    while True:
        line = fin.readline()
        if line == "":
            break
        else:
            line = line.strip().split(sep="||")
            gene_name = line[0]
            binding_inf = line[1]
            gene_id = utils.get_update_dict_index(gene_dict, gene_name)
            binding_pairs = binding_inf.split("\t")
            drug_ids = []
            scores = []
            for pair in binding_pairs:
                drug_name, score = pair.split("|")
                score = float(score)
                if score <= BINDING_THRESHOLD:
                    continue

                drug_id = utils.get_update_dict_index(drug_dict, drug_name)
                drug_ids.append(drug_id)
                scores.append(score)
            gene_id_binding[gene_id] = [drug_ids, scores]
    # print(gene_id_binding)
    return gene_dict, drug_dict, gene_id_binding


def convert_bindinginf_to_matrix(file_in):
    gene_dict, drug_dict, gene_id_binding = extract_genelist_druglist_frombinding(file_in)
    nrows = len(gene_dict)
    ncols = len(drug_dict)
    mat = np.zeros((nrows,ncols))
    for gene_id in range(nrows):
        mat_line = mat[gene_id, :]
        drug_ids, scores = gene_id_binding[gene_id]
        mat_line[drug_ids] = scores
    # print(gene_dict)
    if params.USE_DIFFERENCE:
        _ms = np.mean(mat, axis = 0, keepdims=True)
        _dv = np.std(mat, axis = 0, keepdims=True)
        mat = (mat - _ms) / _dv
    return mat, gene_dict, drug_dict

## Extracting gene pathway assocations from file

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
        _vec = _vec.reshape(len(_vec),1)
        arr = np.append(arr, _vec, axis=1)
    ## UNCOMMENT line below to see full gene list missing from predictions
    # print("Warning: gene not found in binding prediction:", list(missing_genes.keys()))

    return arr, pwys



## Extract pathway categories

def extract_pathway_categories(file_in):
    fin = open(file_in, "r")
    pathway_dict = dict()
    categories_dict = dict()
    pathway_categories_dict = dict()
    cc = 0
    while True:
        line = fin.readline()
        if line == "":
            break
        if cc != 0:
            line = line.strip().split("\t")
            pathway = line[0]
            pathway_id = utils.get_update_dict_index(pathway_dict, pathway)
            categories = line[2:]
            category_ids = []
            membership_vector = []
            for category in categories:
                category_id = utils.get_update_dict_index(categories_dict, category)
                category_ids.append(category_id)
                membership_vector.append(1)
            pathway_categories_dict[pathway_id] = [category_ids, membership_vector]
        cc += 1
    return pathway_dict, categories_dict, pathway_categories_dict

def convert_pathway_categories_to_matrix(file_in):
    pathway_dict, categories_dict, pathway_categories_dict = extract_pathway_categories(file_in)
    nrows = len(pathway_dict)
    # print(nrows)
    ncols = len(categories_dict)
    arr = np.zeros((nrows, ncols))
    for pathway_id in range(nrows):
        arr_line = arr[pathway_id, :]
        category_ids, memberships = pathway_categories_dict[pathway_id]
        # print(category_ids, memberships)
        arr_line[category_ids] = memberships
    # print(arr)
    arr = np.transpose(arr)
    return arr, categories_dict, pathway_dict

if __name__ == "__main__":
    pass