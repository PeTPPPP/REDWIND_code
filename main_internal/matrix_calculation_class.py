import numpy as np

import data_factory.extraction as extraction
import params
import utils


class data_calculations:
    def __init__(self, genexp_S1_file_in, geneexp_S2_file_in, gene_drug_binding_file_in, gene_pathway_file, categories_pathways_file_in):
        self.genexp_S1_file_in, self.geneexp_S2_file_in, self.gene_drug_binding_file_in, self.gene_pathway_file, self.categories_pathways_file_in = genexp_S1_file_in, geneexp_S2_file_in, gene_drug_binding_file_in, gene_pathway_file, categories_pathways_file_in
        ## Initializing a dictionary for gene expression with keys: genes, values: mean/sum of log2 gene expression values
        self.use_mean = params.use_mean_of_expression_values
        self.normalize = params.NORMALIZE
        self.gene_expression_dict = extraction.complement(genexp_S1_file_in, geneexp_S2_file_in, self.use_mean)
        ## Initializing gene - drug binding predictions array with rows: genes, columns: drugs, also dictionaries for index of genes and of drugs in array
        if params.restrict_to_specific_drugs:
            extraction.restrict_drugs_to_selective_ones(gene_drug_binding_file_in, params.PREDICTION_COUNTS_FILE,
                                                        params.GENES_PREDICTION_CORRECTED_FILE)
            gene_drug_binding_file_in = params.GENES_PREDICTION_CORRECTED_FILE
        self.gene_drugbinding_arr, self.gene_id_dict, self.drug_id_dict = extraction.convert_bindinginf_to_matrix(gene_drug_binding_file_in)
        if params.restrict_to_specific_drugs == False:
            print("WARNING: no filtering for multi-binding drugs!")

        ## Initializing gene - pathway association array with rows: genes, columns: pathways and a list of the pathways
        self.gene_pathwy_arr, self.pwys_list = extraction.extract_gene_pathway_associations(gene_pathway_file, self.gene_id_dict)
        ## Initializing pathway - categories array with rows: categories, columns: pathways and dictionaries of category ID and pathway IDs
        self.categories_pwy_arr, self.category_id_dict, self.pathway_id_dict = extraction.convert_pathway_categories_to_matrix(categories_pathways_file_in)
        self.re_order()

        self.category_id_list = [k for (k, _) in utils.sort_dict(self.category_id_dict)[::-1]]
        self.normalize_bygenenr = params.NORMALIZE_RESULTS_BYGENENR
        self.number_of_genes_per_pathway = np.sum(self.gene_pathwy_arr, axis=0)
        self.results = self.generate_gene_expr_pwy_matrix()
        self.nonweighted_results = self.generate_non_weighted_matrix()

    def re_order(self):
        old_pw_ids = []
        # print(self.pwys_list)
        # print(self.pathway_id_dict)
        # exit(-1)
        assert len(self.pwys_list) == len(self.pathway_id_dict)
        new_col_id = 0
        new_pw_name_list = []
        new_pw_dict = dict()
        for pw_name in self.pwys_list:
            col_id = utils.get_dict(self.pathway_id_dict, pw_name, -1)
            if col_id == -1:
                continue
            new_pw_dict[pw_name] = new_col_id
            old_pw_ids.append(col_id)
            new_col_id += 1
            new_pw_name_list.append(pw_name)
        self.categories_pwy_arr = self.categories_pwy_arr[:, old_pw_ids]
        self.pathway_id_dict = new_pw_dict
        self.pwys_list = new_pw_name_list
        file_out = "%s/pathways_cate_tmp.csv" % params.RES_DIR
        utils.write_nparray_to_file(self.categories_pwy_arr, rownames=self.category_id_dict, colnames=self.pwys_list, file_out=file_out)

    def generate_normalization_tensor(self):

        nr_of_drugs = len(self.drug_id_dict.keys())
        # print(nr_of_drugs)
        norm_matrix = np.tile(self.number_of_genes_per_pathway, reps=[nr_of_drugs,1])
        # print(norm_matrix.shape)
        norm_tensor = np.stack([norm_matrix, norm_matrix, norm_matrix], axis=2)
        return norm_tensor


    def __str__(self):
        dims_of_drugbinding_arr = "Nr. of genes: %s\nNr. of drugs: %s\n" % (self.gene_drugbinding_arr.shape[0], self.gene_drugbinding_arr.shape[1])
        dims_of_gene_pathways_arr = "Nr. of genes: %s\nNr. of pathways: %s\n" % (self.gene_pathwy_arr.shape[0], self.gene_pathway_file[1])
        dims_of_categories_pathway_arr = "Nr. of categories: %s\nNr. of pathways: %s\n" % (self.categories_pwy_arr.shape[0], self.categories_pwy_arr.shape[1])
        text_to_print = "DATA CLASS\nInput files: %s\n%s\n%s\n%s\n%s\nDimensions of arrays:\nGene - drug binding array: " \
                        "%s\nGene - pathway membership array: %s\nCategory - pathway membership array: %s\n" % \
                        (self.genexp_S1_file_in, self.geneexp_S2_file_in, self.gene_drug_binding_file_in,
                         self.gene_pathway_file, self.categories_pathways_file_in, dims_of_drugbinding_arr,
                         dims_of_gene_pathways_arr, dims_of_categories_pathway_arr)
        return text_to_print



    def create_gene_expression_pathway_matrix(self):
        gene_expression_dict = self.gene_expression_dict
        gene_pathway_arr = self.gene_pathwy_arr

        gene_id_dict = self.gene_id_dict
        # print(gene_id_dict)
        gene_id_list = [k for (k, _) in utils.sort_dict(gene_id_dict)[::-1]]
        # print(gene_id_list)
        nrows = gene_pathway_arr.shape[0]
        ncols = gene_pathway_arr.shape[1]

        _arr = np.zeros((nrows, ncols))
        for i, gene in enumerate(gene_id_list):
            expression = gene_expression_dict[gene]
            _vec = np.repeat(expression, ncols)
            _arr[i] = _vec


        assert _arr.shape == gene_pathway_arr.shape
        output_arr = np.multiply(_arr, gene_pathway_arr)
        if self.normalize:
            output_arr_max = output_arr.flat[abs(output_arr).argmax()]
            output_arr = output_arr/output_arr_max
        file_out = "%s/tmp_gene_path.txt" % params.TMP_DIR
        colnames = list(self.pwys_list)
        utils.write_nparray_to_file(output_arr, rownames=gene_id_list, colnames=colnames,
                                    file_out=file_out)
        return output_arr

    @staticmethod
    def get_mean_of_values(arr):
        means = np.mean(arr, axis=1)
        return means

    @staticmethod
    def get_L1_norms(arr):
        L1_norm = np.linalg.norm(arr, ord=1, axis=1)
        return L1_norm

    @staticmethod
    def get_L2_norms(arr):
        L2_norm = np.linalg.norm(arr, axis=1)
        return L2_norm

    def generate_non_weighted_matrix(self):
        ## NOTE: gene_expression_pathway_arr should be/is organized: rows: genes, cols: expressions
        ## NOTE: gene_drug_binding_arr should be/is organized: rows: genes, cols: drugs
        ## transposing these into required format
        gene_pathway_arr = self.gene_pathwy_arr
        # print(gene_pathway_arr[1,:])
        gene_pathway_arr = np.transpose(gene_pathway_arr)
        nr_of_pwys = gene_pathway_arr.shape[0]
        gene_drug_binding_arr = self.gene_drugbinding_arr
        gene_drug_binding_arr = np.transpose(gene_drug_binding_arr)
        # print(gene_drug_binding_arr)
        if params.NORMALIZE:
            # print(gene_drug_binding_arr).argmax()
            gene_drug_binding_arr_max = gene_drug_binding_arr.flat[abs(gene_drug_binding_arr).argmax()]
            gene_drug_binding_arr = gene_drug_binding_arr/gene_drug_binding_arr_max
        nr_of_drugs = gene_drug_binding_arr.shape[0]

        self.nonweighted_results = np.empty((nr_of_drugs, nr_of_pwys, 3))
        for i in range(nr_of_drugs):
            _tmp_arr_per_drug = np.tile(gene_drug_binding_arr[i], reps=[nr_of_pwys, 1])
            results_for_one_drug = np.multiply(_tmp_arr_per_drug, gene_pathway_arr)
            mean_v = self.get_mean_of_values(results_for_one_drug)
            L1_norm = self.get_L1_norms(results_for_one_drug)
            L2_norm = self.get_L2_norms(results_for_one_drug)
            output_for_one_drug = np.column_stack([mean_v, L1_norm, L2_norm])
            self.nonweighted_results[i] = output_for_one_drug

        if self.normalize_bygenenr == True:
            # print(self.normalize_bygenenr)
            norm_tensor = self.generate_normalization_tensor()
            self.nonweighted_results = np.divide(self.nonweighted_results, norm_tensor)
        return self.nonweighted_results


    def generate_gene_expr_pwy_matrix(self):
        ## NOTE: gene_expression_pathway_arr should be/is organized: rows: genes, cols: expressions
        ## NOTE: gene_drug_binding_arr should be/is organized: rows: genes, cols: drugs
        ## transposing these into required format
        gene_expression_pathway_arr = self.create_gene_expression_pathway_matrix()
        gene_expression_pathway_arr = np.transpose(gene_expression_pathway_arr)
        nr_of_pwys = gene_expression_pathway_arr.shape[0]
        gene_drug_binding_arr = self.gene_drugbinding_arr
        gene_drug_binding_arr = np.transpose(gene_drug_binding_arr)
        # print(gene_drug_binding_arr)
        if params.NORMALIZE:
            # print(gene_drug_binding_arr).argmax()
            gene_drug_binding_arr_max = gene_drug_binding_arr.flat[abs(gene_drug_binding_arr).argmax()]
            gene_drug_binding_arr = gene_drug_binding_arr/gene_drug_binding_arr_max
        nr_of_drugs = gene_drug_binding_arr.shape[0]

        self.results = np.empty((nr_of_drugs, nr_of_pwys, 3))
        for i in range(nr_of_drugs):
            _tmp_arr_per_drug = np.tile(gene_drug_binding_arr[i], reps=[nr_of_pwys, 1])
            results_for_one_drug = np.multiply(_tmp_arr_per_drug, gene_expression_pathway_arr)
            # print(results_for_one_drug)
            mean_v = self.get_mean_of_values(results_for_one_drug)
            L1_norm = self.get_L1_norms(results_for_one_drug)
            L2_norm = self.get_L2_norms(results_for_one_drug)
            output_for_one_drug = np.column_stack([mean_v, L1_norm, L2_norm])
            self.results[i] = output_for_one_drug

        if self.normalize_bygenenr == True:
            # print(self.normalize_bygenenr)
            norm_tensor = self.generate_normalization_tensor()
            self.results = np.divide(self.results, norm_tensor)
        return self.results




