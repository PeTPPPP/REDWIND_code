import params
import utils
import main_internal.matrix_calculation_class as internal
import numpy as np
import main_internal.stats_int as st
import matplotlib.pyplot as plt
from optparse import OptionParser

def parse_config(opts):
    params.GEN_SUFFIX = opts.s
    params.I_SEG = opts.sid
    params.DEVICE = opts.device

def get_drug_sign(dat):
    dat = dat
    res = dat.results
    typ = ["L2"]
    nrw = len(dat.pwys_list) * res.shape[0]
    o = np.zeros((nrw, 3))
    for i, ty in enumerate(typ):
        pred_list = []
        for j, pwy in enumerate(dat.pwys_list):
            _res = res[:, j, i]
            pred_list = np.concatenate((pred_list, _res))

        r = st.fastCalPValue(pred_list)
        # print("SS", np.sum(r < 4.478681476173415e-06))
        o[:, i] = r
        # xvals = range(len(r))[:100]
        # yvals = np.log10(np.sort(r))[:100]
        # plt.scatter(xvals, yvals, color="orange")
        # plt.show()

    pwy_len = len(dat.pwys_list)
    drg_len = res.shape[0]
    drg_names = dat.drug_id_dict
    drg_names = sorted(drg_names.items(), key=lambda x: x[1])
    return o, drg_len, pwy_len, drg_names, dat.pwys_list


def get_sig_pwy_drug_values(sig, dat):
    p_df, drg_len, pwy_len, drg_names, pws_list = get_drug_sign(dat)
    typ = ["L2"]
    if sig is not None:
        if params.SIGN_MC_TYPE.lower() == "bonferroni":
            sig_threshold = params.SIGN / params.NR_TOPKDRUGS
            print("Significance threshold from params:", params.SIGN)
            print("Drug length:", drg_len)
            print("Significance threshold:", sig_threshold)
        # if params.SIGN_MC_TYPE.lower() == "bh":
        else:
            sig_threshold = params.SIGN
    res = np.zeros((drg_len, pwy_len, len(typ)))
    for i, ty in enumerate(typ):
        _int = p_df[:, i]
        _int = np.reshape(_int, (pwy_len, drg_len))
        _int = _int.transpose()
        res[:, :, i] = _int

    thresh_arr = np.full(res.shape, sig_threshold)
    sign_bool = np.less_equal(res, thresh_arr)
    return res, sign_bool, drg_names, pws_list


def extracting_sign_pathway_drugs(sig, dat):
    res, boo, drg_names, pws_list = get_sig_pwy_drug_values(sig, dat)
    res_sign = np.asarray(boo).nonzero()
    # print(res_sign)
    drg_output = []
    typ = ["L2"]
    o = []
    if res_sign != None:
        for i in range(len(res_sign[0])):
            triplets = [drg_names[res_sign[0][i]], pws_list[res_sign[1][i]], typ[res_sign[2][i]], res[res_sign[0][i], res_sign[1][i],res_sign[2][i]]]
            o.append(triplets)
            if typ[res_sign[2][i]] == "mean":
                drg_output.append(drg_names[res_sign[0][i]])
    return o

def run():
    dat_full = internal.data_calculations(params.DEG_GENES_RESULTS_ALL_S1_FILE, params.DEG_GENES_RESULTS_ALL_S2_FILE,
                                          params.GENES_PREDICTION_ALL_FILE, params.PATHWAY_LIST_ALL_FILE,
                                          params.PATHWAY_CATEGORIES_ALL_FILE)

    ###### FULL

    res = dat_full.results[:, :, 2]
    drug_names = utils.sort_dict(dat_full.drug_id_dict, onlykey=True, reverse=False)
    drug_names = [key for key in drug_names]
    print("DrugNames: ", drug_names[:10])
    file_out = "%s/full_data_array_L2.txt" % params.RESULTS_DIR
    print("Shape: ", res.shape, len(drug_names), len(dat_full.pwys_list))
    utils.write_nparray_to_file(res, rownames=drug_names, colnames=dat_full.pwys_list, file_out=file_out)

    triplets_full = extracting_sign_pathway_drugs(sig="", dat=dat_full)

    fout = open("%s/significant_results_%s.txt" % (params.RESULTS_DIR, round(params.SIGN / params.NR_TOPKDRUGS, 6)), "w")
    for ls in triplets_full:
        for i, l1 in enumerate(ls):
            if i == 0:
                fout.write("%s\t" % l1[0])
            else:
                fout.write("%s\t" % l1)
        fout.write("\n")
    fout.close()


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-d", "--device", dest="device", type='str', default='cuda:0', help="Device")

    parser.add_option("-s", "--suffix", dest="s", type='str', default='',
                      help="Gene list type: '' for all, 'M' for male, 'F' for female")
    parser.add_option("-c", "--create", dest="c", action='store_true', help="Flag for creating a test file")
    parser.add_option("-p", "--predict", dest="p", action='store_true', help="Flag for predicting")
    parser.add_option("-i", "--sid", dest="sid", type='int', default=-1, help="Segment id")

    parser.add_option("-m", "--merge", dest="m", action='store_true', help="Merge prediction")

    parser.add_option("-x", "--export", dest="x", action='store_true', help="Flag for exporting scores")
    parser.add_option("-w", "--ways", dest="w", action='store_true', help="Flag for running pathways")

    (options, args) = parser.parse_args()
    parse_config(options)
    if options.c:
        from exportSelectedGeneProteins import match
        from createNewTestFile import create_test

        match()
        create_test(params.N_SEG)
    elif options.p:
        from HyperAttentionDTI.predicting import run_predict

        print("Run prediction...")
        run_predict()
    elif options.m:
        from HyperAttentionDTI.predicting import merge

        print("Run merging...")
        merge()

    elif options.x:
        from exportSelectedGeneProteins import export_gene_drug

        print("Exporting gene-drug scores...")
        export_gene_drug()

    elif options.w:
        print("Calculating pathway level binding predictions...")
        run()
