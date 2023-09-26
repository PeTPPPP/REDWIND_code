import numpy as np
from scipy.stats import norm
import pandas as pd



def fastCalPValue(selectedList):
    r"""

    Args:
        selectedList: All scores

    Returns:
        corresponding empirical p-values
    """
    if type(selectedList) == list:
        selectedList = np.asarray(selectedList)
    sortedIndices = np.argsort(selectedList)[::-1]
    nS = len(selectedList)
    # print(sortedIndices)

    trueIndices = np.zeros(nS)
    for i in range(nS-1, 0, -1):
        # print(i)
        trueIndices[sortedIndices[i]] = i
        # print(trueIndices)
        if i < nS-1:
            if selectedList[sortedIndices[i]] == selectedList[sortedIndices[i+1]]:
                trueIndices[sortedIndices[i]] = trueIndices[sortedIndices[i+1]]
    # print("trueIndices", trueIndices)
    trueIndices += 1.0
    # print("DB: ")
    # print(sortedIndices)
    # sortedIndices = np.arange(0, nS)
    # random.shuffle(sortedIndices)
    return trueIndices / nS



def retrieve_binding_drugs_for_genes(file_in, binding_file, nr_of_drugs_to_retrieve, file_out):
    fin = open(file_in, "r")
    fin_binding = open(binding_file, "r")

    dc = dict()
    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip().split("\t")
        dc.update({line[0]: "NONE"})


    while True:
        line = fin_binding.readline()
        if line == "":
            break
        line = line.strip().split(sep="||")
        gene_name = line[0]
        binding_inf = line[1]
        if gene_name in dc.keys():
            binding_pairs = binding_inf.split("\t")
            ls = []
            for i in range(nr_of_drugs_to_retrieve):
                binding_pair = binding_pairs[i].split("|")
                for elem in binding_pair:
                    ls.append(elem)
            dc[gene_name] = ls
    fin.close()
    fin_binding.close()
    fout = open(file_out, "w")
    for k,v in dc.items():
        name = k
        # print(name)
        # print(name)
        # print(line[0])
        # print(i[1])
        values = list(v)
        values = "\t".join(values)
        line = name + "\t" + values
        fout.writelines(line + "\n")
    # print(dc)
    return dc


def restrict_to_LE(file_in, fin_delim, file_out):
    df = pd.read_csv("/home/petschnerp/PycharmProjects/REDWIND/results/ALL_gene_list_with_values_onlyLE.txt", header = None, sep="\t")
    ls_LE = df.iloc[:,0]
    # print(list(ls_LE))
    df_binding = pd.read_csv(file_in, header = 0, delimiter=fin_delim)
    # print(list(df_binding.columns))
    na = [n.split("_")[0] for n in list(df_binding.columns)]
    bl = [n in list(ls_LE) for n in na[1:]]
    bl.insert(0, True)
    df_out = df_binding.loc[:,bl]
    diff = set(ls_LE) - set(na)
    print("Not found in Alphafold predictions: ", diff)
    df_out.to_csv(file_out, sep="\t", index = False, header=True)

## USED for Z.tests
def calc_testZx(num):
    return norm.cdf(num)

def calculate_Zstats_func(file_in, file_out):
    df = pd.read_csv(file_in, header = 0, sep = "\t", index_col=0)
    _z = np.vectorize(calc_testZx)
    df_out = _z(df)
    df_out = pd.DataFrame(df_out)
    df_out.columns = df.columns
    df_out.index = df.index
    # print(df_out)
    df_out.to_csv(file_out, index = True, header = True, sep = "\t")



if __name__ == "__main__":
    pass
