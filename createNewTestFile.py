import os

import numpy as np
from utils import ensure_dir
import params

CHARISOSMISET = {"#": 29, "%": 30, ")": 31, "(": 1, "+": 32, "-": 33, "/": 34, ".": 2,
                 "1": 35, "0": 3, "3": 36, "2": 4, "5": 37, "4": 5, "7": 38, "6": 6,
                 "9": 39, "8": 7, "=": 40, "A": 41, "@": 8, "C": 42, "B": 9, "E": 43,
                 "D": 10, "G": 44, "F": 11, "I": 45, "H": 12, "K": 46, "M": 47, "L": 13,
                 "O": 48, "N": 14, "P": 15, "S": 49, "R": 16, "U": 50, "T": 17, "W": 51,
                 "V": 18, "Y": 52, "[": 53, "Z": 19, "]": 54, "\\": 20, "a": 55, "c": 56,
                 "b": 21, "e": 57, "d": 22, "g": 58, "f": 23, "i": 59, "h": 24, "m": 60,
                 "l": 25, "o": 61, "n": 26, "s": 62, "r": 27, "u": 63, "t": 28, "y": 64}

ensure_dir("HyperAttentionDTI/tmp")
TEST_PATH = "HyperAttentionDTI/tmp/NewTest.txt"


def label_smiles(line, smi_ch_ind, MAX_SMI_LEN=100):
    X = np.zeros(MAX_SMI_LEN, dtype=np.int64())
    for i, ch in enumerate(line[:MAX_SMI_LEN]):
        X[i] = smi_ch_ind[ch]
    return X


def ex_smiles():
    fin = open("wdir/DrugBankNameID_SMILE.txt")
    fout = open("wdir/DrugBankNameSMILE_Filtered.txt", "w")
    while True:

        line = fin.readline()
        if line == "":
            break

        parts = line.strip().split("||")
        ss = parts[-1]
        print(ss)
        try:
            label_smiles(ss, CHARISOSMISET, 100)
            fout.write("%s" % line)
        except:
            print("Error", line)

    fin.close()
    fout.close()


def create_test(N_SEG=1):
    ex_smiles()

    from HyperAttentionDTI import dataset
    drug_names = sorted(list(dataset.DRUGNAME_2_SMILE.keys()))
    protein_ids = sorted(list(dataset.UNIPROT_2_SEQUENCE.keys()))
    # drug_names = ["n-[(3z)-5-tert-butyl-2-phenyl-1,2-dihydro-3h-pyrazol-3-ylidene]-n'-(4-chlorophenyl)urea","purvalanol","mgb-bp-3", "(4r)-n-[4-({[2-(dimethylamino)ethyl]amino}carbonyl)-1,3-thiazol-2-yl]-4-methyl-1-oxo-2,3,4,9-tetrahydro-1h-beta-carboline-6-carboxamide"]
    # protein_ids = protein_ids[:10]
    if N_SEG == 1:
        fout = open("%s%s" % (TEST_PATH, params.GEN_SUFFIX), "w")
        for protein in protein_ids:
            ss = ["%s||%s\n" % (drug, protein) for drug in drug_names]
            fout.write("%s" % "".join(ss))
        fout.close()
    else:
        SEG_SIZE = len(protein_ids) // N_SEG
        fx = open("%s%s___.Info" % (TEST_PATH, params.GEN_SUFFIX), "w")
        fx.write("%s" % SEG_SIZE)
        fx.close()
        for segId in range(N_SEG):
            startId = SEG_SIZE * segId
            endID = SEG_SIZE * (segId + 1)
            if segId == N_SEG - 1:
                endID = len(protein_ids)
            fout = open("%s%s___%s_%s" % (TEST_PATH, params.GEN_SUFFIX, segId, startId), "w")

            for protein in protein_ids[startId:endID]:
                ss = ["%s||%s\n" % (drug, protein) for drug in drug_names]
                fout.write("%s" % "".join(ss))
            fout.close()

        nHafl = N_SEG // 2
        for i in range(2):
            path = "%s/rd%s.sh" % (params.C_DIR, i)
            fo = open(path, "w")
            if i == 0:
                ar = ["%s" % j for j in range(0, nHafl)]
            else:
                ar = ["%s" % j for j in range(nHafl, N_SEG)]
            fo.write("for i in %s\n" % " ".join(ar))
            fo.write("do\n")
            fo.write("    python main.py -p -i $i -d cuda:%s\n" % i)
            fo.write("done\n")
            fo.close()
            os.system("chmod +x \"%s\"" % path)


if __name__ == "__main__":
    create_test(params.N_SEG)
