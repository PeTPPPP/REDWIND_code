import torch
from torch.utils.data import Dataset, DataLoader
import random
import numpy as np
import os
import sys
import subprocess

C_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(C_DIR)
import sys
[sys.path.append(i) for i in ['.', '..']]
import conf

CHARISOSMISET = {"#": 29, "%": 30, ")": 31, "(": 1, "+": 32, "-": 33, "/": 34, ".": 2,
                 "1": 35, "0": 3, "3": 36, "2": 4, "5": 37, "4": 5, "7": 38, "6": 6,
                 "9": 39, "8": 7, "=": 40, "A": 41, "@": 8, "C": 42, "B": 9, "E": 43,
                 "D": 10, "G": 44, "F": 11, "I": 45, "H": 12, "K": 46, "M": 47, "L": 13,
                 "O": 48, "N": 14, "P": 15, "S": 49, "R": 16, "U": 50, "T": 17, "W": 51,
                 "V": 18, "Y": 52, "[": 53, "Z": 19, "]": 54, "\\": 20, "a": 55, "c": 56,
                 "b": 21, "e": 57, "d": 22, "g": 58, "f": 23, "i": 59, "h": 24, "m": 60,
                 "l": 25, "o": 61, "n": 26, "s": 62, "r": 27, "u": 63, "t": 28, "y": 64}

CHARISOSMILEN = 64

CHARPROTSET = {"A": 1, "C": 2, "B": 3, "E": 4, "D": 5, "G": 6,
               "F": 7, "I": 8, "H": 9, "K": 10, "M": 11, "L": 12,
               "O": 13, "N": 14, "Q": 15, "P": 16, "S": 17, "R": 18,
               "U": 19, "T": 20, "W": 21, "V": 22, "Y": 23, "X": 24, "Z": 25}

CHARPROTLEN = 25


def load_map(path):
    d = {}
    fin = open(path)
    ic = 0
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("||")
        d[parts[0]] = parts[-1]
        if ic == 0:
            print("DB Map: ", parts[0], parts[-1])
        ic += 1
    fin.close()
    return d


def label_smiles(line, smi_ch_ind, MAX_SMI_LEN=100):
    X = np.zeros(MAX_SMI_LEN, dtype=np.int64())
    for i, ch in enumerate(line[:MAX_SMI_LEN]):
        X[i] = smi_ch_ind[ch]
    return X


def label_sequence(line, smi_ch_ind, MAX_SEQ_LEN=1000):
    X = np.zeros(MAX_SEQ_LEN, np.int64())
    for i, ch in enumerate(line[:MAX_SEQ_LEN]):
        X[i] = smi_ch_ind[ch]
    return X


DRUGNAME_2_SMILE = load_map("%s/../wdir/DrugBankNameSMILE_Filtered.txt" % C_DIR)
UNIPROT_2_SEQUENCE = load_map("%s/../wdir/ProteinTarget%s.txt" % (C_DIR, conf.GENE_SUFFIX))


class CustomDataSet(Dataset):
    def __init__(self, pairs):
        self.pairs = pairs

    def __getitem__(self, item):
        return self.pairs[item]

    def __len__(self):
        return len(self.pairs)


def line_count(filename):
    return int(subprocess.check_output(['wc', '-l', filename]).split()[0])


class CustomFileDataSet(Dataset):
    def __init__(self, path, bufferSize=10000):
        self.path = path
        self.ln = line_count(path)
        self.buffer = []
        self.currentBufferId = -1
        self.bufferSize = bufferSize
        self.fin = open(self.path)

    def __getNextBuffer__(self):
        self.buffer.clear()
        self.currentBufferId += 1
        for i in range(self.bufferSize):
            line = self.fin.readline()
            if line == "":
                break
            self.buffer.append(line)
    def __getitem__(self, item):
        if self.currentBufferId == -1:
            self.__getNextBuffer__()
        bufferId = item // self.bufferSize
        offset = item % self.bufferSize
        if bufferId == self.currentBufferId:
            return self.buffer[offset]
        elif bufferId == self.currentBufferId + 1:
            self.__getNextBuffer__()
            return self.buffer[offset]
        else:
            print("Fatal error: CustomFileDataSet needed to be access in order.")
            exit(-1)

    def __len__(self):
        return self.ln


def collate_fn(batch_data):
    N = len(batch_data)
    drug_ids, protein_ids = [], []
    compound_max = 100
    protein_max = 1000
    compound_new = torch.zeros((N, compound_max), dtype=torch.long)
    protein_new = torch.zeros((N, protein_max), dtype=torch.long)
    labels_new = torch.zeros(N, dtype=torch.long)
    for i, pair in enumerate(batch_data):
        # _, _, compoundstr, proteinstr, label = pair.strip().split()
        pair = pair.strip().split()
        drug_id, protein_id, compoundstr, proteinstr, label = pair[-5], pair[-4], pair[-3], pair[-2], pair[-1]
        drug_ids.append(drug_id)
        protein_ids.append(protein_id)
        compoundint = torch.from_numpy(label_smiles(compoundstr, CHARISOSMISET, compound_max))
        compound_new[i] = compoundint
        proteinint = torch.from_numpy(label_sequence(proteinstr, CHARPROTSET, protein_max))
        protein_new[i] = proteinint
        label = float(label)
        labels_new[i] = int(label)
    # return (drug_ids,protein_ids, compound_new, protein_new, labels_new)
    return (compound_new, protein_new, labels_new)


def collate_fn2(batch_data):
    N = len(batch_data)
    compound_max = 100
    protein_max = 1000
    compound_new = torch.zeros((N, compound_max), dtype=torch.long)
    protein_new = torch.zeros((N, protein_max), dtype=torch.long)
    for i, pair in enumerate(batch_data):
        pair = pair.strip().split("||")
        drug_id, protein_id = pair[0], pair[1]
        compoundstr = DRUGNAME_2_SMILE[drug_id]
        proteinstr = UNIPROT_2_SEQUENCE[protein_id]
        compoundint = torch.from_numpy(label_smiles(compoundstr, CHARISOSMISET, compound_max))
        compound_new[i] = compoundint
        proteinint = torch.from_numpy(label_sequence(proteinstr, CHARPROTSET, protein_max))
        protein_new[i] = proteinint

    return (compound_new, protein_new)
