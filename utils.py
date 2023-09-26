import os
import csv
import joblib
import numpy as np
from io import open
import time

def write_results(res, file_out):
    with open(file_out, "w") as f:
        writer = csv.writer(f)
        writer.writerows(res)


def read_all_lines(path):
    fin = open(path)
    lines = [line.strip() for line in fin.readlines()]
    fin.close()
    return lines

def get_key(val, d):
    for k, v in d.items():
        if val == v:
            return k
        return "Key not found!"


def get_insert_dict(d, k, v):
    try:
        v = d[k]
    except:
        d[k] = v

    return v

def expand_dict_by_complement_values(dc1, dc2, val):
    for key in dc2.keys():
        if key not in dc1.keys():
            dc1.update({key: val})
    return dc1

def write_list_to_file(lst, file_out):
    with open(file_out, "w") as f:
        for elem in lst:

            f.write(elem + "\n")

def getCurrentTimeString():
    t = time.localtime()
    currentTime = time.strftime("%Y-%m-%d %H:%M:%S", t)
    return currentTime

def write_dict_to_file(dc, file_out):
    with open(file_out, "w") as f:
        writer = csv.writer(f)
        for k, v in dc.items():
            writer.writerow([k,v])


def write_nparray_to_file(arr, rownames, colnames, file_out):
    if rownames == "" and colnames == "":
        np.savetxt(file_out, arr, delimiter="\t")
    else:
        import pandas as pd
        arr = pd.DataFrame(arr, index=rownames, columns=colnames)
        arr.to_csv(file_out, index=True, header=True, sep="\t")


def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def convertHexToBinString888(hexString):
    # scale = 16  ## equals to hexadecimal
    # num_of_bits = 888
    return bin(int(hexString, 16))[2:].zfill(888)


def convertBinString888ToArray(binString888):
    ar = np.ndarray(888, dtype=float)
    ar.fill(0)
    for i in range(887, -1, -1):
        if binString888[i] == "1":
            ar[i] = 1
    return ar


def convertHex888ToArray(hex888):
    return convertBinString888ToArray(convertHexToBinString888(hex888))


def get_dict(d, k, v=-1):
    try:
        v = d[k]
    except:
        pass
    return v


def get_insert_key_dict(d, k, v=0):
    try:
        v = d[k]
    except:
        d[k] = v
    return v


def add_dict_counter(d, k, v=1):
    try:
        v0 = d[k]
    except:
        v0 = 0
    d[k] = v0 + v


def sort_dict(dd, onlykey=False, reverse=True):
    # Reverse=True: Descending

    kvs = []
    for key, value in sorted(dd.items(), key=lambda p: (p[1], p[0])):
        if onlykey:
            kvs.append(key)
        else:
            kvs.append([key, value])
    if reverse:
        kvs = kvs[::-1]
    return kvs

# def sort_dict_by_value(dd, reverse = True):
#     kvs = sorted(dd.items(), key=lambda x:x[1], reverse=reverse)
#     return kvs
#
#

def sum_sort_dict_counter(dd):
    cc = 0
    for p in dd:
        cc += p[1]
    return cc


def get_update_dict_index(d, k):
    try:
        current_index = d[k]
    except:
        current_index = len(d)
        d[k] = current_index
    return current_index


def get_dict_index_only(d, k):
    try:
        current_index = d[k]
    except:
        current_index = -1

    return current_index


def load_list_from_file(path):
    list = []
    fin = open(path)
    while True:
        line = fin.readline()
        if line == "":
            break
        list.append(line.strip())
    fin.close()
    return list


def reverse_dict(d):
    d2 = dict()
    for k, v in d.items():
        d2[v] = k
    return d2


def save_obj(obj, path):
    joblib.dump(obj, path)


def load_obj(path):
    return joblib.load(path)


def loadMapFromFile(path, sep="\t", keyPos=0, valuePos=1):
    fin = open(path)
    d = dict()
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split(sep)
        d[parts[keyPos]] = parts[valuePos]
    fin.close()
    return d


def loadMapSetFromFile(path, sep="\t", keyPos=0, valuePos=1, sepValue="", isStop=""):
    fin = open(path)
    dTrain = dict()

    if isStop != "":
        dTest = dict()

    d = dTrain

    while True:
        line = fin.readline()
        if line == "":
            break
        if isStop != "":
            if line.startswith(isStop):
                d = dTest
                continue
        parts = line.strip().split(sep)
        v = get_insert_key_dict(d, parts[keyPos], set())
        if sepValue == "":
            v.add(parts[valuePos])
        else:
            values = parts[valuePos]
            values = values.split(sepValue)
            for value in values:
                v.add(value)
    fin.close()
    if isStop != "":
        return dTrain, dTest
    return dTrain

def getTanimotoScore(ar1, ar2):
    c1 = np.sum(ar1)
    c2 = np.sum(ar2)
    bm = np.dot(ar1, ar2)
    return bm * 1.0 / (c1 + c2 - bm + 1e-10)

def getCosin(ar1, ar2):

    return np.dot(ar1, ar2)

def get3WJaccardOnSets(set1, set2):
    len1 = len(set1)
    len2 = len(set2)
    nMatch = 0
    for s in set1:
        if s in set2:
            nMatch += 1
    return 3.0 * nMatch / (len1 + len2 + nMatch + 0.01)