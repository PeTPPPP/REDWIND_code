import random
import os

import params
from .model import AttentionDTI
from .dataset import CustomDataSet, CustomFileDataSet, collate_fn2
from torch.utils.data import DataLoader
from prefetch_generator import BackgroundGenerator
from tqdm import tqdm
from hyperparameter import hyperparameter
from pytorchtools import EarlyStopping
import timeit
from tensorboardX import SummaryWriter
import numpy as np
import torch
import torch.nn as nn2
import torch.optim as optim
import torch.nn.functional as F
import os
from . import conf

C_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_PATH = "%s/KIBA/%s/valid_best_checkpoint.pth" % (C_DIR, conf.IFOLD)


def loadModel():
    SEED = 1234
    random.seed(SEED)
    torch.manual_seed(SEED)
    torch.cuda.manual_seed_all(SEED)
    # torch.backends.cudnn.deterministic = True

    """init hyperparameters"""
    hp = hyperparameter()

    """Load preprocessed data."""

    device = torch.device(params.DEVICE)

    model = AttentionDTI(hp).cuda(device=device)
    model.load_state_dict(torch.load(MODEL_PATH))
    return model, hp
    # test_dataset = CustomDataSet(test_dataset)
    # test_dataset_load = DataLoader(test_dataset, batch_size=hp.Batch_size, shuffle=False, num_workers=0,
    #                               collate_fn=collate_fn)


def loadTestData(path):
    with open(path, "r") as f:
        data_list = f.read().strip().split('\n')

        return data_list


def predict(model, pbar):
    model.eval()
    P, S = [], []
    with torch.no_grad():
        for i, data in pbar:
            '''data preparation '''
            compounds, proteins = data
            compounds = compounds.cuda(params.DEVICE)
            proteins = proteins.cuda(params.DEVICE)

            predicted_scores = model(compounds, proteins)
            predicted_scores = F.softmax(predicted_scores, 1).to('cpu').data.numpy()
            predicted_labels = np.argmax(predicted_scores, axis=1)
            predicted_scores = predicted_scores[:, 1]

            P.extend(predicted_labels)
            S.extend(predicted_scores)
    return P, S


def getSegSuffix(p):
    ss = p.split("___")
    sg = ""
    offset = -1
    if len(ss) > 1:
        suffix = ss[-1]
        sg = "___%s" % suffix
        offset = int(suffix.split("_")[-1])
    return sg, offset


def run_test(testPath="./tmp/NN.txt"):
    model, hp = loadModel()
    # test_data = loadTestData("./data/Testing.txt")
    # test_data = loadTestData(test_path)
    # test_dataset = CustomDataSet(test_data)
    test_dataset = CustomFileDataSet(testPath)
    test_dataset_load = DataLoader(test_dataset, batch_size=70, shuffle=False, num_workers=0,
                                   collate_fn=collate_fn2)
    test_pbar = tqdm(
        enumerate(
            BackgroundGenerator(test_dataset_load)),
        total=len(test_dataset_load))
    P, S = predict(model, test_pbar)
    sg, _ = getSegSuffix(testPath)
    np.save("%s/tmp/TestPredictionLabel_%s%s%s.npy" % (C_DIR, conf.IFOLD, conf.GENE_SUFFIX, sg), np.asarray(P))
    np.save("%s/tmp/TestPredictionScore_%s%s%s.npy" % (C_DIR, conf.IFOLD, conf.GENE_SUFFIX, sg), np.asarray(S))

    return len(P)


def getPredictionIndices(testPath="", n=-1):
    print("getPredictionIndices...")
    sg, _ = getSegSuffix(testPath)

    if n == -1:
        ar = np.load("%s/tmp/TestPredictionLabel_%s%s%s.npy" % (C_DIR, conf.IFOLD, conf.GENE_SUFFIX, sg))
        n = len(ar)
        # ss = np.nonzero(ar)[0]
    indices = np.arange(0, n)
    np.save("%s/tmp/BindingIndices_%s%s%s.npy" % (C_DIR, conf.IFOLD, conf.GENE_SUFFIX, sg), indices)

    print("Finish getPredictionIndices")


def get_insert_dict(d, k, v):
    try:
        v = d[k]
    except:
        d[k] = v
    return v


def protein2Drug(testPath="", all=True):
    import dataset
    import joblib
    sg, proteinOffset = getSegSuffix(testPath)

    drug_names = sorted(list(dataset.DRUGNAME_2_SMILE.keys()))
    protein_ids = sorted(list(dataset.UNIPROT_2_SEQUENCE.keys()))

    N_DRUGS = len(drug_names)
    print("Loading prediction scores")
    scores = np.load("%s/tmp/TestPredictionScore_%s%s%s.npy" % (C_DIR, conf.IFOLD, conf.GENE_SUFFIX, sg))

    if len(scores) == 0:
        print("No Drug-Protein binding found on the given test file.")
        print("Exit.")
        exit(-1)
    d_protein2drug = {}
    if all:

        assert len(scores) % N_DRUGS == 0
        assert proteinOffset != -1
        nProtein = len(scores) // N_DRUGS
        # print("N protein: ", nProtein, proteinOffset, N_DRUGS, testPath)
        proteinNameList = []
        drugNameList = []
        scores = []
        for ii in range(nProtein):
            proteinId = ii + proteinOffset
            startId = ii * N_DRUGS
            endId = (ii + 1) * N_DRUGS
            # print(startId, endId, len(scores))
            xscores = scores[startId:endId]
            proteinName = protein_ids[proteinId]
            drugScores = get_insert_dict(d_protein2drug, proteinName, [])
            for p in zip(drug_names, xscores):
                drugScores.append(p)
            if len(drugScores) == 0:
                print("????")
                exit(-1)
        joblib.dump(d_protein2drug, "%s/tmp/protein2Drugs_%s%s%s.npy" % (C_DIR, conf.IFOLD, conf.GENE_SUFFIX, sg))

    else:
        indices = np.load("%s/tmp/BindingIndices_%s%s%s.npy" % (C_DIR, conf.IFOLD, conf.GENE_SUFFIX, sg))
        for ind in indices:
            ind = int(ind)
            drug_id = ind % N_DRUGS
            protein_id = ind // N_DRUGS
            sc = scores[ind]
            # if sc > 1:
            #     print(ind, sc)

            # assert sc < 1.01

            drugs = get_insert_dict(d_protein2drug, protein_ids[protein_id], [])
            drugs.append((drug_names[drug_id], sc))
        joblib.dump(d_protein2drug, "%s/tmp/protein2Drugs_%s%s%s.npyx" % (C_DIR, conf.IFOLD, conf.GENE_SUFFIX, sg))


def getSuffix(iSeg):
    if params.N_SEG == 1:
        return ""
    finfo = open("%s/tmp/NewTest.txt%s___.Info" % (C_DIR, conf.GENE_SUFFIX))
    segSize = int(finfo.readline())
    finfo.close()
    return "___%s_%s" % (iSeg, segSize * iSeg)


def getTestPath(iSeg=-1):
    if iSeg == -1:
        test_path = "%s/tmp/NewTest.txt%s" % (C_DIR, conf.GENE_SUFFIX)
    else:
        suffix = getSuffix(iSeg)
        test_path = "%s/tmp/NewTest.txt%s%s" % (C_DIR, conf.GENE_SUFFIX, suffix)
    return test_path


def run_predict():
    if not os.path.exists("%s/KIBA/%s/valid_best_checkpoint.pth" % (C_DIR, conf.IFOLD)):
        from HyperAttentionDTI_main import run
        run()
    testPath = getTestPath(params.I_SEG)
    nSample = -1
    print("Testing path: ", testPath)
    nSample = run_test(testPath=testPath)
    if not params.FULL_PREDICTION:
        getPredictionIndices(testPath=testPath, n=nSample)
        protein2Drug(testPath=testPath)


def merge():
    if params.FULL_PREDICTION:
        merge_full()
    else:
        merge_sub()


def merge_sub():
    import joblib
    d = {}
    if params.N_SEG <= 1:
        return

    for sid in range(params.N_SEG):
        sg = getSuffix(sid)
        path = "%s/tmp/protein2Drugs_%s%s%s.npy" % (C_DIR, conf.IFOLD, conf.GENE_SUFFIX, sg)
        print("Loading %s" % path)
        p = joblib.load(path)

        for k, v in p.items():
            d[k] = v
    joblib.dump(d, "%s/tmp/protein2Drugs_%s%s%s.npy" % (C_DIR, conf.IFOLD, conf.GENE_SUFFIX, ""))


def merge_full():
    import joblib
    import dataset
    # if params.N_SEG <= 0:
    #    return
    drug_names = sorted(list(dataset.DRUGNAME_2_SMILE.keys()))
    # drug_names = ["n-[(3z)-5-tert-butyl-2-phenyl-1,2-dihydro-3h-pyrazol-3-ylidene]-n'-(4-chlorophenyl)urea","purvalanol","mgb-bp-3", "(4r)-n-[4-({[2-(dimethylamino)ethyl]amino}carbonyl)-1,3-thiazol-2-yl]-4-methyl-1-oxo-2,3,4,9-tetrahydro-1h-beta-carboline-6-carboxamide"]

    protein_ids = sorted(list(dataset.UNIPROT_2_SEQUENCE.keys()))
    # protein_ids = protein_ids[:10]
    protein_names = [protein_ids[i] for i in range(len(protein_ids))]
    all_scores = []
    for sid in range(params.N_SEG):
        # testPath = getTestPath(params.I_SEG)

        sg = getSuffix(sid)
        path_score_i = "%s/tmp/TestPredictionScore_%s%s%s.npy" % (C_DIR, conf.IFOLD, conf.GENE_SUFFIX, sg)
        print("Loading scores at: %s" % path_score_i)
        scores_i = np.load(path_score_i)
        all_scores.append(scores_i)

    all_scores = np.concatenate(all_scores).reshape((-1, len(drug_names)))
    predictions = (drug_names, protein_names, all_scores)

    # joblib.dump(d, "%s/tmp/protein2Drugs_%s%s%s.npy" % (C_DIR, conf.IFOLD, conf.GENE_SUFFIX, ""))
    joblib.dump(predictions, "%s/tmp/protein2DrugsXF_%s%s%s.npy" % (C_DIR, conf.IFOLD, conf.GENE_SUFFIX, ""))


if __name__ == "__main__":
    # cmd = "head -n 100000 ./tmp/NewTest.txt >> ./tmp/NN.txt"
    # os.system(cmd)
    # if not os.path.exists("./KIBA/%s/valid_best_checkpoint.pth" % conf.IFOLD):
    #    from HyperAttentionDTI_main import run

    #    run()
    # run_test(test_path="./tmp/NewTest.txt%s" % conf.GENE_SUFFIX)
    # vef()
    # protein2Drug()
    print(C_DIR)
