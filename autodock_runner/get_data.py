import glob
import os
import random

import params
import utils
from autodock_runner import atd_conf as conf
import time

import re


def load_except_drug():
    except_drugs = set()
    fin = open(conf.EXCEPT_DRUGS)
    for l in fin.readlines():
        except_drugs.add(l.strip())
    fin.close()
    return except_drugs


def load_drugbank_map():
    p = "%s/../data/DrugBankNameX.txt" % conf.C_DIR
    fin = open(p)
    drugNameMap = dict()
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("||")
        dname = parts[0].lower()
        drugbankId = parts[2]
        syns = parts[-2].lower().split(",")
        salts = parts[-1].lower().split(",")
        drugNameMap[dname] = drugbankId
        for syn in syns:
            if len(syn) > 1:
                drugNameMap[syn] = drugbankId
        for salt in salts:
            if len(salt) > 1:
                drugNameMap[salt] = drugbankId

    fin.close()
    return drugNameMap


def load_dti_prediction():
    fin = open(conf.DTI_PREDICTED_PATH)
    except_drugs = load_except_drug()
    drugbankNameMap = load_drugbank_map()
    genList = []
    predictedDrugs = []
    drugSet = set()
    drugRemap = dict()
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("\t")
        assert len(parts) > 2
        geneName = parts[0]
        genList.append(geneName)
        drugList = []
        for i in range((len(parts) - 1) // 2):
            dname = parts[2 * i + 1]
            if dname.startswith("DB"):
                if dname in except_drugs:
                    continue
                drugList.append(dname)
                drugSet.add(dname)
                drugRemap[dname] = dname
            else:
                dbid = drugbankNameMap[dname.lower()]
                if dbid in except_drugs:
                    continue
                drugList.append(dbid)
                drugSet.add(dbid)
                drugRemap[dbid] = dname.lower()
        predictedDrugs.append(drugList)

    return genList, predictedDrugs, drugSet, drugRemap


def get_drug_pdb(drugSet=None, drugRemap=None):
    if drugSet is None:
        _, _, drugSet, drugRemap = load_dti_prediction()
    import requests
    headers = {
        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/39.0.2171.95 Safari/537.36'}

    DB_PDB_PATTERN = "https://go.drugbank.com/structures/small_molecule_drugs/%s.pdb"
    dbIds = sorted(list(drugSet))
    for dbId in dbIds:
        URL = DB_PDB_PATTERN % dbId
        print(URL)

        contents = requests.get(URL, headers=headers).text
        spath = "%s/%s.pdb" % (conf.ATD_DIR, drugRemap[dbId])
        print(spath)
        fout = open(spath, "w")
        fout.write("%s" % contents)
        fout.close()
        time.sleep(3)
        # break


def gen_drug_pdbqt(drugSet=None, drugRemap={}):
    if drugSet is None:
        _, _, drugSet, drugRemap = load_dti_prediction()
    bash_path = "%s/gen_drug_pdbqt.sh" % conf.ATD_DIR
    f_bash = open(bash_path, "w")
    f_bash.write("#!/bin/bash\n")
    f_bash.write("source %s\n" % conf.CONDA_SH)
    f_bash.write("conda activate %s\n" % conf.ENV_NAME)
    f_bash.write("echo $SHELL\n")
    f_bash.write("conda info\n")
    f_bash.write("%s\n" % conf.EXPORT_BIN_PATH_1)
    for d in drugSet:
        dname = utils.get_dict(drugRemap, d, d)
        cmd = "prepare_ligand -l \"%s/%s.pdb\" -o \"%s/%s.pdbqt\"" % (
            conf.ATD_DIR, dname.replace(" ", "_"), conf.ATD_DIR, dname.replace(" ", "_"))
        f_bash.write("%s\n" % cmd)
    f_bash.close()

    cmd = "chmod +x %s" % bash_path
    os.system(cmd)
    os.system("%s" % bash_path)


def get_protein_gene_pdbqt(gene_name):
    pattern = "%s/%s*.pdbqt" % (conf.ATD_DIR, gene_name)

    paths = glob.glob(pattern)
    # print(gene_name,paths)
    assert len(paths) == 1
    return paths[0]


def get_gene_grid_conf_file(gene_name):
    path = "%s/conf_%s.txt" % (conf.ATD_DIR, gene_name)
    if os.path.exists(path):
        return path
    else:
        return "%s/conf_%s.txt" % (conf.ATD_DIR, "default")


def gen_receptor_pdbqt_sh(genList=None):
    if genList is None:
        genList, predictedDrugs, drugSet, drugRemap = load_dti_prediction()
    bash_path = "%s/gen_receptor_pdbqt.sh" % conf.ATD_DIR
    f_bash = open(bash_path, "w")
    f_bash.write("#!/bin/bash\n")
    f_bash.write("source %s\n" % conf.CONDA_SH)
    f_bash.write("conda activate %s\n" % conf.ENV_NAME)
    f_bash.write("echo $SHELL\n")
    f_bash.write("conda info\n")
    f_bash.write("%s\n" % conf.EXPORT_BIN_PATH_1)

    for gene in genList:
        pdb_file_pattern = "%s/%s*.pdb" % (conf.ATD_DIR, gene)
        paths = glob.glob(pdb_file_pattern)
        # print(pdb_file_pattern, paths)
        assert len(paths) == 1
        path = paths[0]
        cmd = "prepare_receptor  -A \"hydrogens\" -r \"%s\" -o \"%s\"" % (path, path[:-4] + ".pdbqt")
        f_bash.write("%s\n" % cmd)
    f_bash.close()
    cmd = "chmod +x %s" % bash_path
    os.system(cmd)
    os.system(bash_path)


def prepare_restricted_protein_pdbqt_bash3(RESTRICTED_PROTEIN_DIR="/home/gpux1/Downloads/RestrictedGeneProteins",
                                           W_DIR="/home/gpux1/Downloads/",
                                           bash_path="/home/gpux1/Downloads/bash_gen_proteins2.sh"):
    pdbqt_dir = "%s/GeneProH_PDBQT" % W_DIR
    utils.ensure_dir(pdbqt_dir)
    f_bash = open(bash_path, "w")
    f_bash.write("#!/bin/bash\n")
    f_bash.write("source %s\n" % conf.CONDA_SH)
    f_bash.write("conda activate %s\n" % conf.ENV_NAME)
    f_bash.write("echo $SHELL\n")
    f_bash.write("conda info\n")
    f_bash.write("%s\n" % conf.EXPORT_BIN_PATH_1)

    genProList = utils.load_list_from_file(conf.CUSTOM_GENE_LIST_PATH)

    for genePro in genProList:
        pdb_file = "%s/%s.pdb" % (RESTRICTED_PROTEIN_DIR, genePro)
        assert os.path.exists(pdb_file), pdb_file
        targetPDBQT = "%s/%s.pdbqt" % (pdbqt_dir, genePro)
        cmd = "prepare_receptor -A \"hydrogens\" -r \"%s\" -o \"%s\"\n" % (pdb_file, targetPDBQT)
        f_bash.write("%s\n" % cmd)
    f_bash.close()
    cmd = "chmod +x %s" % bash_path
    os.system(cmd)
    os.system(bash_path)


def gen_autodock_vina_script(genList=None, predictedDrugs=None, drugSet=None, drugRemap={}):
    print(predictedDrugs)
    if genList is None:
        genList, predictedDrugs, drugSet, drugRemap = load_dti_prediction()
    bash_path = "%s/vina_scripts.sh" % conf.ATD_DIR
    f_bash = open(bash_path, "w")
    f_bash.write("#!/bin/bash\n")
    f_bash.write("source %s\n" % conf.CONDA_SH)
    f_bash.write("conda activate %s\n" % conf.ENV_NAME)
    f_bash.write("echo $SHELL\n")
    f_bash.write("conda info\n")
    f_bash.write("%s\n" % conf.EXPORT_BIN_PATH_1)
    for i, geneName in enumerate(genList):
        gene_grid_path = get_gene_grid_conf_file(geneName)
        if gene_grid_path.__contains__("default"):
            receptor_path = get_protein_gene_pdbqt(geneName)
            pref = "vina --config \"%s\" --exhaustiveness 32 --receptor \"%s\"" % (gene_grid_path, receptor_path)
        else:
            pref = "vina --config \"%s\" --exhaustiveness 32" % gene_grid_path
        # if len(predictedDrugs) > 1:
        #     drugs = predictedDrugs[i]
        # else:
        drugs = predictedDrugs
        print(drugs)
        for drugId in drugs:
            drugName = utils.get_dict(drugRemap, drugId, drugId)
            cmd = "%s --ligand \"%s\" --out binding_%s_%s.pdbqt >> logs__%s_%s" % (
                pref, "%s/%s.pdbqt" % (conf.ATD_DIR, drugName), geneName, drugName, geneName, drugName)
            f_bash.write("%s\n" % cmd)
    f_bash.close()
    cmd = "chmod +x %s" % bash_path
    os.system(cmd)


def get_core_name(p):
    p = p.split("/")[-1].split(".")[0]
    return p


def gen_autodock_vina_script_input_folder(proteinFolder="/home/gpux1/Downloads/GeneProH_PDBQT",
                                          drugFolder="/home/gpux1/Downloads/LigandXAll_PDBQT", N_SEG=1):
    conf_default = "%s/conf_%s.txt" % (conf.ATD_DIR, "default")

    gene_paths = glob.glob("%s/*.pdbqt" % proteinFolder)
    drug_paths = glob.glob("%s/*.pdbqt" % drugFolder)
    SEG_SIZE = len(drug_paths) // N_SEG
    curSegID = -1
    f_bash = None
    bash_path = None
    drug_paths = sorted(drug_paths)
    gene_paths = sorted(gene_paths)
    for i, drug_path in enumerate(drug_paths):
        segId = i // SEG_SIZE
        if segId == curSegID + 1:
            bash_path = "%s/vina_scripts_batch_%s.sh" % (conf.ATD_DIR, segId)
            curSegID = segId
            f_bash = open(bash_path, "w")
            f_bash.write("#!/bin/bash\n")
            f_bash.write("source %s\n" % conf.CONDA_SH)
            f_bash.write("conda activate %s\n" % conf.ENV_NAME)
            f_bash.write("echo $SHELL\n")
            f_bash.write("conda info\n")
            f_bash.write("%s\n" % conf.EXPORT_BIN_PATH_1)
        drugName = get_core_name(drug_path)

        pref = "vina --config \"%s\" --exhaustiveness 32 --ligand \"%s\"" % (conf_default, drug_path)

        for j, gene_path in enumerate(gene_paths):
            rep_name = get_core_name(gene_path)
            cmd = "%s --receptor \"%s\" --out binding_%s_%s.pdbqt >> logs__%s_%s" % (
                pref, gene_path, drugName, rep_name, drugName, rep_name)
            f_bash.write("%s\n" % cmd)

        if i % SEG_SIZE == SEG_SIZE - 1 or i == len(drug_paths) - 1:
            f_bash.close()
            cmd = "chmod +x %s" % bash_path
            os.system(cmd)


def gen_autodock_vina_script_from_sample_file(sample_file="%s/GenePro-DrugRandomSample.txt" % params.DATA_DIR,
                                              proteinFolder="/home/gpux1/Downloads/GeneProX100HD_PDBQT",
                                              drugFolder="/home/gpux1/Downloads/LigandXAll_PDBQT", N_SEG=1):
    from get3DData.getProtein3D import loadGeneProDrugSample

    drug2GenePro = loadGeneProDrugSample(sample_file)
    allGeneProPDBQT = glob.glob("%s/*.pdbqt" % proteinFolder)
    allgenePros = list()
    for p in allGeneProPDBQT:
        allgenePros.append(p.split("/")[-1].split(".")[0])

    def sampleGenePro(targetList, excludeList, nSample):
        r2 = []
        # print(len(targetList), len(excludeList), nSample)
        for r in targetList:
            if r not in excludeList:
                r2.append(r)
        return random.sample(r2, k=nSample)

    conf_default = "%s/conf_%s.txt" % (conf.ATD_DIR, "default")
    bash_path = "%s/vina_scripts_batch_%s.sh" % (conf.ATD_DIR, 0)

    f_bash = open(bash_path, "w")
    f_bash.write("#!/bin/bash\n")
    f_bash.write("source %s\n" % conf.CONDA_SH)
    f_bash.write("conda activate %s\n" % conf.ENV_NAME)
    f_bash.write("echo $SHELL\n")
    f_bash.write("conda info\n")
    f_bash.write("%s\n" % conf.EXPORT_BIN_PATH_1)
    for drug, genePros in drug2GenePro.items():
        genProPathList = []
        nc = 0
        selectedGene = []
        for genePro in genePros:
            path = "%s/%s.pdbqt" % (proteinFolder, genePro)
            if not os.path.exists(path):
                nc += 1
            else:
                selectedGene.append(genePro)
        if nc > 0:
            ss = sampleGenePro(allgenePros, set(genePros), nc)

        else:
            ss = []
        selectedGene = selectedGene + ss

        drug_path = "%s/%s.pdbqt" % (drugFolder, drug)
        pref = "vina --config \"%s\" --exhaustiveness 32 --ligand \"%s\"" % (conf_default, drug_path)

        for j, genePro in enumerate(selectedGene):
            gene_path = "%s/%s.pdbqt" % (proteinFolder, genePro)

            rep_name = genePro
            cmd = "%s --receptor \"%s\" --out binding_%s_%s.pdbqt >> logs__%s_%s" % (
                pref, gene_path, drug, rep_name, drug, rep_name)
            f_bash.write("%s\n" % cmd)

    f_bash.close()
    cmd = "chmod +x %s" % bash_path
    os.system(cmd)


def get_aff(log_path):
    token = "-----+"
    fin = open(log_path)
    while True:
        line = fin.readline()

        if line == "":
            break
        if not line.startswith(token):
            continue
        else:
            break

    assert line.startswith(token)
    line = fin.readline().strip()
    line = re.sub("\s\s+", " ", line)
    aff = line.split(" ")[1]
    return aff


def extract_aff():
    img_paths = glob.glob("%s/*.png" % conf.IMG_FOLDER)
    IMG_AFF_FOLDER = "%s/IMG_AFF" % conf.IMG_FOLDER
    utils.ensure_dir(IMG_AFF_FOLDER)
    for img_path in img_paths:
        img_name = os.path.basename(img_path)
        binding_name = img_name.split(".")[0]
        log_path = "%s/logs__%s" % (conf.ATD_DIR, binding_name)
        aff = get_aff(log_path)
        cmd = "cp \"%s\" \"%s\"" % (img_path, "%s/%s_%s.png" % (IMG_AFF_FOLDER, img_name[:-4], aff.replace(".", "_")))
        print("CMD: ", cmd)
        os.system(cmd)


def append_aff_binding():
    import re
    def get_aff(log_path):
        token = "-----+"
        fin = open(log_path)
        while True:
            line = fin.readline()

            if line == "":
                break
            if not line.startswith(token):
                continue
            else:
                break

        assert line.startswith(token)
        line = fin.readline().strip()
        line = re.sub("\s\s+", " ", line)
        aff = line.split(" ")[1]
        return aff

    img_paths = glob.glob("%s/*.png" % conf.IMG_FOLDER)
    IMG_AFF_FOLDER = "%s/IMG_AFF" % conf.IMG_FOLDER
    utils.ensure_dir(IMG_AFF_FOLDER)
    for img_path in img_paths:
        img_name = os.path.basename(img_path)
        binding_name = img_name.split(".")[0]
        log_path = "%s/logs__%s" % (conf.ATD_DIR, binding_name)
        aff = get_aff(log_path)
        cmd = "cp \"%s\" \"%s\"" % (img_path, "%s/%s_%s.png" % (IMG_AFF_FOLDER, img_name[:-4], aff.replace(".", "_")))
        print("CMD: ", cmd)
        os.system(cmd)


def extract_all_logs():
    all_logs = glob.glob("%s/logs__*" % conf.ATD_DIR)
    drug2Gene = dict()
    for log_path in all_logs:
        log_name = log_path.split("/")[-1].split("__")[-1]
        ss = log_name.split("_")
        drugName = ss[0]
        geneName = "_".join(ss[1:])
        try:
            aff = get_aff(log_path)
        except:
            continue
        dGeneScore = utils.get_insert_key_dict(drug2Gene, drugName, {})
        dGeneScore[geneName] = aff

    drugList = sorted(list(drug2Gene.keys()))
    geneSet = drug2Gene[drugList[0]]

    priorPref = ["RA", "RX"]
    priorList = []
    suffList = []
    for gene in geneSet:
        isInPrior = False
        for prior in priorPref:
            if gene.startswith(prior):
                priorList.append(gene)
                isInPrior = True
                break
        if not isInPrior:
            suffList.append(gene)
    priorList = sorted(priorList)
    suffList = sorted(suffList)
    geneList = priorList + suffList

    fout = open("%s/../out/binding_aff_all5.csv" % conf.C_DIR, "w")
    for gene in geneList:
        fout.write(",%s" % gene)
    fout.write("\n")
    for drug in drugList:
        affs = drug2Gene[drug]
        fout.write(drug)
        for i, gene in enumerate(geneList):
            try:
                aff = affs[gene]
            except:
                continue
            fout.write(",%s" % aff)
        fout.write("\n")
    fout.close()


def extract_all_logs2():
    all_logs = glob.glob("%s/logs__*" % conf.ATD_DIR)
    drug2Gene = dict()
    for log_path in all_logs:
        log_name = log_path.split("/")[-1].split("__")[-1]
        ss = log_name.split("_")
        drugName = ss[0]
        geneName = "_".join(ss[1:])
        aff = get_aff(log_path)
        dGeneScore = utils.get_insert_key_dict(drug2Gene, drugName, {})
        dGeneScore[geneName] = aff

    fout = open("%s/sampleDrugGeneBinding_new.txt" % params.W_DIR, "w")
    for drugName, geneScores in drug2Gene.items():
        fout.write("%s\n" % drugName)
        gss = utils.sort_dict(geneScores)
        for g, s in gss:
            fout.write("%s,%s\n" % (g, s))
        fout.write("\n")
    fout.close()


if __name__ == "__main__":
    # prepare_restricted_protein_pdbqt_bash3()
    # gen_receptor_pdbqt_sh(genList=utils.load_list_from_file(conf.CUSTOM_GENE_LIST_PATH))
    # get_drug_pdb()
    # gen_drug_pdbqt(drugSet=utils.load_list_from_file(conf.CUSTOM_DRUG_LIST_PATH))
    # gen_drug_pdbqt(utils.load_list_from_file(conf.CUSTOM_DRUG_LIST_PATH))
    # gen_autodock_vina_script()
    # gen_autodock_vina_script(genList=utils.load_list_from_file(conf.CUSTOM_GENE_LIST_PATH), predictedDrugs=utils.load_list_from_file(conf.CUSTOM_DRUG_LIST_PATH))
    # extract_aff()
    # gen_autodock_vina_script_input_folder()
    # extract_all_logs()
    gen_autodock_vina_script_from_sample_file()
    # extract_all_logs2()
    pass
