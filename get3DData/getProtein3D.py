import os.path
import random

from selenium import webdriver
import time
import utils
from bs4 import BeautifulSoup
import params
from selenium.webdriver.common.by import By
import requests
import shutil
import glob
from autodock_runner import atd_conf as conf
RAW_REPS = "%s/out/RawProteinURLReps.dat" % params.C_DIR
MISSING_SET = "%s/out/MissingAlphaFoldURL.dat" % params.C_DIR

Protein_REMAP = "%s/out/ProteinRemap.dat" % params.C_DIR
RESTRICTED_PROTEIN_DIR = "/home/gpux1/Downloads/GeneH"
def loadDGeneProtein(path="%s/out/SelectedGene2Protein.txt" % params.C_DIR):
    d = dict()
    fin = open(path)
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("||")
        geneName = parts[0]
        protein1 = parts[1].split(",")[0]
        d[geneName] = protein1
    fin.close()
    return d


def getProteinUniprotURL(proteinID):
    return "https://www.uniprot.org/uniprotkb/%s/entry" % proteinID


def getAlphaFoldURL(remappedProteinID):
    return "https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_v4.pdb" % remappedProteinID


def downloadProteinURL():
    r"""
    Getting raw responses from drugs.com for mapping from drug names to drug ids
    The result is saved in RAW_DRUG_TEXT
    """
    browser = webdriver.Chrome()

    gen2Protein = loadDGeneProtein()
    try:
        dProtein2URLP = utils.load_obj(RAW_REPS)
    except:
        dProtein2URLP = dict()

    try:
        missingIds = utils.load_obj(MISSING_SET)
    except:
        missingIds = set()
    proteinIds = sorted(list(gen2Protein.values()))
    for proteinId in proteinIds:
        if proteinId in dProtein2URLP:
            continue
        if proteinId in missingIds:
            continue
        urlx = getProteinUniprotURL(proteinId)
        browser.get(urlx)
        time.sleep(6)
        try:
            html = browser.find_element(By.XPATH, "//a[starts-with(@href,'https://alphafold.ebi.ac.uk/')]")
            url = html.get_attribute('href')
            remapProteinID = url.split("/")[-1]
            alphaFoldURL = getAlphaFoldURL(remapProteinID)
            dProtein2URLP[proteinId] = alphaFoldURL
            print(alphaFoldURL)
        except Exception as e:
            print(e)
            missingIds.add(proteinId)
            continue

        if len(dProtein2URLP) % 10 == 0:
            utils.save_obj(dProtein2URLP, RAW_REPS)
            utils.save_obj(missingIds, MISSING_SET)
        # exit(-1)
    utils.save_obj(dProtein2URLP, RAW_REPS)
    utils.save_obj(missingIds, MISSING_SET)

def download3DFromURL(url,targetPath):
    response = requests.get(url)
    fout = open(targetPath, "w")
    fout.write(response.text)
    fout.close()


def extractProteinFromAlphaFoldURL(url):
    parts = url.split("-")
    return parts[1]

def downloadAlphaFold3D(outDir = "/home/gpux1/Downloads/ProteinPDB"):
    utils.ensure_dir(outDir)
    dProtein2AlphaFoldURL = utils.load_obj(RAW_REPS)
    try:
        dProteinRemap = utils.load_obj(Protein_REMAP)
    except:
        dProteinRemap = dict()
    protein_rds = list(dProtein2AlphaFoldURL.keys())
    random.shuffle(protein_rds)
    for _ , proteinId in enumerate(protein_rds):
        url = dProtein2AlphaFoldURL[proteinId]
        remappedProteinId = extractProteinFromAlphaFoldURL(url)
        targetPath = "%s/%s.pdb" % (outDir, remappedProteinId)
        if proteinId in dProteinRemap and os.path.exists(targetPath):
            print("Skip ", proteinId)
            continue
        targetPath = "%s/%s.pdb" % (outDir, remappedProteinId)
        try:
            download3DFromURL(url, targetPath)
            time.sleep(4)
        except Exception as e:
            print(e)
            continue
        print("Downloaded: ", url)

        dProteinRemap[proteinId] = remappedProteinId
        if len(dProteinRemap) % 10 == 0:
            utils.save_obj(dProteinRemap, Protein_REMAP)
    utils.save_obj(dProteinRemap, Protein_REMAP)

def loadRestrictedGene():
    p = "%s/data/ALL_gene_list_with_values_for_drug prediction.txt" % params.C_DIR
    restrictedGenes = set()
    f = open(p)
    while True:
        line = f.readline()
        if line == "":
            break

        parts = line.strip().split("\t")
        geneName = parts[0]
        if len(geneName) > 1:
            restrictedGenes.add(geneName)
    return sorted(list(restrictedGenes))



def loadRestrictedGene2():
    p = "%s/data/gen_list_full" % params.C_DIR
    restrictedGenes = set()
    f = open(p)
    while True:
        line = f.readline()
        if line == "":
            break

        parts = line.strip().split("\t")
        geneName = parts[0]
        if len(geneName) > 1:
            restrictedGenes.add(geneName)
    x = list(restrictedGenes)
    random.shuffle(x)
    return x

def downloadAlphaFold3DWithRestrictedGenes(outDir1 = "/home/gpux1/Downloads/ProteinPDB", outDir2 = RESTRICTED_PROTEIN_DIR):
    utils.ensure_dir(outDir1)
    utils.ensure_dir(outDir2)
    dProtein2AlphaFoldURL = utils.load_obj(RAW_REPS)
    # restrictedGenes = loadRestrictedGene()
    restrictedGenes = loadRestrictedGene2()
    gene2Protein = loadDGeneProtein()
    try:
        dProteinRemap = utils.load_obj(Protein_REMAP)
    except:
        dProteinRemap = dict()

    for geneName in restrictedGenes:
        proteinId = utils.get_dict(gene2Protein, geneName, -1)
        if proteinId == -1:
            print("Missing gene: ", geneName)
            continue
        alphaFoldURL = utils.get_dict(dProtein2AlphaFoldURL, proteinId, -1)
        if alphaFoldURL == -1:
            print("Missing alphaFoldURL for: ", proteinId, geneName)
            continue
        remappedProteinId = extractProteinFromAlphaFoldURL(alphaFoldURL)
        targetPath = "%s/%s.pdb" % (outDir1, remappedProteinId)
        if proteinId in dProteinRemap and os.path.exists(targetPath):
            cmd = "cp %s %s/%s_%s.pdb" % (targetPath, outDir2, geneName, remappedProteinId)
            os.system(cmd)
            continue
        targetPath = "%s/%s.pdb" % (outDir1, remappedProteinId)
        try:
            download3DFromURL(alphaFoldURL, targetPath)
            time.sleep(5)
            cmd = "cp %s %s/%s_%s.pdb" % (targetPath, outDir2, geneName, remappedProteinId)
            os.system(cmd)
        except Exception as e:
            print(e)
            continue
        print("Downloaded: ", alphaFoldURL)

        dProteinRemap[proteinId] = remappedProteinId
        if len(dProteinRemap) % 10 == 0:
            utils.save_obj(dProteinRemap, Protein_REMAP)
    utils.save_obj(dProteinRemap, Protein_REMAP)




def prepare_protein_pdbqt_bash(inpDir = RESTRICTED_PROTEIN_DIR, bash_path = "/home/gpux1/Downloads/bash_gen_proteins.sh"):
    pdb_paths = glob.glob("%s/*.pdb" % RESTRICTED_PROTEIN_DIR)
    pdbqt_dir = "%s_PDBQT" % inpDir
    utils.ensure_dir(pdbqt_dir)
    f_bash = open(bash_path, "w")
    f_bash.write("#!/bin/bash\n")
    f_bash.write("source %s\n" % conf.CONDA_SH)
    f_bash.write("conda activate %s\n" % conf.ENV_NAME)
    f_bash.write("echo $SHELL\n")
    f_bash.write("conda info\n")
    f_bash.write("%s\n" % conf.EXPORT_BIN_PATH_1)

    for pdb_path in pdb_paths:
        proteinName = pdb_path.split("/")[-1].split(".")[0]
        targetPDBQT = "%s/%s.pdbqt" % (pdbqt_dir, proteinName)
        # f_bash.write("cp \"%s\" \"%s\"\n"% (pdb_path, targetPDBQT))
        f_bash.write("prepare_receptor  -A \"hydrogens\" -r \"%s\" -o \"%s\"\n" % (pdb_path, targetPDBQT))
    f_bash.close()
    cmd = "chmod +x %s" % bash_path
    os.system(cmd)
    os.system(bash_path)


def getMapGene2ProteinPDBPath():
    gene2Protein = loadDGeneProtein()
    dProteinReMap = utils.load_obj(Protein_REMAP)

    dGenePro2Path = dict()
    PDB_DIR = "/home/gpux1/Downloads/ProteinPDB"
    for gene, protein in gene2Protein.items():
        try:
            genePro  = "%s_%s" % (gene, protein)
            proteinRemapped = dProteinReMap[protein]
            path = "%s/%s.pdb" % (PDB_DIR, proteinRemapped)
            if os.path.exists(path):
                dGenePro2Path[genePro] = path
        except:
            pass
    print(len(dGenePro2Path))
    return dGenePro2Path





def getSubList(l, indices):
    return [l[i] for i in indices]

def genRandomGeneTarget():
    dGenePro2Path = getMapGene2ProteinPDBPath()

    exceptGenePro = set(utils.load_list_from_file("%s/excludeGeneList" % params.DATA_DIR))
    geneProCandidates = list()
    for genePro in dGenePro2Path.keys():
        if genePro in exceptGenePro:
            continue
        if len(genePro) >= 13:
            continue
        geneProCandidates.append(genePro)
    print("GenePro Len: ", len(geneProCandidates))

    import numpy as np
    np.random.seed(1)
    drugList = utils.load_list_from_file("%s/selectedDrugList" % params.DATA_DIR)
    fout = open("%s/GenePro-DrugRandomSample.txt" % params.DATA_DIR, "w")
    ids = np.arange(0, len(geneProCandidates))

    rdIndices = []
    SZ = 100
    for i in range(len(drugList)):
        rdIndices.append(np.random.choice(ids,replace=False, size=SZ))
    for i, drugName in enumerate(drugList):
        genePros = getSubList(geneProCandidates, rdIndices[i])
        fout.write("%s,%s\n" % (drugName, ",".join(genePros)))
    fout.close()



def loadGeneProDrugSample(p = "%s/GenePro-DrugRandomSample.txt" % params.DATA_DIR):

    d = dict()

    with open(p) as f:
        lines = f.readlines()
        for line in lines:
            parts = line.strip().split(",")
            d[parts[0]] = parts[1:]
    return d

def prepare_protein_pdbqt_bash2(W_DIR = "/home/gpux1/Downloads/", bash_path = "/home/gpux1/Downloads/bash_gen_proteins2.sh"):
    dGenePro2Path = getMapGene2ProteinPDBPath()
    dGenePro2Drug = loadGeneProDrugSample()


    pdbqt_dir = "%s/GeneProX100HD_PDBQT" % W_DIR
    utils.ensure_dir(pdbqt_dir)
    f_bash = open(bash_path, "w")
    f_bash.write("#!/bin/bash\n")
    f_bash.write("source %s\n" % conf.CONDA_SH)
    f_bash.write("conda activate %s\n" % conf.ENV_NAME)
    f_bash.write("echo $SHELL\n")
    f_bash.write("conda info\n")
    f_bash.write("%s\n" % conf.EXPORT_BIN_PATH_1)
    genePros = set()

    for geneProsx in dGenePro2Drug.values():
        for v in geneProsx:
            genePros.add(v)
    for genePro in genePros:
        pdb_path = dGenePro2Path[genePro]
        # proteinName = pdb_path.split("/")[-1].split(".")[0]
        targetPDBQT = "%s/%s.pdbqt" % (pdbqt_dir, genePro)
        # f_bash.write("cp \"%s\" \"%s\"\n"% (pdb_path, targetPDBQT))
        f_bash.write("prepare_receptor -A \"hydrogens\" -r \"%s\" -o \"%s\"\n" % (pdb_path, targetPDBQT))
    f_bash.close()
    cmd = "chmod +x %s" % bash_path
    os.system(cmd)
    # os.system(bash_path)



def prepare_protein_pdbqt_bash3(W_DIR = "/home/gpux1/Downloads/", bash_path = "/home/gpux1/Downloads/bash_gen_proteins2.sh"):
    dGenePro2Path = getMapGene2ProteinPDBPath()



    pdbqt_dir = "%s/GeneProH_PDBQT" % W_DIR
    utils.ensure_dir(pdbqt_dir)
    f_bash = open(bash_path, "w")
    f_bash.write("#!/bin/bash\n")
    f_bash.write("source %s\n" % conf.CONDA_SH)
    f_bash.write("conda activate %s\n" % conf.ENV_NAME)
    f_bash.write("echo $SHELL\n")
    f_bash.write("conda info\n")
    f_bash.write("%s\n" % conf.EXPORT_BIN_PATH_1)

    genePros = set()
    genList=utils.load_list_from_file(conf.CUSTOM_GENE_LIST_PATH)

    print(len(genList))
    for geneProsx in genList:
        if geneProsx in dGenePro2Path:
            genePros.add(geneProsx)
        else:
            print(geneProsx)
    print(len(genePros))
    exit(-1)
    for genePro in genePros:
        pdb_path = dGenePro2Path[genePro]
        # proteinName = pdb_path.split("/")[-1].split(".")[0]
        targetPDBQT = "%s/%s.pdbqt" % (pdbqt_dir, genePro)
        # f_bash.write("cp \"%s\" \"%s\"\n"% (pdb_path, targetPDBQT))
        f_bash.write("prepare_receptor -A \"hydrogens\" -r \"%s\" -o \"%s\"\n" % (pdb_path, targetPDBQT))
    f_bash.close()
    cmd = "chmod +x %s" % bash_path
    os.system(cmd)
    os.system(bash_path)



if __name__ == "__main__":
    # downloadProteinURL()
    # downloadAlphaFold3D()
    # downloadAlphaFold3DWithRestrictedGenes()
    # prepare_protein_pdbqt_bash()
    # genRandomGeneTarget()
    prepare_protein_pdbqt_bash2()
    # prepare_protein_pdbqt_bash3()
