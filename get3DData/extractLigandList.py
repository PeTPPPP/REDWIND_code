import glob

from params import C_DIR
import os.path
import urllib.parse
from selenium import webdriver
import time
import utils
from bs4 import BeautifulSoup
import params
from selenium.webdriver.common.by import By
import requests
import shutil
import autodock_runner.atd_conf as conf
Ligand2CID_PATH = "%s/out/LigandCID.map" % C_DIR
Ligand_SDF3D_DIR = "/home/gpux1/Downloads/Ligand3"
def extract_restricted_ligand():
    ligands = set()
    fin = open("%s/data/ligand_dup.txt" % C_DIR)
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split(",")
        ligands.add(parts[0])
    fout = open("%s/data/restricted_ligands.txt" % C_DIR, "w")
    for l in sorted(list(ligands)):
        fout.write("%s\n" % l)
    fout.close()


def getCIDIDURL():
    browser = webdriver.Chrome()
    try:
        dLigand2CidURL = utils.load_obj(Ligand2CID_PATH)
    except:
        dLigand2CidURL = dict()
    ligandNames = utils.load_list_from_file("%s/data/restricted_ligandsx.txt" % C_DIR)
    for ligand in ligandNames:
        if ligand in dLigand2CidURL:
            print("Skip: ", ligand)
            continue
        pubchemURL = "https://pubchem.ncbi.nlm.nih.gov/#query="+urllib.parse.quote(ligand)
        print(pubchemURL)
        # exit(-1)
        browser.get(pubchemURL)
        time.sleep(5)
        try:
            featureDiv = browser.find_element(By.ID, "featured-results")
            url = featureDiv.find_element(By.XPATH, "//a[starts-with(@href,'https://pubchem.ncbi.nlm.nih.gov/compound/')]")
            url = url.get_attribute("href")
            dLigand2CidURL[ligand] = url
        except:
            print("Not found CID URL of: ", ligand)
            continue

        print(ligand, url)
        if len(dLigand2CidURL) % 10 == 0:
            utils.save_obj(dLigand2CidURL, Ligand2CID_PATH)

    utils.save_obj(dLigand2CidURL, Ligand2CID_PATH)
def download3DFromURL(url, targetPath):
    response = requests.get(url)
    fout = open(targetPath, "w")
    # print(response.text)
    fout.write(response.text)
    fout.close()
    # exit(-1)
def downloadPubChem3D(outDir = Ligand_SDF3D_DIR):
    utils.ensure_dir(outDir)
    dLigand2CidURL = utils.load_obj(Ligand2CID_PATH)
    for ligandName, cidURL in dLigand2CidURL.items():
        cid = cidURL.split("/")[-1]
        url3DS = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/%s/record/SDF/?record_type=3d&response_type=text&response_basename=Conformer3D_CID_%s" % (cid, cid)
        targetPath = "%s/%s.sdf" % (outDir, ligandName.replace(" ", "-"))
        if os.path.exists(targetPath):
            print("Skip: ", targetPath)
            continue
        try:
            print(url3DS)
            download3DFromURL(url3DS, targetPath)

            print("Downloaded ", ligandName, cid)
            time.sleep(4)
        except:
            print("Error: Unable to download: ", ligandName, cid)


def prepare_ligand_pdbqt_bash(inpDir = Ligand_SDF3D_DIR, bash_path = "/home/gpux1/Downloads/bash_gene_ligands.sh"):
    sdf_paths = glob.glob("%s/*.sdf" % Ligand_SDF3D_DIR)
    pdb_dir = "%s_PDBQT" % inpDir
    utils.ensure_dir(pdb_dir)
    f_bash = open(bash_path, "w")
    f_bash.write("#!/bin/bash\n")
    f_bash.write("source %s\n" % conf.CONDA_SH)
    f_bash.write("conda activate %s\n" % conf.ENV_NAME)
    f_bash.write("echo $SHELL\n")
    f_bash.write("conda info\n")
    f_bash.write("%s\n" % conf.EXPORT_BIN_PATH_1)

    for sdf_path in sdf_paths:
        ligandName = sdf_path.split("/")[-1].split(".")[0]
        targetPDB = "%s/%s.pdb" % (pdb_dir, ligandName)
        f_bash.write("obabel -isdf \"%s\" -o pdb >> \"%s\"\n" % (sdf_path, targetPDB))
        f_bash.write("prepare_ligand -l \"%s\" -o \"%sqt\"\n" % (targetPDB, targetPDB))
        f_bash.write("if [ ! -f \"%sqt\" ]\n" % targetPDB)
        f_bash.write("then\n")
        f_bash.write("\techo \"Remove: %s\"\n" % targetPDB)
        f_bash.write("\trm \"%s\"\n"%targetPDB)
        f_bash.write("fi\n")
    f_bash.close()
    cmd = "chmod +x %s" % bash_path
    os.system(cmd)
    os.system(bash_path)


def cp_ligand_pdbqt():
    import shutil
    ligand_list = utils.load_list_from_file("%s/selectedDrugList" % params.DATA_DIR)
    ligand_list = sorted(list(set(ligand_list)))
    targetDir = "/home/gpux1/Downloads/LigandXAll_PDBQT"
    utils.ensure_dir(targetDir)
    srcDir = "/home/gpux1/Downloads/3DLigands_PDBQT"
    for ligand in ligand_list:
        ligand = ligand.replace(" ", "-")
        srcPath = "%s/%s.pdbqt" % (srcDir, ligand)
        if os.path.exists(srcPath):
            shutil.copy(srcPath, targetDir)
        else:
            print("mISSING", srcPath)


if __name__ == "__main__":
    # extract_restricted_ligand()
    # getCIDIDURL()
    # downloadPubChem3D()
    # prepare_ligand_pdbqt_bash()
    cp_ligand_pdbqt()