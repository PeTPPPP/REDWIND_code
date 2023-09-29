import params
import utils


def load_uniprot_gene2protein():
    fin = open("%s/UniprotReviewedGene2Proteins.txt" % params.W_DIR)
    d = {}
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("||")
        d[parts[0]] = parts[1].split(",")
    fin.close()
    return d


def load_candidate_gen_list(path):
    gene_lines = open(path).readlines()
    genes = [gene.strip() for gene in gene_lines]
    return genes


def loadProtein2Sequence():
    fin = open("%s/UniprotReviewedMapping2.txt" % params.W_DIR)
    d = {}
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("||")
        pp = parts[1].split(";")
        for p in pp:
            p = p.strip()
            if len(p) > 1:
                d[p] = parts[-1]
    return d


def match():

    d_geneprotein = load_uniprot_gene2protein()
    gene_list = load_candidate_gen_list(params.GENE_PATH)
    print("Total Input Gene: ", len(gene_list))
    missing_list = []
    matching_list = []
    proteinsets = set()
    for gene in gene_list:
        if gene not in d_geneprotein:
            missing_list.append(gene)
        else:
            matching_list.append(gene)
            for protein in d_geneprotein[gene]:
                proteinsets.add(protein)
    print("Missing: ", len(missing_list))
    print(missing_list)
    fout = open("%s/SelectedGene2Protein%s.txt" % (params.W_DIR, params.GEN_SUFFIX), "w")
    for gene in matching_list:
        fout.write("%s||%s\n" % (gene, ",".join(d_geneprotein[gene])))
    fout.close()

    # Exporting protein-sequence
    d = loadProtein2Sequence()
    fProteinTarget = open("%s/ProteinTarget%s.txt" % (params.W_DIR, params.GEN_SUFFIX), "w")
    for p in sorted(list(proteinsets)):
        fProteinTarget.write("%s||%s\n" % (p, d[p]))
    fProteinTarget.close()


def get_dict(d, k, v=-1):
    try:
        v = d[k]
    except:
        pass
    return v


def loaddrugname2id():
    fin = open("wdir/DrugBankNameID_SMILE.txt")
    dname2id = {}
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("||")
        dname2id[parts[0]] = parts[1]
    return dname2id

def export_gene_drug():
    if params.FULL_PREDICTION:
        export_gene_drug_full()
    else:
        export_gene_drug_sub()

def export_gene_drug_full():
    import joblib
    print("Loading full scores...")
    drugNames, proteinNames, scores = joblib.load("./HyperAttentionDTI/tmp/protein2DrugsXF_0%s.npy" % params.GEN_SUFFIX)
    print("Score shape: ", scores.shape)
    dProtein2Id = dict()
    for i, name in enumerate(proteinNames):
        dProtein2Id[name] = i
    def getScores(proteinName):
        proteinid = utils.get_dict(dProtein2Id, proteinName, -1)
        if proteinid != -1:
            return scores[proteinid]
        else:
            return None

    ddrug2id = loaddrugname2id()
    fin = open("%s/SelectedGene2Protein%s.txt" % (params.W_DIR, params.GEN_SUFFIX))
    fout = open("%s/Gene2Drugs%s.txt" % (params.W_DIR, params.GEN_SUFFIX), "w")
    drugId2Chem = {}
    cc = 0
    n = 0
    print("Exporting...")
    ic = 0
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("||")
        protins = parts[1].split(",")
        ds = dict()
        for p in protins:
            iscores = getScores(p)
            if iscores is not None:
                for idrug, di in enumerate(drugNames):
                    s = iscores[idrug]
                    if len(di) > 15:
                        dx = ddrug2id[di]
                        drugId2Chem[di] = dx
                        di = dx
                    s0 = utils.get_dict(ds, di, 0)
                    s = max(s0, s)
                    ds[di] = s
            else:
                print("Missing protein binding: ", p)
        ic += 1
        print(ic, len(ds))

        if len(ds) > 0:
            pairs_str = []
            pairs_sorted = utils.sort_dict(ds)
            for nm, sc in pairs_sorted:
                pairs_str.append("%s|%.2f" % (nm, sc))
            fout.write("%s||%s\n" % (parts[0], "\t".join(pairs_str)))
            cc += len(ds)
            n += 1
        else:
            print(line)
    fout.close()
    fin.close()
    fmap = open("%s/DBID2Chem%s.txt" % (params.W_DIR, params.GEN_SUFFIX), "w")
    for k, v in drugId2Chem.items():
        fmap.write("%s||%s\n" % (v, k))
    fmap.close()
    print("Avg: ", cc / n)


def export_gene_drug_sub():
    import joblib
    dProtein2Drug = joblib.load("./HpyerAttentionDTI/tmp/protein2Drugs_0%s.npy" % params.GEN_SUFFIX)
    print("Num d protein 2 drug: ", len(dProtein2Drug))
    ddrug2id = loaddrugname2id()
    fin = open("%s/SelectedGene2Protein%s.txt" % (params.W_DIR, params.GEN_SUFFIX))
    fout = open("%s/Gene2Drugs%s.txt" % (params.W_DIR, params.GEN_SUFFIX), "w")
    drugId2Chem = {}
    cc = 0
    n = 0
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("||")
        protins = parts[1].split(",")
        ds = dict()
        for p in protins:
            dd = get_dict(dProtein2Drug, p, -1)

            if dd != -1:
                for di, s in dd:

                    if len(di) > 15:
                        dx = ddrug2id[di]
                        drugId2Chem[di] = dx
                        di = dx
                    s0 = utils.get_dict(ds, di, 0)
                    s = max(s0, s)
                    ds[di] = s
            else:
                print("Missing protein binding: ", p)

        print(len(ds))

        if len(ds) > 0:
            pairs_str = []
            pairs_sorted = utils.sort_dict(ds)
            for nm, sc in pairs_sorted:
                pairs_str.append("%s|%.2f" % (nm, sc))
            fout.write("%s||%s\n" % (parts[0], "\t".join(pairs_str)))
            cc += len(ds)
            n += 1
        else:
            print(line)
    fout.close()
    fin.close()
    fmap = open("%s/DBID2Chem%s.txt" % (params.W_DIR, params.GEN_SUFFIX), "w")
    for k, v in drugId2Chem.items():
        fmap.write("%s||%s\n" % (v, k))
    fmap.close()
    print("Avg: ", cc / n)


if __name__ == "__main__":
    match()
    # export_gene_drug()
