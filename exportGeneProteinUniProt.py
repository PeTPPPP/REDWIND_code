import params
import utils
import re

RE_ALPHABET = re.compile('[^a-zA-Z]')


class CounterX:
    def __init__(self, c=0):
        self.c = c

    def inc(self):
        self.c += 1


def get_line(f, c):
    line = f.readline()
    c.inc()
    return line


def parse_genename(s):
    parts = s.split(";")
    genes = []
    for part in parts:

        part = part.strip()
        if len(part) < 2:
            continue
        # print("P0: ", part)
        try:
            space_idx = part.index(' {')
        except:
            space_idx = -1

        if space_idx > 0:
            part = part[:space_idx]

        # print("P1:", part)
        if not part.__contains__("="):
            continue
        gene = part.split("=")[1]
        genes.append(gene)
    return genes

def extract2():
    fin = open(params.UNIPORT_REVIEWED_PATH)
    fout = open("%s/UniprotReviewedMapping2.txt" % params.W_DIR, "w")
    cc = 0
    c = CounterX()

    def read_next_sequence(f, cLine=""):
        while not cLine.startswith("SQ"):
            cLine = f.readline()
        ss = []
        cLine = f.readline()
        while not cLine.startswith("//"):
            ss.append(cLine)
            cLine = f.readline()
        ss = "".join(ss)
        ss = RE_ALPHABET.sub("", ss)
        return ss

    while True:
        line = get_line(fin, c)
        if line == "":
            break
        if line.startswith("ID "):

            cc += 1

            if cc % 10 == 0:
                print("\r%s" % cc, end="")
            proid1 = line[5:line.index(" ", 6)]
            line = get_line(fin, c)  # fin.readline()
            acid = line.strip()[5:]
            all_gene_names = []
            while not line.startswith("GN"):
                line = get_line(fin, c)
            while line.startswith("GN"):
                genename = line.strip()[5:]
                all_gene_names.append(genename)
                line = get_line(fin, c)
            genename = "".join(all_gene_names)
            # print(genename)
            osname = line.strip()[5:]
            ss = read_next_sequence(fin, line)
            fout.write("%s||%s||%s||%s||%s\n" % (proid1, acid, "#".join(parse_genename(genename)), osname, ss))

        else:
            continue

    fin.close()
    fout.close()
    print("Finished")


def print_seg(lid=46934642, ws=10):
    fin = open(params.UNIPORT_REVIEWED_PATH)
    c = 0
    while c < lid - ws:
        fin.readline()
        c += 1
    for i in range(2 * ws):
        print(fin.readline())


def export_gene_2_protein():
    fin = open("%s/UniprotReviewedMapping2.txt" % params.W_DIR)
    fout = open("%s/UniprotReviewedGene2Proteins.txt" % params.W_DIR, "w")
    d_gene_2proteins = dict()
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("||")
        proteins = parts[1].split(";")
        genes = parts[2]
        species = parts[3]
        gens = genes.split("#")
        all_gens = []
        if not species.__contains__("Human"):
            continue
        for gen in gens:
            gs = gen.split(",")
            for g in gs:
                all_gens.append(g.strip())

        for g in all_gens:
            protein_sets = utils.get_insert_dict(d_gene_2proteins, g, set())
            for protein in proteins:
                protein = protein.strip()
                if len(protein) > 2:
                    protein_sets.add(protein)
            # print(protein_sets)
    for k, v in d_gene_2proteins.items():
        fout.write("%s||%s\n" % (k, ",".join(list(v))))
    fout.close()


if __name__ == "__main__":
    utils.ensure_dir(params.W_DIR)
    extract2()
    export_gene_2_protein()
