import params

def export_drugname_smile():
    fin = open("%s/DrugBankNameX.txt" % params.DATA_DIR)
    fout = open("%s/DrugBankNameID_SMILE.txt" % params.W_DIR, "w")
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("||")
        if parts[1] == "small molecule":
            nn = parts[5].strip()
            if len(nn) > 0:
                fout.write("%s||%s||%s\n" % (parts[0], parts[2], nn))
    fin.close()
    fout.close()

if __name__ == "__main__":
    export_drugname_smile()