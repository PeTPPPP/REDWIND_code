import os
C_DIR = os.path.dirname(os.path.abspath(__file__))

ENV_NAME = "p39"

CONDA_SH = "/home/gpux1/anaconda3/etc/profile.d/conda.sh"
ATD_DIR = "/home/gpux1/Downloads/ADTX"

DTI_PREDICTED_PATH = "%s/GENEBINDINGFORTOPDRG_ALL20_filtered.csv" % ATD_DIR
EXCEPT_DRUGS = "%s/except_drugs.txt" % ATD_DIR
CON_DIR = "%s/dat" % C_DIR
CUSTOM_DRUG_LIST_PATH = "%s/custom_drug_list.txt" % CON_DIR
CUSTOM_GENE_LIST_PATH = "%s/custom_gene.txt" % CON_DIR
IMG_FOLDER = "/home/gpux1/Pictures"


EXPORT_BIN_PATH_1 = "export PATH=/home/gpux1/Downloads/ADFRsuite_x86_64Linux_1.0/bin:$PATH"
