import os

## FOLDERS
C_DIR = os.path.dirname(os.path.abspath(__file__))
RES_DIR = "%s/resource" % C_DIR #
TMP_DIR = "%s/tmp" % C_DIR #
RESULTS_DIR = "%s/results" % C_DIR #

W_DIR = "%s/wdir" % C_DIR
DATA_DIR = "%s/data" % C_DIR

## PARAMS
GEN_SUFFIX = ""
DEVICE = "cuda:0"
N_SEG = 1
I_SEG = -1
FULL_PREDICTION = True
use_mean_of_expression_values = True ## will use sum if False
restrict_to_specific_drugs = True
number_of_genes_to_bind_in_tops = 10

NORMALIZE = False
USE_DIFFERENCE = True
USE_ABS_EXPR = False # not used

## Ordering results
use_category = ""
order_by = "drugs"
GENEXP_LOG2 = True
NORMALIZE_RESULTS_BYGENENR = False
SIGN = 0.05
SIGN_MC_TYPE = "bonferroni"
NR_TOPKDRUGS = 100
## not tested:
BINDING_THRESHOLD = 0


## INPUT/OUTPUT FILES
#
GENE_PATH = "%s/gene_list_females.txt" % DATA_DIR
if GEN_SUFFIX == "M":
    GENE_PATH = "%s/gene_list_male.txt" % DATA_DIR
elif GEN_SUFFIX == "":
    GENE_PATH = "%s/gene_list.txt" % DATA_DIR
GENES_PREDICTION_ALL_FILE = "%s/Gene2Drug/9657Genes2Drugs/Gene2Drugs.txt" % RES_DIR #
PREDICTION_COUNTS_FILE = "%s/Top10Stats.txt" % RES_DIR #

PATHWAY_LIST_ALL_FILE = "%s/ALL_pathway_gene_list_total.txt" % RES_DIR #

DEG_GENES_RESULTS_ALL_S1_FILE = "%s/DEG_results_S1_all_SAAS.csv" % RES_DIR
DEG_GENES_RESULTS_ALL_S2_FILE = "%s/DEG_results_S2_all_SAAS.csv" % RES_DIR

PATHWAY_CATEGORIES_ALL_FILE = "%s/ALLMIGvsALLCONTR_SAAS_categories.csv" % RES_DIR
DrugBankDrugGenePathO = "%s/data/DrugBankProteinGeneTarget.txt" % C_DIR
DrugBankDrugGenePath = "%s/data/DrugBankProteinGeneTargetNN.txt" % C_DIR


GENE_DRUG_JSON = "%s/GeneDrug.json" % W_DIR
DRUGBANK_PROTEINIDNOMATPPING_PATH = "%s/data/BENoUniProt_filter.txt" % C_DIR

UNIPORT_REVIEWED_PATH = "/home/gpux1/Data/UniProtReviewed/uniprot_sprot.dat"

GENES_PREDICTION_CORRECTED_FILE = "%s/Genes2Drugs_modded.txt" % TMP_DIR


