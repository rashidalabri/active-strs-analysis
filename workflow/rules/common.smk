import pandas as pd

SAMPLES = pd.read_table(config["samples"]).set_index("Sample", drop=False)
# CATALOGS = config['catalogs']
CATALOGS = ['ActiveSTRs', 'TRF', 'InactiveSTRs']

def create_encode_histone_df(file):
    f = open(file, "r")
    lines = f.readlines()
    metadata = pd.read_table(str(lines[0]).replace("\"", ""))
    metadata = metadata.set_index("File accession")
    metadata = metadata[metadata["File analysis status"] == "released"]
    return metadata

ENCODE_HISTONE_FILES = create_encode_histone_df(config["encode_histone_file"])

def get_histone_pdf_file_names(wildcards):
    names = []
    for i, row in ENCODE_HISTONE_FILES.iterrows():
        target = row["Experiment target"]
        biosample = row["Biosample term name"]
        # accession = row["File accession"]
        # name = "results/histone/{}_{}_{}.pdf".format(target, biosample, i)
        name = "results/histone/{}.pdf".format(i)
        names.append(name)
    return names


