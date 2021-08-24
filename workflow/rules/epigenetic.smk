rule all_epigenetic:
    input:
        get_histone_pdf_file_names

rule download_histone_file:
    output:
        gz=temp("resources/encode/histone/{accession}.bed.gz"),
        unsorted=temp("resources/encode/histone/{accession}.bed.temp"),
        sorted="resources/encode/histone/{accession}.sorted.bed"
    params:
        url=lambda wildcards: ENCODE_HISTONE_FILES.loc[wildcards.accession, 'File download URL']
    shell:
        "wget -O {output.gz} {params.url} && "
        "gzip -dc {output.gz} > {output.unsorted} && "
        "sort -k1,1 -k2,2n {output.unsorted} > {output.sorted}"

# rule sort_histone_file:
#     input:
#         "resources/encode/histone/{accession}.bed"
#     output:
#         "resources/encode/histone/{accession}.sorted.bed"
#     shell:
#         "sort -k1,1 -k2,2n {input} > {output}"

rule histone_bedtools_closest:
    input:
        catalog="resources/catalog_bed_standard/{catalog}.bed",
        histone="resources/encode/histone/{accession}.sorted.bed"
    output:
        "resources/distance/histone/{catalog}_{accession}.bed"
    conda:
        "../envs/bedtools.yaml"
    envmodules:
        "bedtools/2.27.1"
    shell:
        "bedtools closest -D ref -a {input.catalog} -b {input.histone} > {output}"

rule plot_distance_histone:
    input:
        expand("resources/distance/histone/{catalog}_{accession}.bed", catalog=CATALOGS, allow_missing=True)
    params:
        assay=lambda wildcards: ENCODE_HISTONE_FILES.loc[wildcards.accession, 'Assay'],
        target=lambda wildcards: ENCODE_HISTONE_FILES.loc[wildcards.accession, 'Experiment target'],
        biosample=lambda wildcards: ENCODE_HISTONE_FILES.loc[wildcards.accession, 'Biosample term name'],
        audit_warning=lambda wildcards: ENCODE_HISTONE_FILES.loc[wildcards.accession, 'Audit WARNING'],
        audit_not_compliant=lambda wildcards: ENCODE_HISTONE_FILES.loc[wildcards.accession, 'Audit NOT_COMPLIANT'],
        audit_error=lambda wildcards: ENCODE_HISTONE_FILES.loc[wildcards.accession, 'Audit ERROR']
    output:
        "results/histone/{accession}.pdf"
    conda:
        "../envs/plot.yaml"
    notebook:
        "../notebooks/encode.py.ipynb"

