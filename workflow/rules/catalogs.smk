localrules:
    collect_catalogs,
    download_estrs,
    download_gangstr,
    download_trf,
    convert_estrs_xlsx_to_tsv,
    unzip_raw_catalogs,
    convert_catalog_vcf_to_tsv,
    convert_catalog_txt_to_tsv,
    preprocess_novel_json,
    standardize_estrs,
    standardize_gangstr,
    standardize_trf,
    standardize_novel,
    convert_tsv_to_bed,

rule collect_catalogs:
    input:
        expand("resources/catalog_tsv_standard/{catalog}.tsv", catalog=config['catalogs']),
        expand("resources/catalog_bed_standard/{catalog}.bed", catalog=config['catalogs']),

# rule download_chain_file:
#     output:
#         "resources/misc/hg19ToHg38.over.chain"
#     shell:
#         "wget -O - {config[hg19_to_hg38_chain_url]} | gzip -dc > {output}"

rule download_estrs:
    output:
        temp("resources/catalog_raw/estrs.xlsx")
    shell:
        "wget --user-agent=\"Mozilla/5.0 (X11; Linux x86_64; rv:60.0) Gecko/20100101 Firefox/60.0\" "
        "-O {output} {config[catalog_urls][estrs]}"

rule download_trf:
    output:
        temp("resources/catalog_raw/trf.txt.gz")
    shell:
        "wget {config[catalog_urls][trf]} -O {output}"

rule download_gangstr:
    output:
        temp("resources/catalog_raw/gangstr.vcf.gz")
    shell:
        "wget {config[catalog_urls][gangstr]} -O {output}"

rule convert_estrs_xlsx_to_tsv:
    input:
        "resources/catalog_raw/estrs.xlsx"
    params:
        sheet="Dataset S1 - All eSTRs"
    output:
        "resources/catalog_tsv/estrs.tsv"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/xlsx_to_tsv.py"

# rule liftover_estrs:
#     input:
#         bed="resources/catalog_bed/estrs_hg19.bed",
#         chain="resources/misc/hg19ToHg38.over.chain",
#     output:
#         lifted="resources/catalog_bed/estrs.bed",
#         unmapped="resources/misc/estrs_unmapped.bed",
#     conda:
#         "../envs/liftover.yaml"
#     shell:
#         "liftOver {input[bed]} {input[chain]} {output[lifted]} {output[unmapped]}"

# rule liftover_estrs_bed_to_tsv:
#     input:
#         old_tsv="resources/catalog_tsv/estrs_hg19.tsv",
#         bed="resources/catalog_bed/estrs.bed"
#     output:
#         "resources/catalog_tsv/estrs.tsv"
#     conda:
#         "../envs/pandas.yaml"
#     script:
#         "../scripts/liftover_estrs_bed_to_tsv.py"

rule unzip_raw_catalogs:
    input:
        "resources/catalog_raw/{catalog}.{ext}.gz"
    output:
        temp("resources/catalog_raw/{catalog}.{ext}")
    shell:
        "gzip -dc {input} > {output}"

rule convert_catalog_vcf_to_tsv:
    input:
        "resources/catalog_raw/{catalog}.vcf"
    output:
        "resources/catalog_tsv/{catalog}.tsv"
    conda:
        "../envs/vcf-kit.yaml"
    shell:
        "vk vcf2tsv wide --print-header {input} > {output}"

rule convert_catalog_txt_to_tsv:
    input:
        "resources/catalog_raw/{catalog}.txt"
    output:
        "resources/catalog_tsv/{catalog}.tsv"
    shell:
        "cat {input} > {output}"

rule preprocess_novel_json:
    input:
        "resources/catalog_raw/novel_{version}.json"
    output:
        "resources/catalog_tsv/novel_{version}.tsv",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/preprocess_novel_json.py"

'''
Standardization rules
---
Rearrange columns of catalogs such that first four
columns are: chr, start, stop, motif in addition to
applying global filters specified in config file
'''

rule standardize_estrs:
    input:
        "resources/catalog_tsv/estrs.tsv"
    output:
        hg38="resources/catalog_tsv_standard/estrs.tsv",
        hg19="resources/catalog_tsv_standard/estrs_hg19.tsv"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/standardize_catalog/estrs.py"

rule standardize_gangstr:
    input:
        "resources/catalog_tsv/gangstr.tsv"
    output:
        "resources/catalog_tsv_standard/gangstr.tsv"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/standardize_catalog/gangstr.py"

rule standardize_trf:
    input:
        "resources/catalog_tsv/trf.tsv"
    output:
        "resources/catalog_tsv_standard/trf.tsv"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/standardize_catalog/trf.py"

rule standardize_novel:
    input:
        "resources/catalog_tsv/novel_{version}.tsv"
    output:
        "resources/catalog_tsv_standard/novel_{version}.tsv"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/standardize_catalog/novel.py"


rule convert_standard_tsv_to_bed:
    input:
        "resources/catalog_tsv_standard/{catalog}.tsv"
    output:
        "resources/catalog_bed_standard/{catalog}.bed"
    shell:
        "tail -n +2 {input} | cut -f1-4 | sort -k1,1 -k2,2n > {output}"