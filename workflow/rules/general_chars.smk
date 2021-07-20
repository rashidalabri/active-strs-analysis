rule all_general_chars:
    input:
        expand("results/general_chars/motif_length_hist_{catalog}.pdf", catalog=config['catalogs']),
        expand("results/general_chars/ref_length_hist_{catalog}.pdf", catalog=config['catalogs']),
        expand("results/general_chars/gc_content_hist_{catalog}.pdf", catalog=config['catalogs']),
        expand("results/general_chars/histplot_ref_motif_len_{catalog}.pdf", catalog=config['catalogs']),

rule all_histplot_motif_length:
    input:
        expand("results/general_chars/motif_length_hist_{catalog}.pdf", catalog=config['catalogs'])
    
rule histplot_motif_length:
    input:
        "resources/catalog_tsv_standard/{catalog}.tsv"
    output:
        "results/general_chars/motif_length_hist_{catalog}.pdf"
    conda:
        "../envs/plot.yaml"
    notebook:
        "../notebooks/histplot_motif_length.py.ipynb"

rule all_histplot_ref_length:
    input:
        expand("results/general_chars/ref_length_hist_{catalog}.pdf", catalog=config['catalogs'])
    
rule histplot_ref_length:
    input:
        "resources/catalog_tsv_standard/{catalog}.tsv"
    output:
        "results/general_chars/ref_length_hist_{catalog}.pdf"
    conda:
        "../envs/plot.yaml"
    notebook:
        "../notebooks/histplot_ref_length.py.ipynb"

rule all_histplot_gc_content:
    input:
        expand("results/general_chars/gc_content_hist_{catalog}.pdf", catalog=config['catalogs'])
    
rule histplot_gc_content:
    input:
        "resources/catalog_tsv_standard/{catalog}.tsv"
    output:
        "results/general_chars/gc_content_hist_{catalog}.pdf"
    conda:
        "../envs/plot.yaml"
    notebook:
        "../notebooks/histplot_gc_content.py.ipynb"

rule all_histplot_ref_motif_len:
    input:
        expand("results/general_chars/histplot_ref_motif_len_{catalog}.pdf", catalog=config['catalogs'])
    
rule histplot_ref_motif_len:
    input:
        "resources/catalog_tsv_standard/{catalog}.tsv"
    output:
        "results/general_chars/histplot_ref_motif_len_{catalog}.pdf"
    conda:
        "../envs/plot.yaml"
    notebook:
        "../notebooks/histplot_ref_motif_len.py.ipynb"
