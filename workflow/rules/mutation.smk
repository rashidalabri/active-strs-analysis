rule parse_sample_mismatches_target:
    input:
        expand("resources/mutations/{variant_catalog}/{sample}.tsv", sample=SAMPLES.index, variant_catalog=['active'])


rule parse_sample_mismatches:
    input:
        bam="resources/realigned_bam/{variant_catalog}/{sample}/{sample}_realigned.bam",
        bai="resources/realigned_bam/{variant_catalog}/{sample}/{sample}_realigned.bam.bai"
    output:
        "resources/mutations/{variant_catalog}/{sample}.tsv"
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/parse_relative_mismatches.py"

rule plot_mismatches:
    input:
        expand("resources/mutations/{variant_catalog}/{sample}.tsv", sample=SAMPLES.index, variant_catalog=['active'])
    output:
        "results/mismatches.pdf"
    conda:
        "../envs/plot.yaml"
    notebook:
        "../notebooks/plot_mismatch.py.ipynb"