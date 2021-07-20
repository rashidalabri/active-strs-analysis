'''
exon/intron/UTR/CDS
cCRE
eQTL
Overlap with eSTRs (Gymrek et al. Nature Genetics 2015; Fotsing et al. Nature Genetics 2019)
Early/late replication regions
Pathogenic regions
Cancer hotspots
TAD boundaries (as mentioned in Sun et al. Cell 2018)
'''

rule download_ccre:
    output:
        "resources/annotations/ccre.bed"
    shell:
        "wget {config[annotation_urls][ccre]} -O {output}"

rule distance_catalog_annotation:
    input:
        ann="resources/annotations/{annotation}.bed",
        cat="resources/catalog_bed_standard/{catalog}.bed"
    output:
        "resources/distance/{catalog}_{annotation}.bed"
    conda:
        
    shell:
        "bedtools closest -d -a {input.ann} -b {input.cat} > {output}"

rule distance_catalogs_ccre:
    input:
        expand("resources/distance/{catalog}_ccre.bed", catalog=config['catalogs'])

# rule plot_distance_catalog_annotation:
#     input:
#         expand("resources/distance/{catalog}_{annotation}.bed", catalog=config['catalogs'])
#     output:
#         "results/regulation/boxplot_distance_closest_{annotation}.pdf"

rule intersect_catalogs:
    input:
        a="resources/catalog_standard_bed/{catalog_a}.bed",
        b="resources/catalog_standard_bed/{catalog_b}.bed",
    output:
        "resources/overlap/{catalog_a}_{catalog_b}.bed"
    shell:
        "bedtools intersect -wa -a {a} -b {b} > {output}"

# rule plot_pct_overlap_with_catalogs:
#     input:
#         expand("resources/overlap/{a}_{b}.bed", b=config['catalogs'])
#     output:
#         "results/catalog_overlap.pdf"

    
rule annovar_setup:
    input:
        "resources/tools/annovar/table_annovar.pl"
    shell:
        ""
# [kaiwang@biocluster ~/]$ annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/

# [kaiwang@biocluster ~/]$ annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/

# [kaiwang@biocluster ~/]$ annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/ 

# [kaiwang@biocluster ~/]$ annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp147 humandb/ 

# [kaiwang@biocluster ~/]$ annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp30a humandb/