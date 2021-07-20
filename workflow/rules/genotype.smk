# rule download_genome:
#     input:
#         "resources/misc/creds.json"
#     output:
#         ""
#     container:
#         "docker://biodepot/gen3-client:2020.09-10-g0e21292__alpine_3.12__b0d3587c"
#     shell:
#         "gen3-client configure --profile=rashidalabri --cred={input} "
#         "--apiendpoint=https://icgc.bionimbus.org/ && "
#         "gen3-client download-single --profile=rashidalabri --guid={guid} "
#         "--no-prompt --skip-completed --download-path={tmpdir}"

# rule download_reference:
#     output:
#         "resources/reference/hg38.fa"
#     shell:
#         "wget -O - {config[reference_fasta_url]} | gunzip -c > {output}"

# rule genotype:
#     input:
#         bam="resources/bam/{donor}.bam",
#         bai="resources/bam/{donor}.bam.bai",
#         ref="resources/reference/hg38.fa",
#         catalog="resources/catalog_genotype_json/{catalog}.json"
#     output:
#         "resources/genotype/"
#     conda:
#         "../envs/expansionhunter.yaml"
#     shell:
#         "ExpansionHunter --reads {input[bam]} "
#          "--reference {input[ref]} "
#          "--variant-catalog {catalog} "
#          "--output-prefix <Prefix for the output files>"

rule extract_al