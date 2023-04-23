configfile: "config.yaml"


def get_raw_fastqc_input(wildcards):
    if wildcards.sampleType == "Collibri_standard_protocol":
        return expand("qc/raw/Collibri_standard_protocol-{tissueType}-Collibri-100_{repeats1}_L001_R{read}_001_fastqc.html", zip, sampleType=wildcards.sampleType, tissueType=wildcards.tissueType, k=wildcards.k, repeats1=wildcards.repeats1, read=[1,2])
    else:
        return expand("qc/raw/KAPA_mRNA_HyperPrep_-{tissueType}-KAPA-100_ng_total_{rna}_L001_R{read}_001_fastqc.html", zip, sampleType=wildcards.sampleType, tissueType=wildcards.tissueType, k=wildcards.k, rna=wildcards.rna, read=[1,2])

def get_fastqc_input(wildcards):
    if wildcards.sampleType == "Collibri_standard_protocol":
        if wildcards.tissueType == "HBR":
            repeats1_list = ["ng-2_S1", "ng-3_S2"]
        else:
            repeats1_list = ["ng-2_S3", "ng-3_S4"]
        k_str = "Collibri"
    else:
        if wildcards.tissueType == "HBR":
            repeats1_list = ["ng_total_RNA-2_S5", "ng_total_RNA-3_S6"]
        else:
            repeats1_list = ["ng_total_RNA-2_S7", "ng_total_RNA-3_S8"]
        #k_str = "Collibri"
        #repeats1_list = [ "ng_total_RNA-2_S7", "ng_total_RNA-3_S8"]
        k_str = "KAPA"

    #return "data/RNA-seq-sample-data/{wildcards.sampleType}-{wildcards.tissueType}-{k_str}-100-{repeats1_list}_L001_R{wildcards.read}_001.fastq.gz"
    return expand("data/RNA-seq-sample-data/{sampleType}-{tissueType}-{k}-100_{repeats1}_L001_R{read}_001.fastq.gz", sampleType=wildcards.sampleType, tissueType=wildcards.tissueType, k=k_str, repeats1=repeats1_list, read=[1,2], allow_missing=True)

def test(sampleType, tissueType):
    if sampleType == "Collibri_standard_protocol":
        sample = "Collibri_standard_protocol"
        if tissueType == "HBR":
            repeats1_list = ["ng-2_S1", "ng-3_S2"]
        else:
            repeats1_list = ["ng-2_S3", "ng-3_S4"]
        k_str = "Collibri"
    else:
        sample = "KAPA_mRNA_HyperPrep_"
        if tissueType == "HBR":
            repeats1_list = ["ng_total_RNA-2_S5", "ng_total_RNA-3_S6"]
        else:
            repeats1_list = ["ng_total_RNA-2_S7", "ng_total_RNA-3_S8"]
        #k_str = "Collibri"
        #repeats1_list = [ "ng_total_RNA-2_S7", "ng_total_RNA-3_S8"]
        k_str = "KAPA"

    if tissueType == "HBR":
        tissue = "HBR"
    else:
        tissue = "UHRR"

    #return "data/RNA-seq-sample-data/{wildcards.sampleType}-{wildcards.tissueType}-{k_str}-100-{repeats1_list}_L001_R{wildcards.read}_001.fastq.gz"
    return expand("qc/raw/{sampleType}-{tissueType}-{k}-100_{repeats1}_L001_R{read}_001_fastqc.html", sampleType=sample, tissueType=tissue, k=k_str, repeats1=repeats1_list, read=[1,2], allow_missing=True)


rule extract_data:
    input:
        "data/play_data_ref_annot.tar.gz"
    output:
        "genome/chr19_20Mb.bed",
        "genome/chr19_20Mb.fa",
        "genome/chr19_20Mb.gtf"
    shell:
        "tar -xf {input} -C genome"

rule fastqc_raw:
    input:
        #"data/RNA-seq-sample-data/{sample}_R{read}_001.fastq.gz"
        #"data/RNA-seq-sample-data/{sampleType}-{tissueType}-{k}-100_{repeats}_{rna}_{s}_L001_R{read}_001.fastqc.gz"
        #fastq_input = "if config["sampleType"] == "Collibri_standard_protocol" \
        #     data/RNA-seq-sample-data/{sampleType}-{tissueType}-{k}-100_{repeats}_{s}_L001_R{read}_001.fastq.gz" \ 
        #      else "data/RNA-seq-sample-data/{sampleType}-{tissueType}-{k}-100_ng_total_{rna}_{s}_L001_R{read}_001.fastq.gz"
        #fastq_input = lambda wildcards: f"data/RNA-seq-sample-data/{wildcards.sampleType}-{wildcards.tissueType}-Collibri-100_{wildcards.repeats}_L001_R{wildcards.read}_001.fastq.gz" \
        #                                 if wildcards.sampleType == "Collibri_standard_protocol" \ 
        #                                 else f"data/RNA-seq-sample-data/{wildcards.sampleType}-{wildcards.tissueType}-KAPA-100_ng_total_RNA-{wildcards.rna}_L001_R{wildcards.read}_001.fastq.gz"
        get_fastqc_input
        #"data/RNA-seq-sample-data/{sampleType}-{tissueType}-{k}-100_{repeats1}_L001_R{read}_001.fastq.gz"
    output:
        #"qc/raw/{sample}_R{read}_001_fastqc.html"
        #"qc/raw/{sampleType}-{tissueType}-{k}-100_{repeats}_{rna}_{s}_L001_R{read}_001_fastqc.html"
        #"qc/raw/{sampleType}-{tissueType}-Collibri-100_{sample_type}_L001_R{read}_001_fastqc.html",
        #"qc/raw/{sampleType}-{tissueType}-KAPA-100_{sample_type}_L001_R{read}_001_fastqc.html"
        #test("{sampleType}", "{tissueType}")
        "qc/raw/{sampleType}-{tissueType}-{k}-100_{repeats1}_L001_R{read}_001_fastqc.html",
    #params:
        #sample_type=lambda wildcards: {wildcards.repeates} if wildcards.sampleType == "Collibri_standard_protocol" else ng_total_RNA-{wildcards.rna}
    shell:
        "fastqc {input} -o qc/raw/"

rule multiqc_raw:
    input:
        #expand("qc/raw/{sampleType}-{tissueType}-Collibri-100_{repeats}_L001_R{read}_001_fastqc.html", zip, sampleType=config["sampleType"], tissueType=config["tissueType"], repeats=config["repeats"], read=[1,2],  allow_missing=True),
        #expand("qc/raw/{sampleType}-{tissueType}-KAPA-100_ng_total_{rna}_L001_R{read}_001_fastqc.html", zip, sampleType=config["sampleType"], tissueType=config["tissueType"], rna=config["rna"], read=[1,2],  allow_missing=True)
        test(config["sampleType"], config["tissueType"])
        #$expand("qc/raw/{sampleType}-{tissueType}-Collibri-100_{repeats1}_L001_R{read}_001_fastqc.html", sampleType=config["sampleType"], tissueType=config["tissueType"], k=config["k"], repeats1=config["repeats1"], read=[1,2],  allow_missing=True),
        #expand("qc/raw/KAPA_mRNA_HyperPrep_-{tissueType}-KAPA-100_ng_total_{rna}_L001_R{read}_001_fastqc.html", sampleType=config["sampleType"], tissueType=config["tissueType"], k=config["k"], rna=config["rna"], read=[1,2],  allow_missing=True),
        #get_raw_fastqc_input
    output:
        "qc/multiqc_report.html"
    shell:
        "multiqc qc/raw -o qc/"

        #Collibri_standard_protocol-HBR-Collibri-100_ng-2_S1_L001_R1_001.fastq.gz