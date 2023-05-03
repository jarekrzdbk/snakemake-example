configfile: "config.yaml"

(READS,) = glob_wildcards("data/RNA-seq-sample-data/{read}_R1_001.fastq.gz")

def test_input(sampleType, tissueType, directory="RNA-seq-sample-data"):
    print("sample type")
    print(sampleType)
    print("tissueType")
    print(tissueType)
    results = []
    for key, valueSampleType in sampleType.items():
        for tissueKey, valueTissue in tissueType.items():
                if valueSampleType == "Collibri_standard_protocol":
                    sample = "Collibri_standard_protocol"
                    if valueTissue == "HBR":
                        repeats1_list = ["ng-2_S1", "ng-3_S2"]
                    else:
                        repeats1_list = ["ng-2_S3", "ng-3_S4"]
                    k_str = "Collibri"
                else:
                    sample = "KAPA_mRNA_HyperPrep_"
                    if valueTissue == "HBR":
                        repeats1_list = ["ng_total_RNA-2_S5", "ng_total_RNA-3_S6"]
                    else:
                        repeats1_list = ["ng_total_RNA-2_S7", "ng_total_RNA-3_S8"]
                    #k_str = "Collibri"
                    #repeats1_list = [ "ng_total_RNA-2_S7", "ng_total_RNA-3_S8"]
                    k_str = "KAPA"

                if valueTissue == "HBR":
                    tissue = "HBR"
                else:
                    tissue = "UHRR"

                result= expand("data/{directory}/{sampleType}-{tissueType}-{k}-100_{repeats1}_L001_R{read}_001.fastq.gz", sampleType=sample, tissueType=tissue, k=k_str, repeats1=repeats1_list, read=[1,2], directory=directory, allow_missing=True)     
                results.extend(result)       
    print(results)
    return results

rule all:
    input:
        #expand("star/{read}/aligned.bam", read=READS),
        #expand("mapped/{read}.bam.bai", read=READS),
        expand("results/{read}.featureCounts", read=READS)
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
        test_input(config["sampleType"], config["tissueType"])
    output:
        directory("qc/raw")
    shell:
        r"""
        mkdir -p qc/raw
        fastqc {input} -o qc/raw
        """
        #"mkdir -p qc/raw" && "fastqc {input} -o qc/raw/"

rule multiqc_raw:
    input:
        directory("qc/raw")
    output:
        "qc/multiqc_report.html"
    shell:
        "multiqc qc/raw -o qc/"

rule trim_reads:
    input:
        test_input(config["sampleType"], config["tissueType"])
    output:
        directory("data/trimmed")
    shell:
        """
        mkdir -p qc/trimmed
        for i in {input}
        do
        echo $i
        filepath=$i
        drn=$(dirname "$filepath")
        basename=$(basename "$filepath")
        prefix="$drn/"
        suffix="$basename#$prefix"
        echo "$suffix"
        bbduk.sh in=$i out="data/trimmed/$basename" ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10
        done
        """
        
rule fastqc_trimmed:
    input:
        #test_input(config["sampleType"], config["tissueType"], "trimmed")
        directory("data/trimmed")
    output:
        directory("qc/trimmed_qc")
    shell:
        """
        mkdir -p qc/trimmed_qc
        dir=data/trimmed/
        for i in data/trimmed/*;
        do
        echo "$i"
        fastqc "$i"  -o qc/trimmed_qc
        done
        """

rule multiqc_trimmed:
    input:
        directory("qc/trimmed_qc")
    output:
        "qc/multiqc_trimmed_report_multiqc_report.html"
    shell:
        "multiqc qc/trimmed_qc -o qc/ --title multiqc_trimmed_report"

rule star_index:
    input:
        "genome/chr19_20Mb.fa"
    output:
        directory("genome/chr19_20Mb")
    params:
        sjdbOverhang = 100
    shell:
        """
        echo {params.sjdbOverhang}
        STAR --runThreadN 4 --runMode genomeGenerate --genomeDir {output} \
        --genomeFastaFiles {input}
        """

#rule star_align:
#    input:
#        reads = test_input(config["sampleType"], config["tissueType"]),
#        index = directory("genome/chr19_20Mb")
#    output:
#        "star_aligned/Aligned.sortedByCoord.out.bam"
#    shell:
#        """
#        mkdir -p star_aligned
#        STAR --genomeDir {input.index} --readFilesIn {input.reads} \
#        --outFileNamePrefix star_aligned/ --outSAMtype BAM SortedByCoordinate
#        """

rule star_pe:
    input:
        fq1=["data/trimmed/{read}_R1_001.fastq.gz", "data/trimmed/{read}_R2_001.fastq.gz"],
        #fq2=["data/RNA-seq-sample-data/{sample}_R2_001.fastaq.gz"],
        index= directory("genome/chr19_20Mb")
    output:
        aln="star/{read}.bam"
    shell:
        """
        STAR --genomeDir {input.index} --readFilesIn {input.fq1} --outFileNamePrefix star/ --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c --outStd BAM_SortedByCoordinate > {output.aln}
        """

rule samtools_index:
    input:
        bam="star/{read}.bam",
    output:
        bai="star/{read}.bam.bai",
    shell:
        """
        echo {input.bam}
        samtools index {input.bam} {output.bai}
        """

rule featureCounts:
    input:
        samples="star/{read}.bam",
        bai="star/{read}.bam.bai",
        annotation="genome/chr19_20Mb.gtf"
    output:
        multiext(
            "results/{read}",
            ".featureCounts",
            ".featureCounts.summary"
        ),
    params:
        strand = "1" #if input(" ").startswith("Collibri") else "2"# or "2", depending on the sample preparation method
    shell:
        """
        if [[ "$test" == *"Collibri"* ]]; then
            strand="1"
        else
            strand="2"
        fi
        echo "strand"
        echo $strand

        samtools sort -n {input.samples} | featureCounts -p -t exon -g gene_id -a {input.annotation} -o {output[0]} -s $strand
        """
