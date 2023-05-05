configfile: "config.yaml"

(READS,) = glob_wildcards("data/RNA-seq-sample-data/{read}_R1_001.fastq.gz")

rule all:
    input:
        "qc/multiqc_report.html",
        "qc/multiqc_trimmed_report_multiqc_report.html",
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
        fq1=["data/RNA-seq-sample-data/{read}_R{r}_001.fastq.gz"],
        #test_input(config["sampleType"], config["tissueType"])
    output:
        fq1=["qc/raw/{read}_R{r}_001_fastqc.html"],
        #directory("qc/raw")
    shell:
        r"""
        fastqc {input.fq1} -o qc/raw
        """
        #"mkdir -p qc/raw" && "fastqc {input} -o qc/raw/"

rule multiqc_dir:
    input:
        expand("qc/raw/{read}_R{r}_001_fastqc.html", read=READS, r=["1","2"])
    output:
        "qc/multiqc_report.html"
    params:
        extra=""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log"
    wrapper:
        "v1.28.0/bio/multiqc"

rule bbduk_pe:
    input:
        sample=["data/RNA-seq-sample-data/{read}_R1_001.fastq.gz", "data/RNA-seq-sample-data/{read}_R2_001.fastq.gz"],
        adapters="adapters.fa",
    output:
        trimmed=["data/trimmed/{read}_R1_001.fastq.gz", "data/trimmed/{read}_R2_001.fastq.gz"]
    params:
        extra = lambda w, input: "ref={},adapters,artifacts ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10".format(input.adapters),
    resources:
        mem_mb=4000,
    threads: 7
    wrapper:
        "v1.28.0/bio/bbtools/bbduk"
        
rule fastqc_trimmed:
    input:
        #test_input(config["sampleType"], config["tissueType"], "trimmed")
        #directory("data/trimmed")
        fq1=["data/trimmed/{read}_R1_001.fastq.gz", "data/trimmed/{read}_R2_001.fastq.gz"],
    output:
        fq1=["qc/trimmed_qc/{read}_R1_001_fastqc.html", "qc/trimmed_qc/{read}_R2_001_fastqc.html"],
        #directory("qc/trimmed_qc")
    shell:
        r"""
        mkdir -p qc/trimmed_qc
        fastqc {input.fq1} -o qc/trimmed_qc
        """

rule multiqc_dir1:
    input:
        expand("qc/trimmed_qc/{read}_R{r}_001_fastqc.html", read=READS, r=["1","2"])
    output:
        "qc/multiqc_trimmed_report_multiqc_report.html"
    params:
        extra=""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log"
    wrapper:
        "v1.28.0/bio/multiqc"

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

rule star_pe:
    input:
        fq1=["data/trimmed/{read}_R1_001.fastq.gz", "data/trimmed/{read}_R2_001.fastq.gz"],
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
        export test={input.samples}
        if [[ "$test" == *"Collibri"* ]]; then
            strand="1"
        else
            strand="2"
        fi
        echo "strand"
        echo $strand

        samtools sort -n {input.samples} | featureCounts -p -t exon -g gene_id -a {input.annotation} -o {output[0]} -s $strand
        """
