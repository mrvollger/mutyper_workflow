#
# make an alignment chain if it doesn't exist
#
if "chain" not in config:

    rule make_bam:
        input:
            ref=REF,
            outgroup=OUTGROUP,
        output:
            bam=temp("results/chain/out-to-ref.bam"),
        log:
            "logs/chain/bam.log",
        conda:
            "../envs/env.yml"
        threads: 8
        shell:
            """
            minimap2 -ax asm20 -Y --eqx -t {threads} \
                {input.ref} {input.outgroup} \
                    | samtools view -u - \
                    | samtools sort -@ {threads} -m 8G - \
                > {output.bam}
            """

    rule make_psl:
        input:
            bam=rules.make_bam.output.bam,
        output:
            psl=temp("results/chain/out-to-ref.psl"),
        log:
            "logs/chain/psl.log",
        conda:
            "../envs/env.yml"
        shell:
            """
            bamToPsl {input.bam} {output.psl}  
            """

    rule make_chain:
        input:
            psl=rules.make_psl.output.psl,
        output:
            chain=CHAIN,
        log:
            "logs/chain/chain.log",
        conda:
            "../envs/env.yml"
        shell:
            """
            pslToChain {input.psl} {output.chain}
            """


#
#
#
rule prep_vcf:
    input:
        vcf=VCF,
    output:
        bcf=temp("results/vcf/input/{chrm}.bcf"),
        csi=temp("results/vcf/input/{chrm}.bcf.csi"),
    log:
        "logs/vcf/input/{chrm}.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        bcftools view -v snps {input.vcf} {wildcards.chrm} \
            | bcftools sort -m 8G - \
            | bcftools +fill-tags \
            | bcftools filter -i 'AN>0' \
            -Ob -o {output.bcf}

        bcftools index -f {output.bcf}
        """


rule prep_ref:
    input:
        reference=REF,
    output:
        ref=temp("results/ref-fasta/{chrm}.fa"),
        fai_ref=temp("results/ref-fasta/{chrm}.fa.fai"),
    log:
        "results/ref-fasta/{chrm}.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        samtools faidx {input.reference} {wildcards.chrm} \
            | seqtk seq -l 60 > {output.ref}
        samtools faidx {output.ref}
        """


if "ancestor" in config:

    rule prep_ancestor:
        input:
            ancestor=ANCESTOR_FA,
        output:
            fasta=temp("results/ancestral-fasta/{chrm}.fa"),
            fai=temp("results/ancestral-fasta/{chrm}.fa.fai"),
        log:
            "results/ancestral-fasta/{chrm}.log",
        conda:
            "../envs/env.yml"
        shell:
            """
            samtools faidx {input.ancestor} {wildcards.chrm} \
                | seqtk seq -l 60 > {output.fasta}
            samtools faidx {output.fasta}
            """


else:

    rule make_ancestor_per_chr:
        input:
            bcf=rules.prep_vcf.output.bcf,
            ref=rules.prep_ref.output.ref,
            out=OUTGROUP,
            chain=CHAIN,
        output:
            fasta=temp("results/ancestral-fasta/{chrm}.fa"),
            fai=temp("results/ancestral-fasta/{chrm}.fa.fai"),
        log:
            "results/ancestral-fasta/{chrm}.log",
        conda:
            "../envs/env.yml"
        shell:
            """
            mutyper ancestor \
                {input.bcf} \
                {input.ref} \
                {input.out} \
                {input.chain} \
            {output.fasta}
            samtools faidx {output.fasta}
            """


rule annotate_vcf:
    input:
        fasta="results/ancestral-fasta/{chrm}.fa",
        bcf=rules.prep_vcf.output.bcf,
    output:
        bcf=temp("results/vcf/mutyper/{chrm}.bcf"),
        csi=temp("results/vcf/mutyper/{chrm}.bcf.csi"),
    log:
        "logs/vcf/mutyper/{chrm}.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        mutyper variants {input.fasta} {input.bcf} \
            | bcftools sort -Ob -m 8G - \
            > {output.bcf}
        bcftools index -f {output.bcf}
        """


rule mutyper_vcf:
    input:
        bcf=expand(rules.annotate_vcf.output.bcf, chrm=CHRS),
        csi=expand(rules.annotate_vcf.output.csi, chrm=CHRS),
    output:
        bcf="results/vcf/mutyper.bcf",
        csi="results/vcf/mutyper.bcf.csi",
    log:
        "logs/mutyper_vcf.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        bcftools concat -Ob -a \
            {input.bcf} > {output.bcf}
        bcftools index -f {output.bcf}
        """


rule mutyper_spectra:
    input:
        bcf=rules.mutyper_vcf.output.bcf,
    output:
        spectra="results/spectra/spectra.txt",
    log:
        "logs/spectra.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        mutyper spectra {input.bcf} > {output.spectra}
        """


rule filter_stratify_bed:
    input:
        filter_bed=config["include"],
        bed=lambda wc: config["stratify"][wc.rgn],
    output:
        bed=temp("temp/spectra/stratify/{rgn}_filtered.bed"),
    log:
        "logs/spectra.{rgn}.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        bedtools intersect -a {input.bed} -b {input.filter_bed} \
            | bedtools sort -i - \
            | bedtools merge -i - \
        > {output.bed}
        """


rule seq_content_stratify:
    input:
        ref=REF,
        bed=rules.filter_stratify_bed.output.bed,
    output:
        tbl="results/spectra/stratify/{rgn}_seq_content.tbl",
    log:
        "logs/seq_content_{rgn}.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        bedtools nuc -fi {input.ref} -bed {input.bed} \
            | sed 's/$/\\t{wildcards.rgn}/' \
        > {output.tbl}
        """


rule seq_content:
    input:
        tbl=expand(rules.seq_content_stratify.output.tbl, rgn=RGNS),
    output:
        tbl="results/spectra/stratify/seq_content.tbl",
    log:
        "logs/seq_content.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        head -n1 {input.tbl[0]} > {output.tbl} 
        cat {input.tbl} | grep -v "#" >> {output.tbl}
        """


rule mutyper_spectra_stratify:
    input:
        bcf=rules.mutyper_vcf.output.bcf,
        #bed=lambda wc: config["stratify"][wc.rgn],
        bed=rules.filter_stratify_bed.output.bed,
    output:
        spectra="results/spectra/stratify/{rgn}_spectra.txt",
    log:
        "logs/spectra.{rgn}.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        bcftools view --regions-file {input.bed} {input.bcf} \
            | mutyper spectra - \
            > {output.spectra}
        """


rule mutyper_spectra_stratify_population:
    input:
        bcf=rules.mutyper_vcf.output.bcf,
        #bed=lambda wc: config["stratify"][wc.rgn],
        bed=rules.filter_stratify_bed.output.bed,
    output:
        spectra="results/spectra/stratify/{rgn}_spectra_population.txt",
    log:
        "logs/spectra.{rgn}.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        bcftools view --regions-file {input.bed} {input.bcf} \
            | mutyper spectra --population - \
            > {output.spectra}
        """


rule mutyper_spectra_targets:
    input:
        fasta=expand("results/ancestral-fasta/{chrm}.fa", chrm=CHRS),
        #bed=lambda wc: config["stratify"][wc.rgn],
        bed=rules.filter_stratify_bed.output.bed,
    output:
        targets="results/spectra/stratify/{rgn}_targets.txt",
        fasta=temp("temp/spectra/stratify/fasta/{rgn}_targets.fa"),
        fai=temp("temp/spectra/stratify/fasta/{rgn}_targets.fa.fai"),
    log:
        "logs/targets.{rgn}.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        cat {input.fasta} \
            | seqtk subseq - {input.bed} \
            | sed 's/:\|-/_/g' \
            > {output.fasta}
        samtools faidx {output.fasta}

        mutyper targets \
            {output.fasta} \
            > {output.targets}
        """


# | seqtk seq -M {input.bed} -c \


rule mutyper_spectra_ksfs:
    input:
        bcf=rules.mutyper_vcf.output.bcf,
        #bed=lambda wc: config["stratify"][wc.rgn],
        bed=rules.filter_stratify_bed.output.bed,
    output:
        ksfs="results/ksfs/{rgn}/ksfs.txt",
    log:
        "logs/ksfs/{rgn}.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        bcftools view \
            --regions-file {input.bed} \
                {input.bcf} \
            | bcftools filter --include 'AN=2*N_SAMPLES' \
            | mutyper ksfs -  \
        > {output.ksfs}
        """


rule mutyper_spectra_correction:
    input:
        targets=rules.mutyper_spectra_targets.output.targets,
        spectra=rules.mutyper_spectra_stratify.output.spectra,
    output:
        spectra="results/spectra/stratify/spectra_corrected_{rgn}.txt",
    log:
        "logs/spectra/correct.{rgn}.log",
    conda:
        "../envs/R.yml"
    script:
        "../scripts/correct.R"
