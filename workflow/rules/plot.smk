include: "mutyper.smk"


from itertools import combinations


rule plot_spectra:
    input:
        stratify=expand(rules.mutyper_spectra_stratify.output.spectra, rgn=RGNS),
        full=rules.mutyper_spectra.output.spectra,
    output:
        violin="results/plots/violin.pdf",
        heatmap="results/plots/heatmap.pdf",
        pca="results/plots/pca.pdf",
    log:
        "logs/plots/spectra.log",
    conda:
        "../envs/R.yml"
    script:
        "../scripts/spectra.R"


rule plot_spectra_stratify:
    input:
        full=rules.mutyper_spectra_stratify.output.spectra,
    output:
        violin="results/plots/stratify/{rgn}/violin.pdf",
        heatmap="results/plots/stratify/{rgn}/heatmap.pdf",
        pca="results/plots/stratify/{rgn}/pca.pdf",
    log:
        "logs/plots/spectra_{rgn}.log",
    conda:
        "../envs/R.yml"
    script:
        "../scripts/spectra.R"


rule plot_comparison:
    input:
        lambda wc: (rules.mutyper_spectra_stratify.output.spectra).format(rgn=wc.name1),
        lambda wc: (rules.mutyper_spectra_stratify.output.spectra).format(rgn=wc.name2),
        lambda wc: (rules.mutyper_spectra_targets.output.targets).format(rgn=wc.name1),
        lambda wc: (rules.mutyper_spectra_targets.output.targets).format(rgn=wc.name2),
    output:
        plot="results/plots/stratify/compare/{name1}_{name2}/spectra.pdf",
        fold="results/plots/stratify/compare/{name1}_{name2}/log_fold.pdf",
        targets="results/plots/stratify/compare/{name1}_{name2}/targets.pdf",
        pca="results/plots/stratify/compare/{name1}_{name2}/pca.pdf",
        heatmap="results/plots/stratify/compare/{name1}_{name2}/heatmap_{name1}_vs_{name2}.pdf",
        heatmap2="results/plots/stratify/compare/{name1}_{name2}/heatmap_{name1}_divide_{name2}.pdf",
        spectra_table="results/tables/stratify/spectra_{name1}_{name2}.tbl",
    log:
        "logs/plots/compare/{name1}_{name2}.log",
    conda:
        "../envs/R.yml"
    script:
        "../scripts/compare.R"


rule plot_comparison_corrected:
    input:
        lambda wc: (rules.mutyper_spectra_correction.output.spectra).format(
            rgn=wc.name1
        ),
        lambda wc: (rules.mutyper_spectra_correction.output.spectra).format(
            rgn=wc.name2
        ),
        lambda wc: (rules.mutyper_spectra_targets.output.targets).format(rgn=wc.name1),
        lambda wc: (rules.mutyper_spectra_targets.output.targets).format(rgn=wc.name2),
    output:
        plot="results/plots/stratify/compare_corrected/{name1}_{name2}/spectra.pdf",
        fold="results/plots/stratify/compare_corrected/{name1}_{name2}/log_fold.pdf",
        targets="results/plots/stratify/compare_corrected/{name1}_{name2}/targets.pdf",
        pca="results/plots/stratify/compare_corrected/{name1}_{name2}/pca.pdf",
        heatmap="results/plots/stratify/compare_corrected/{name1}_{name2}/heatmap_{name1}_vs_{name2}.pdf",
        heatmap2="results/plots/stratify/compare_corrected/{name1}_{name2}/heatmap_{name1}_divide_{name2}.pdf",
        spectra_table="results/tables/stratify/compare_corrected/spectra_{name1}_{name2}.tbl",
    log:
        "logs/plots/compare_corrected/{name1}_{name2}.log",
    conda:
        "../envs/R.yml"
    script:
        "../scripts/compare.R"


def make_plot_comparison_outputs(wc):
    for name1, name2 in combinations(RGNS, 2):
        rtn = expand(rules.plot_comparison.output, name1=name1, name2=name2)
        rtn += expand(rules.plot_comparison_corrected.output, name1=name1, name2=name2)
        for f in rtn:
            yield f


rule make_comparisons:
    input:
        make_plot_comparison_outputs,
