library(vroom)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(scales)
library(cowplot)
library(ggrepel)
library(data.table)
library(tidyverse)
library(ggforce)
library(glue)
# library(tidylog)

spectra_f_1 <- "results/spectra/stratify/SD_spectra.txt"
spectra_f_2 <- "results/spectra/stratify/Unique_spectra.txt"
spectra_f_1 <- snakemake@input[[1]]
spectra_f_2 <- snakemake@input[[2]]

name1 <- "SD"
name2 <- "Unique"
name1 <- snakemake@wildcards$name1
name2 <- snakemake@wildcards$name2

out_plot <- "~/Desktop/1_2.pdf"
out_fold <- "~/Desktop/log_fold.pdf"
out_fold <- snakemake@output$fold
out_plot <- snakemake@output$plot



make_spectra_matrix <- function(spectra) {
    spec.m <- as.matrix(spectra[, -1])
    row.names(spec.m) <- t(spectra[, 1])
    spec.m
}

make_spectra_long <- function(spectra) {
    long <- colnames(spectra)[grepl(".*>.*", colnames(spectra))]
    spectra %>%
        pivot_longer(
            cols = long,
            names_to = "spectra",
            values_to = "count"
        ) %>%
        mutate(first_two_bases = substr(spectra, 1, 2)) %>%
        data.table()
}


read_spectra <- function(f) {
    spectra <- fread(f, sep = "\t")
    spectra.m <- make_spectra_matrix(spectra)
    pca_res <- prcomp(spectra.m, center = TRUE, scale. = TRUE)
    spectra$PC1 <- pca_res$x[, 1]
    spectra$PC2 <- pca_res$x[, 2]

    list(
        df = spectra,
        long = make_spectra_long(spectra),
        pca = pca_res,
        m = spectra.m
    )
}
spec1 <- read_spectra(spectra_f_1)
spec2 <- read_spectra(spectra_f_2)
l <- list(spec1$long, spec2$long)
names(l) <- c(name1, name2)
spec <- bind_rows(l, .id = "stratify") %>%
    group_by(stratify) %>%
    mutate(percent = 100 * count / sum(count)) %>%
    data.table()

pdf(out_plot, height = 11, width = 8)
for (two in unique(spec$first_two_bases)) {
    p <- spec %>%
        filter(first_two_bases == two) %>%
        ggplot(
            aes(
                x = spectra,
                y = percent,
                color = spectra,
            )
        ) +
        geom_violin() +
        geom_jitter() +
        scale_x_discrete(guide = guide_axis(n.dodge = 6)) +
        cowplot::theme_minimal_vgrid() +
        facet_grid(~stratify) +
        theme(legend.position = "none")
    print(p)
}
dev.off()


fold_change.df <- spec %>%
    group_by(stratify, spectra) %>%
    summarise(count = sum(count)) %>%
    ungroup() %>%
    group_by(stratify) %>%
    mutate(percent = 100 * count / sum(count)) %>%
    ungroup() %>%
    pivot_wider(
        id_cols = spectra,
        names_from = stratify,
        values_from = "percent"
    ) %>%
    mutate(
        log_change = log2(
            (!!as.name(name1)) / (!!as.name(name2))
        )
    ) %>%
    drop_na() %>%
    arrange(log_change) %>%
    mutate(spectra = factor(spectra, levels = unique(spectra))) %>%
    data.table()


pval.df <- spec %>%
    mutate(
        spectra = factor(spectra, levels = levels(fold_change.df$spectra))
    ) %>%
    group_by(stratify, spectra) %>%
    summarise(count = sum(count)) %>%
    group_by(stratify) %>%
    mutate(catagory_count = sum(count)) %>%
    group_by(spectra) %>%
    pivot_wider(
        id_cols = spectra,
        names_from = stratify,
        values_from = c("catagory_count", "count")
    ) %>%
    drop_na() %>%
    rowwise() %>%
    mutate(
        mat = list(matrix(as.numeric(cur_data()[1:4]), nrow = 2)),
        pval = fisher.test(mat)$p.value,
    ) %>%
    ungroup() %>%
    mutate(
        pval.corrected = pval * n()
    ) %>%
    data.table()

p.fold <- fold_change.df %>%
    merge(pval.df, by = "spectra") %>%
    ggplot(
        aes(
            x = spectra,
            y = log_change,
            fill = pval.corrected < 0.05,
        )
    ) +
    ylab(glue("log2({name1}/{name2})")) +
    geom_bar(stat = "identity") +
    scale_x_discrete(guide = guide_axis(n.dodge = 6)) +
    cowplot::theme_minimal_grid() +
    theme(legend.position = "top")
scale <- 1
ggsave(out_fold, height = 12 * scale, width = 24 * scale, plot = p.fold)




mutate(
    mat = list(matrix(
        c(count_SD, catagory_count_SD - count_SD, count_Unique, catagory_count_Unique - count_Unique),
        ncol = 2, nrow = 2
    )),
) %>%
    ungroup() %>%
    mutate(
        pval.corrected = pval * n()
    ) %>%
    data.table()
pval.df[pval.corrected < 0.05]

pval.df$mat[[1]]
fisher.test(pval.df$mat[[1]])