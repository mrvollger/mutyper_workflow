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
library(RColorBrewer)
# library(tidylog)

pop <- fread("https://eichlerlab.gs.washington.edu/help/mvollger/share/mutyper/hprc_year1_sample_metadata.txt", fill = TRUE)
spectra_f_1 <- "https://eichlerlab.gs.washington.edu/help/mvollger/share/mutyper/SD_spectra.txt"
target_f_1 <- "https://eichlerlab.gs.washington.edu/help/mvollger/share/mutyper/SD_targets.txt"
spectra_f_2 <- "https://eichlerlab.gs.washington.edu/help/mvollger/share/mutyper/Unique_spectra.txt"
target_f_2 <- "https://eichlerlab.gs.washington.edu/help/mvollger/share/mutyper/Unique_targets.txt"
spectra_f_1 <- snakemake@input[[1]]
spectra_f_2 <- snakemake@input[[2]]
target_f_1 <- snakemake@input[[3]]
target_f_2 <- snakemake@input[[4]]




name1 <- "SD"
name2 <- "Unique"
name1 <- snakemake@wildcards$name1
name2 <- snakemake@wildcards$name2

out_plot <- "~/Desktop/1_2.pdf"
out_fold <- "~/Desktop/log_fold.pdf"
out_targets <- "~/Desktop/targets_fold.pdf"
out_pca <- "~/Desktop/pca.pdf"
out_heatmap <- "~/Desktop/pca_heatmap.pdf"
spectra_table_f <- "~/Desktop/spectra_table.tbl"

out_fold <- snakemake@output$fold
out_plot <- snakemake@output$plot
out_targets <- snakemake@output$targets
out_pca <- snakemake@output$pca
out_heatmap <- snakemake@output$heatmap
out_heatmap2 <- snakemake@output$heatmap2
spectra_table_f <- snakemake@output$spectra_table



read_targets <- function(f, tag) {
    fread(f, col.names = c("spectra", "count")) %>%
        mutate(
            freq = count / sum(count),
            id = tag
        )
}

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
        mutate(first_base = substr(spectra, 1, 1)) %>%
        mutate(mid_base = substr(spectra, 2, 2)) %>%
        mutate(last_base = substr(spectra, 3, 3)) %>%
        mutate(first_two_bases = substr(spectra, 1, 2)) %>%
        mutate(first_three_bases = substr(spectra, 1, 3)) %>%
        data.table()
}


read_spectra <- function(f) {
    spectra <- fread(f, sep = "\t") %>%
        filter(sample != "CHM1_2") %>%
        filter(sample != "GRCh38_1") %>%
        filter(sample != "GRCh38_2") %>%
        data.table()

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


t1 <- read_targets(target_f_1, name1)
t2 <- read_targets(target_f_2, name2)
ts <- list(t1, t2)
names(ts) <- c(name1, name2)
targets <- bind_cols(ts) %>%
    data.table()
targets$fold_change <- targets[, 3] / targets[, 7]
targets$spectra <- targets[, 1]

spec1 <- read_spectra(spectra_f_1)
spec2 <- read_spectra(spectra_f_2)
l <- list(spec1$long, spec2$long)
names(l) <- c(name1, name2)
spec <- bind_rows(l, .id = "stratify") %>%
    group_by(stratify) %>%
    mutate(percent = 100 * count / sum(count)) %>%
    data.table()

pdf(out_plot, height = 8, width = 16)
for (bases in unique(spec$first_base)) {
    p <- spec %>%
        filter(first_base == bases) %>%
        ggplot(
            aes(
                x = spectra,
                y = percent,
                color = stratify,
            )
        ) +
        geom_jitter(alpha = 0.3, width = 0.25) +
        geom_violin() +
        scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
        cowplot::theme_minimal_vgrid() +
        # facet_grid(~stratify) +
        theme(legend.position = "top")
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
    # scale_y_continuous(trans = "log2") +
    # annotation_logticks(base = 2, sides = "l") +
    scale_x_discrete(guide = guide_axis(n.dodge = 6)) +
    cowplot::theme_minimal_grid() +
    theme(legend.position = "top")
scale <- 1
ggsave(out_fold, height = 12 * scale, width = 24 * scale, plot = p.fold)





targets.p <- targets %>%
    ggplot(
        aes(
            x = spectra,
            y = fold_change,
            fill = fold_change > 1
        )
    ) +
    ylab(glue("{name1}/{name2}")) +
    geom_bar(stat = "identity") +
    # scale_y_continuous(trans = "log2") +
    # annotation_logticks(base = 2, sides = "l") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    cowplot::theme_minimal_grid() +
    theme(legend.position = "top")
ggsave(out_targets, height = 12 * scale, width = 12 * scale, plot = targets.p)





#
#
# compare PCA
#
#
row.names(spec1$m) <- paste(row.names(spec1$m), name1, sep = "_")
row.names(spec2$m) <- paste(row.names(spec2$m), name2, sep = "_")
print(dim(spec1$m))
print(dim(spec2$m))
shared <- intersect(colnames(spec1$m), colnames(spec2$m))
spec_matrix <- rbind(spec1$m[, shared], spec2$m[, shared], fill = TRUE)
spec_matrix <- spec_matrix[row.names(spec_matrix) != "fill", ]
pca_res <- prcomp(spec_matrix, center = TRUE, scale. = TRUE)


spec_df <- data.table(spec_matrix) %>%
    mutate(sample = row.names(spec_matrix)) %>%
    separate(sample, sep = "_", into = c("Sample", "stratify"), remove = FALSE) %>%
    merge(pop, by = "Sample", all.x = T) %>%
    data.table()
spec_df$PC1 <- pca_res$x[, 1]
spec_df$PC2 <- pca_res$x[, 2]

pdf(out_pca, height = 12, width = 12)
autoplot(pca_res, data = spec_df, size = 0.001) +
    geom_point(aes(x = PC1, y = PC2, shape = stratify, color = Superpopulation)) +
    # geom_text_repel(aes(label = sample, color = Superpopulation)) +
    theme_cowplot() +
    theme(legend.position = "top")
dev.off()

pdf(out_heatmap, height = 8, width = 8)
heatmap(spec1$m / sum(spec1$m) - spec2$m / sum(spec2$m),
    col = colorRampPalette(rev(brewer.pal(8, "Spectral")))(25)
)
legend(
    x = "topright", legend = c("min", "ave", "max"),
    fill = colorRampPalette(rev(brewer.pal(8, "Spectral")))(3)
)
dev.off()

pdf(out_heatmap2, height = 8, width = 8)
heatmap(
    (spec1$m / sum(spec1$m)) / (spec2$m / sum(spec2$m)),
    col = colorRampPalette(rev(brewer.pal(8, "Spectral")))(25)
)
legend(
    x = "topright", legend = c("min", "ave", "max"),
    fill = colorRampPalette(rev(brewer.pal(8, "Spectral")))(3)
)
dev.off()


heatmap.df <- fold_change.df %>%
    separate(spectra, ">", into = c("ancestral", "right")) %>%
    mutate(
        derived = substr(right, 2, 2),
        first = substr(ancestral, 1, 1),
        abs_change = (
            (!!as.name(name1)) - (!!as.name(name2))
        ),
    )

p <- ggplot(heatmap.df, aes(y = ancestral, x = derived)) +
    theme_cowplot() +
    facet_grid(first ~ ., scales = "free", space = "free") +
    scale_fill_distiller(palette = "Spectral") +
    theme(legend.position = "bottom")
fig <- cowplot::plot_grid(p + geom_tile(aes(fill = log_change)), p + geom_tile(aes(fill = abs_change)))
ggsave(
    file = out_heatmap,
    plot = fig,
    height = 8, width = 12
)


#
# save long form table
#
fwrite(spec,
    file = spectra_table_f,
    sep = "\t",
    row.names = F,
    quote = F,
)