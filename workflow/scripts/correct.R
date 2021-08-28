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

spectra_f <- "https://eichlerlab.gs.washington.edu/help/mvollger/share/mutyper/SD_spectra.txt"
target_f <- "https://eichlerlab.gs.washington.edu/help/mvollger/share/mutyper/SD_targets.txt"
spectra_f <- snakemake@input$spectra
target_f <- snakemake@input$targets


spectra <- fread(spectra_f)
targets <- fread(target_f, col.names = c("three", "total_count"))



make_spectra_long <- function(spectra) {
    long <- colnames(spectra)[grepl(".*>.*", colnames(spectra))]
    spectra %>%
        pivot_longer(
            cols = long,
            names_to = "spectra",
            values_to = "count"
        ) %>%
        mutate(three = substr(spectra, 1, 3)) %>%
        data.table()
}


df <- make_spectra_long(spectra) %>%
    merge(targets, by = "three") %>%
    mutate(count = 100 * count / total_count) %>%
    select(-three, -total_count) %>%
    pivot_wider(id_cols = sample, names_from = spectra, values_from = count) %>%
    data.table()

fwrite(df,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    file = snakemake@output[[1]]
)