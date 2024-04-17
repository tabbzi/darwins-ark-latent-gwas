#!/usr/bin/env Rscript

#' Truncated singular value decomposition of Darwin's Ark GWAS data
#' 
#' @author Kathleen Morrill Pirovich
#' @date 2023-09-05
#' @version 1.0.0
#' @license MIT

## LIBRARIES ----
library(tidyverse)
library(data.table)
library(devtools)
library(GenomicSEM)
library(ggrepel)
library(ggpubr)
library(plotly)
library(patchwork)


## INPUT DATA ----

# Read in z-score matrix file
mat = fread("~/Documents/GitHub/darwins-ark-latent-gwas/final_z_matrix.csv")

# Read in phenotype index file
phe = read_csv("~/Documents/GitHub/darwins-ark-latent-gwas/final_phe_info.csv")

# Re-index phenotypes and extract information from file names
phe$index = seq(1,nrow(phe))-1
phe$phe_code = sub('.*phe-(.*?)_.*', '\\1', basename(phe$file))

# Re-name phenotypes
items = read_csv("~/Documents/GitHub/darwins-ark-latent-gwas/data_release/dat/DarwinsArk_20221120_questions.csv")
factors = data.frame(factor = seq(1,25),
                     factor_name = c("Human Sociability",
                                     "Food Selectivity",
                                     "Motor Patterns",
                                     "Repetition Frequency",
                                     "Social Strategy",
                                     "Biddability",
                                     "Food Sensitivity",
                                     "Proximity Seeking",
                                     "Shared Diet",
                                     "Remote Reliability",
                                     "Controlled Diet",
                                     "Healthy Weight",
                                     "Inappropriate Eating",
                                     "Dog Sociability",
                                     "Repetition Severity",
                                     "Food Solicitation",
                                     "Stimulus Sensitivity",
                                     "Arousal Level",
                                     "Food Appetite",
                                     "Attention Seeking",
                                     "Physical Mobility",
                                     "Agonistic Threshold",
                                     "Engagement Outdoors",
                                     "Separation Anxiety",
                                     "Engagement at Home"))
phe = phe %>%
  mutate(item = if_else(phe_code %like% "bq.",
                        as.integer(round(parse_number(phe_code) * 1000)),
                        NA_integer_,
                        NA_integer_)) %>%
  mutate(factor = if_else(phe_code %like% "fa.",
                          as.integer(round(parse_number(phe_code) * 100)),
                          NA_integer_,
                          NA_integer_)) %>%
  merge((items %>% select(id,string)), by.x = "item", by.y = "id", all.x = T) %>%
  merge(factors, by = "factor", all.x = T)

# Index factor GWA with only age covariates
id_rm = phe %>% filter(file %like% "qcov-age.fa") %>% pull(index)

# Index item GWA with quantitative scores
id_rm = c(id_rm, (phe %>% filter(phe_code %like% "bq." & !phe_code %like% "mean-binary") %>% pull(index)))

# Remove indexed GWA results
phe = phe %>% arrange(index) %>% filter(!index %in% id_rm)

phe = phe %>%
  mutate(phe = if_else(!is.na(string),
                       paste(string," (Q#",item,")",sep=""),
                       if_else(!is.na(factor_name),
                               paste(factor_name," (F#",factor,")", sep = ""),
                               phe_code)))

# Prepare variant, phenotype, and matrix data
var = mat %>% select(chr=Chr, pos=bp, snp=SNP, ref=A2, alt=A1)
diff = setdiff(as.character(phe$index), colnames(mat))
mat = mat %>% select(any_of(as.character(phe$index)))
phe = phe %>% filter(index != diff)

mat = as.matrix(mat)
mat[is.na(mat)] <- 0  # replacing NAs with 0
mat[is.infinite(mat)] <- 0  # replacing Inf with 0

n_phe = ncol(mat)
n_var = nrow(mat)

## Genomic SEM ----
