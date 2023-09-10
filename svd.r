#!/usr/bin/env Rscript

#' Truncated singular value decomposition of Darwin's Ark GWAS data
#' 
#' @author Kathleen Morrill Pirovich
#' @date 2023-09-05
#' @version 1.0.0
#' @license MIT

library(tidyverse)
library(data.table)
library(irlba)
library(ggrepel)
library(ggpubr)
library(plotly)
library(patchwork)

## GWAS z-scores
# Read in z-score matrix file
mat = fread("~/Documents/GitHub/darwins-ark-latent-gwas/final_z_matrix.csv")

## Phenotypes
# Read in phenotype index file
phe = read_csv("~/Documents/GitHub/darwins-ark-latent-gwas/final_phe_info.csv")
phe$index = seq(1,nrow(phe))-1
phe$phe_code = sub('.*phe-(.*?)_.*', '\\1', basename(phe$file))

# Name phenotypes
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

# Remove those GWA
phe = phe %>% arrange(index) %>% filter(!index %in% id_rm)

phe = phe %>%
  mutate(phe = if_else(!is.na(string),
                       paste(string," (Q#",item,")",sep=""),
                       if_else(!is.na(factor_name),
                               paste(factor_name," (F#",factor,")", sep = ""),
                               phe_code)))

var = mat %>% select(chr=Chr, pos=bp, snp=SNP, ref=A2, alt=A1)
diff = setdiff(as.character(phe$index), colnames(mat))
mat = mat %>% select(any_of(as.character(phe$index)))
phe = phe %>% filter(index != diff)

#mat = as.matrix(mat[,-c(1,2,3,4,5)])
mat = as.matrix(mat)
mat[is.na(mat)] <- 0  # replacing NAs with 0
mat[is.infinite(mat)] <- 0  # replacing Inf with 0

# Run tSVD
res_tsvd <- irlba(t(mat), nv = 40)

# Save outputs
saveRDS(res_tsvd$d, "~/Documents/GitHub/darwins-ark-latent-gwas/d_values.rds")
saveRDS(res_tsvd$u, "~/Documents/GitHub/darwins-ark-latent-gwas/u_matrix.rds")
saveRDS(res_tsvd$v, "~/Documents/GitHub/darwins-ark-latent-gwas/v_matrix.rds")
write_csv(as.data.frame(res_tsvd$d), "~/Documents/GitHub/darwins-ark-latent-gwas/d_values.csv")
write_csv(as.data.frame(res_tsvd$u), "~/Documents/GitHub/darwins-ark-latent-gwas/u_values.csv")
write_csv(as.data.frame(res_tsvd$v), "~/Documents/GitHub/darwins-ark-latent-gwas/v_values.csv")

# Scree Plot: relative variance explained by each of the components
p.scree = data.frame(S = res_tsvd$d,
                     Index = 1:length(res_tsvd$d)) %>%
  ggplot(aes(x = Index, 
             y = S)) +
  geom_point() +
  geom_line() +
  xlab("component index") +
  ylab("singular value") +
  theme_pubr()
p.scree

# Load PHATE output
phate_output = read_csv("~/Documents/GitHub/darwins-ark-latent-gwas/phate_output.csv", col_names = F)

p.phate = phate_output %>%
  mutate(index = row_number()-1) %>%
  merge(phe, by = "index") %>%
  ggplot(aes(x = X1,
             y = X2,
             label = phe)) +
  geom_point() +
  geom_text_repel() +
  theme_pubr()
ggplotly(p.phate)


# Calculate Fv = V * S for the first two components
Fv <- res_tsvd$v[, 1:2] %*% diag(res_tsvd$d[1:2])

# Create a dataframe for Fv
df_variants <- data.frame(PC1 = Fv[, 1], PC2 = Fv[, 2])

# Create a dataframe for U (phenotypes) for the first two components
df_phenotypes <- data.frame(PC1 = res_tsvd$u[, 1], 
                            PC2 = res_tsvd$u[, 2])

# Scatter plot for variants
p <- ggplot(df_variants, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = 'Variant'), alpha = 0.5) +
  scale_color_manual(values = c("Variant" = "blue"))

# Add phenotype vectors as arrows
p <- p + geom_segment(data = df_phenotypes, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(type = "closed"), alpha = 0.5, color = "red") +
  geom_text(data = df_phenotypes, aes(label = Phenotype, x = PC1, y = PC2), vjust = 1, hjust = 1)

# Additional plot settings
p <- p + labs(title = "Biplot of Variants and Phenotypes",
              x = "Principal Component 1",
              y = "Principal Component 2") +
  theme_minimal()

## Squared cosine scores

# Extract U matrix and singular values
U <- res_tsvd$u
S <- diag(res_tsvd$d)

# Compute squared cosine scores (a unit-normalized row of US)
US <- U %*% S

# Compute squared cosine scores for each trait
squared_cosine_scores <- rowSums(US^2)
normalized_US <- sweep(US, 1, sqrt(squared_cosine_scores), FUN="/")

# Prepare data for plotting
data_long <- as.data.frame(normalized_US^2)
colnames(data_long) <- paste0("PC", 1:20)

data_long <- gather(data_long, "PC", "Score")

# Plot
ggplot(data_long, aes(x=PC, y=Score, fill=PC)) +
  geom_bar(stat="identity") +
  ggtitle("Squared Cosine Scores for All Traits Across 20 Components") +
  ylab("Score") +
  xlab("Principal Component")

## Phenotype projections in variant space

# Extract matrices from tSVD results
U <- res_tsvd$u # phenotypes
V <- res_tsvd$v # variants

# Generate phenotype vectors


as.data.frame(V) %>%


# Create data for phenotype vectors (selected phenotypes: e.g., 1, 2, 3)
selected_phenotypes <- c(54, 22, 34)
phenotype_data <- as.data.frame(U[selected_phenotypes, ])
colnames(phenotype_data) <- c("PC1", "PC2")
phenotype_data$Phenotype <- as.factor(selected_phenotypes)

# Plot
ggplot() +
  geom_point(data = var_pca_data, aes(x = PC1, y = PC2), color = "blue") +
  geom_segment(data = phenotype_data, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(type = "closed", length = unit(0.2, "inches")), color = "red") +
  geom_text(data = phenotype_data, aes(x = PC1, y = PC2, label = Phenotype), nudge_x = 0.1) +
  ggtitle("First Two Principal Components of Variants with Phenotype Vectors") +
  xlab("PC1") +
  ylab("PC2")



# 30 - ear shape / HMGA2

# Component 02 - food sensitivity
## lncRNA near Prdm8
## DSCAML1
## C1QTNF3


# Plot variant and phenotype contributions
calculate_contribution_scores <- function(U_matrix, V_matrix, i) {
  # For phenotypes
  cntrphe <- U_matrix[, i]^2
  
  # For variants
  cntrvar <- V_matrix[, i]^2
  
  # Normalize scores so they sum to 1
  cntrphe <- cntrphe / sum(cntrphe)
  cntrvar <- cntrvar / sum(cntrvar)
  
  list(
    phenotype_contribution = cntrphe,
    variant_contribution = cntrvar
  )
}

i = 1

result = calculate_contribution_scores(res_tsvd$u, res_tsvd$v, i)

con.phe = bind_cols((phe %>% select(phe)),
          data.frame(score=result$phenotype_contribution)) %>%
  arrange(desc(score)) %>%
  head(5) %>%
  mutate(neg_log_score = -log10(score))

con.var = bind_cols(var,
                    data.frame(score=result$variant_contribution)) %>%
  mutate(
    neg_log_score = -log10(score),
    chr = as.factor(chr)
  )

p.con.var <- ggplot(con.var, 
                    aes(x = pos, 
                        y = neg_log_score, 
                        color = chr)) +
  geom_point() +
  scale_color_manual(values = rep(c("blue", "red"), length(unique(con.var$chr)) / 2)) +
  facet_wrap(~ chr, scales = "free_x") +
  theme_minimal() +
  labs(
    title = "Genome-wide Scores",
    x = "Position",
    y = "-log10(Score)"
  )

p.con.phe <- ggplot(con.phe, 
                    aes(x = phe, 
                        y = neg_log_score)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "Top 5 Phenotype Contributions",
    x = "Phenotype",
    y = "-log10(Score)"
  )

# Combine using patchwork
combined_plot <- p1 / p2

# Show the plot
combined_plot

