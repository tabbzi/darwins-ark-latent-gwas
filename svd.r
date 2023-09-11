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
library(irlba)
library(lavaan)
library(ggrepel)
library(ggpubr)
library(plotly)
library(patchwork)

## FUNCTIONS ----
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
    phenotype_eigenvalue = U_matrix[, i],
    variant_contribution = cntrvar,
    variant_eigenvalue = V_matrix[, i]
  )
}

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

# Index and separate morphology and behavior item GWA
id_mor = phe %>%
  filter(phe_code %like% "mq." | phe_code %like% "mp.") %>%
  pull(index)

mat_mor = mat %>% select(any_of(as.character(id_mor)))

id_beh = phe %>%
  filter(phe_code %like% "bq.") %>%
  pull(index)

mat_beh = mat %>% select(any_of(as.character(id_beh)))

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

mat_mor = as.matrix(mat_mor)
mat_mor[is.na(mat_mor)] <- 0
mat_mor[is.infinite(mat_mor)] <- 0
tsvd_result_mor <- irlba(t(mat_mor), nv = 5)

mat_beh = as.matrix(mat_beh)
mat_beh[is.na(mat_beh)] <- 0
mat_beh[is.infinite(mat_beh)] <- 0
tsvd_result_beh <- irlba(t(mat_beh), nv = 10)

# Run truncated singular value decomposition
n_svd = 20
tsvd_result <- irlba(t(mat), nv = n_svd)

# Save tSVD results to file
write_csv(as.data.frame(tsvd_result$d), "~/Documents/GitHub/darwins-ark-latent-gwas/d_values.csv")
write_csv(as.data.frame(tsvd_result$u), "~/Documents/GitHub/darwins-ark-latent-gwas/u_values.csv")
write_csv(as.data.frame(tsvd_result$v), "~/Documents/GitHub/darwins-ark-latent-gwas/v_values.csv")

# Scree Plot: relative variance explained by each of the components
p.scree = data.frame(S = tsvd_result$d,
                     Index = 1:length(tsvd_result$d)) %>%
  ggplot(aes(x = Index, 
             y = S)) +
  geom_line(color = "#2B555D") +
  geom_point(color = "#2B555D", shape = 21, fill = "#5FA29B", stroke = 1, size = 3, aes(color = NA)) +
  xlab("component index") +
  ylab("singular value") +
  theme_pubr() +
  theme(legend.position = "none",
        panel.background = element_blank(),
        plot.background = element_blank(),
        text = element_text(color = "#2B555D"),
        axis.title = element_text(color = "#2B555D"),
        axis.text = element_text(color = "#2B555D"),
        axis.ticks = element_line(color = "#2B555D")) 

ggsave(filename = "~/Documents/GitHub/darwins-ark-latent-gwas/scree-plot.png",
       plot = p.scree,
       device = "png",
       bg = "transparent", 
       width = 6,
       height = 5,
       units = "in")

# Summarize results
svd_gwa <- data.frame()

for (i in 1:20) {
  
  # Extract eigenvalues and component loadings
  eigenvalue_i <- tsvd_result$d[i]
  loading_i <- tsvd_result$u[, i]
  
  # Calculate weights for the weighted random effect model
  weights_i <- loading_i * eigenvalue_i
  
  # Calculate meta z-score for all variants for component i
  met_z <- rowSums(mat * weights_i) / sqrt(sum(weights_i^2))
  
  # Calculate meta p-value for all variants for component i
  met_p <- 2 * (1 - pnorm(abs(met_z)))
  
  tmp_df <- data.frame(
    svd = rep(i, length(met_z)),
    chr = var$chr,
    pos = var$pos,
    ref = var$ref,
    alt = var$alt,
    snp = var$snp,
    svd_var_con = tsvd_result$v[, i],
    met_z = met_z,
    met_p = met_p
  )
  
  svd_gwa <- rbind(svd_gwa, tmp_df)
}

svd_gwa = svd_gwa %>%
  group_by(svd) %>%
  mutate(met_p_cor_BH = p.adjust(met_p, method = "BH"),
         met_p_cor_BY = p.adjust(met_p, method = "BY"))

## Exploratory factor analysis ----
efa  = factanal(covmat = Ssmooth, factors = 2, rotation = "promax")

## Structural equation modeling ----
i = 1

# Extract important phenotypes and variants based on the i-th component
phe_ids <- order(abs(tsvd_result$u[,i]), decreasing = TRUE)[1:5]
var_ids <- order(abs(tsvd_result$v[,i]), decreasing = TRUE)[1:5]

# Generate variable names
var_names <- paste("V", var_ids, sep = "")
phe_names <- paste("P", phe_ids, sep = "")

# Construct the SEM model string
sem_model <- paste(
  paste(var_names, collapse = " + "),
  "~",
  paste(phe_names, collapse = " + "),
  sep = " "
)

# Convert matrix to a data frame and rename columns and rows
mat_df <- as.data.frame(mat)
colnames(mat_df) <- paste("P", 1:ncol(mat), sep = "")
rownames(mat_df) <- paste("V", 1:nrow(mat), sep = "")

# Fit the SEM model
fit <- sem(sem_model, data=mat_df)

# Summary of fit
summary(fit, fit.measures=TRUE)

## Unorganized code ----
# Plot meta_p values
i = 1

# Calculate maximum position for each chromosome and cumulative position offset
data_cum <- svd_gwa %>%
  filter(svd==i)
  group_by(chr) %>%
  summarise(max.pos = max(pos)) %>%
  mutate(pos.add = lag(cumsum(as.numeric(max.pos)), default = 0))

# Add cumulative positions to original data
con.var <- inner_join(con.var, data_cum, by = "chr") %>%
  mutate(pos.cum = pos + pos.add)

# Add color for each chromosome based on a rainbow color palette
chrom_colors <- rainbow(length(unique(con.var$chr)))
con.var <- con.var %>%
  mutate(color = chrom_colors[as.integer(chr)])

# Add color for each chromosome based on a bicolor palette
chrom_colors <- rep(x = c("#5FA29B","#2B555D"), length(unique(con.var$chr))/2)
con.var <- con.var %>%
  mutate(color = chrom_colors[as.integer(chr)])

# Calculate axis centers
axis.set <- con.var %>%
  group_by(chr) %>%
  summarize(center = mean(pos.cum))

# Label top SNPs (uppermost quartile per chromosome)
con.var <- con.var %>%
  group_by(chr) %>%
  mutate(top_snp = score >= quantile(score, 1-1e-8)) %>%
  ungroup()

# Generate the variant plot
p <- ggplot(con.var, aes(x = pos.cum, 
                                 y = eigen)) +
  geom_point(aes(color = color),
             alpha = 0.5) +
  geom_text_repel(
    data = subset(con.var, top_snp),
    aes(label = snp)
  ) +
  scale_x_continuous(label = axis.set$chr, breaks = axis.set$center) +
  scale_color_identity() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Genome-wide Scores",
       x = "Position",
       y = "-log10(Score)") +
  theme_pubr()


# Plot variant and phenotype contributions
con.phe.list = list()
for (i in 1:n_comp){
  con.phe.list[[i]] = 
    bind_cols((phe %>% select(phe)),
              data.frame(score=calculate_contribution_scores(tsvd_result$u, tsvd_result$v, i)$phenotype_contribution)) %>%
    as.data.frame() %>%
    arrange(desc(score)) %>%
    head(5) %>%
    mutate(component = i)
    
}
con.phe.list = rbindlist(con.phe.list)

i = 19

result = calculate_contribution_scores(tsvd_result$u, tsvd_result$v, i)

# Prepare phenotype data
con.phe = bind_cols((phe %>% select(phe)),
                    data.frame(score=result$phenotype_contribution,
                               eigen=result$phenotype_eigenvalue)) %>%
  arrange(desc(score)) %>%
  head(5)

# Prepare variant data
con.var = bind_cols(var,
                    data.frame(score=result$variant_contribution,
                               eigen=result$variant_eigenvalue)) %>%
  mutate(
    chr = as.factor(chr)
  )

# Calculate maximum position for each chromosome and cumulative position offset
data_cum <- con.var %>%
  group_by(chr) %>%
  summarise(max.pos = max(pos)) %>%
  mutate(pos.add = lag(cumsum(as.numeric(max.pos)), default = 0))

# Add cumulative positions to original data
con.var <- inner_join(con.var, data_cum, by = "chr") %>%
  mutate(pos.cum = pos + pos.add)

# Add color for each chromosome based on a rainbow color palette
chrom_colors <- rainbow(length(unique(con.var$chr)))
con.var <- con.var %>%
  mutate(color = chrom_colors[as.integer(chr)])

# Add color for each chromosome based on a bicolor palette
chrom_colors <- rep(x = c("#5FA29B","#2B555D"), length(unique(con.var$chr))/2)
con.var <- con.var %>%
  mutate(color = chrom_colors[as.integer(chr)])

# Calculate axis centers
axis.set <- con.var %>%
  group_by(chr) %>%
  summarize(center = mean(pos.cum))

# Label top SNPs (uppermost quartile per chromosome)
con.var <- con.var %>%
  group_by(chr) %>%
  mutate(top_snp = score >= quantile(score, 1-1e-8)) %>%
  ungroup()

# Generate the variant plot
p.con.var <- ggplot(con.var, aes(x = pos.cum, 
                                 y = eigen)) +
  geom_point(aes(color = color),
             alpha = 0.5) +
  geom_text_repel(
    data = subset(con.var, top_snp),
    aes(label = snp)
  ) +
  scale_x_continuous(label = axis.set$chr, breaks = axis.set$center) +
  scale_color_identity() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Genome-wide Scores",
       x = "Position",
       y = "-log10(Score)") +
  theme_pubr()

# Generate the phenotype plot
p.con.phe <- ggplot(con.phe, aes(x = score, y = phe)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Top 5 Phenotype Contributions",
       x = "-log10(Score)",
       y = "Phenotype") +
  theme_pubr()

# Combine the two plots using patchwork
combined_plot <- p.con.var / p.con.phe

# Show the plot
combined_plot




# Calculate all_scores_df using your previous functions
U_matrix <- tsvd_result$u
V_matrix <- tsvd_result$v
n_comp <- dim(U_matrix)[2]

all_scores_df <- NULL

for (i in 1:n_comp) {
  result <- calculate_contribution_scores(U_matrix, V_matrix, i)
  
  if (is.null(all_scores_df)) {
    all_scores_df <- as.data.frame(result$variant_contribution)
  } else {
    all_scores_df <- bind_cols(all_scores_df, result$variant_contribution)
  }
}

colnames(all_scores_df) <- paste("Comp", 1:n_comp, sep="_")

# Run analyze_components
analyze_components(all_scores_df, mat, n_comp)
