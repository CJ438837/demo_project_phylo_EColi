library(tidyverse)
library(rstatix)
library(ggpubr)

df <- read.csv("D:/Programation/phylo E.Coli/genome_features.csv")

df$genus <- as.factor(df$genus)

str(df)
summary(df)

desc_stats <- df %>%
  group_by(genus) %>%
  summarise(
    n_genomes = n(),
    mean_gc = mean(gc_content_percent),
    sd_gc = sd(gc_content_percent),
    mean_size_mb = mean(genome_size_bp) / 1e6
  )

desc_stats

ggplot(df, aes(gc_content_percent)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(
    title = "Distribution of GC content",
    x = "GC content (%)",
    y = "Count"
  )


normality <- df %>%
  group_by(genus) %>%
  shapiro_test(gc_content_percent)

normality


ggplot(df, aes(genome_size_bp / 1e6, gc_content_percent)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(
    title = "GC content vs genome size",
    x = "Genome size (Mb)",
    y = "GC content (%)"
  )

cor_test <- cor.test(
  df$genome_size_bp,
  df$gc_content_percent,
  method = "spearman"
)

cor_test

ggplot(df, aes(genus, gc_content_percent, fill = genus)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "GC content by genus",
    x = "Genus",
    y = "GC content (%)"
  ) +
  theme(legend.position = "none")

library(tidyverse)
library(ggdendro)
library(factoextra)

df <- read.csv("D:/Programation/phylo E.Coli/genome_features.csv")
row.names(df) <- df$genome_id

features <- df %>%
  select(gc_content_percent, genome_size_bp, n_contigs) %>%
  scale()

dist_matrix <- dist(features, method = "euclidean")

hc <- hclust(dist_matrix, method = "ward.D2")

fviz_dend(
  hc,
  k = 3,
  rect = TRUE,
  cex = 0.6,
  main = "Genome similarity clustering"
)

df$cluster <- cutree(hc, k = 3)

ggplot(df, aes(cluster, gc_content_percent)) +
  geom_boxplot(fill = "steelblue") +
  theme_minimal() +
  labs(
    title = "GC content by cluster",
    x = "Cluster",
    y = "GC content (%)"
  )


