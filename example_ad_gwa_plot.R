ad_gwa = read_tsv("~/Documents/GitHub/darwins-ark-latent-gwas/example_ad_gwa.tsv")

pal = c("#EEE1CB", "#5FA29B", "#2B555D", "#FEAE03", "#E06545", "#769762", "#B95E79")

ad_gwa = ad_gwa %>%
  arrange(P) %>%
  mutate(Plog = -log10(P)) %>%
  mutate(gene_name = fct_reorder(Gene, Plog))

p = ad_gwa %>%
  ggplot(aes(x = gene_name, y = Plog)) +
  geom_segment(
    aes(x=gene_name, xend=gene_name, y=0, yend=Plog), 
    color=ifelse(ad_gwa$gene_name %in% c("APOE","LILRB2"), pal[5], pal[1]), 
    size=ifelse(ad_gwa$gene_name %in% c("APOE","LILRB2"), 1.3, 0.7)
  ) +
  geom_point(
    aes(x=gene_name, y=Plog),
    color=ifelse(ad_gwa$gene_name %in% c("APOE","LILRB2"), pal[5], pal[1]), 
    size=ifelse(ad_gwa$gene_name %in% c("APOE","LILRB2"), 5, 2)
  ) +
  coord_flip(clip = "off") +
  theme(
    legend.position="none"
  ) + 
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "gene",
       y = "-log<sub>10</sub>(p)") +
  theme_pubr() +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_line(color = pal[2]),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_markdown(color = pal[3], size = 32),
    axis.text.x = element_markdown(angle = 60, size = 20, vjust = 0.5),
    text = element_text(color = pal[3], size = 32),
    axis.line = element_line(color = pal[3]),
    axis.title.x = element_markdown(color = pal[3], size = 32),
    axis.title.y = element_markdown(color = pal[3], size = 32),
    axis.ticks = element_line(color = pal[3]),
    plot.margin = margin(1, 44, 1, 0)
  )

ggsave(filename = "~/Documents/GitHub/darwins-ark-latent-gwas/p.ad.gwa.png",
       plot = p,
       device = "png",
       bg = "transparent", 
       width = 12/1.75,
       height = 8.5/1.75,
       units = "in")
