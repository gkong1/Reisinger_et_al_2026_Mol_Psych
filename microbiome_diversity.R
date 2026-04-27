# Import libraries
library(tidyverse)

# Import metadata
meta <- read.csv("sample_metadata/NF1_MB1_all_data_forGK_updated.csv")#, row.names = 1)

# Format metadata and remove sample with poor sequencing data
meta <- meta %>% 
  filter(!is.na(X16Ssampleid), X16Ssampleid != 'DMG2400369', X16Ssampleid != "") %>% 
  arrange(X16Ssampleid) %>% 
  mutate(Genotype = factor(Genotype, levels = c("WT","HET")),
         Sex = factor(Sex, levels = c("m","f")),
         sex_geno = paste0(Sex, "_", Genotype)) %>% 
  tibble::column_to_rownames(var = 'X16Ssampleid')
# meta <- meta[order(rownames(meta)),]
# meta <- meta[which(rownames(meta) != "DMG2400369"),]

#########
# Beta diversity
#########
## Import data
data <- read.delim("data/emu-combined-abundance-species.tsv", row.names = 1)

## Replace NA if any
data[is.na(data)] <- 0

## Subset and order
common <- intersect(colnames(data), rownames(meta))
data <- data[, sort(common), drop = FALSE]

colSums(data)

## Filter out taxa < 0.1% (data currently in proportions 0-1)
data.filt <- t(data[rowMeans(data) >= 0.001, ])
dim(data.filt)

## Robust CLR-transformation (for samples as rows x taxa as columns)
data.clr = vegan::decostand(data.filt, MARGIN = 1, method = 'rclr')

### Attach column and rownames
rownames(data.clr) <- rownames(data.filt)
colnames(data.clr) <- colnames(data.filt)

## PCA
pca.res <- prcomp(data.clr, scale = T, center = T)
plot(pca.res)
summ <- summary(pca.res)
comp_var <- round(summ$importance["Proportion of Variance",]*100, digits = 1)

pca.df <- pca.res$x |> as.data.frame()
pca.df <- pca.df[order(rownames(pca.df)),]
table(rownames(pca.df) == rownames(meta))

## Merge with metadata
pca.df <- cbind(pca.df, meta)|> 
  as.data.frame() %>% 
  mutate(sex_geno = factor(sex_geno, levels = c("male_WT","male_HET","female_WT","female_HET")))

## Save data
pca.data <- pca.df %>% 
  select(PC1, PC2, mouse, reads, genotype, sex, label, analysis, sex_geno)
pca_data <- list(pca_data = pca.data,
                 var_exp = comp_var)
saveRDS(pca_data, file = "analysis/beta_div/pca_data.RDS")


## Plot
p_aitch_facet <- pca.df %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(size = 2, aes(color = sex_geno, shape = sex), alpha = 0.85) + # med_hist
  scale_shape_manual(values = c(16,15,17,8)) +  
  xlab(paste0("Component1: ", comp_var[1],"%")) +
  ylab(paste0("Component2: ", comp_var[2],"%")) +
  ggtitle("Beta diversity: Aitchison") +
  facet_wrap(~sex) +
  theme_bw(base_size = 8.5) +
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        legend.position = "right",
        #legend.title = element_blank()
  )

## Save plot
jpeg("analysis/beta_div/aitchison_faceted_by_sex.jpg", res= 300, units = 'cm', 
     width = 10, height = 7)
p_aitch_facet
dev.off()

## PERMANOVA 
dist_aitch <- vegan::vegdist(data.clr, method = 'euclidean')
table(rownames(data.clr) == rownames(meta))
perm_res <- vegan::adonis2(dist_aitch ~ genotype * sex, data = meta, permutations = 999, by = 'terms')

vegan::adonis2(dist_aitch ~ genotype * sex, data = meta, permutations = 999, by = 'margin')

vegan::adonis2(dist_aitch ~ genotype, data = meta, permutations = 999, by = 'margin')
vegan::adonis2(dist_aitch ~ sex, data = meta, permutations = 999, by = 'margin')

perm_res <- as.data.frame(perm_res) %>% 
  tibble::rownames_to_column(var = "variables") %>% 
  filter(!variables %in% c("Residual", "Total")) %>% 
  mutate(sig = case_when(`Pr(>F)` < 0.01 ~ "**",
                         `Pr(>F)` < 0.05 ~ "*", 
                         .default = ""))
write.table(perm_res, file = "analysis/beta_div/aitch_permanova_res_filt_0.1.txt", sep = "\t",quote = F, row.names = F)

## Pairwise PERMANOVA
pairwise_perm <- RVAideMemoire::pairwise.perm.manova(dist_aitch, fact = meta$sex_geno)
pairwise_perm$p.value %>% 
  as.data.frame()

write.table(pairwise_perm$p.value, file = "analysis/beta_div/aitch_pairwise_permanova_res_filt_0.1.txt", sep = "\t", row.names = T, quote = F)

############
# mixOmics DIABLO integration of microbiome with behavioural data
############
## Import behavioural data
meta_behav <- read.csv("sample_metadata/NF1_MB1_all_data_forGK_updated.csv")#, row.names = 1)

## Samples wtih missing values in gut params: f HET: 5/8, f WT:1/8, m WT:2/9, m HET: 0/7 
### Pick samples with 16S sequences + selected behav columns
params <- c("BW_final", "brain", "GTT", "stool_score", "gut_permeability")
meta_behav_filt <- meta_behav_filt %>% 
  filter(!is.na(X16Ssampleid), X16Ssampleid != 'DMG2400369', X16Ssampleid != "") %>% 
  arrange(X16Ssampleid) %>% 
  dplyr::select(Sex, Genotype, ID, cage, sex_geno, X16Ssampleid, all_of(params)) %>% 
  dplyr::select(-calprotectin,-cecum_g, -colon_cm, ) %>% # 
  filter(!is.na(gut_permeability), !is.na(brain)) %>% 
  mutate(sex_geno = paste0(Sex, "_", Genotype)) %>% 
  mutate(sex_geno = factor(sex_geno, levels = c('m_WT',"m_HET","f_WT","f_HET"))) %>% 
  tibble::column_to_rownames(var = 'X16Ssampleid')

meta_behav_filt %>% count(sex_geno)

## Standardize behavioural data
meta_behav_filt_std <- vegan::decostand(meta_behav_filt[,-which(colnames(meta_behav_filt) %in% c('Sex', 'Genotype', 'ID', 'cage','sex_geno'))], MARGIN = 2, method = 'standardize') 

## Subset species data to those with metadata
data.clr.meta <- data.clr[which(rownames(data.clr) %in% rownames(meta_behav_filt_std)),]
data.clr.meta <- data.clr.meta[order(rownames(data.clr.meta)),]
dim(data.clr.meta)

meta_behav_filt_std <- meta_behav_filt_std[which(rownames(meta_behav_filt_std) %in% rownames(data.clr.meta)),]
meta_behav_filt_std <- meta_behav_filt_std[order(rownames(meta_behav_filt_std)),]
dim(meta_behav_filt_std)

## Subset to male samples 
meta_behav_filt.male <- meta_behav_filt %>% 
  filter(Sex == 'm')
table(meta_behav_filt.male$Genotype)

meta_behav_filt_std.male <- meta_behav_filt_std[which(rownames(meta_behav_filt_std) %in% rownames(meta_behav_filt.male)),]
meta_behav_filt_std.male <- meta_behav_filt_std.male[order(rownames(meta_behav_filt_std.male)),]

length(which(rownames(data.clr.meta) %in% rownames(meta_behav_filt.male))) # missing 1

rownames(meta_behav_filt.male)[which(!rownames(meta_behav_filt.male ) %in% rownames(data.clr.meta))]

data.clr.meta.male <- data.clr.meta[which(rownames(data.clr.meta) %in% rownames(meta_behav_filt.male)),]
data.clr.meta.male <- data.clr.meta.male[order(rownames(data.clr.meta.male)),]
table(rownames(meta_behav_filt_std.male) == rownames(data.clr.meta.male))

## Remove features with close to no variance
feat_var <- apply(data.clr.meta.male, MARGIN = 2, var)
feat_var <- feat_var[order(feat_var, decreasing = T)]
summary(feat_var)
feat_var <-feat_var[feat_var > 0.1]
data.clr.meta.male <- data.clr.meta.male[,which(colnames(data.clr.meta.male) %in% names(feat_var))]

## Set DIABLO data
data.diablo.male = list(microb = data.clr.meta.male,
                        meta = meta_behav_filt_std.male)

## Check that samples for all datasets match
lapply(data.diablo.male, dim)

## Set up grid for tuning
test.keepX = list (microb = c(5:9, seq(10, 18, 2), seq(20,30,5)), 
                   meta = c(2:5))

library(mixOmics)
tune.diablo = tune.block.splsda(X = data.diablo.male, 
                                Y = meta_behav_filt.male$Genotype, ncomp = 2, 
                                test.keepX = test.keepX, design = design,
                                validation = 'loo', #folds = 10, nrepeat = 1,
                                dist = "centroids.dist")
list.keepX = tune.diablo$choice.keepX # set the optimal values of features to retain
list.keepX

## Final number of features to keep
list.keepX$microb <- c(7,5)
list.keepX$meta <- c(2,1)

## Final DIABLO model
final.diablo.male = block.splsda(X = data.diablo.male, 
                                  Y = meta_behav_filt.male$Genotype, 
                                  ncomp = 2, 
                                  keepX = list.keepX, design = design)

### Sample plot
jpeg('analysis/diablo_omics_int/DIABLO_sample_plot_male.jpg', res = 300, units = 'cm', width = 19, height = 11)
plotIndiv(final.diablo.male, ind.names = FALSE, legend = TRUE, 
          group = meta_behav_filt.male$Genotype, abline = TRUE,
          title = 'DIABLO Sample Plot: Male')
dev.off()

diablo_sample_plot <- plotIndiv(final.diablo.male, ind.names = FALSE, legend = TRUE, 
                                group = meta_behav_filt.male$Genotype, abline = TRUE,
                                title = 'DIABLO Sample Plot: Male')
diablo_sample_plot_microb_coords <- diablo_sample_plot[['df']] |> as.data.frame() %>% 
  tibble::rownames_to_column(var = 'sample')
write.table(diablo_sample_plot_microb_coords, file = 'analysis/diablo_omics_int/DIABLO_sample_plot_male_coords.tsv', sep = '\t', quote = F, row.names = F)

### Variable plot
jpeg('analysis/diablo_omics_int/DIABLO_var_plot_male.jpg', res = 300, units = 'cm', width = 14, height = 14)
plotVar(final.diablo.male, var.names = c(FALSE, TRUE), 
        style = 'graphics', legend = TRUE,
        title = 'Correlation circle: Male',
        pch = c(16, 17), cex = c(1,0.8), 
        col = c('darkorchid', 'lightgreen'))
dev.off()

### Loadings 
loadings_comp1_male <- rbind(selectVar(final.diablo.male)$microb$value, selectVar(final.diablo.male)$meta$value) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = 'features') %>% 
  mutate(sex = 'm',
         block = case_when(features %in% colnames(meta_behav) ~ 'Metadata',
                           .default = 'Microbiome'))

loadings_plot_male <- loadings_comp1_male %>% # loadings_comp1_male
  # filter(features != 'Clostridium sp. SY8519') %>% 
  mutate(Outcome = case_when(value.var > 0 ~ "WT",
                             .default = 'HET'),
         block = factor(block, levels = c('Microbiome','Metadata'))) %>% 
  ggplot(aes(x = value.var, y = reorder(features, value.var, decreasing = F))) +
  geom_col(aes(fill = Outcome)) +
  labs(title = 'Comp 1 Loadings: Male', x = 'Loadings', y = 'Features') + 
  facet_wrap(~block, scales = 'free_y') +
  scale_fill_manual(values = c('#388ECC', '#F68B33')) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
  )


jpeg('analysis/diablo_omics_int/Loadings_comp1_male.jpg', res = 300, units = 'cm', width = 18, height = 8)
loadings_plot_male
dev.off()

#### Export as table
write.table(loadings_comp1_male, file = 'analysis/diablo_omics_int/Loadings_comp1_Male.tsv', sep = '\t', row.names = F, quote = F)

### Circos plot
jpeg('analysis/diablo_omics_int/DIABLO_circos_comp1_male_blank.jpg', res = 300, units = 'cm', width = 17, height = 16)
circosPlot(final.diablo.male, cutoff = 0.7, line = TRUE, comp = 1, color.Y = c("#00cc99","#0000ff"),legend.title = "", var.names = var_names_hide, 
           color.blocks= c('darkorchid', 'lightgreen'),
           color.cor = c("chocolate3","grey20"), 
           size.labels = 0, 
           size.variables = 0)
dev.off()

############
# Check if the subset of samples are representative of the whole cohort
############
meta_boxplot <- meta_behav %>% 
  rename('sample' = 'X16Ssampleid') %>% 
  pivot_longer(cols = c('BW_final', 'brain', 'GTT', 'stool_score', 'gut_permeability','colon_cm'), values_to = 'measure', names_to = 'readout') %>% # params
  mutate(sex_geno = paste0(Sex, "_", Genotype)) %>% 
  mutate(sex_geno = factor(sex_geno, levels = c('m_WT',"m_HET","f_WT","f_HET"))) %>% 
  ggplot(aes(x = sex_geno, y = measure)) +
  geom_boxplot(outlier.shape = NA, aes(fill = sex_geno)) +
  geom_jitter(width = 0.25, alpha = 0.9) +
  facet_wrap(~ readout, scales = 'free_y') +
  labs(title = 'All samples') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5))

meta_boxplot_sub <-  meta_behav %>% 
  filter(!is.na(X16Ssampleid), X16Ssampleid != 'DMG2400369', X16Ssampleid != "") %>% 
  dplyr::select(Sex, Genotype, ID, cage, X16Ssampleid, all_of(params)) %>% 
  dplyr::select(-calprotectin, -cecum_g) %>% 
  filter(!is.na(gut_permeability), !is.na(brain)) %>% 
  mutate(sex_geno = paste0(Sex, "_", Genotype)) %>% 
  mutate(sex_geno = factor(sex_geno, levels = c('m_WT',"m_HET","f_WT","f_HET"))) %>% 
  pivot_longer(cols = any_of(params), values_to = 'measure', names_to = 'readout') %>% 
  ggplot(aes(x = sex_geno, y = measure)) +
  geom_boxplot(outlier.shape = NA, aes(fill = sex_geno)) +
  geom_jitter(width = 0.25, alpha = 0.9) +
  facet_wrap(~ readout, scales = 'free_y') +
  labs(title = 'Subset of samples for DIABLO') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5))

jpeg('analysis/diablo_omics_int/meta_boxplot.jpg', res = 300, width = 18, height = 20, units = 'cm')
cowplot::plot_grid(meta_boxplot, meta_boxplot_sub, nrow = 2)
dev.off()

meta_boxplot_sub_data <-  meta_behav %>% 
  filter(!is.na(X16Ssampleid), X16Ssampleid != 'DMG2400369', X16Ssampleid != "") %>% 
  dplyr::select(Sex, Genotype, ID, cage, X16Ssampleid, all_of(params)) %>% 
  dplyr::select(-calprotectin, -cecum_g) %>% 
  filter(!is.na(gut_permeability), !is.na(brain)) %>% 
  mutate(sex_geno = paste0(Sex, "_", Genotype)) %>% 
  mutate(sex_geno = factor(sex_geno, levels = c('m_WT',"m_HET","f_WT","f_HET")))

### Export as table
write.table(meta_boxplot_sub_data, file = 'analysis/diablo_omics_int/meta_boxplot_sub_data.tsv', sep = "\t", row.names = F,quote = F)

