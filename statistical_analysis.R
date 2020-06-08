# load packages

library(tidyverse)
library(pheatmap)

library(ggrepel)
library(VennDiagram)
library(viridis)

library(edgeR)


RdBu_pal2 <- c(RColorBrewer::brewer.pal(3, "RdBu"))[-2]

# Load DATA -----------------------------------------------------------
MQ.wide <- read.csv("data/MaxQuant_proteingroups.csv")

# list of Sfps and corresponding IDs
dmon_Sfps <- read.csv("data/DmonSfps.csv")

# % Proteins without Dvir ortholog

MQ.wide %>% dplyr::select(FBgn) %>% 
  na.omit() %>% 
  nrow/MQ.wide %>% 
  dplyr::select(PROTEIN) %>% 
  na.omit() %>% nrow * 100


# Subset data - protein must be IDd by 2+ unique peptides
MQ.ion <- MQ.wide %>% 
  select(PROTEIN, FBgn, FBgn_mel, SYMBOL, Sec, Unique.peptides, starts_with("LFQ.intensity.")) %>% 
  filter(Unique.peptides >=2) %>% 
  select(-Unique.peptides)

#remove LFQ.intensity from colnames
colnames(MQ.ion) <- gsub(pattern = "LFQ.intensity.", replacement = "", x = colnames(MQ.ion))


# percent proteins found after filtering
nrow(MQ.ion)/nrow(MQ.wide) * 100


# Number of secretome/Sfps/other proteins
MQ.ion %>% 
  group_by(Sec) %>% summarise(n_distinct(PROTEIN))


# heatmap of all proteins
pheatmap(MQ.ion[, -c(1:5)], 
         scale = "row", # scale by protein
         show_colnames = T, 
         show_rownames = F, 
         color = viridis::viridis(n = 100),
         clustering_distance_cols = "euclidean", 
         clustering_method = "average") 



# Find unique proteins in each tissue/population

long.dat <- MQ.ion %>% 
  reshape2::melt(id.vars = c("PROTEIN", "FBgn", "Sec", "SYMBOL", "FBgn_mel")) %>% 
  filter(value > 0) %>% 
  mutate(pop = factor(substr(variable, 1, 1)), # Population = Colorado or Vancouver
         tissue = factor(substr(variable, 2, 3)), # Tissue = AG accessory gland
         p_t = factor(paste0(pop, tissue)), # Pop + Tissue
         BIO = factor(paste0(pop, str_extract(variable, "[0-9]+"))), # Biological replicate C1:C3 or V1:V3
         TEC = factor(if_else(grepl('_1A|2_1', variable), "A", 
                              if_else(grepl("_1B|2_2", variable), "B", 
                                      "C")))) # Technical replicate


# Table S1 ----------------------------------------------------------------
# n proteins IDd in X/5 replicates
cbind(long.dat %>% 
        filter(tissue == 'AG') %>% 
        count(pop, PROTEIN) %>% 
        count(pop, n) %>% 
        mutate(tissue = 'AG') %>% 
        pivot_wider(names_from = c(pop, tissue), 
                    values_from = nn),
      long.dat %>% 
        filter(tissue == 'EB') %>% 
        count(pop, PROTEIN) %>% 
        count(pop, n) %>% 
        mutate(tissue = 'EB') %>% 
        pivot_wider(names_from = c(pop, tissue), 
                    values_from = nn) %>% 
        select(-n))




# Figure S4 ---------------------------------------------------------------

# Overlap of proteins ID'd with known SFPs

# Proteins unique to each tissue
split.tiss <- split(long.dat$PROTEIN, long.dat$tissue, drop = T)

n_distinct(split.tiss[[1]]) # Accessory gland proteins
n_distinct(split.tiss[[2]]) # Ejacylatory bulb proteins

intersect(split.tiss[[1]], split.tiss[[2]]) %>% length/nrow(MQ.ion) * 100

setdiff(split.tiss[[1]], split.tiss[[2]]) %>% length/ nrow(MQ.ion) * 100
setdiff(split.tiss[[2]], split.tiss[[1]]) %>% length/ nrow(MQ.ion) * 100


# venn.diagram(x = list(split.tiss[[1]],
#                       split.tiss[[2]],
#                       MQ.ion %>% 
#                         filter(Sec == 'Sfp', !duplicated(PROTEIN)) %>% 
#                         pull(PROTEIN)),
#              category.names = c("Accessory gland\nproteome", "Ejaculatory bulb\nproteome", "Sfps"),
#              filename = "figs/raw/Fig_S4a.tiff",
#              compression = "lzw", output = TRUE, imagetype = "tiff",
#              height = 1700, euler.d = TRUE, width = 1700, resolution = 600,
#              lty = 'solid',
#              fill = c('turquoise', 'orange', 'red'),
#              cat.cex = 0.6, cex = .8,
#              cat.pos = c(-150, 150, 0),
#              cat.dist = c(0.05, 0.05, -0.05),
#              fontfamily = "sans", cat.fontface = "bold",
#              cat.default.pos = "outer", cat.fontfamily = "sans")


# Split proteins found in each population
split.pop <- split(long.dat$PROTEIN, long.dat$pop, drop = T)

n_distinct(split.pop[[1]]) # Colorado proteins
n_distinct(split.pop[[2]]) # Vancouver proteins

# percent proteins identified in both populations
intersect(split.pop[[1]], split.pop[[2]]) %>% length/nrow(MQ.ion) * 100

# percent proteins only identified in each population
setdiff(split.pop[[1]], split.pop[[2]]) %>% length/ nrow(MQ.ion) * 100
setdiff(split.pop[[2]], split.pop[[1]]) %>% length/ nrow(MQ.ion) * 100

# venn.diagram(x = list(split.pop[[1]],
#                       split.pop[[2]]),
#              category.names = c("Colorado", "Vancouver"),
#              filename = paste0("figs/raw/Fig_S4b.tiff"),
#              compression = "lzw", output = TRUE , imagetype = "tiff",
#              height = 1700, width = 1700, resolution = 600,
#              lty = 'solid', euler.d = TRUE, fill = RdBu_pal2,
#              cat.cex = .8, cex = .8,
#              fontfamily = "sans", cat.fontface = "bold",
#              cat.default.pos = "outer", cat.fontfamily = "sans",
#              cat.pos = c(-35, 145), cat.dist = c(0.05, 0.05))



# calculate mean abundance for each protein across all replicates/tissues

mn_dat <- data.frame(PROTEIN = MQ.ion$PROTEIN, 
                     mn = rowMeans(MQ.ion %>% 
                                     select(starts_with("C"), starts_with('V')))) %>% 
  mutate(log_mn = log2(mn + 1), # add log2 value
         tissue = factor(
           if_else(PROTEIN %in% unique(split.tiss[[1]]) == TRUE & PROTEIN %in% unique(split.tiss[[2]]) == FALSE, "AG", 
                   if_else(PROTEIN %in% unique(split.tiss[[2]]) == TRUE & PROTEIN %in% unique(split.tiss[[1]]) == FALSE, "EB", 
                           "Shared"))),
         pop.sp = factor(
           if_else(PROTEIN %in% unique(split.pop[[2]]) == FALSE, "Colorado",
                   if_else(PROTEIN %in% unique(split.pop[[1]]) == FALSE, "Vancouver", "Shared"))),
         tissue = relevel(tissue, ref = 'Shared'),
         pop.sp = relevel(pop.sp, ref = 'Shared'))


# log2 ion intensity for proteins shared/unique to tissue
ggplot(mn_dat, aes(x = tissue, y = log_mn, fill = tissue)) +
  geom_boxplot() +
  scale_fill_manual(values = c('grey', 'turquoise', 'orange')) +
  scale_x_discrete(labels = c("Shared", "Accessory glands", "Ejaculatory bulbs")) +
  labs(y = 'log2(mean ion intensity)') +
  theme_bw() +
  theme(legend.position = 'none',
        legend.title = element_blank(),       
        legend.background = element_blank(), 
        legend.spacing.x = unit(0.25, 'cm'),
        axis.title.x = element_blank(),
        strip.text = element_text(face = "bold", size = 15),
        strip.background = element_rect(fill = "grey"),
        plot.title = element_text(hjust = 0.5),
        plot.background = element_rect(colour = NA)) +
  geom_text(data = mn_dat2 %>% 
              group_by(tissue) %>% 
              summarise(N = n()),
            aes(y = 12, label = paste0("n = ", N)), 
            size = 5, colour = "black") 

# stats
mn_dat %>% 
  na.omit() %>% 
  group_by(tissue) %>% 
  summarise(mn_intensity = mean(mn),
            log2_intensity = mean(log_mn)) %>% 
  mutate(fc = mn_intensity[tissue == 'Shared']/mn_intensity,
         lfc = log2_intensity[tissue == 'Shared']/log2_intensity)


# log2 ion intensity for proteins shared/unique to population
ggplot(mn_dat, aes(x = pop.sp, y = log_mn, fill = pop.sp)) +
  geom_boxplot() +
  scale_fill_manual(values = c('grey', RdBu_pal2)) +
  labs(y = 'log2(mean ion intensity)') +
  theme_bw() +
  theme(legend.position = 'none',
        legend.title = element_blank(),       
        legend.background = element_blank(), 
        legend.spacing.x = unit(0.25, 'cm'),
        axis.title.x = element_blank(),
        strip.text = element_text(face = "bold", size = 15),
        strip.background = element_rect(fill = "grey"),
        plot.title = element_text(hjust = 0.5),
        plot.background = element_rect(colour = NA)) +
  geom_text(data = mn_dat2 %>% 
              group_by(pop.sp) %>% 
              summarise(N = n()),
            aes(y = 12, label = paste0("n = ", N)), 
            size = 5, colour = "black") 

# stats
mn_dat %>% 
  na.omit() %>% 
  group_by(pop.sp) %>% 
  summarise(mn_intensity = mean(mn),
            log2_intensity = mean(log_mn)) %>% 
  mutate(fc = mn_intensity[pop.sp == 'Shared']/mn_intensity,
         lfc = log2_intensity[pop.sp == 'Shared']/log2_intensity)


# Figure 1. MDS plot ---------------------------------------------------------------

pr_dat <- as.matrix(MQ.ion[, -c(1:5)])
rownames(pr_dat) <- MQ.ion$PROTEIN

pop <- factor(rep(c("C", "V"), each = 10))

# create DGElist
pr_dat2 <- DGEList(na.omit(pr_dat), group = pop)
dim(pr_dat2)

# run the TMM normalization
pr_dat2 <- calcNormFactors(pr_dat2)

# estiamte common dispersion
pr_dat2 <- estimateDisp(pr_dat2)

mdsObj <- plotMDS(pr_dat2, plot = F, dim.plot = c(1,3))$cmdscale.out
mdsObj <- as.data.frame(as.matrix(mdsObj))
mdsObj$replicate <- rownames(mdsObj)
colnames(mdsObj) = c("dim1", "dim2", "dim3", "replicate")

mdsObj %>% 
  mutate(pop = substr(replicate, 1, 1),
         tissue = substr(replicate, 2, 3)) %>% 
  mutate(pop = gsub("C", "Colorado", x = pop),
         pop = gsub("V", "Vancouver", x = pop),
         tissue = gsub("AG", "AcgP", x = tissue),
         tissue = gsub("EB", "EbP", x = tissue),
         condition = paste(pop, tissue, sep = "_")) %>% 
  ggplot(aes(x = dim1, y = dim2, colour = condition)) +
  geom_point(size = 5, alpha = .7) +
  labs( x = "Dimension 1", y = "Dimension 2") +
  scale_colour_brewer(palette = "RdBu") +
  theme_bw() +
  theme(legend.position = c(0.9, .9), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        legend.background = element_blank(),
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12)) +
  NULL



# Differential abundance analysis between populations -----------------------------------------

# Grouping variables
pop <- factor(rep(c("C", "V"), each = 5))
batch <- factor(c(1, 1, 2, 2, 3, 1, 1, 2, 2, 3))
samp <- factor(paste0(pop, batch))

# model matrix to test differences between populations with replicate (batch) effect
mod1 <- model.matrix(~0 + pop + batch)

# Split data by tissue

# Accessory glands only
MQ.ags <- MQ.ion[, c(1:5, grep(pattern = ".AG_.*", x = colnames(MQ.ion)))]

# keep only proteins found in all 5 replicates
MQ.ags2 <- MQ.ags[apply(MQ.ags[, 6:15]!=0, 1, all), ]

# make matrix for DGE list
AG_dat2 <- as.matrix(MQ.ags2[, -c(1:5)])
rownames(AG_dat2) <- MQ.ags2$PROTEIN

# make DGElist
acgs_dat <- DGEList(AG_dat2, group = pop)
dim(acgs_dat)

# run the TMM normalization
acgs_dat2 <- calcNormFactors(acgs_dat)

# estimate dispersion
acgs_dat2 <- estimateDisp(acgs_dat2)

# transform data for linear modelling
acgs_dat3 <- voom(acgs_dat2, mod1, plot = FALSE)

# Duplicate correlations - see 'https://support.bioconductor.org/p/77093/'
acgs_corfit <- duplicateCorrelation(acgs_dat3, mod1, block = samp)

# fit the linear models
acgs_fit <- lmFit(acgs_dat3, mod1, correlation = acgs_corfit$consensus)

# Make contrasts
acgs_contr <- makeContrasts(popC - popV, levels = colnames(coef(acgs_fit)))
tmp_acgs <- contrasts.fit(acgs_fit, acgs_contr)

# Compute empirical bayes stats
tmp_acgs <- eBayes(tmp_acgs)

# Diagnostic plots
par(mfrow=c(2,2))
# Biological coefficient of variation
plotBCV(acgs_dat2, main = "Dispersion trends")
# mean-variance trend
virginoom = voom(acgs_dat, mod1, plot = TRUE)
# QQ-plot
qqt(tmp_acgs$t,df=tmp_acgs$df.prior+tmp_acgs$df.residual,pch=16,cex=0.2)
abline(0,1)
# log2 transformed and normalize boxplot of counts across samples
boxplot(virginoom$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
abline(h=median(virginoom$E),col="blue")
par(mfrow=c(1,1))


# Results table
acg_results <- topTable(tmp_acgs, sort.by = "P", n = Inf)
acg_results$id <- rownames(acg_results)

acg_results2 <- left_join(acg_results %>% 
                            select(-AveExpr, -t, -P.Value, -B) %>% 
                            mutate(ID = rownames(acg_results)),
                          MQ.ags2 %>% 
                            select(PROTEIN, FBgn, FBgn_mel, Sec, SYMBOL),
                          by = c("ID" = "PROTEIN")) %>% 
  left_join(dmon_Sfps %>% select(PROTEIN, Gene.symbol),
            by = c("ID" = "PROTEIN")) %>%  
  mutate(threshold = factor(if_else(adj.P.Val < 0.05, "SD", "NS")),
         Acp = factor(if_else(is.na(Gene.symbol) == TRUE, "N", "Acp")),
         ID = factor(ID),
         Sec = factor(Sec))

acg_results2$Sec <- as.character(acg_results2$Sec)
acg_results2$Sec[is.na(acg_results2$Sec)] <- "All"

# Volcano plot
acg_results2 %>% 
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), colour = Sec)) +
  geom_point(alpha = .7) +
  scale_colour_manual(values = c("black", "blueviolet", "hotpink")) +
  labs(x = 'logFC', y = expression(-log[10](P))) +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        plot.background = element_rect(colour = NA)) +  
  annotate("text", x = 1.9, y = 0, size = 5, hjust = 0.5,
           label = "Colorado", colour = 'grey60') +  
  annotate("text", x = -2.5, y = 0, size = 5, hjust = 0.5,
           label = "Vancouver", colour = 'grey60') +
  geom_label_repel(
    data = acg_results2 %>% filter(threshold == 'SD' & Acp == "Acp"),
    aes(label = Gene.symbol),
    size = 4,
    colour = 'black',
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )


# heatmap of differentially abundant proteins
daa <- acg_results2 %>% 
  filter(threshold == 'SD') %>% 
  pull(id)

pheatmap(MQ.ags2 %>%
           filter(PROTEIN %in% daa) %>%
           select(-c(1:5)) %>% log2(),
         scale = "row", # scale by protein
         show_colnames = T,
         show_rownames = F,
         color = viridis::viridis(n = 150),
         border_color = NA,
         clustering_distance_cols = "euclidean",
         clustering_method = "average")




# Ejaculatory bulbs only
MQ.ebs <- MQ.ion[, c(1:5, grep(pattern = ".EB_.*", x = colnames(MQ.ion)))]

# keep only those found in all 5 replicates
MQ.ebs2 <- MQ.ebs[apply(MQ.ebs[, 6:15]!=0, 1, all), ]

# Make DGElist
EB_dat <- as.matrix(MQ.ebs2[, -c(1:5)])
rownames(EB_dat) <- MQ.ebs2$PROTEIN

dEB <- DGEList(EB_dat, group = pop)

# run the TMM normalization
dEB1 <- calcNormFactors(dEB)

dEB1 <- estimateDisp(dEB1)

yEB <- voom(dEB1, mod1, plot = FALSE)

corfitEB <- duplicateCorrelation(yEB, mod1, block = samp)

fitEB <- lmFit(yEB, mod1, correlation = corfitEB$consensus)

contrEB <- makeContrasts(popC - popV, levels = colnames(coef(fitEB)))

tmpEB <- contrasts.fit(fitEB, contrEB)

tmpEB <- eBayes(tmpEB)


# Diagnostic plots
par(mfrow=c(2,2))
# Biological coefficient of variation
plotBCV(dEB1, main = "Dispersion trends")
# mean-variance trend
virginoom = voom(dEB, mod1, plot = TRUE)
# QQ-plot
qqt(tmpEB$t,df=tmpEB$df.prior+tmpEB$df.residual,pch=16,cex=0.2)
abline(0,1)
# log2 transformed and normalize boxplot of counts across samples
boxplot(virginoom$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
abline(h=median(virginoom$E),col="blue")
par(mfrow=c(1,1))


# results table
top.tableEB <- topTable(tmpEB, sort.by = "P", n = Inf)
top.tableEB$id <- rownames(top.tableEB)

ddatEB <- left_join(top.tableEB %>% 
                      select(-AveExpr, -t, -P.Value, -B) %>% 
                      mutate(ID = rownames(top.tableEB)),
                    MQ.ebs2 %>% 
                      select(PROTEIN, FBgn, FBgn_mel, Sec),
                    by = c("ID" = "PROTEIN")) %>% 
  left_join(dmon_Sfps %>% select(PROTEIN, Gene.symbol),
            by = c("ID" = "PROTEIN")) %>%  
  mutate(threshold = factor(if_else(adj.P.Val < 0.05, "SD", "NS")),
         Acp = factor(if_else(is.na(Gene.symbol) == TRUE, "N", "Acp")),
         ID = factor(ID),
         Sec = factor(Sec))

ddatEB$Sec <- as.character(ddatEB$Sec)
ddatEB$Sec[is.na(ddatEB$Sec)] <- "All"

# Volcano plot
ddatEB %>% 
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), colour = Sec)) +
  geom_point(alpha = .7) +
  scale_colour_manual(values = c("black", "blueviolet", "hotpink")) +
  labs(x = 'logFC', y = expression(-log[10](P))) +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        plot.background = element_rect(colour = NA)) +
  annotate("text", x = 2, y = 7, size = 5, hjust = 1,
           label = "n = 929") +
  annotate("text", x = 2, y = 0, size = 5, hjust = 0.5,
           label = "Colorado", colour = 'grey60') +  
  annotate("text", x = -2.75, y = 0, size = 5, hjust = 0.5,
           label = "Vancouver", colour = 'grey60') +
  geom_label_repel(
    data = ddatEB %>% filter(threshold == 'SD' & Acp == "Acp"),
    aes(label = Gene.symbol),
    size = 4,
    colour = 'black',
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )


# heatmap for differentially abundant proteins
dae <- ddatEB %>% 
  filter(threshold == 'SD') %>% 
  pull(id)

pheatmap(MQ.ebs2 %>% 
           filter(PROTEIN %in% dae) %>% 
           select(-c(1:5)) %>% 
           log2(), 
         scale = "row", # scale by protein
         show_colnames = T, 
         show_rownames = F, 
         color = viridis::viridis(n = 100),
         border_color = NA,
         clustering_distance_cols = "euclidean", 
         clustering_method = "average") 



# Figure 2. Combined Volcano plots ----------------------------------------

rbind(ddat2 %>% 
        select(logFC, adj.P.Val, ID, Gene.symbol, Sec, threshold) %>% 
        mutate(Tissue = "a) Accessory gland proteome"), 
      ddatEB %>% 
        select(logFC, adj.P.Val, ID, Gene.symbol, Sec, threshold) %>% 
        mutate(Tissue = "b) Ejaculatory bulb proteome")) %>% 
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), colour = Sec)) +
  geom_hline(yintercept = -log10(0.05), lty = 2, colour = 'grey') +
  geom_point(alpha = .7) +
  scale_colour_manual(values = c("black", "blueviolet", "hotpink")) +
  labs(x = 'logFC', y = expression(-log[10](P))) +
  facet_wrap(~ Tissue) +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        strip.text = element_text(face = "bold", size = 15),
        plot.background = element_rect(colour = NA)) +
  annotate("text", x = 2, y = 0, size = 5, hjust = 0.5,
           label = "Colorado", colour = 'grey60') +  
  annotate("text", x = -2.75, y = 0, size = 5, hjust = 0.5,
           label = "Vancouver", colour = 'grey60') +
  geom_label_repel(
    data = rbind(ddat2 %>% 
                   select(logFC, adj.P.Val, ID, Gene.symbol, Acp, Sec, threshold) %>% 
                   mutate(Tissue = "a) Accessory gland proteome") %>% 
                   filter(threshold == 'SD' & Acp == "Acp"), 
                 ddatEB %>% 
                   select(logFC, adj.P.Val, ID, Gene.symbol, Acp, Sec, threshold) %>% 
                   mutate(Tissue = "b) Ejaculatory bulb proteome") %>% 
                   filter(threshold == 'SD' & Acp == "Acp")),
    aes(label = Gene.symbol),
    size = 4,
    colour = 'black',
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )


# Representation of secreted proteins in DA proteins -------

# Accessory gland proteome
table(ddat2$Sec, ddat2$threshold)

chisq.test(table(ddat2$Sec, ddat2$threshold))

# Ejaculatory bulb proteome

table(ddatEB$Sec, ddatEB$threshold)

chisq.test(table(ddatEB$Sec, ddatEB$threshold))


# Figure S5 ---------------------------------------------------------------
rep_dat <- rbind(ddat2 %>% 
                   select(id, Sec, threshold) %>% 
                   mutate(Tissue = "a) Accessory gland proteome"), 
                 ddatEB %>% 
                   select(id, Sec, threshold) %>% 
                   mutate(Tissue = "b) Ejaculatory bulb proteome")) %>% 
  mutate(threshold = fct_relevel(threshold, "SD", "NS")) %>% 
  mutate(Sec = fct_recode(Sec, "Secretome" = "Sec"))

rep_dat %>%
  ggplot(aes(x = Sec, fill = threshold)) +
  geom_bar(position = 'fill') +
  scale_fill_brewer(palette = "Set1", labels = c("Differentially abundant", "Not significant")) +
  facet_wrap(~Tissue) +
  theme_bw() +
  theme(legend.position = 'top',
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(face = "bold", size = 15),
        strip.background = element_rect(fill = "grey"),
        plot.background = element_rect(colour = NA)) +
  geom_text(data = rep_dat %>% 
              group_by(Tissue, Sec) %>% 
              summarise(N = n()) %>% 
              mutate(threshold = 'NS'),
            aes(y = -0.05, label = paste0("N = ", N)),
            size = 5)


# Figure S6 concordant log2fold change between populations ----------------------------

concord_dat <- inner_join(top.table,
                          top.tableEB,
                          by = 'id') %>% 
  mutate(sig.x = if_else(adj.P.Val.x < 0.05, "sig", "ns"),
         sig.y = if_else(adj.P.Val.y < 0.05, "sig", "ns"),
         s.x = sign(logFC.x),
         s.y = sign(logFC.y),
         concord = if_else(s.y == s.x, "same", "diff"),
         direction = if_else(sig.x == 'sig' & sig.y == 'sig' & concord == 'same' & logFC.y > 0, 'up.col', 
                             if_else(sig.x == 'sig' & sig.y == 'sig' & concord == 'same' & logFC.y < 0, 'up.van',
                                     if_else(sig.x == 'sig' & sig.y == 'sig' & concord == 'diff', 'discord',
                                             'ns')))) %>% 
  rename(lfc_ag = logFC.x,
         lfc_eb = logFC.y)

concord_dat %>% 
  ggplot(aes(x = lfc_ag, y = lfc_eb, colour = direction)) +
  geom_hline(yintercept = 0, lty = 2, colour = 'grey') +
  geom_vline(xintercept = 0, lty = 2, colour = 'grey') +
  geom_point(alpha = .5, size = 3) +
  scale_colour_manual(values = c('black', 'grey', RdBu_pal2[1], RdBu_pal2[2])) +
  labs(x= "Accessory gland proteome logFC", y = "Ejaculatory bulb proteome logFC") +
  theme_bw() +
  theme(legend.position = 'none') +
  annotate("text", x = 1.7, y = 2, size = 5, hjust = 1, colour = 'grey60',
           label = "Up Colorado") +
  annotate("text", x = -0.4, y = -2.5, size = 5, hjust = 1, colour = 'grey60',
           label = "Up Vancouver")



# Differential abundance analysis between tissues -----------------------------------------

# Grouping variables
tissue <- factor(rep(c("AG", "EB"), each = 5))
batch <- factor(c(1, 1, 2, 2, 3, 1, 1, 2, 2, 3))
sampleID <- factor(paste0(tissue, batch))

# model 
mod2 <- model.matrix(~0 + tissue + batch)

# subset data by population

# Colorado 
MQ.Col <- MQ.ion[, c(1:5, grep(pattern = "C.*", x = colnames(MQ.ion)))] 

# keep only those found in all 5 replicates
MQ.Col2 <- MQ.Col[apply(MQ.Col[, 6:15]!=0, 1, all), ]

# Make DGElist
Col_dat <- as.matrix(MQ.Col2[, -c(1:5)])
rownames(Col_dat) <- MQ.Col2$PROTEIN

dCol <- DGEList(Col_dat, group = tissue)
dim(dCol)

# run the TMM normalization
dCol1 <- calcNormFactors(dCol)

dCol1 <- estimateDisp(dCol1)

yCol <- voom(dCol1, mod2, plot = FALSE)

corfitCol <- duplicateCorrelation(yCol, mod2, block = sampleID)

fitCol <- lmFit(yCol, mod2, correlation = corfitCol$consensus)

contrCol <- makeContrasts(tissueAG - tissueEB, levels = colnames(coef(fitCol)))

tmpCol <- contrasts.fit(fitCol, contrCol)
tmpCol <- eBayes(tmpCol)

# Diagnostic plots
par(mfrow=c(2,2))
# Biological coefficient of variation
plotBCV(dCol1, main = "Dispersion trends")
# mean-variance trend
virginoom = voom(dCol, mod2, plot = TRUE)
# QQ-plot
qqt(tmpCol$t,df=tmpCol$df.prior+tmpCol$df.residual,pch=16,cex=0.2)
abline(0,1)
# log2 transformed and normalize boxplot of counts across samples
boxplot(virginoom$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
abline(h=median(virginoom$E),col="blue")
par(mfrow=c(1,1))


# Results table
top.tableCol <- topTable(tmpCol, sort.by = "P", n = Inf)
top.tableCol$id <- rownames(top.tableCol)

ddatCOL <- left_join(top.tableCol %>% 
                       select(-AveExpr, -t, -P.Value, -B),
                     MQ.Col %>% 
                       select(PROTEIN, FBgn, FBgn_mel, Sec, SYMBOL),
                     by = c("id" = "PROTEIN")) %>% 
  left_join(dmon_Sfps %>% select(PROTEIN, Gene.symbol),
            by = c("id" = "PROTEIN")) %>%  
  mutate(threshold = factor(if_else(adj.P.Val < 0.05, "SD", "NS")),
         Acp = factor(if_else(is.na(Gene.symbol) == TRUE, "N", "Acp")),
         id = factor(id),
         Sec = factor(Sec))

ddatCOL$Sec[is.na(ddatCOL$Sec)] <- "NS"


# Volcano plot
ddatCOL %>% 
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), colour = Sec)) +
  geom_hline(yintercept = -log10(0.05), lty = 2, colour = 'grey') +
  geom_point(alpha = .5) +
  labs(x = 'logFC', y = expression(-log[10](p-value))) +
  theme_bw() +
  theme(legend.title = element_blank(),
        plot.background = element_rect(colour = NA)) +
  scale_shape_manual(values = c(1, 16)) +
  geom_text_repel(
    data = ddatCOL %>% filter(threshold == 'SD' & Acp == "Acp"),
    aes(label = Gene.symbol),
    size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    colour = 'blue'
  )



# Vancouver
MQ.Van <- MQ.ion[, c(1:5, grep(pattern = "V.*", x = colnames(MQ.ion)))] # look at how the tech reps compare

# keep only those found in all 5 replicates
MQ.Van2 <- MQ.Van[apply(MQ.Van[, 6:15]!=0, 1, all), ]

# Make DGElist
Van_dat <- as.matrix(MQ.Van2[, -c(1:5)])
rownames(Van_dat) <- MQ.Van2$PROTEIN

dVan <- DGEList(Van_dat, group = tissue)
dim(dVan)

# run the TMM normalization
dVan1 <- calcNormFactors(dVan)
dVan1 <- estimateDisp(dVan1)

yVan <- voom(dVan1, mod2, plot = FALSE)

corfitVan <- duplicateCorrelation(yVan, mod2, block = sampleID)

fitVan <- lmFit(yVan, mod2, correlation = corfitVan$consensus)

contrVan <- makeContrasts(tissueAG - tissueEB, levels = colnames(coef(fitVan)))

tmpVan <- contrasts.fit(fitVan, contrVan)
tmpVan <- eBayes(tmpVan)

# Diagnostic plots
par(mfrow=c(2,2))
# Biological coefficient of variation
plotBCV(dVan1, main = "Dispersion trends")
# mean-variance trend
virginoom = voom(dVan, mod2, plot = TRUE)
# QQ-plot
qqt(tmpVan$t,df=tmpVan$df.prior+tmpVan$df.residual,pch=16,cex=0.2)
abline(0,1)
# log2 transformed and normalize boxplot of counts across samples
boxplot(virginoom$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
abline(h=median(virginoom$E),col="blue")
par(mfrow=c(1,1))


# Results table
top.tableVan <- topTable(tmpVan, sort.by = "P", n = Inf)
top.tableVan$id <- rownames(top.tableVan)

ddatVan <- left_join(top.tableVan %>% 
                       select(-AveExpr, -t, -P.Value, -B),
                     MQ.Col %>% 
                       select(PROTEIN, FBgn, FBgn_mel, Sec, SYMBOL),
                     by = c("id" = "PROTEIN")) %>% 
  left_join(dmon_Sfps %>% select(PROTEIN, Gene.symbol),
            by = c("id" = "PROTEIN")) %>%  
  mutate(threshold = factor(if_else(adj.P.Val < 0.05, "SD", "NS")),
         Acp = factor(if_else(is.na(Gene.symbol) == TRUE, "N", "Acp")),
         id = factor(id),
         Sec = factor(Sec))

ddatVan$Sec[is.na(ddatVan$Sec)] <- "NS"


# Volcano plot
ddatVan %>% 
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), colour = Acp)) +
  geom_hline(yintercept = -log10(0.05), lty = 2, colour = 'grey') +
  geom_point(alpha = .5) +
  scale_color_manual(values = c("black", "blue")) +
  labs(x = 'logFC', y = expression(-log[10](p-value))) +
  ggtitle("Vacouver: AGs vs. EBs") + 
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        plot.background = element_rect(colour = NA)) +
  scale_shape_manual(values = c(1, 16)) +
  geom_text_repel(
    data = ddatVan %>% filter(threshold == 'SD' & Acp == "Acp"),
    aes(label = Gene.symbol),
    size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    colour = 'blue'
  )



# Overlap of differentially abundant proteins between tissues in each population
# venn.diagram(x = list(top.tableCol$id, top.tableVan$id),
#              category.names = c("Colorado", "Vancouver"),
#              filename = "figs/raw/Fig_3_inset.tiff",
#              compression = "lzw", output = TRUE , imagetype = "tiff",
#              height = 1700, width = 1700, resolution = 600,
#              lty = 'solid', euler.d = TRUE, fill = RdBu_pal2,
#              cat.cex = .8, cex = .8,
#              fontfamily = "sans", cat.fontface = "bold",
#              cat.default.pos = "outer", cat.fontfamily = "sans",
#              cat.pos = c(-35, 145), cat.dist = c(0.05, 0.05))



# rank correlation of DA proteins in each tissue shared between pops

t_dat <- inner_join(top.tableCol %>% 
                      select(logFC, adj.P.Val, id) %>% 
                      mutate(pop = "Col",
                             s = sign(logFC)), 
                    top.tableVan %>% 
                      select(logFC, adj.P.Val, id) %>% 
                      mutate(pop = "Van",
                             s = sign(logFC)),
                    by = 'id') %>% 
  mutate(concord = if_else(s.y == s.x, "same", "diff"),
         up_in = if_else(s.x > 0, "Ag", "Eb"),
         sig = if_else(adj.P.Val.x < 0.05 & adj.P.Val.y < 0.05, "sig", "not")) %>% 
  rename(Col_logFC = logFC.x,
         Van_logFC = logFC.y) 

# Numbers of proteins in each category
t_dat %>% 
  group_by(up_in, sig, concord) %>% 
  summarise(N = n())

# Correlation of log2foldchange in abundance between populations
cor.test(t_dat$Col_logFC, t_dat$Van_logFC, method = 's')

# Figure 3a
t_dat %>% 
  ggplot(aes(x = Col_logFC, y = Van_logFC)) +
  geom_hline(yintercept = 0, lty = 2, colour = 'grey') +
  geom_vline(xintercept = 0, lty = 2, colour = 'grey') +
  geom_abline(slope = 1, lty = 2) + 
  geom_point(data = t_dat %>% filter(sig == 'not'), colour = 'grey', size = 3, alpha = .5) +
  geom_point(data = t_dat %>% filter(sig == 'sig', concord != 'same'), colour = 'black', size = 3, alpha = .5) +
  geom_point(data = t_dat %>% filter(sig == 'sig', concord == 'same'), aes(fill = up_in), size = 3, pch = 21) +
  scale_colour_manual(values = c('grey', 'black')) +
  scale_fill_manual(values = c('turquoise', 'orange')) +
  labs(x = "Colorado logFC", y = "Vancouver logFC") +
  theme_bw() +
  theme(legend.position = 'none',
        legend.text = element_text(size = 20),
        legend.key.height = unit(3, "line"),
        legend.title = element_blank(),
        legend.background = element_blank(),
        axis.title = element_text(size = 20)) +
  NULL


# Venn for figure 3b.
ag_bias <- t_dat %>% 
  filter(up_in == 'Ag', sig == 'sig', concord == 'same')

eb_bias <- t_dat %>% 
  filter(up_in == 'Eb', sig == 'sig', concord == 'same')

the_rest <- t_dat %>% 
  filter(id %in% c(ag_bias$id, eb_bias$id) == FALSE)

# venn.diagram(x = list(c(ag_bias$id, the_rest$id),
#                       c(eb_bias$id, the_rest$id)),
#              category.names = c("", ""),
#              filename = "figs/raw/Fig_3b_inset.tiff",
#              compression = "lzw", output = FALSE, imagetype = "tiff",
#              height = 1700, euler.d = TRUE, scaled = TRUE, width = 1700, resolution = 600,
#              lty = 'solid',
#              fill = c('turquoise', 'orange'),
#              cex = 1, cat.cex = 0.6,
#              cat.pos = c(-150, 25), cat.dist = c(0.05, 0.05),
#              fontfamily = "sans", cat.fontface = "bold",
#              cat.default.pos = "outer", cat.fontfamily = "sans")


# Evolutionary rates analysis ---------------------------------------------

omega_vals <- read.csv("data/omega_data.csv")

# omega for all proteins we identified by LC-MS/MS
omega_vals %>% 
  filter(PROTEIN %in% MQ.ion$PROTEIN == TRUE) %>% 
  ggplot(aes(x = dN.dS)) +
  geom_histogram()


# IDs for assigning proteins for evolutionary rates analysis
up_in_AG <- t_dat %>% filter(up_in == 'Ag', sig == 'sig', concord == 'same') %>% pull(id)
up_in_EB <- t_dat %>% filter(up_in == 'Eb', sig == 'sig', concord == 'same') %>% pull(id)


omega_dat <- omega_vals %>% 
  filter(PROTEIN %in% MQ.ion$PROTEIN == TRUE) %>% 
  mutate(Tissue = factor(if_else(PROTEIN %in% MQ.ion$PROTEIN[MQ.ion$Sec == 'Sfp'] == TRUE, "Sfps", 
                                 if_else(PROTEIN %in% MQ.ion$PROTEIN[MQ.ion$Sec == 'Sec'] == TRUE, "Secretome",
                                         if_else(PROTEIN %in% up_in_AG == TRUE, "Acc. gland proteome",
                                                 if_else(PROTEIN %in% up_in_EB == TRUE, "Ejac. bulb proteome",
                                                         "Background")))))) %>% na.omit()


omega_dat$Tissue <- ordered(omega_dat$Tissue,
                            levels = c("Background", "Acc. gland proteome", 
                                       "Ejac. bulb proteome", "Secretome", 
                                       "Sfps"))

# Figure 4.
omega_dat %>% 
  group_by(Tissue) %>% 
  summarise(mn_w = mean(dN.dS),
            se_w = sd(dN.dS)/sqrt(n()),
            N = n()) %>% 
  ggplot(aes(x = Tissue, y = mn_w, fill = Tissue)) + 
  geom_errorbar(aes(ymin = mn_w - se_w, ymax = mn_w + se_w), width = .1) + 
  geom_point(pch = 21, size = 5) +
  scale_fill_manual(values = c('grey', 'turquoise', 'orange', 'blueviolet', 'hotpink')) +
  labs(y = expression(omega)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15, colour = 'black', angle = 30, hjust = 1),
        axis.title.y = element_text(size = 30),
        strip.text = element_text(face = "bold", size = 25),
        strip.background = element_rect(fill = "grey"),
        plot.background = element_rect(colour = NA)) + 
  geom_text(aes(y = 0.03, label = paste0("n = ", N)), 
            size = 5, colour = "black") +
  NULL

# Figure S7
omega_dat %>% 
  dplyr::select(PROTEIN, Tissue, dN, dS) %>% 
  reshape2::melt() %>% 
  ggplot(aes(x = Tissue, y = value, colour = Tissue)) + 
  geom_boxplot() + 
  scale_colour_manual(values = c('grey', 'turquoise', 'orange', 'blueviolet', 'hotpink')) +
  facet_wrap(~variable, scales = 'free_y') +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15, colour = 'black', angle = 30, hjust = 1),
        axis.title.y = element_blank(),
        strip.text = element_text(face = "bold", size = 25),
        strip.background = element_rect(fill = "grey"),
        plot.background = element_rect(colour = NA)) + 
  NULL

# Kruskal-Wallis test
kruskal.test(dN.dS ~ as.factor(Tissue), data = omega_dat)

# Wilcoxon rank sum test
pairwise.wilcox.test(omega_dat$dN.dS, omega_dat$Tissue,
                     p.adjust.method = "BH")
