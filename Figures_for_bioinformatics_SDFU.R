## Figures for bioinformatics submission (produced from SDFU)

## Load Packages ----------------------------------------------------------------
library(reshape2)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
source("utilities.R")

## Import and format data -------------------------------------------------------
data <- read.table("utax.tsv", sep = "\t", header = TRUE)

data <- data %>% mutate(amplicon = ifelse(amplicon == "V4V4", "V4", as.character(amplicon)))

data$nb_OTU <- as.numeric(sub("sp", "", data$nb_OTU))
if ("OTU.lost" %in% names(data)) {
  data <- data %>% mutate(OTU.lost = -OTU.lost) %>%
    melt(id.vars = c("databank", "nb_OTU", "dataset", 
                     "set_number", "amplicon", "abundance_law", "method"), 
         value.name = "divergence", 
         variable.name = "rank")
  levels(data$rank) <- c(levels(data$rank)[1:(length(levels(data$rank))-2)], "FN", "FP")
} else {
  data <- data %>% 
    melt(id.vars = c("databank", "nb_OTU", "dataset", 
                     "set_number", "amplicon", "abundance_law", "method"), 
         value.name = "divergence", 
         variable.name = "rank")
}
data.2 <- data %>% dcast(databank + nb_OTU + dataset + set_number + amplicon + abundance_law + rank ~ method, value.var = "divergence", fun.aggregate = mean) %>% 
  melt(id.vars = c("databank", "nb_OTU", "dataset", 
                   "set_number", "amplicon", "abundance_law", "rank", "frogs"), 
       value.name = "divergence", 
       variable.name = "method")


## Order methods and keep only frogs, uparse, mothur and qiime
data.otus <- filter(data, (rank %in% c("FP", "FN"))) %>% mutate(error_type = rank) %>% 
  mutate(method = factor(method, levels = c("frogs", "uparse_sop", "uparse", "mothur_sop", "mothur", "qiime_sop", "qiime"), 
                         labels = c("frogs", "uparse", "uparse (MA)", "mothur", "mothur (MA)", "qiime", "qiime (MA)")))
data.otus.2 <- filter(data.2, (rank %in% c("FP", "FN"))) %>% mutate(error_type = factor(rank)) %>% 
  mutate(method = factor(method, levels = c("uparse_sop", "uparse", "mothur_sop", "mothur", "qiime_sop", "qiime"),
                         labels = c("uparse (SOP)", "uparse (MA)", "mothur (SOP)", "mothur (MA)", "qiime (SOP)", "qiime (MA)")))
data.2 <- filter(data.2, !(rank %in% c("FP", "FN"))) %>% mutate(rank = factor(rank)) %>% 
  mutate(method = factor(method, levels = c("uparse_sop", "uparse", "mothur_sop", "mothur", "qiime_sop", "qiime"),
                         labels = c("uparse (SOP)", "uparse (MA)", "mothur (SOP)", "mothur (MA)", "qiime (SOP)", "qiime (MA)")))

## Custom palette --------------------------------------------------------------------------
manual.palette <- c("frogs"      = "olivedrab3", 
                    "tied"       = "grey60", 
                    "competitor" = "steelblue1", 
                    "mothur (MA)"     = rgb(1, 0.2, 0, alpha = 0.6, maxColorValue = 1), 
                    "mothur" = rgb(1, 0.5, 0, alpha = 0.6, maxColorValue = 1),                   
                    "uparse (MA)"     = rgb(0, 0.2, 1, alpha = 0.6, maxColorValue = 1), 
                    "uparse" = rgb(0, 0.5, 1, alpha = 0.6, maxColorValue = 1), 
                    "qiime (MA)"      = "deeppink", 
                    "qiime"  = "hotpink")

## Compute statistics -----------------------------------------------------------------------

data.wilcox.test <- data.2 %>% group_by(databank, nb_OTU, dataset, amplicon, abundance_law, rank, method) %>% 
  summarize(pval = wilcox.test(frogs, divergence, paired = TRUE)$p.value, 
            measure = wilcox.test(frogs, divergence, paired = TRUE, conf.int = TRUE)$estimate) %>% 
  mutate(best = ifelse(measure >0, "competitor", "frogs")) %>% 
  mutate(best = ifelse(pval < 0.05, best, "tied")) %>%
  mutate(best = factor(best, levels = c("competitor", "tied", "frogs")))

data.otus.wilcox.test <- data.otus.2 %>% group_by(databank, nb_OTU, dataset, amplicon, abundance_law, error_type, method) %>% 
  summarize(pval = my.paired.wilcox.test(frogs, divergence)$p.value, 
            measure = sum(sign(frogs - divergence) * rank(abs(frogs - divergence)))) %>% 
  mutate(best = ifelse(measure >0, "competitor", "frogs")) %>% 
  mutate(best = ifelse(pval < 0.05, best, "tied")) %>%
  mutate(best = factor(best, levels = c("frogs", "tied", "competitor")))

## Plot graphic ----------------------------------------------------------------------------

## Figure 3 (main text, Divergence, SOP methods, SFDU)
p <- ggplot(mapping = aes(x = rank, fill = best)) + facet_grid(abundance_law + amplicon ~ nb_OTU) + theme_bw() + 
  scale_fill_manual(values = manual.palette, name = "Best:", guide = guide_legend(reverse = TRUE)) + 
  theme(axis.text.x = element_text(angle = 90), 
        axis.text = element_text(size = 6), 
        axis.title = element_text(size = 8), 
        strip.text = element_text(size = 6), 
        plot.title = element_text(size = 10, hjust = 0.5), 
        legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 8))
plot.list <- generate_barplot(p, 
                              data = data.wilcox.test %>% 
                                filter(method %in% c("uparse (SOP)", "mothur (SOP)", "qiime (SOP)")) %>% 
                                mutate(method = droplevels(method)), 
                              vbl = "method", ylab = "Number of communities")
g <- grid_arrange_shared_legend(plot.list, ncol = 1, plot = FALSE)
## export to tiff
ggsave(plot = g, file = "Escudie_FROGS_Fig3.tiff", width = 86, height = 231, units = "mm", dpi = 350)

## Figure 4 (main text, FP/FN, SFDU)
# Keep only frogs, uparse, mothur and qiime
data.otus.shortname <- data.otus %>% filter(method %in% c("frogs", "uparse", "mothur", "qiime"))
p <- ggplot(data.otus.shortname, mapping = aes(x = method, fill = method, y = divergence, alpha = amplicon)) + 
  facet_grid(nb_OTU ~ abundance_law, scales = "free_y") + 
  scale_y_continuous(limits = c(0, NA)) + 
  scale_fill_manual(values = manual.palette, name = "Method:") + 
  scale_color_manual(values = manual.palette) + 
  scale_alpha_discrete(guide = guide_legend(override.aes = list(fill = "grey10")), range = c(0.2, 0.6), name = "Amplicon:") + 
  theme_bw() + 
  theme(axis.text = element_text(color = "black"), 
        axis.text.x = element_text(angle = 90)) + 
  ylab(NULL)
p <- p + theme(axis.text.x = element_text(size = 12), 
               axis.text.y = element_text(size = 10), 
               strip.text = element_text(size = 12), 
               legend.title = element_text(size = 8, face = "bold"), 
               legend.text = element_text(size = 8))
p.fn <- p + geom_boxplot(data = filter(data.otus.shortname, error_type == "FN")) + ggtitle("False Negative OTUs")
p.fp <- p + geom_boxplot(data = filter(data.otus.shortname, error_type == "FP")) + ggtitle("False Positive OTUs")
## Export to tiff
g <- grid_arrange_shared_legend(p.fp, p.fn, ncol = 2, plot = FALSE)
ggsave(plot = g, file = "Escudie_FROGS_Fig4.tiff", width = 178, height = 147, units = "mm", dpi = 350)

## Figure S4A (SOM, Divergence, MA methods, SFDU)
p <- ggplot(mapping = aes(x = rank, fill = best)) + facet_grid(abundance_law + amplicon ~ nb_OTU) + theme_bw() + 
  scale_fill_manual(values = manual.palette, name = "Best:", guide = guide_legend(reverse = TRUE)) + 
  theme(axis.text.x = element_text(angle = 90), 
        axis.text = element_text(size = 6), 
        axis.title = element_text(size = 8), 
        strip.text = element_text(size = 6), 
        plot.title = element_text(size = 10, hjust = 0.5), 
        legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 8))
plot.list <- generate_barplot(p, 
                              data = data.wilcox.test %>% 
                                filter(method %in% c("uparse (MA)", "mothur (MA)", "qiime (MA)")) %>% 
                                mutate(method = droplevels(method)), 
                              vbl = "method", ylab = "Number of communities")
g <- grid_arrange_shared_legend(plot.list, ncol = 2, plot = FALSE)
## export to tiff
ggsave(plot = g, file = "Escudie_FROGS_FigS4A.tiff", width = 178, height = 178, units = "mm", dpi = 350)

## Figure S4B (SOM, Divergence, MA methods, SFDU, focus on frogs and qiime)
p <- ggplot(data %>% filter(method %in% c("frogs", "qiime"), ! (rank %in% c("FP", "FN"))) %>% 
              mutate(method = ifelse(method == "qiime", "qiime (MA)", as.character(method))), 
             mapping = aes(x = rank, y = (divergence), fill = method)) + 
  facet_grid(amplicon + abundance_law ~ nb_OTU, scales = "free_y") + 
  theme_bw() + scale_y_continuous(limits = c(0, 10), breaks = c(0, 3, 6, 9)) + 
  scale_fill_manual(values = manual.palette, name = "Method:", guide = guide_legend(byrow = TRUE)) + 
  scale_color_manual(values = manual.palette, guide = "none") + 
  xlab("Number of OTUs") + ylab("Divergence") + 
  geom_boxplot(aes(group = interaction(rank, method), color = method), outlier.size = 0.8) + 
  geom_boxplot(aes(group = interaction(rank, method)), outlier.color = "transparent") + 
  theme(axis.text.x = element_text(angle = 90), 
        legend.position = "bottom", 
        legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 8))
## export to tiff
ggsave(plot = p, file = "Escudie_FROGS_FigS4B.tiff", width = 178, height = 178, units = "mm", dpi = 350)

## Figure S5 (SOM, FP/FN, SFDU, exclude mothur and qiime)
p <- ggplot(mapping = aes(x = method, fill = method, y = divergence, alpha = amplicon)) + 
  facet_grid(nb_OTU ~ abundance_law, scales = "free_y") + 
  scale_y_continuous(limits = c(0, NA)) + 
  scale_fill_manual(values = manual.palette, name = "Method:") + 
  scale_color_manual(values = manual.palette) + 
  scale_alpha_discrete(guide = guide_legend(override.aes = list(fill = "grey10")), range = c(0.2, 0.6), name = "Amplicon:") + 
  theme_bw() + 
  theme(axis.text = element_text(color = "black"), 
        axis.text.x = element_text(angle = 90), 
        legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 8)
        ) + 
  ylab(NULL)
p.fn <- p + geom_boxplot(data = data.otus %>% filter(error_type == "FN", method %in% c("frogs", "uparse")) %>% 
                           mutate(method = droplevels(method))) + ggtitle("False Negative OTUs")
p.fp <- p + geom_boxplot(data = data.otus %>% filter(error_type == "FP", method %in% c("frogs", "uparse")) %>% 
                           mutate(method = droplevels(method))) + ggtitle("False Positive OTUs")
## Export to tiff
g <- grid_arrange_shared_legend(p.fp, p.fn, ncol = 2, plot = FALSE)
ggsave(plot = g, file = "Escudie_FROGS_FigS5.tiff", width = 178, height = 127, units = "mm", dpi = 350)

## Figure S6 (SOM, FP/FN, SFDU, keep only MA variants)
p <- ggplot(mapping = aes(x = error_type, fill = best)) + 
  facet_grid(abundance_law + amplicon ~ nb_OTU) + 
  theme_bw() + 
  scale_fill_manual(values = manual.palette, name = "Best:", guide = guide_legend(reverse = TRUE)) + 
  theme(axis.text.x = element_text(angle = 90), 
        axis.text = element_text(size = 6), 
        axis.title = element_text(size = 8), 
        strip.text = element_text(size = 6), 
        plot.title = element_text(size = 10, hjust = 0.5), 
        legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 8)) + 
  xlab("Error Type")
plot.list <- generate_barplot(p, 
                              data = data.otus.wilcox.test %>% 
                                filter(method %in% c("uparse (MA)", "mothur (MA)", "qiime (MA)")) %>% 
                                mutate(method = droplevels(method)), 
                              vbl = "method", ylab = "Number of communities")
g <- grid_arrange_shared_legend(plot.list, ncol = 2, plot = FALSE)
## export to tiff
ggsave(plot = g, file = "Escudie_FROGS_FigS6.tiff", width = 178, height = 178, units = "mm", dpi = 350)
