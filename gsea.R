rm(list=ls())
library(fgsea)
library(data.table)
library(ggplot2)
library(tidyverse)
set.seed(42)

res <- read.csv("/home/jaro/projects/transcriptomics/homework/results/res.csv")
res
geneEnsbl=gsub("\\..*","",res$X)
res$X=gsub("\\..*","",res$X)

library(org.Hs.eg.db)
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=geneEnsbl, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
ens2symbol <- as_tibble(ens2symbol)
ens2symbol
ens2symbol <- na.omit(ens2symbol)

res <- inner_join(res, ens2symbol, by=c("X"="ENSEMBL"))
res

res2 <- res %>% 
  dplyr::select(SYMBOL, log2FoldChange) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(log2FoldChange))
res2

ranks <- deframe(res2)
head(ranks, 20)
library(reactome.db)

library(org.Hs.eg.db)
entrez <- select(org.Hs.eg.db, keys=geneEnsbl, columns="ENTREZID", keytype="ENSEMBL")


#pathways <- reactomePathways(entrez$ENTREZID)
#fgseaRes <- fgsea(pathways, ranks)
#head(fgseaRes)

pathways.hallmark <- gmtPathways("/home/jaro/projects/transcriptomics/homework/pathways/c2.cp.kegg.v6.2.symbols.gmt")

# Look at them all if you want (uncomment)
 pathways.hallmark

pathways.hallmark %>% 
  head() %>% 
  lapply(head)

fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
library(DT)


# Show in a nice table:
#fgseaResTidy %>% 
#  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
#  arrange(padj) %>% 
#  DT::datatable()


ggplot(fgseaRes, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.12)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="KEGG pathways NES from GSEA") + 
  theme_minimal()

pathways.hallmark %>% 
  enframe("pathway", "SYMBOL") %>% 
  unnest() %>% 
  inner_join(res, by="SYMBOL")

plotEnrichment(pathways.hallmark[["KEGG_PATHWAYS_IN_CANCER"]],
               ranks) + labs(title="Programmed Cell Death")

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways.hallmark[topPathways], ranks, fgseaRes, 
              gseaParam=0.5)

collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.15], 
                                      pathways.hallmark, ranks)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(pathways.hallmark[mainPathways], ranks, fgseaRes, 
              gseaParam = 0.5)


resKAPA <- read.csv("/home/jaro/projects/transcriptomics/homework/results/resLFSKAPA.csv")
resKAPA
geneEnsblKAPA=gsub("\\..*","",resKAPA$X)
resKAPA$X=gsub("\\..*","",resKAPA$X)

library(org.Hs.eg.db)
ens2symbolKAPA <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=geneEnsblKAPA, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
ens2symbolKAPA <- as_tibble(ens2symbolKAPA)
ens2symbolKAPA

resKAPA <- inner_join(resKAPA, ens2symbolKAPA, by=c("X"="ENSEMBL"))
resKAPA

res2KAPA <- resKAPA %>% 
  dplyr::select(SYMBOL, log2FoldChange) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(log2FoldChange))
res2KAPA

ranksKAPA <- deframe(res2KAPA)
head(ranksKAPA, 20)

pathways.hallmark <- gmtPathways("/home/jaro/projects/transcriptomics/homework/pathways/c2.cp.kegg.v6.2.symbols.gmt")

# Look at them all if you want (uncomment)
pathways.hallmark

pathways.hallmark %>% 
  head() %>% 
  lapply(head)

fgseaResKAPA <- fgsea(pathways=pathways.hallmark, stats=ranksKAPA)

fgseaResTidyKAPA <- fgseaResKAPA %>%
  as_tibble() %>%
  arrange(desc(NES))
library(DT)


# Show in a nice table:
#fgseaResTidyKAPA %>% 
#  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
#  arrange(padj) %>% 
#  DT::datatable()


ggplot(fgseaResKAPA, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.15)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="KEGG pathways NES from GSEA") + 
  theme_minimal()

pathways.hallmark %>% 
  enframe("pathway", "SYMBOL") %>% 
  unnest() %>% 
  inner_join(res, by="SYMBOL")

plotEnrichment(pathways.hallmark[["KEGG_PATHWAYS_IN_CANCER"]],
               ranksKAPA) + labs(title="Programmed Cell Death")
