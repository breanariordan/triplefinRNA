# Creating visual tables of DE genes

To create all tables used as visual aid for determining differentially expressed genes for all pairwise models, R Studio was utilised.
``` r
# Upload relevant files and checking for identical heading
annotation <- read.csv("annotationblast.csv")
DEgenes <- read.csv("DE_results_lapillum.csv")

head(anno)
head(blastx)

# Merge files and download new file
DEtable <- merge(annotation, DEgenes, by = "GeneID", all = TRUE)

write.csv(DEtable, "lbrainfullDEgnes.csv")
```
You should now have a .csv file which you can edit to display all important data for differentially expressed genes, see [Lap_DE_gene](https://github.com/breanariordan/triplefinRNA/blob/main/DE_results/lap_DE_gene.csv) and others.
