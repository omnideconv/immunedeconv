---
title: "Preparing the single cell dataset by Schelker et al. "
output: html_notebook
---

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code. 

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*. 

```{r, eval=FALSE, include=TRUE}
sc_normalized = read_csv("../data/single_cell_benchmark_publication_data/export/sc_data_normalized.csv")

sample_id = read_csv("../data/single_cell_benchmark_publication_data/export/sample_id.csv") %>% 
  separate("Var1", into=c("source", "donor"))

tsne = read_csv("../data/single_cell_benchmark_publication_data/export/tsne_merged.csv")

cell_names = read_csv("../data/single_cell_benchmark_publication_data/export/cellnames.txt") %>%
  rename(cell_type=Var1) %>%
  mutate(cell_id=as.integer(row_number()-1))

cell_types = read_csv("../data/single_cell_benchmark_publication_data/export/classified_data_donor_abc_12k_merged_tcell_celltype.csv") %>% 
  rename(cell_id=Var1) %>%
  inner_join(cell_names)

expr = sc_normalized %>% select(-Row) %>% as.matrix()
pdata = cbind(sample_id, tsne, cell_types) %>%
  mutate(sample=colnames(expr)) %>%
  column_to_rownames(var="sample")
fdata = sc_normalized %>% select(Row) %>% rename(gene_symbol=Row)
rownames(fdata) = fdata$gene_symbol
rownames(expr) = fdata$gene_symbol

sc_eset = ExpressionSet(expr,
                        phenoData = new("AnnotatedDataFrame", data=pdata),
                        featureData = new("AnnotatedDataFrame", data=fdata))

save(sc_eset, file="/storage/home/sturm/projects/immune_deconvolution_methods/single_cell_schelker.rda", compress = TRUE)
```

Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Ctrl+Alt+I*.

When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Ctrl+Shift+K* to preview the HTML file).

The preview shows you a rendered HTML copy of the contents of the editor. Consequently, unlike *Knit*, *Preview* does not run any R code chunks. Instead, the output of the chunk when it was last run in the editor is displayed.
