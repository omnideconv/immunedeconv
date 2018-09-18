test_that("rescaled sample sums to 1,000,000", {
  sample = abs(rnorm(10000))
  assert("sums to milltion", round(sum(scale_to_million(sample))) == 1e6)
})

test_that("random bulk is correct", {
  gene_expression = diag(5)
  gene_names = c("marker_A", "marker_B", "marker_C", "marker_D", "marker_E")
  cell_types = c("A", "B", "C", "D", "E")
  pdata_df = data.frame(cell_type=cell_types)
  rownames(pdata_df) = cell_types
  colnames(gene_expression) = cell_types
  rownames(gene_expression) = gene_names
  eset = ExpressionSet(gene_expression, 
                       phenoData = AnnotatedDataFrame(pdata_df))
  assert(make_random_bulk(eset, cell_fractions = c(A=1), n_cells = 100)/1e6 == c(1,0,0,0,0))
  assert(make_random_bulk(eset, cell_fractions = c(A=0, B=1, C=0, D=0, E=0), n_cells = 100)/1e6 == c(0,1,0,0,0))
  assert(make_random_bulk(eset, cell_fractions = c(A=0, B=2, C=2, D=0, E=0), n_cells = 100)/1e6 == c(0,.5,.5,0,0))
  assert(make_random_bulk(eset, cell_fractions = c(A=.1, B=.2, C=.3, D=.13, E=.27), n_cells = 100)/1e6 == c(.1, .2, .3, .13, .27))
})

test_that("random bulk eset is correct", {
  gene_expression = diag(5)
  gene_names = c("marker_A", "marker_B", "marker_C", "marker_D", "marker_E")
  cell_types = c("A", "B", "C", "D", "E")
  pdata_df = data.frame(cell_type=cell_types)
  rownames(pdata_df) = cell_types
  colnames(gene_expression) = cell_types
  rownames(gene_expression) = gene_names
  eset = ExpressionSet(gene_expression, 
                       phenoData = AnnotatedDataFrame(pdata_df))
  cell_fractions = rbind(diag(5), diag(5))
  colnames(cell_fractions) = cell_types
  simulated_eset = make_bulk_eset(eset, cell_fractions, n_cells = 500)
  assert("pdata is consistent with input", all(pData(simulated_eset) == cell_fractions))
  assert("cell fractions are correct", all(t(exprs(simulated_eset))/1e6 == cell_fractions))
})