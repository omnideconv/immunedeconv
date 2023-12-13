test_mat <- read_tsv("bulk_mat.tsv") %>%
  as.data.frame() %>%
  tibble::column_to_rownames("gene_symbol")


test_that("estimate works", {
  res <- deconvolute_estimate(test_mat)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
  assert("matrix dimensions consistent", nrow(res) == 4)
})

test_that("remapping with estimate works", {
  res <- deconvolute(test_mat, "estimate")
  res <- map_result_to_celltypes(res, c(
    "stroma score", "immune score",
    "estimate score", "tumor purity"
  ))
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
  assert("matrix dimensions consistent", nrow(res) == 4)
})

test_that("timer works", {
  res <- deconvolute_timer(test_mat, indications = rep("brca", ncol(test_mat)))
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
})

test_that("timer with multiple indications", {
  res <- deconvolute_timer(test_mat, indications = c("brca", "brca", "brca", "chol", "chol", "chol", "chol", "chol"))
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
  assert("column order consistent", all(colnames(res) == colnames(test_mat)))
})

test_that("timer with multiple indications, unordered indication vector", {
  res <- deconvolute_timer(test_mat, indications = c("brca", "brca", "brca", "chol", "chol", "blca", "brca", "brca"))
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
  assert("column order consistent", all(colnames(res) == colnames(test_mat)))
})

test_that("mcp_counter works", {
  res <- deconvolute_mcp_counter(test_mat)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
})

test_that("epic works", {
  res <- deconvolute_epic(test_mat, tumor = TRUE, scale_mrna = TRUE)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
  res <- deconvolute_epic(test_mat, tumor = FALSE, scale_mrna = FALSE)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
})

test_that("quantiseq works", {
  res <- deconvolute_quantiseq(test_mat, tumor = TRUE, arrays = FALSE, scale_mrna = TRUE)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
  res <- deconvolute_quantiseq(test_mat, tumor = FALSE, arrays = TRUE, scale_mrna = FALSE)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
})

test_that("xcell works", {
  res <- deconvolute_xcell(test_mat, arrays = FALSE)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
  res <- deconvolute_xcell(test_mat, arrays = TRUE)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
})

test_that("xcell works with reduced set of expected cell types", {
  expected_cell_types <- c("T cell CD4+", "T cell CD8+", "Myeloid dendritic cell", "Macrophage M1", "Macrophage M2")
  res <- deconvolute_xcell(test_mat, arrays = FALSE, expected_cell_types)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
  res <- deconvolute_xcell(test_mat, arrays = TRUE, expected_cell_types)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
})

test_that("abis works", {
  res <- deconvolute_abis(test_mat, arrays = FALSE)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
  res <- deconvolute_abis(test_mat, arrays = TRUE)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
})

test_that("consensus_tme works", {
  res <- deconvolute_consensus_tme(test_mat, indications = rep("brca", ncol(test_mat)))
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
})


test_that("consensus_tme with multiple indications, ordered and unordered", {
  indications_1 <- c("blca", "blca", "brca", "brca", "brca", "brca", "chol", "dlbc")
  indications_2 <- c("brca", "brca", "blca", "chol", "dlbc", "brca", "brca", "blca")
  res_1 <- deconvolute_consensus_tme(test_mat, indications = indications_1)
  res_2 <- deconvolute_consensus_tme(test_mat, indications = indications_2)
  assert("matrix dimensions consistent", ncol(res_1) == ncol(test_mat))
  assert("matrix dimensions consistent", ncol(res_2) == ncol(test_mat))

  order.samples <- order(toupper(indications_2))

  assert("proper results", all(colnames(res_1)[order.samples] == colnames(res_2)))
})


test_that("consensus_tme with multiple indications, including 'Unfiltered'", {
  indications_1 <- c("blca", "blca", "brca", "brca", "chol", "chol", "Unfiltered", "Unfiltered")
  indications_2 <- c("brca", "Unfiltered", "brca", "blca", "Unfiltered", "chol", "chol", "blca")
  res_1 <- deconvolute_consensus_tme(test_mat, indications = indications_1)
  res_2 <- deconvolute_consensus_tme(test_mat, indications = indications_2)
  assert("matrix dimensions consistent", ncol(res_1) == ncol(test_mat))
  assert("matrix dimensions consistent", ncol(res_2) == ncol(test_mat))

  order.samples <- order(toupper(indications_2))

  assert("proper results", all(colnames(res_1)[order.samples] == colnames(res_2)))
})

test_that("generic deconvolution works for all methods", {
  lapply(deconvolution_methods, function(method) {
    # cibersort requires the binary path to be set, n/a in unittest.
    if (!method %in% c("cibersort", "cibersort_abs")) {
      print(paste0("method is ", method))
      res <- deconvolute(test_mat, method,
        indications = rep("brca", ncol(test_mat)),
        tumor = TRUE, arrays = FALSE, rmgenes = c("ALB", "ERBB2"),
        expected_cell_types = c("T cell CD4+", "T cell CD8+", "Macrophage", "NK cell"),
        scale_mrna = FALSE
      )
      # matrix has the 'cell type' column -> +1
      assert("matrix dimensions consistent", ncol(res) == ncol(test_mat) + 1)
      assert("cell type column exists", colnames(res)[1] == "cell_type")
      assert("sample names consistent with input", colnames(res)[-1] == colnames(test_mat))
    }
  })
})

test_that("generic deconvolution works for all methods, without extra arguments", {
  lapply(deconvolution_methods, function(method) {
    # cibersort requires the binary path to be set, n/a in unittest.
    if (!method %in% c("cibersort", "cibersort_abs")) {
      print(paste0("method is ", method))
      res <- deconvolute(test_mat, method, indications = rep("brca", ncol(test_mat)))
      # matrix has the 'cell type' column -> +1
      assert("matrix dimensions consistent", ncol(res) == ncol(test_mat) + 1)
      assert("cell type column exists", colnames(res)[1] == "cell_type")
      assert("sample names consistent with input", colnames(res)[-1] == colnames(test_mat))
    }
  })
})

test_that("cell-types are mapped correctely in a method that has different cell-types in different expression matrices", {
  res_tref <- deconvolute(test_mat, "epic", tumor = TRUE)
  res_bref <- deconvolute(test_mat, "epic", tumor = FALSE)
  assert("CAF in TRef", "Cancer associated fibroblast" %in% res_tref$cell_type)
  assert("CAF not in BRef", !("Cancer associated fibroblast" %in% res_bref$cell_type))
  assert("Neutrophil in BRef", "Neutrophil" %in% res_bref$cell_type)
  assert("Neutrophil not in TRef", !("Neutrophil" %in% res_tref$cell_type))
})

test_that("additional arguments are properly passed to original method", {
  tmp_file <- tempfile()
  res <- deconvolute_xcell(test_mat, arrays = TRUE, file.name = tmp_file)
  assert("File is not created by xcell although specified. ", file.exists(tmp_file))
})

test_that("additional (native) arguments take precedence over the corresponding immuedeconv parameters", {
  # this test is only for xCell.
  # For the other methods it's not so easy to judge by the result if it worked.
  # they use exactely the same pattern, though, so it should be ok.
  res <- deconvolute_xcell(test_mat,
    arrays = TRUE,
    expected_cell_types = c("T cell CD8+", "T cell CD4+"),
    cell.types.use = c("Basophils", "Astrocytes")
  )
  assert("native argument does not take precedent", rownames(res) == c("Basophils", "Astrocytes"))
  res <- deconvolute_xcell(test_mat,
    arrays = TRUE,
    cell.types.use = c("Basophils", "Astrocytes"),
    expected_cell_types = c("T cell CD8+", "T cell CD4+")
  )
  assert("native argument does not take precedent when in different order", rownames(res) == c("Basophils", "Astrocytes"))
})
