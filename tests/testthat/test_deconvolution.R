test_mat = read_tsv("bulk_mat.tsv") %>% as.data.frame() %>% tibble::column_to_rownames("gene_symbol")

test_that("timer works", {
  res = deconvolute_timer(test_mat, indications=rep("brca", ncol(test_mat)))
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
})

test_that("timer with multiple indications", {
  res = deconvolute_timer(test_mat, indications=c("brca", "brca", "brca", "chol", "chol", "chol", "chol", "chol"))
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
  assert("column order consistent", all(colnames(res) == colnames(test_mat)))
})

test_that("timer with multiple indications, unordered indication vector", {
  res = deconvolute_timer(test_mat, indications=c("brca", "brca", "brca", "chol", "chol", "blca", "brca", "brca"))
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
  assert("column order consistent", all(colnames(res) == colnames(test_mat)))
})

test_that("mcp_counter works", {
  res = deconvolute_mcp_counter(test_mat)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
})

test_that("epic works", {
  res = deconvolute_epic(test_mat, tumor=TRUE, scale_mrna=TRUE)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
  res = deconvolute_epic(test_mat, tumor=FALSE, scale_mrna=FALSE)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
})

test_that("quantiseq works", {
  res = deconvolute_quantiseq(test_mat, tumor=TRUE, arrays=FALSE, scale_mrna=TRUE)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
  res = deconvolute_quantiseq(test_mat, tumor=FALSE, arrays=TRUE, scale_mrna=FALSE)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
})

test_that("xcell works", {
  res = deconvolute_xcell(test_mat, arrays=FALSE)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
  res = deconvolute_xcell(test_mat, arrays=TRUE)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
})

test_that("xcell works with reduced set of expected cell types", {
  expected_cell_types = c("T cell CD4+", "T cell CD8+", "Myeloid dendritic cell", "Macrophage M1", "Macrophage M2")
  res = deconvolute_xcell(test_mat, arrays=FALSE, expected_cell_types)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
  res = deconvolute_xcell(test_mat, arrays=TRUE, expected_cell_types)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
})

test_that("generic deconvolution works for all methods", {
  lapply(deconvolution_methods, function(method) {
    # cibersort requires the binary path to be set, n/a in unittest.
    if(!method %in% c("cibersort", "cibersort_abs")) {
      print(paste0("method is ", method))
      res = deconvolute(test_mat, method, indications=rep("brca", ncol(test_mat)),
                        tumor=TRUE, arrays=FALSE, rmgenes=c("ALB", "ERBB2"),
                        expected_cell_types=c("T cell CD4+", "T cell CD8+", "Macrophage", "NK cell"),
                        scale_mrna=FALSE)
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
    if(!method %in% c("cibersort", "cibersort_abs")) {
      print(paste0("method is ", method))
      res = deconvolute(test_mat, method, indications=rep("brca", ncol(test_mat)))
      # matrix has the 'cell type' column -> +1
      assert("matrix dimensions consistent", ncol(res) == ncol(test_mat) + 1)
      assert("cell type column exists", colnames(res)[1] == "cell_type")
      assert("sample names consistent with input", colnames(res)[-1] == colnames(test_mat))
    }
  })
})
