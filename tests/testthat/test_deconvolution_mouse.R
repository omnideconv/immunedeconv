test_mat = read_tsv("bulk_mat_mouse.tsv") %>% as.data.frame() %>% tibble::column_to_rownames("gene_symbol")

test_that("mmcp_counter works", {
  res = deconvolute_mmcp_counter(test_mat)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
})

test_that("seqimmucc works", {
  res = deconvolute_seqimmucc(test_mat, algorithm='LLSR')
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
})

test_that("DCQ works", {
  res = deconvolute_dcq(test_mat)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
})

test_that("BASE works", {
  res = deconvolute_base_algorithm(test_mat)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
})

test_that("DCQ works without reducing cell types", {
  res = deconvolute_dcq(test_mat, combine_cells=FALSE)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
})


test_that("generic deconvolution works for all methods", {
  lapply(deconvolution_methods_mouse, function(method) {
    print(paste0("method is ", method))
    res = deconvolute_mouse(test_mat, method, algorithm='LLSR')
    # matrix has the 'cell type' column -> +1
    assert("matrix dimensions consistent", ncol(res) == ncol(test_mat) + 1)
    assert("cell type column exists", colnames(res)[1] == "cell_type")
    assert("sample names consistent with input", colnames(res)[-1] == colnames(test_mat))
  })
})


test_that("mouse gene names can be converted into their human orthologues", {
  test_mat_newGenes = mouse_genes_to_human(test_mat)
  assert("matrix dimensions consistent", ncol(test_mat_newGenes) == ncol(test_mat))
})