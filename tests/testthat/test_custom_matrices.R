test_mat <- read_tsv("bulk_mat.tsv") %>%
  as.data.frame() %>%
  tibble::column_to_rownames("gene_symbol")
test_mat <- as.matrix(test_mat)



test_that("seqimmucc works with a custom signature matrix", {
  sign_mat <- matrix(120 * runif(1500), ncol = 10)
  colnames(sign_mat) <- c(
    "A", "B", "C", "D",
    "E", "F", "G", "H",
    "I", "J"
  )

  rownames(sign_mat) <- sample(rownames(test_mat), nrow(sign_mat))

  res <- deconvolute_seqimmucc_custom(test_mat, sign_mat)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
  assert("matrix dimensions consistent", nrow(res) == ncol(sign_mat))
})



test_that("EPIC works with a custom signature matrix", {
  sign_mat <- matrix(120 * runif(10 * nrow(test_mat)), ncol = 10)
  colnames(sign_mat) <- c(
    "A", "B", "C", "D",
    "E", "F", "G", "H",
    "I", "J"
  )

  rownames(sign_mat) <- rownames(test_mat)
  marker_genes_epic <- sample(rownames(sign_mat), 50)


  res <- deconvolute_epic_custom(
    gene_expression_matrix = test_mat,
    signature_matrix = sign_mat,
    signature_genes = marker_genes_epic
  )
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
  assert("matrix dimensions consistent", nrow(res) == (ncol(sign_mat) + 1))
})

test_that("EPIC works with a custom signature matrix and variances", {
  sign_mat <- matrix(120 * runif(10 * nrow(test_mat)), ncol = 10)
  colnames(sign_mat) <- c(
    "A", "B", "C", "D",
    "E", "F", "G", "H",
    "I", "J"
  )

  rownames(sign_mat) <- rownames(test_mat)
  marker_genes_epic <- sample(rownames(sign_mat), 50)

  variance_genes_epic <- matrix(runif(10 * nrow(test_mat)), ncol = 10)
  dimnames(variance_genes_epic) <- list(
    rownames(sign_mat),
    colnames(sign_mat)
  )

  res <- deconvolute_epic_custom(
    gene_expression_matrix = test_mat,
    signature_matrix = sign_mat,
    signature_genes = marker_genes_epic,
    genes_var = variance_genes_epic
  )
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
  assert("matrix dimensions consistent", nrow(res) == (ncol(sign_mat) + 1))
})




test_that("ConsensusTME works with a custom signature matrix", {
  sign_mat <- matrix(120 * runif(1500), ncol = 10)
  colnames(sign_mat) <- c(
    "A", "B", "C", "D",
    "E", "F", "G", "H",
    "I", "J"
  )

  rownames(sign_mat) <- sample(rownames(test_mat), nrow(sign_mat))

  sign_genes <- list()
  for (i in 1:ncol(sign_mat)) {
    sign_genes[[i]] <- sample(rownames(sign_mat), 10)
  }
  names(sign_genes) <- colnames(sign_mat)

  res <- deconvolute_consensus_tme_custom(test_mat, sign_genes)
  assert("matrix dimensions consistent", ncol(res) == ncol(test_mat))
  assert("matrix dimensions consistent", nrow(res) == ncol(sign_mat))
})
