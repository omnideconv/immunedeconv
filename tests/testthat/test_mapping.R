test_that("cell_type mapping to a single method_cell_type are calculated correctly", {
  fractions <- data.frame(sample1 = c("B cell" = .3, "T cell CD4+" = .5, "Myeloid dendritic cell" = .2))
  cell_type <- "T cell CD4+"
  node <- node_by_name[[cell_type]]
  result <- find_children(node, fractions, "timer")
  assert("Correct result", result[cell_type, "sample1"] == fractions[cell_type, "sample1"])
})

test_that("cell_type mapping to multiple sub method_cell_type are calculated correctly", {
  fractions <- data.frame(sample1 = c("Macrophage" = .3, "Monocyte" = .2, "NK cell" = .5))
  cell_type <- "Macrophage/Monocyte"
  node <- node_by_name[[cell_type]]
  result <- find_children(node, fractions, "xcell")
  assert("Correct result", result[cell_type, "sample1"] == fractions["Macrophage", "sample1"] + fractions["Monocyte", "sample1"])
})

test_that("cell_type mapping is NA if no corresponding method_cell_type exists", {
  fractions <- data.frame(sample1 = c("B cell" = .3, "T cell CD4+" = .5, "Myeloid dendritic cell" = .2))
  cell_type <- "NK cell activated"
  node <- node_by_name[[cell_type]]
  result <- find_children(node, fractions, "timer")
  assert("Result is NA", is.na(result[cell_type, "sample1"]))
})

test_that("cell_type mapping is NA if one of the children is not available.", {
  # test_method1 does not predict M2 Macrophages. While M0 are optional, M2 are not, the result should therefore be NA
  fractions <- data.frame(sample1 = c("Macrophage M1" = .3, "Macrophage M0" = .5))
  cell_type <- "Macrophage"
  node <- node_by_name[[cell_type]]
  result <- find_children(node, fractions, "test_method1")
  assert("Result is NA", is.na(result[cell_type, "sample1"]))
})

test_that("cell_type mapping to multiple sub method_cell_type are calculated correctly (on 2nd level)", {
  fractions <- data.frame(sample1 = c("Monocyte" = .3, "Macrophage M0" = .1, "Macrophage M1" = .2, "Macrophage M2" = .4, "B cell plasma" = .125))
  cell_type <- "Macrophage/Monocyte"
  node <- node_by_name[[cell_type]]
  result <- find_children(node, fractions, "cibersort")
  assert("Correct result", result[cell_type, "sample1"] == 1.)
})

test_that("cell type mapping also works with multiple columns", {
  fractions <- data.frame(
    sample1 = c("Macrophage" = .3, "Monocyte" = .2, "NK cell" = .5),
    sample2 = c("Macrophage" = .15, "Monocyte" = .7, "NK cell" = .15)
  )
  cell_type <- "Macrophage/Monocyte"
  node <- node_by_name[[cell_type]]
  result <- find_children(node, fractions, "xcell")
  assert("Correct result", result[cell_type, "sample1"] == .5)
  assert("Correct result", result[cell_type, "sample2"] == .85)
})

test_that("cell type mapping also works with a named vector", {
  fractions <- c("Macrophage" = .3, "Monocyte" = .2, "NK cell" = .5)
  result <- map_cell_types(c("Macrophage/Monocyte", "NK cell"), fractions, "xcell")
  assert(result[, 1] == c(.5, .5))
})

test_that("children marked as 'optional' are not neccessary. ", {
  fractions <- c("Macrophage M1" = .2, "Macrophage M2" = .5)
  result <- map_cell_types(c("Macrophage"), fractions, "quantiseq")
  assert(result[, 1] == .7)
})

test_that("all children are taken into account despite being optional ", {
  fractions <- c("Macrophage M0" = .3, "Macrophage M1" = .2, "Macrophage M2" = .5)
  result <- map_cell_types(c("Macrophage"), fractions, "cibersort")
  assert(result[, 1] == 1.)
})

test_that("if no children are available, NA is returned instead of 0", {
  fractions <- c("T cell" = .5)
  result <- map_cell_types(c("Macrophage/Monocyte"), fractions, "racle")
  assert(is.na(result[, 1]))
})

test_that("get_all_children throws an error on unknown cell type", {
  err <- expect_error(get_all_children("foo cell type"))
  assert(startsWith(err$message, "unknown cell type"))
})

test_that("get_all_children works without specifying a method", {
  cell_type <- "Myeloid dendritic cell"
  assert(all(sort(get_all_children(cell_type)) ==
    sort(c("Myeloid dendritic cell resting", "Myeloid dendritic cell activated", "Myeloid dendritic cell"))))

  cell_type <- "T cell CD4+"
  assert(all(sort(get_all_children(cell_type)) ==
    sort(c(
      "T cell CD4+", "T cell CD4+ (non-regulatory)", "T cell CD4+ central memory",
      "T cell CD4+ effector memory", "T cell CD4+ memory", "T cell CD4+ memory activated",
      "T cell CD4+ memory resting", "T cell CD4+ naive", "T cell CD4+ Th1",
      "T cell CD4+ Th2", "T cell follicular helper", "T cell regulatory (Tregs)"
    ))))
})

test_that("get_all_children works with a method specified", {
  # unlike before, subsequent children should not be returned if the parent
  # is already mapped to a cell type.
  cell_type <- "Myeloid dendritic cell"
  method <- "xcell"
  assert(all(sort(get_all_children(cell_type, method)) == c("DC")))

  cell_type <- "Macrophage"
  method <- "xcell"
  assert(all(sort(get_all_children(cell_type, method)) == c("Macrophages")))

  cell_type <- "Macrophage/Monocyte"
  method <- "quantiseq"
  assert(all(sort(get_all_children(cell_type, method)) == c("Macrophages.M1", "Macrophages.M2", "Monocytes")))
})
