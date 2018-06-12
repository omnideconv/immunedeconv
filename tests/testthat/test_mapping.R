
test_that("cell_type mapping to a single method_cell_type are calculated correctly", {
  fractions = data.frame(sample1=c("B cell"=.3, "T cell CD4+"=.5, "Dendritic cell"=.2))
  cell_type = "T cell CD4+"
  node = node_by_name[[cell_type]]
  result = find_children(node, fractions, "TIMER")
  assert("Correct result", result[cell_type, "sample1"] == fractions[cell_type, "sample1"])
})

test_that("cell_type mapping to multiple sub method_cell_type are calculated correctly", {
  fractions = data.frame(sample1=c("Macrophage"=.3, "Monocyte"=.2, "NK cell"=.5))
  cell_type = "Macrophage/Monocyte"
  node = node_by_name[[cell_type]]
  result = find_children(node, fractions, "xCell")
  assert("Correct result", result[cell_type, "sample1"] == fractions["Macrophage", "sample1"] + fractions["Monocyte", "sample1"])
})

test_that("cell_type mapping is NA if no corresponding method_cell_type exists", {
  fractions = data.frame(sample1=c("B cell"=.3, "T cell CD4+"=.5, "Dendritic cell"=.2))
  cell_type = "NK cell activated"
  node = node_by_name[[cell_type]]
  result = find_children(node, fractions, "TIMER")
  assert("Result is NA", is.na(result[cell_type, "sample1"]))
})

test_that("cell_type mapping is NA if one of the children is not available.", {
  # gamma delta t cells are not predicted by TIMER
  fractions = data.frame(sample1=c("T cell CD8+"=.3, "T cell CD4+"=.5))
  cell_type = "T cell"
  node = node_by_name[[cell_type]]
  result = find_children(node, fractions, "TIMER")
  assert("Result is NA", is.na(result[cell_type, "sample1"]))
})

test_that("cell_type mapping to multiple sub method_cell_type are calculated correctly (on 2nd level)", {
  fractions = data.frame(sample1=c("Monocyte"=.3, "Macrophage M0"=.1, "Macrophage M1"=.2, "Macrophage M2"=.4, "B cell plasma"=.125))
  cell_type = "Macrophage/Monocyte"
  node = node_by_name[[cell_type]]
  result = find_children(node, fractions, "CIBERSORT")
  assert("Correct result", result[cell_type, "sample1"] == 1.)
})

test_that("cell type mapping also works with multiple columns", {
  fractions = data.frame(sample1=c("Macrophage"=.3, "Monocyte"=.2, "NK cell"=.5),
                         sample2=c("Macrophage"=.15, "Monocyte"=.7, "NK cell"=.15))
  cell_type = "Macrophage/Monocyte"
  node = node_by_name[[cell_type]]
  result = find_children(node, fractions, "xCell")
  assert("Correct result", result[cell_type, "sample1"] == .5)
  assert("Correct result", result[cell_type, "sample2"] == .85)
})

test_that("cell type mapping also works with a named vector", {
  fractions = c("Macrophage"=.3, "Monocyte"=.2, "NK cell"=.5)
  result = map_cell_types(c("Macrophage/Monocyte", "NK cell"), fractions, "xCell")
  assert(result[,1] == c(.5, .5))
})

