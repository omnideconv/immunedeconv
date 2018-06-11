test_that("cell_type mapping to a single method_cell_type are calculated correctly", {
  fractions = list("B_cell"=.3, "T_cell.CD4"=.5, "DC"=.2)
  cell_type = "T cell CD4+"
  node = node_by_name[[cell_type]]
  result = find_children(node, fractions, "TIMER")
  assert("Correct result", result == fractions$T_cell.CD4)
})

test_that("cell_type mapping to multiple sub method_cell_type are calculated correctly", {
  fractions = list("Macrophages"=.3, "Monocytes"=.2, "NK cells"=.5)
  cell_type = "Macrophage/Monocyte"
  node = node_by_name[[cell_type]]
  result = find_children(node, fractions, "xCell")
  assert("Correct result", result == fractions$Macrophages + fractions$Monocytes)
})

test_that("cell_type mapping is NA if no corresponding method_cell_type exists", {
  fractions = list("B_cell"=.3, "T_cell.CD4"=.5, "DC"=.2)
  cell_type = "NK cell activated"
  node = node_by_name[[cell_type]]
  result = find_children(node, fractions, "TIMER")
  assert("Result is NA", is.na(result))
})

test_that("cell_type mapping is NA if one of the children is not available.", {
  # gamma delta t cells are not predicted by TIMER
  fractions = list("T_cell.CD8"=.3, "T_cell.CD4"=.5)
  cell_type = "T cell"
  node = node_by_name[[cell_type]]
  result = find_children(node, fractions, "TIMER")
  assert("Result is NA", is.na(result))
})

test_that("cell_type mapping to multiple sub method_cell_type are calculated correctly (on 2nd level)", {
  fractions = list("Monocytes"=.3, "Macrophages M0"=.1, "Macrophages M1"=.2, "Macrophages M2"=.4, "Plasma cells"=.125)
  cell_type = "Macrophage/Monocyte"
  node = node_by_name[[cell_type]]
  result = find_children(node, fractions, "CIBERSORT")
  assert("Correct result", result == fractions$`Macrophages M0` + fractions$Monocytes + fractions$`Macrophages M1` + fractions$`Macrophages M2`)
})
