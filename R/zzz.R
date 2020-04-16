#' assign the following global variables
#' .onLoad because of 'staged-install'
#'
#' See https://developer.r-project.org/Blog/public/2019/02/14/staged-install/index.html
#' for more details.
#'
#' @importFrom utils assignInMyNamespace

.onLoad = function(libname, pkgname) {
  assignInMyNamespace('cell_type_map', .get_cell_type_map())
  assignInMyNamespace('available_datasets', .get_available_datasets())
  assignInMyNamespace('cell_type_list', .get_cell_type_list())
  assignInMyNamespace('cell_type_tree', .get_cell_type_tree())
  assignInMyNamespace('node_by_name', .get_node_by_name())
}
