# assign the following global variables
# .onLoad because of 'staged-install'
#
# See https://developer.r-project.org/Blog/public/2019/02/14/staged-install/index.html
# for more details.
#

.onLoad = function(libname, pkgname) {
  assign('cell_type_map', .get_cell_type_map(), envir=.GlobalEnv)
  assign('available_datasets', .get_available_datasets(), envir=.GlobalEnv)
  assign('cell_type_list', .get_cell_type_list(), envir=.GlobalEnv)
  assign('cell_type_tree', .get_cell_type_tree(), envir=.GlobalEnv)
  assign('node_by_name', .get_node_by_name(), envir=.GlobalEnv)

}
