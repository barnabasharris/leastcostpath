#' Calculates an offset value from a given target cell number
#' 
#' Calculates an offset value from a given target cell number so that the target cell and its neighbouring cells can be transformed to values relative to the neighbourhood size. Not be called by the end user. 
#' 
#' @param target_cell \code{integer} Cell number from parent raster
#' 
#' @param target_cell_rn \code{integer} Row number of target cell from parent raster
#' 
#' @param cs_ncol \code{integer} Number of columns in parent raster
#' 
#' @param num_neigh_cells_from_target \code{integer} The number of cells between the target cell and the edge of the neighbourhood. For a neighbourhood with dimensions of 5*5, this would equal 2 (i.e. \code{floor(neigh_dim/2)}).
#' 
#' @author Barney Harris
#' 
#' @return \code{integer} Offset value to be summed with parent cell values to transform them into values relative to the neighbourhood size.
#' 
#' @export 

get_origin_offset <- 
  function(target_cell, target_cell_rn, cs_ncol, num_neigh_cells_from_target) {
    
    target_cell_origin <- 
      (target_cell - num_neigh_cells_from_target) - 
      (target_cell_rn - 1) * cs_ncol + 
      (target_cell_rn - (num_neigh_cells_from_target + 1) ) * cs_ncol
    
    target_cell_origin_offset <- 
      1 - target_cell_origin
    
    return(target_cell_origin_offset)
  }
