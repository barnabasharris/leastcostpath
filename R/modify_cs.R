#' Adjust cost values of a pre-calculated cost surface, using a `SpatVector` object.
#' 
#' Updates nodes within the cppRouting graph, based on whether they spatially overlap with a supplied spatial object. New costs can be supplied when moving from, to and between overlapping nodes.
#' 
#' @param cs \code{conductanceMatrix} 
#' 
#' @param y \code{Spat*} object.
#' 
#' @param cost_moving_to \code{numeric}. Cost of moving from a non-spatially overlapping pixel onto a spatially overlapping pixel. Can be a function where \code{x} is the 'from' node and \code{y} is the 'to' node.
#' 
#' @param cost_moving_from \code{numeric}. Cost of moving from a spatially overlapping pixel onto a non-spatially overlapping pixel. Can be a function where \code{x} is the 'from' node and \code{y} is the 'to' node.
#' 
#' @param cost_moving_between \code{numeric}. Cost of moving between spatially overlapping pixels. Can be a function where \code{x} is the 'from' node and \code{y} is the 'to' node.
#' 
#' @author Barney Harris
#' 
#' @return returns object of class \code{conductanceMatrixProfile}
#' 
#' @export
#' 
#' @examples 
#'
#' 
#' r <- terra::rast(system.file("extdata/SICILY_1000m.tif", package="leastcostpath"))
#' 
#' slope_cs <- create_slope_cs(x = r, cost_function = "tobler", neighbours = 4)
#'  
#' cs_profile <- profile_cs(cs = slope_cs)
#' 
#' r <- terra::rast(system.file("extdata/SICILY_1000m.tif", package="leastcostpath"))
#' 
#' slope_cs <- create_slope_cs(x = r, cost_function = "tobler", neighbours = 4)
#' 
#' locs <- sf::st_sf(geometry = sf::st_sfc(
#'   sf::st_point(c(839769, 4199443)),
#'   sf::st_point(c(1038608, 4100024)),
#'   sf::st_point(c(907695, 4145478)),
#'   crs = terra::crs(r)))
#'
#' ln <- sf::st_combine(locs[1:2,]) %>% 
#'   st_cast('LINESTRING')
#' y <- vect(ln)
#' 
#' cs_mod <- modify_cs(cs = cs, y = vect(ln), cost_moving_to=200, cost_moving_from=200, cost_moving_between=200)
#' 
#' lcp <- create_lcp(cs = cs_mod, origin = locs, destination = locs, permute = TRUE, plot_paths = FALSE)

modify_cs <- 
  function(cs, y, cost_moving_to=NULL, cost_moving_from=NULL, cost_moving_between=NULL) {
    # get cells / node names from spatial object
    new_nodes <- extract(unwrap(cs$cs_dem),y,cells=T)
    
    # use cells to lookup graph node IDs
    new_nodes_with_refs <- left_join(
      data.frame(cell = as.character(new_nodes$cell)),
      cs$cppRoutingGraph$dict,
      by = c('cell' = 'ref')
    ) 
    #TODO: if (class(cost_moving_from)=='function') {}
    #TODO: if (class(y) =='SpatRaster') {}
    
    if (cs$neighbours <= 8) {
      
      if (is.numeric(cost_moving_from)) {
        cs$cppRoutingGraph$data$dist[cs$cppRoutingGraph$data$from %in% new_nodes_with_refs$id] <- 
          cost_moving_from
      }
      
      if (is.numeric(cost_moving_to)) {
        cs$cppRoutingGraph$data$dist[cs$cppRoutingGraph$data$to %in% new_nodes_with_refs$id] <- 
          cost_moving_to
      }
      
      if (is.numeric(cost_moving_from) & is.numeric(cost_moving_to)) {
        cs$cppRoutingGraph$data$dist[
          cs$cppRoutingGraph$data$to %in% new_nodes_with_refs$id & 
            cs$cppRoutingGraph$data$from %in% new_nodes_with_refs$id] <- 
          cost_moving_between
      }
      
    }
    
    if (cs$neighbours >= 16) {
      
      
      
      if (is.numeric(cost_moving_from)) {
        from_moves <- cs$cppRoutingGraph$data$dist[cs$cppRoutingGraph$data$from %in% new_nodes_with_refs$id]
        
      }
      
      if (is.numeric(cost_moving_to)) {
        cs$cppRoutingGraph$data$dist[cs$cppRoutingGraph$data$to %in% new_nodes_with_refs$id] <- 
          cost_moving_to
      }
      
      if (is.numeric(cost_moving_from) & is.numeric(cost_moving_to)) {
        cs$cppRoutingGraph$data$dist[
          cs$cppRoutingGraph$data$to %in% new_nodes_with_refs$id & 
            cs$cppRoutingGraph$data$from %in% new_nodes_with_refs$id] <- 
          cost_moving_between
      }
      
    }
    
    
    
    
    
    

    return(cs)
    
  }
