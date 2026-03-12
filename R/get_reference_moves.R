#' Generates reference tables of moves available and adjacent cell numbers
#' 
#' Generates tables containing all possible moves available to the least-cost path routing algorithm, including interstitial cells, from the central cell, according to the supplied neighborhood value. Move cell numbers are given in values relative to the neighbourhood window value and source \code{cs_dem}.
#' 
#' @param cs_dem \code{spatRaster} 
#' 
#' @param neighbours \code{numeric} value. Number of directions used in the conductance matrix calculation. Expected numeric values are 4, 8, 16, 32, 48. 16 (default)
#' 
#' @author Barney Harris
#' 
#' @return \code{list} containing `reference_tables` of class \code{data.table} containing all possible moves from central given as parent cell number \code{cell} and relative cell number \code{rel_next_node}, including interstitial cells \code{rel_next_node_path}, according to the cost-surface \code{neighbourhood} value and supplied \code{cs_dem} dimensions, and \code{vector} of cells immediately adajacent to the central cell. `distance_multipliers` contains the euclidean distance of each move and a multiplier normalized to the distance of a straight move to an adjacent cell. Not be called by the end user. 

get_reference_moves <- 
  function(cs_dem, neighbours, plot_move_paths=F) {
    
    if (neighbours==4) neigh_dim <- 3; dir_val <- 4
    if (neighbours==8) neigh_dim <- 3; dir_val <- 8
    if (neighbours==16) neigh_dim <- 5; dir_val <- 16
    if (neighbours==32) neigh_dim <- 7; dir_val <- neighbourhood(32)
    if (neighbours==48) neigh_dim <- 9; dir_val <- neighbourhood(48)
    
    neigh_cells_from_ref <- floor(neigh_dim/2)
    
    cs_centroid <-
      terra::centroids(terra::as.polygons(terra::ext(cs_dem)))
    
    cs_central_cell <- terra::cells(cs_dem,cs_centroid)[1,2]
    
    os <- get_origin_offset(
      target_cell = cs_central_cell,
      target_cell_rn = terra::rowFromCell(cs_dem,cs_central_cell),
      cs_ncol = terra::ncol(cs_dem),
      num_neigh_cells_from_target = neigh_cells_from_ref
    )
    
    # first get relative cell numbers for those cells immediately
    # adjacent to the central cell
    adjacent_rel_cells <- 
      (terra::adjacent(cs_dem,cs_central_cell,directions=8) + os)[1,]

    # now get relative cell numbers for those cells which represent
    # possible moves, as defined by the neighbourhood
    poss_moves <- 
      terra::adjacent(cs_dem,cs_central_cell,directions=dir_val)
    # make cell number relative using offset
    poss_moves_rel <- poss_moves + os
    # construct vector lines using cell centroids
    poss_moves_xy <- terra::xyFromCell(cs_dem,poss_moves[1,])
    
    poss_moves_xy_mat <- 
      cbind(
        poss_moves_xy,
        object = 1:nrow(poss_moves_xy),
        part = rep(1,nrow(poss_moves_xy)),
        hole = rep(0,nrow(poss_moves_xy))
      )
    
    central_xy <- 
      terra::xyFromCell(cs_dem,rep(cs_central_cell,nrow(poss_moves_xy)))
    
    central_xy_mat <- 
      cbind(
        central_xy,
        object = 1:nrow(poss_moves_xy),
        part = rep(1,nrow(poss_moves_xy)),
        hole = rep(0,nrow(poss_moves_xy))
      )
    
    all_xy <- 
      rbind(central_xy_mat,poss_moves_xy_mat)
    
    all_xy_ordered <- 
      all_xy[order(all_xy[,'object']),
             c('object','part','x','y','hole')]
    
    all_xy_ordered_df <- as.data.frame(all_xy_ordered)
    
    move_paths <- 
      terra::vect(
        all_xy_ordered,
        type='lines',
        atts = 
          data.frame(
            target_cell = poss_moves_rel[1,]
          ),
        crs = crs(cs_dem)
      )
    
    move_points <- 
      terra::vect(
        as.data.frame(all_xy_ordered) |>
          dplyr::group_by(object) |>
          slice_tail() |> 
          as.matrix(),
        type='points',
        atts = 
          data.frame(
            target_cell = poss_moves_rel[1,]
          )
        # ,crs = crs(cs_dem)
        )
    
    distance_multipliers <- 
      data.frame(
      target_cell = move_paths$target_cell,
      euclidean_distance = terra::perim(move_paths)
    ) |>
      mutate(
        distance_multiplier = 
          euclidean_distance / min(euclidean_distance)
      )
    # plot_move_paths <- T
    if (plot_move_paths) {
      plot(move_paths)
      plot(cs_dem,add=T)
      lines(move_paths)
      points(move_points)
      text(move_points,
           labels = move_points$target_cell,
           halo=T,
           hw=.05,
           hc='white'
           )
    }
    
    ex_cells <- 
      terra::cells(cs_dem, move_paths, touches = F) |>
      as.data.frame()
    
    ex_cells$rel_next_node_path <- 
      ex_cells$cell + os
    
    ex_cells$rel_next_node <- 
      move_paths[ex_cells$ID,]$target_cell
    
    ex_cells <- 
      ex_cells[ex_cells$rel_next_node_path != ex_cells$rel_next_node,]
    
    ex_cells_dt <- 
      data.table::as.data.table(ex_cells)
    
    l <- 
      list(adjacent_rel_cells = adjacent_rel_cells,
           reference_tables = ex_cells_dt,
           distance_multipliers = distance_multipliers)
    
    return(l)
  }
