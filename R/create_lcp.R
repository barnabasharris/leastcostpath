#' Calculate Least-cost Path from Origin to Destinations
#'
#' Calculates the Least-cost path from origin location(s) to destination location(s). Applies Dijkstra's algorithm as implemented in the igraph / cppRouting R package.
#'
#' @param cs \code{costSurface}
#'
#' @param origin \code{sf} 'POINT' or 'MULTIPOINT', \code{SpatVector}, \code{data.frame} or \code{matrix} containing the origin coordinate(s).
#'
#' @param destination \code{sf} 'POINT' or 'MULTIPOINT', \code{SpatVector}, \code{data.frame}, \code{matrix} or \code{NULL}, containing the destination coordinates. If the object contains multiple coordinates then least-cost paths will be calculated from the origin to all destinations. If NULL then FETE paths will be calculated for all origin points.
#'
#' @param permute \code{logical} if TRUE least-cost paths are calculated from all origin locations to all destination locations. If FALSE, then a single least-cost path is calculated for each origin and location pair. Can only be FALSE if \code{origin} and \code{destination} have the same number of rows. TRUE (default)
#'
#' @param cost_distance \code{logical} if TRUE computes total accumulated cost from origin to the destinations. FALSE (default)
#'
#' @param paths_vect \code{logical} if TRUE least-cost paths are converted into line geometries and returned, unless \code{paths_vect_filename} not NULL. TRUE (default)
#'
#' @param paths_vect_filename \code{character} File path to write least-cost paths line geometries to. No line geometries will be returned. NULL (default)
#'
#' @param density_rast \code{logical} if TRUE a raster counting the number times a least-cost path touches a cell is produced and returned, unless \code{density_rast_filename} not NULL. FALSE (default)
#'
#' @param density_rast_filename \code{character} File path to write density raster to. No density raster will be returned. NULL (default)
#'
#' @param plot_paths \code{logical} if TRUE least-cost paths and origins / destinations are plotted after calculations. Can take a long time with many paths. FALSE (default)
#'
#' @author Joseph Lewis
#' @author Barney Harris
#'
#' @return \code{lcpObject} Least-cost path list object containing least-cost path line geometries, least-cost path density raster and or least-cost path distances, based on the supplied \code{costSurface}.
#'
#' @importFrom data.table := .N .GRP
#'
#' @export
#'
#' @examples
#'
#' r <- terra::rast(system.file("extdata/SICILY_1000m.tif", package="leastcostpath"))
#'
#' slope_cs <- create_slope_cs(x = r, cost_function = "tobler", neighbours = 4)
#'
#' locs <- sf::st_sf(geometry = sf::st_sfc(
#' sf::st_point(c(839769, 4199443)),
#' sf::st_point(c(1038608, 4100024)),
#' sf::st_point(c(907695, 4145478)),
#' crs = terra::crs(r)))
#'
#' lcp <- create_lcp(cs = slope_cs, origin = locs, destination = locs, permute = TRUE, plot_paths = FALSE)
#' lcp <- create_lcp(cs = slope_cs, origin = locs, destination = NULL, permute = TRUE, plot_paths = FALSE)
#' lcp <- create_lcp(cs = slope_cs, origin = locs, destination = rev(locs), permute = FALSE, plot_paths = FALSE)
create_lcp <-
  function(
    cs,
    origin,
    destination = NULL,
    permute = TRUE,
    paths_vect = TRUE,
    paths_vect_filename = NULL,
    density_rast = TRUE,
    density_rast_filename = NULL,
    cost_distance = FALSE,
    plot_paths = TRUE
  ) {
    `%!in%` <- Negate(`%in%`)

    if (all(c(paths_vect, density_rast, cost_distance) == F)) {
      message('All calculation flags set to FALSE. Nothing to do...')
      message(
        'Re-run with any or all `paths_vect = TRUE`, `density_rast = TRUE`, `cost_distance = TRUE`'
      )
      return(NULL)
    }

    if (permute == FALSE) {
      if (is.null(destination)) {
        message('No destination locations supplied. Setting `permute = TRUE`')
        permute = TRUE
      }

      if (length(nrow) != length(nrow)) {
        message(
          'Different numbers of origin and destination locations. Setting `permute = TRUE`'
        )
        permute = TRUE
      }
    }

    # add locations to list
    locations <- list(origin = origin)
    if (!is.null(destination)) {
      locations$destination <- destination
    }

    # extract raster
    cs_rast <- terra::unwrap(cs$cs_dem)

    # process locations ----
    # following routine checks for locations with shared pixels, warning where these are
    # detected and returning a vector of unique cells numbers.
    locations_names <- names(locations)
    names(locations_names) <- locations_names
    # x <- locations_names[[2]]
    if (permute == T) {
      locations_cells <-
        lapply(
          locations_names,
          FUN = \(x) {
            check_duplicate_locations(
              locations_data = locations[[x]],
              locations_name = x,
              cs_rast = cs_rast
            )
          }
        )
    }

    if (permute == FALSE) {
      # routes between origin and destination are caluclated 'as is' including
      # repeats of routes from the same origin/destination pixels.
      locations_names <- names(locations)
      names(locations_names) <- locations_names

      locations_cells <-
        lapply(locations_names, FUN = \(x) {
          if (any(class(locations[[x]]) != "SpatVector")) {
            v <- vect(locations[[x]])
          } else {
            v <- locations[[x]]
          }

          locations_cells <-
            terra::extract(cs_rast, v, cells = T)$cell
        })
    }

    if (any(locations_cells$origin %in% locations_cells$destination)) {
      message(
        'Some destination locations fall within the same pixels as origin locations...'
      )
      message(
        'No least-cost paths will returned between these origins and destinations.'
      )
    }

    # calculate least-cost paths ----
    if (any(c(paths_vect, density_rast))) {
      message(
        paste0(
          'Using cppRouting engine with ',
          RcppParallel::defaultNumThreads(),
          ' threads...'
        )
      )

      ## origin and destination provided (permute) ----
      if (length(locations_cells) == 2) {
        if (permute == TRUE) {
          routes <-
            expand.grid(
              from = locations_cells$origin,
              to = locations_cells$destination,
              stringsAsFactors = F
            )

          routes <-
            routes[routes$from != routes$to, ]

          message(
            paste0(
              'Attempting ',
              nrow(routes),
              ' least-cost paths from ',
              length(locations_cells$origin),
              ' origin pixels and ',
              length(locations_cells$destination),
              ' destination pixels...'
            )
          )

          st <- Sys.time()
          suppressWarnings({
            lcp_graph <-
              cppRouting::get_multi_paths(
                cs$cppRoutingGraph,
                from = locations_cells$origin,
                to = locations_cells$destination,
                long = T
              )
          })
          d <- Sys.time() - st

          ## cppRouting::get_* can produce duplicates!!
          message(paste0('done calculating LCPs!'))
          print(d)
          theoretical_total_paths <- nrow(routes)
        }

        ## origin and destination provided (no permute) ----
        if (permute == FALSE) {
          routes <-
            data.frame(
              from = locations_cells$origin,
              to = locations_cells$destination
            )
          routes <-
            routes[routes$from != routes$to, ]

          message(
            paste0(
              'Attempting ',
              nrow(routes),
              ' least-cost paths from ',
              length(locations_cells$origin),
              ' origin pixels and ',
              length(locations_cells$destination),
              ' destination pixels...'
            )
          )

          lcp_graph <-
            cppRouting::get_path_pair(
              cs$cppRoutingGraph,
              from = locations_cells$origin,
              to = locations_cells$destination,
              # to = rev(locations_cells$origin),
              long = T,
              algorithm = 'Dijkstra'
            )

          theoretical_total_paths <- nrow(routes)

          message('done!')
        }
      }

      ## just origin provided (permute) -----
      if (length(locations_cells) == 1) {
        routes <-
          expand.grid(
            from = locations_cells$origin,
            to = locations_cells$origin
          )
        routes <-
          routes[routes$from != routes$to, ]

        message(
          paste0(
            'Attempting ',
            nrow(routes),
            ' least-cost paths from ',
            length(locations_cells$origin),
            ' origin / destination pixels ...'
          )
        )

        # get_multi_paths almost always faster than get_pairs, so
        # easier to run and discard unwanted paths / self-paths etc.
        # Warning! Can produce duplicate paths where nodes repeated across
        # from and to !! easier to remove as list !! see deduplicateRouteList
        # function above.
        lcp_graph <-
          cppRouting::get_multi_paths(
            cs$cppRoutingGraph,
            from = locations_cells$origin,
            to = locations_cells$origin,
            long = T
          )

        theoretical_total_paths <- nrow(routes)
      }

      if (nrow(lcp_graph) == 0) {
        message('No least-cost paths could be calculated')
        message(
          'Are some locations inaccessible i.e. surrounded by all NaN cells?'
        )
        return(NULL)
      }

      message('Converting to data.table...')
      data.table::setDTthreads(
        threads = RcppParallel::defaultNumThreads()
      )

      message(paste0(
        'data.table using ',
        data.table::getDTthreads(),
        ' threads...'
      ))

      lcp_graph_dt <-
        data.table::as.data.table(lcp_graph)
      rm(lcp_graph)
      gc()
      message('Adding grouping variables...')
      lcp_graph_dt[, from_to := paste0(from, '_', to)]

      output_list <- list()

      # construct least-cost path vector lines -----
      # get dt of coords
      if (paths_vect) {
        message('building vector data from route nodes...')
        message('retrieving coordinates...')
        cell_dt <-
          as.data.frame(cs_rast, cells = T, xy = T, row.names = F)[, 1:3] |>
          dplyr::rename('node' = cell) |>
          dplyr::mutate(node = as.character(node)) |>
          data.table::as.data.table()

        lcp_graph_coords <-
          cell_dt[lcp_graph_dt, on = 'node']

        lcp_graph_coords_dist <-
          unique(lcp_graph_coords, by = 'from_to')

        # divide in geoms for terra
        result <-
          lcp_graph_coords[,
            .(object = .GRP, num_nodes = .N, part = 1, hole = 0, x, y),
            by = from_to
          ]

        # Select and convert to matrix
        result_matrix <- as.matrix(result[, .(object, part, x, y, hole)])

        message('constructing geoms...')
        lcps <-
          terra::vect(
            x = result_matrix,
            crs = terra::crs(cs_rast),
            type = 'lines',
            atts = lcp_graph_coords_dist
          )

        # TODO: explore recovery options in this scenario
        if (terra::geomtype(lcps) != 'lines') {
          message('some paths failed to resolve')
          return(NULL)
        }

        actual_total_paths <-
          nrow(lcp_graph_coords_dist)

        if (actual_total_paths < theoretical_total_paths) {
          message(
            paste0(
              theoretical_total_paths - actual_total_paths,
              ' least-cost paths could not be calculated.',
              ' Are some locations inaccessible i.e. surrounded by or near to NaN cells?'
            )
          )
        } else {
          message('All least-cost paths succesfully calculated.')
        }

        if (is.null(paths_vect_filename)) {
          output_list$paths_vect <- lcps
        } else {
          terra::writeVector(lcps, paths_vect_filename)
          output_list$paths_vect <- paths_vect_filename
        }
      }

      # calculate least-cost path density raster -----
      if (density_rast) {
        message('Calculating least-cost path density...')
        if (cs$neighbours >= 16) {
          message('Inferring cells between chosen route nodes...')
          lcp_graph_dt[,
            next_node := data.table::shift(node, -1L, fill = NA),
            by = .(from_to)
          ]
          lcp_graph_dt[, `:=`(
            node = as.integer(node),
            next_node = as.integer(next_node)
          )]
          lcp_graph_dt[, node_rn := terra::rowFromCell(cs_rast, node)]
          lcp_graph_dt[, seq_id := 1:.N, by = .(from_to)]

          lcp_graph_dt[,
            origin_offset := get_origin_offset(
              target_cell = node,
              target_cell_rn = node_rn,
              cs_ncol = terra::ncol(cs_rast),
              num_neigh_cells_from_target = floor(cs$neighDim / 2)
            )
          ]

          lcp_graph_dt[,
            `:=`(
              rel_node = node + origin_offset,
              rel_next_node = next_node + origin_offset
            )
          ]

          # unpack relative moves info
          adjacent_rel_cells <- cs$relativeMoves$adjacent_rel_cells
          possible_moves <- cs$relativeMoves$reference_tables

          # filter lcps graphs into those
          lcp_graph_dt_adjacent <- lcp_graph_dt[
            rel_next_node %in% adjacent_rel_cells,
          ]
          lcp_graph_dt_nonadjacent <- lcp_graph_dt[
            rel_next_node %!in% adjacent_rel_cells,
          ]

          nrow(lcp_graph_dt_adjacent) + nrow(lcp_graph_dt_nonadjacent) ==
            nrow(lcp_graph_dt)
          #> TRUE

          # join to reference table to retrieve interstitial cells
          missing_cells <-
            possible_moves[
              lcp_graph_dt_nonadjacent,
              on = 'rel_next_node',
              allow.cartesian = T
            ]

          # transform relative to parent cell numbers
          missing_cells[, `:=`(
            next_node_path = rel_next_node_path - origin_offset
          )]

          lcp_graph_dt_adjacent_sum <- lcp_graph_dt_adjacent[,
            .(count = .N),
            by = node
          ]
          lcp_graph_dt_nonadjacent_sum <- lcp_graph_dt_nonadjacent[,
            .(count = .N),
            by = node
          ]
          missing_cells_sum <- missing_cells[,
            .(count = .N),
            by = next_node_path
          ]
          data.table::setnames(missing_cells_sum, 'next_node_path', 'node')

          nodes_sum <- rbind(
            lcp_graph_dt_adjacent_sum,
            missing_cells_sum
          )[, .(all_sum = sum(count)), by = node]
        } else {
          nodes_sum <-
            lcp_graph_dt[, .(all_sum = .N), by = node]
        }

        # create density raster
        lcp_density <- terra::rast(cs_rast)
        terra::set.values(
          lcp_density,
          cells = terra::cells(cs_rast),
          values = 0
        )
        terra::set.values(
          lcp_density,
          cells = as.integer(nodes_sum$node),
          values = nodes_sum$all_sum
        )

        if (is.null(density_rast_filename)) {
          output_list$density_rast <- lcp_density
        } else {
          terra::writeRaster(lcp_density, density_rast_filename)
          output_list$density_rast <- density_rast_filename
        }
      }
    }

    # calculate least-cost distances -----
    if (cost_distance) {
      message(
        paste0(
          'Using cppRouting engine with ',
          RcppParallel::defaultNumThreads(),
          ' threads...'
        )
      )

      # origin and destination provided
      if (length(locations_cells) == 2) {
        if (permute == TRUE) {
          message(
            paste0(
              'Attempting distances for ',
              nrow(routes),
              ' least-cost paths from ',
              length(locations_cells$origin),
              ' origin pixels and ',
              length(locations_cells$destination),
              ' destination pixels...'
            )
          )

          lcp_distances <-
            cppRouting::get_distance_matrix(
              cs$cppRoutingGraph,
              from = locations_cells$origin,
              to = locations_cells$destination
            ) |>
            as.data.frame()

          lcp_distances$from <-
            row.names(lcp_distances)

          lcp_distances_out <-
            lcp_distances |>
            tidyr::pivot_longer(
              cols = -from,
              names_to = 'to'
            ) |>
            dplyr::mutate(
              from_to = paste0(from, '_', to)
            )

          message('done!')
        }

        if (permute == FALSE) {
          message(
            paste0(
              'Attempting distances for ',
              nrow(routes),
              ' least-cost paths from ',
              length(locations_cells$origin),
              ' origin pixels and ',
              length(locations_cells$destination),
              ' destination pixels...'
            )
          )

          lcp_distances_out <-
            cppRouting::get_distance_pair(
              cs$cppRoutingGraph,
              from = locations_cells$origin,
              to = locations_cells$destination,
              algorithm = 'Dijkstra'
            ) |>
            as.data.frame() |>
            dplyr::mutate(
              from_to = paste0(from, '_', to)
            )

          message('done!')
        }
      }

      if (length(locations_cells) == 1) {
        message(
          paste0(
            'Attempting distances for ',
            (length(locations_cells$origin)^2) - length(locations_cells$origin),
            ' least-cost paths from ',
            length(locations_cells$origin),
            ' origin / destination pixels ...'
          )
        )

        lcp_distances <-
          cppRouting::get_distance_matrix(
            cs$cppRoutingGraph,
            from = locations_cells$origin,
            to = locations_cells$origin
          ) |>
          as.data.frame()

        lcp_distances$from <-
          row.names(lcp_distances)

        lcp_distances_out <-
          lcp_distances |>
          tidyr::pivot_longer(
            cols = -from,
            names_to = 'to'
          ) |>
          dplyr::mutate(
            from_to = paste0(from, '_', to)
          )

        message('done!')
      }

      output_list$cost_distance <- lcp_distances_out
    } # / if (cost_distance)

    if (!is.function(cs$costFunction)) {
      output_list$costFunction <- cs$costFunction
    } else if (is.function(cs$costFunction)) {
      output_list$costFunction <- deparse(body(cs$costFunction)[[2]])
    }

    return(output_list)
  }
