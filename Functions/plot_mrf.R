library('mgcv')
library('ggplot2')
library('dplyr')

plot_mrf <- function(.model, .terms = NULL, .newdata, .type = 'link',
                     .fun = identity, .limits = c(-1, 1),
                     .full_model = FALSE, .rug = FALSE, ti_surfaces = FALSE,
                     .n = 100) {
  # make a custom NDVI palette that is colorblind-friendly
  .pal <- khroma::color('bukavu')(30) # wrong order
  .pal <- .pal[c(seq(1, 13, by = 2), 30, 23, 21:16)]
  .pal <- colorRampPalette(.pal)(100)
  
  # shapefile of canada
  shp <- st_read('Data/ecodistricts/Canada_Ecodistricts.shp') %>%
    st_geometry() %>%
    st_as_sf() %>%
    st_transform(crs(rast(file_names[1])))

  # find unique cells
  .newdata <- .newdata %>%
    group_by(x, y) %>%
    slice(1) %>%
    ungroup() %>%
    filter(! is.na(ecodistrict))
  
  if(.full_model) {
    # find predictions
    .newdata <- .newdata %>%
      mutate(mu_hat = predict(.model, newdata = .newdata,
                              terms = 's(ecodistrict)', type  = .type,
                              se.fit = FALSE) %>%
               # apply transformation, if necessary
               .fun())
    
    p <-
      ggplot(.newdata, aes(x, y)) +
      geom_sf(data = shp, inherit.aes = FALSE) +
      geom_raster(aes(fill = mu_hat)) +
      geom_contour(aes(z = mu_hat), color = 'black') +
      labs(title = 's(ecodistrict)', caption = 'Basis: MRF') +
      scale_fill_distiller('Partial\neffect', type = 'div', palette = 5,
                           limits = max(abs(.newdata$mu_hat)) * c(-1, 1))
    
    .smooths <- smooths(.model)
    # terms with "ecodistrict" and "ti(" or "te("
    mrf_terms <- .smooths[grepl('ecodistrict', .smooths) & grepl('t.\\(', .smooths)]
    .smooths <- .smooths[! grepl('ecodistrict', .smooths)] # drop mrf smooths
    
    .draws <- purrr::map(
      .smooths, \(.s) draw(.model, select = .s, rug = .rug, n = .n)
    )
    
    .draws[[length(.draws) + 1]] <- p
    .draws[1:length(.draws)] <-
      .draws[c(length(.draws), 1:(length(.draws) - 1))]
    
    if(length(mrf_terms) > 0) {
      ti_plots <- purrr::map(
        mrf_terms, \(.ti) {
          other_term <- gsub('t.\\(', '', .ti) %>%
            gsub('\\)', '', .) %>%
            gsub(',', '', .) %>%
            gsub('ecodistrict', '', .)
          
          if(ti_surfaces) {
          d_ti <-
            tibble(
              facets = gratia:::seq_min_max(.model$model[[other_term]],
                                            n = 4) %>%
                round(),
              nested_d = list(.newdata)) %>%
            tidyr::unnest(nested_d) %>%
            mutate(.,
                   mu_hat = predict(.model, newdata = ., terms = .ti,
                                    type = .type, se.fit = FALSE) %>%
                     .fun())
          
          p_ti <-
            ggplot(d_ti, aes(x, y)) +
            facet_wrap(~ facets) +
            geom_sf(data = shp, inherit.aes = FALSE) +
            geom_raster(aes(fill = mu_hat), na.rm = TRUE) +
            geom_contour(aes(z = mu_hat), color = 'black', na.rm = TRUE) +
            labs(x = NULL, y = NULL, title = .ti,
                 caption = 'Basis: Tensor product int.') +
            scale_fill_distiller('Partial\neffect', type = 'div', palette = 5,
                                 limits = max(abs(d_ti$mu_hat)) * c(-1, 1))
          } else {
            d_ti <-
              tibble(
                x1 = gratia:::seq_min_max(.model$model[[other_term]],
                                              n = .n) %>%
                  round(),
                nested_d = list(.newdata)) %>%
              tidyr::unnest(nested_d) %>%
              mutate(.,
                     mu_hat = predict(.model, newdata = ., terms = .ti,
                                      type = .type, se.fit = FALSE) %>%
                       .fun())
            
            p_ti <-
              ggplot(d_ti) +
              geom_line(aes(x1, mu_hat, color = ecodistrict)) +
              labs(x = other_term, y = 'Partial effect', title = .ti,
                   caption = 'Basis: Tensor product int.') +
              scale_color_manual(values = color('smoothrainbow')(n_distinct(.newdata$ecodistrict))) +
              theme(legend.position = 'none')
          }
          return(p_ti)
        } #' close function for `purrr::map()`
      )
      
      .draws[[length(.draws) + (1:length(ti_plots))]] <-
        ti_plots[[1:length(ti_plots)]]
    } #' close `if(length(mrf_terms) > 0)`
    
    plots <- plot_grid(plotlist = .draws)
    
    return(plots)
  } else {
    # find predictions
    .newdata <- mutate(.newdata,
                       mu_hat = predict(.model, newdata = .newdata,
                                        terms = .terms, type = .type,
                                        se.fit = FALSE) %>%
                         # apply transformation
                         .fun())
    
    p <-
      ggplot(.newdata, aes(x, y)) +
      geom_sf(data = shp, inherit.aes = FALSE) +
      geom_raster(aes(fill = mu_hat)) +
      geom_contour(aes(z = mu_hat), color = 'black') +
      labs(x = NULL, y = NULL, title = 's(ecodistrict)', caption = 'Basis: MRF') +
      scale_fill_gradientn(expression(bold(widehat(NDVI))), colours = .pal,
                           limits = .limits)
    return(p)
  }
}
