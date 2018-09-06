
convert_fragment_list <- function(fragment_list, to = c("cuts","footprints")) {
  out_list <- list()
  if(to == "cuts") {
    for(i in 1:length(fragment_list)) {
      prime5 <- fragment_list[[i]]
      end(prime5) <- start(prime5)
      prime3 <- fragment_list[[i]]
      start(prime3) <- end(prime3)
      out_list[[i]] <- sort(c(prime5,prime3))
      names(out_list)[i] <- names(fragment_list)[i]
    }
  } else if(to == "footprints") {
    for(i in 1:length(fragment_list)) {
      prime5 <- fragment_list[[i]]
      end(prime5) <- start(prime5) + 19
      start(prime5) <- start(prime5) -10
      prime3 <- fragment_list[[i]]
      start(prime3) <- end(prime3) - 18
      end(prime3) <- end(prime3) + 10
      both <- c(prime5,prime3)
      start(both)[start(both) < 1] <- 1
      out_list[[i]] <- sort(both)
      names(out_list)[i] <- names(fragment_list)[i]
    }
  }
  out_list
}

# Pileup reads, fragments, or cuts from a BAM file over target GRanges regions.
pileup_gr_list <- function(gr_list,
                           gr_target,
                           gr_groups = NULL,
                           norm="PM",
                           window_size = NULL,
                           window_mode = c("max","mean","median")) {

  if(is.null(gr_groups)) {
    gr_groups <- 1
    groups <- 1
  } else {
    groups <- unique(gr_groups)

  }

  out_list <- list()

  for(i in 1:length(groups)) {

    group <- groups[i]
    pile <- data.frame(pos = start(gr_target):end(gr_target), val = 0)
    group_gr <- gr_list[gr_groups == group]

    for(j in 1:length(group_gr)) {

      frags <- group_gr[[j]]

      ol <- suppressWarnings(as.data.frame(findOverlaps(gr_target,frags)))
      names(ol) <- c("target_hit","frags_hit")

      if(nrow(ol) > 0) {

        target_start <- start(gr_target)
        target_end <- end(gr_target)
        target_strand <- as.character(strand(gr_target))
        target_width <- width(gr_target)
        frags_start <- start(frags)[ol$frags_hit]
        frags_end <- end(frags)[ol$frags_hit]
        frags_width <- width(frags)[ol$frags_hit]

        for(k in 1:nrow(ol)) {

          if(target_strand %in% c("+","*")) {
            hit_start <- frags_start[k] - target_start + 1
          } else if (target_strand == "-") {
            hit_start <- target_end - frags_end[k] + 1
          }
          hit_end <- hit_start + frags_width[k] - 1

          if(hit_start < 1) { hit_start <- 1 }
          if(hit_end > target_width) { hit_end <- target_width }

          if( hit_start <= target_width & hit_end >= 1) {
            pile$val[hit_start:hit_end] <- pile$val[hit_start:hit_end] + 1
          }

        }
      }

    }

    if(!is.null(window_size)) {
      pile <- pile %>%
        mutate(bin = floor(pos/window_size))
      if(window_mode == "max") {
        pile <- pile %>%
          group_by(bin) %>%
          summarise(pos = bin[1]*window_size,
                    val = max(val))
      } else if(window_mode == "mean") {
        pile <- pile %>%
          group_by(bin) %>%
          summarise(pos = bin[1]*window_size,
                    val = mean(val))
      } else if(window_mode == "max") {
        pile <- pile %>%
          group_by(bin) %>%
          summarise(pos = bin[1]*window_size,
                    val = median(val))
      }

      pile <- pile %>%
        select(pos, val)
    }

    if(norm == "PM") {
      m <- sum(unlist(lapply(group_gr,length)))
      pile$val <- pile$val/m*1e6
    }

    out_list[[i]] <- pile
    names(out_list)[i] <- groups[i]

  }


  if(length(groups) == 1) {
    out_list <- out_list[[1]]
  }

  out_list

}

build_pile_plot <- function(gr_list,
                            ucsc_loc,
                            highlight_loc = NULL,
                            padding = c(1e5,1e5),
                            gr_groups = NULL,
                            group_colors = NULL, #named vector
                            norm = "PM",
                            max_val = NULL,
                            window_size = NULL,
                            window_mode = c("max","mean","median"),
                            target_color = "#B7B7B7",
                            highlight_color = "#F9ED32") {

  gr_target <- ucsc_loc_to_gr(ucsc_loc)
  target_start <- start(gr_target)
  target_end <- end(gr_target)

  start(gr_target) <- start(gr_target) - padding[1]
  end(gr_target) <- end(gr_target) + padding[2]

  piles <- pileup_gr_list(gr_list,
                          gr_target,
                          gr_groups,
                          norm,
                          window_size,
                          window_mode)

  target_rect <- data.frame(xmin = target_start,
                            xmax = target_end,
                            ymin = 1,
                            ymax = length(piles) + 1,
                            fill = target_color)

  chr <- as.character(seqnames(gr_target))


  pile_plot <- ggplot() +
    theme_classic() +
    scale_x_continuous(chr, expand = c(0,0)) +
    scale_y_continuous("", expand = c(0,0)) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank()) +
    scale_color_identity() +
    scale_fill_identity() +
    geom_rect(data = target_rect,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = fill))

  if(!is.null(highlight_loc)) {
    hi_target <- ucsc_loc_to_gr(highlight_loc)
    hi_start <- start(hi_target)
    hi_end <- end(hi_target)

    hi_rect <- data.frame(xmin = hi_start,
                          xmax = hi_end,
                          ymin = 1,
                          ymax = length(piles) + 1,
                          fill = highlight_color)

    pile_plot <- pile_plot +
      geom_rect(data = hi_rect,
                aes(xmin = xmin, xmax = xmax,
                    ymin = ymin, ymax = ymax,
                    fill = fill))

  }

  if(is.null(max_val)) {
    max_val <- max(unlist(lapply(piles, function(x) max(x$val))))
  }

  for(i in 1:length(piles)) {
    pile <- piles[[i]]
    pile_color <- group_colors[names(group_colors) == names(piles)[i]]

    pile$val <- i + pile$val / max_val
    pile$min <- i

    pile$val[pile$val > i + 1] <- i + 1

    baseline <- data.frame(x = min(pile$pos),
                           xend = max(pile$pos),
                           y = i,
                           yend = i,
                           color = pile_color)

    pile_plot <- pile_plot +
      geom_segment(data = baseline,
                   aes(x = x, xend = xend,
                       y = y, yend = y,
                       color = color),
                   size = 0.1) +
      geom_ribbon(data = pile,
                  aes(x = pos, ymin = min, ymax = val),
                  color = NA,
                  fill = pile_color)
  }

  pile_plot
}

build_pile_heatmap <- function(gr_list,
                               ucsc_loc,
                               highlight_loc = NULL,
                               padding = c(1e5,1e5),
                               gr_groups = NULL,
                               colorset = c("white","black"),
                               norm = "PM",
                               max_val = NULL,
                               window_size = NULL,
                               window_mode = c("max","mean","median"),
                               target_color = "#CBF8FF",
                               highlight_color = "#F9ED32",
                               baselines = TRUE) {

  gr_target <- ucsc_loc_to_gr(ucsc_loc)
  target_start <- start(gr_target)
  target_end <- end(gr_target)

  start(gr_target) <- start(gr_target) - padding[1]
  end(gr_target) <- end(gr_target) + padding[2]

  piles <- pileup_gr_list(gr_list,
                          gr_target,
                          gr_groups,
                          norm,
                          window_size,
                          window_mode)

  target_rect <- data.frame(xmin = target_start,
                            xmax = target_end,
                            ymin = 1,
                            ymax = length(piles) + 1,
                            fill = target_color)

  chr <- as.character(seqnames(gr_target))


  pile_plot <- ggplot() +
    theme_classic() +
    scale_x_continuous(chr, expand = c(0,0)) +
    scale_y_continuous("", expand = c(0,0)) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank()) +
    scale_color_identity() +
    scale_fill_identity() +
    geom_rect(data = target_rect,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = fill))

  if(!is.null(highlight_loc)) {
    hi_target <- ucsc_loc_to_gr(highlight_loc)
    hi_start <- start(hi_target)
    hi_end <- end(hi_target)

    hi_rect <- data.frame(xmin = hi_start,
                          xmax = hi_end,
                          ymin = 1,
                          ymax = length(piles) + 1,
                          fill = highlight_color)

    pile_plot <- pile_plot +
      geom_rect(data = hi_rect,
                aes(xmin = xmin, xmax = xmax,
                    ymin = ymin, ymax = ymax,
                    fill = fill))

  }

  if(is.null(max_val)) {
    max_val <- max(unlist(lappy(piles, function(x) max(x$val))))
  }

  for(i in 1:length(piles)) {
    pile <- piles[[i]]
    #pile_color <- group_colors[names(group_colors) == names(piles)[i]]

    pile$val <- pile$val / max_val

    baseline <- data.frame(x = min(pile$pos),
                           xend = max(pile$pos),
                           y = i,
                           yend = i,
                           color = "#000000")
    pile <- pile %>%
      filter(val > 0)

    pile$ypos <- i

    pile$color <- values_to_colors(pile$val,
                                   min_val = 0,
                                   max_val = 1,
                                   colorset = colorset)

    if(baselines) {
      pile_plot <- pile_plot +
        geom_segment(data = baseline,
                     aes(x = x, xend = xend,
                         y = y, yend = y,
                         color = color),
                     size = 0.1)
    }

    pile_plot <- pile_plot +
      geom_rect(data = pile,
                aes(xmin = pos,
                    xmax = pos + window_size,
                    ymin = ypos + 0.05,
                    ymax = ypos + 0.95,
                    fill = color),
                color = NA)
  }

  pile_plot
}
