trim_to_bankfull <- function(xsection,
                             edge_prop = 0.15,
                             tol = 0.0) {  # small tolerance if you want, e.g. 0.05
  
  xsection %>%
    group_by(OBJ_ORDER) %>%
    arrange(distance, .by_group = TRUE) %>%
    group_modify(~ {
      df <- .x
      
      # 1. Define edge zones (just for estimating bankfull)
      width <- max(df$distance, na.rm = TRUE) - min(df$distance, na.rm = TRUE)
      left_edge  <- df$distance <= min(df$distance, na.rm = TRUE) + edge_prop * width
      right_edge <- df$distance >= max(df$distance, na.rm = TRUE) - edge_prop * width
      
      if (!any(left_edge) || !any(right_edge)) return(df)
      
      # 2. Max elevation in each edge zone
      left_max  <- max(df$RASTERVALU[left_edge],  na.rm = TRUE)
      right_max <- max(df$RASTERVALU[right_edge], na.rm = TRUE)
      
      # Common bankfull elevation = smaller of the two
      target <- min(left_max, right_max)
      
      # 3. Find innermost bankfull point on each side -------------------------
      # Left: in left_edge, elevation ≥ target (minus tol), choose largest distance
      left_cand <- which(left_edge & df$RASTERVALU >= (target - tol))
      if (length(left_cand) == 0) return(df)
      left_bank_dist <- max(df$distance[left_cand])
      
      # Right: in right_edge, elevation ≥ target (minus tol), choose smallest distance
      right_cand <- which(right_edge & df$RASTERVALU >= (target - tol))
      if (length(right_cand) == 0) return(df)
      right_bank_dist <- min(df$distance[right_cand])
      
      # 4. Trim section between those two distances ---------------------------
      df %>%
        filter(distance >= left_bank_dist,
               distance <= right_bank_dist)
    }) %>%
    ungroup()
}
