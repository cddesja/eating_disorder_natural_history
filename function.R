#' Calculate prevalence or cumulative incidence
#' @param x the name of eating disorder
#' @param incidence calculate incidence? 
calc_prevalence <- function(x, incidence = FALSE, sum = FALSE) {
  if (incidence) {

    ## subset only participants without the eating disorder (i.e., for incidence)
    elig_parts <- bp |>
      select(id, paste0(x, ".00")) |>
      filter(!!rlang::sym(paste0(x, ".00")) == 0) |>
      pull(id)

    ## determine whether the participants developed the eating disorder during the study
    tmp <- bp |>
      select(id, starts_with(x)) |>
      select(id, paste0(x, ".01"):paste0(x, ".36")) |>
      filter(id %in% elig_parts) |>
      pivot_longer(cols = -id) |>
      drop_na() |>
      group_by(id) |>
      summarize(dev_ed = max(value, na.rm = TRUE)) 
 
  } else {
  
    ## calculate prevalence, this includes baseline eating disorder
    tmp <- bp |>
      select(id, starts_with(x)) |>
      select(id, paste0(x, ".00"):paste0(x, ".36")) |>
      pivot_longer(cols = -id) |>
      drop_na() |>
      group_by(id) |>
      summarize(dev_ed = max(value, na.rm = TRUE)) |>
      mutate(ed = x)
  }
  if (sum) {
    
    ## how many participants developed the eating disorder
    tmp |> pull(dev_ed) |> sum()
    
  } else {
    
    ## return a vector indicating whether the eating disorder was observed (1) or not (0)
    tmp
  }
}

#' Calculate 95% CI for prevalence or cumulative incidence
#' @param x the name of eating disorder
#' @param incidence calculate incidence? 
calc_prevalence_ci <- function(x, incidence = FALSE) {
  if (incidence) {
    
    ## subset only participants without the eating disorder (i.e., for incidence)
    elig_parts <- bp |> 
      select(id, paste0(x, ".00")) |>
      filter(!!rlang::sym(paste0(x, ".00")) == 0) |>
      pull(id)
    
    ## extract a vector indicating whether the eating disorder was observed (1) or not (0)
    y <- bp |> 
      select(id, starts_with(x)) |>
      select(id, paste0(x, ".01"):paste0(x, ".36")) |>
      filter(id %in% elig_parts) |>
      pivot_longer(cols = -id) |>
      drop_na() |>
      group_by(id) |>
      summarize(dev_ed = max(value, na.rm = TRUE)) |>
      pull(dev_ed)

    ## what proportion of participants had the eating disorder?
    p_hat <- mean(y)

    ## calculate the 95% CI using the normal approximation
    se <- sqrt((p_hat * (1 - p_hat)) / length(y))
    ci <- (p_hat + c(-1.96, 1.96) * se)
    p_hat <- round(p_hat * 100, 1)
    ci <- round(ci * 100, 1)

    ## return the proportion and 95% CI
    paste0(p_hat, " (", ci[1], ", ", ci[2], ")")
  } else {

    ## extract a vector indicating whether the eating disorder was observed (1) or not (0)
    y <- bp |> 
      select(id, starts_with(x)) |>
      select(id, paste0(x, ".00"):paste0(x, ".36")) |>
      pivot_longer(cols = -id) |>
      drop_na() |>
      group_by(id) |>
      summarize(dev_ed = max(value, na.rm = TRUE)) |>
      pull(dev_ed) 
    
    ## what proportion of participants had the eating disorder?
    p_hat <- mean(y)
    
    ## calculate the 95% CI using the normal approximation
    se <- sqrt((p_hat * (1 - p_hat)) / length(y))
    ci <- (p_hat + c(-1.96, 1.96) * se)
    p_hat <- round(p_hat * 100, 1)
    ci <- round(ci * 100, 1)
    
    ## return the proportion and 95% CI
    paste0(p_hat, " (", ci[1], ", ", ci[2], ")")
  }
}

#' Calculate incidence per 100,000 person years
#' @param x the name of eating disorder
calc_incidence_person_yrs <- function(x) {
  
  ## create a month/year data frame excluding baseline
  my <- tibble(month = c(paste0("0",1:9), 
    paste0("1",0:9),
    paste0("2",0:9),
    paste0("3",0:9),
    paste0("4",0:8)),
    year = rep(c(1, 2, 3, 4), c(12, 12, 12, 12)))
  
  ## calculate the total number of months and years per person
  tmp <- bp |> 
    select(id, starts_with(x)) |>
    pivot_longer(cols = -id) |>
    mutate(month = str_split_fixed(name, "[.]", n = 2)[, 2]) |>
    inner_join(my, by = "month") |>
    
    ## exclude year 4 because this isn't occurring in most analyses
    filter(year != 4) |>
    
    ## dropping the missing data
    drop_na() |>
    
    group_by(year) |>
    summarize(totl = length(unique(id)),
      totl_mo = n(),
      totl_yr = totl_mo / 12)
  
  ## calculate the total years
  pop_tf <- sum(tmp$totl)
  
  incidence <- calc_prevalence(x, incidence = TRUE, sum = TRUE)
  person_yrs <- 100000
  ((incidence / pop_tf) * person_yrs) |>
    round()
}

#' Calculate the annual prevalence of an eating disorder
#' @param x the name of eating disorder
annual_prevalence <- function(x) {
  
  ## create a month/year data frame 
  my <- tibble(month = c(paste0("0",0:9), 
    paste0("1",0:9),
    paste0("2",0:9),
    paste0("3",0:9),
    paste0("4",0:8)),
    year = rep(c(1, 2, 3, 4), c(13, 12, 12, 12)))
  
  ## exclude year 4 because this isn't occurring in most analyses
  tmp <- bp |> 
    select(id, starts_with(x)) |>
    pivot_longer(cols = -id) |>
    mutate(month = str_split_fixed(name, "[.]", n = 2)[, 2]) |>
    left_join(my, by = "month") |>
    filter(year != 4) 
  
  ## calculate whether a participant had the eating disorder in a given year
  tmp |>
    drop_na() |>
    mutate(ed = x) |>
    group_by(id, ed, year) |>
    summarize(dev_ed = max(value)) |>
    ungroup() |>

    ## summarize over persons to find yearly total nubmers of cases
    group_by(ed, year) |>
    summarize(
      annual_prev = sum(dev_ed)) 
}

#' Calculate the duration of an episode
#' @param x the name of eating disorder
#' @param summarize should the raw data or summary information be provided?
calc_epi_dur <- function(x, summarize = FALSE) {
  
  ## calculate the duration of the episodes
  epi_dur <- bp |> 
    select(id, starts_with(x)) |>
    select(id, paste0(x, ".00"):paste0(x, ".36")) |>
    pivot_longer(cols = -id) |>
    mutate(time = str_split_fixed(name, "[.]", 2)[,2] |> as.numeric()) |>
    group_by(id) |>
    
    ## see calc_len_episodes() below
    mutate(num_months = calc_len_episodes(value)) |>
    filter((value == 1 | is.na(value)) & num_months != 0) |>
    mutate(ed = x,
      wave = find_wave(time)) 
  
  if(summarize) {
    epi_dur |>
      group_by(id, wave) |>
      summarize(epi_dur = max(num_months)) |>
      ungroup() |>
      summarize(
        ed = x,
        ave = mean(epi_dur),
      std_dev = sd(epi_dur),
      minimum = min(epi_dur),
      maximum = max(epi_dur))
  } else {
      epi_dur
    }
}

#' Calculate remission rate
#' @param x the name of eating disorder
calc_remission <- function(x) {
  ## identify the participants that developed an eating disorder 
  parts_w_ed <- bp |> 
    select(id, starts_with(x)) |>
    select(id, paste0(x, ".00"):paste0(x, ".36")) |>
    pivot_longer(cols = -id) |>
    mutate(time = rep(0:36, n_distinct(id))) |>
    filter(value == 1) |>
    pull(id) |>
    unique()
  
  ## identify the timing of the first episode
  first_epi <- bp |> 
    filter(id %in% parts_w_ed) |>
    select(id, starts_with(x)) |>
    select(id, paste0(x, ".00"):paste0(x, ".36")) |>
    pivot_longer(cols = -id) |>
    filter(value == 1) |>
    group_by(id) |>
    slice(1) |>
    mutate(time = str_split_fixed(name, "[.]", 2)[,2] |> as.numeric()) |>
    select(id, time)

  ## create a sequence from the start of the first episode to remission
  df <- create_sequence(first_epi)

  ## identify time when eating disorder went into remission
  bp |> 
    filter(id %in% parts_w_ed) |>
    select(id, starts_with(x)) |>
    select(id, paste0(x, ".00"):paste0(x, ".36")) |>
    pivot_longer(cols = -id) |>
    mutate(time = str_split_fixed(name, "[.]", 2)[,2] |> as.numeric()) |>
    inner_join(df, by = join_by(id, time)) |>
    filter(value == 0) |>
    group_by(id) |>
    slice(1) |>
    rename(remission_time = time) |>

    ## join the data back with the time of remission to calculate time-to-remission
    right_join(first_epi, by = join_by(id)) |>
    mutate(time_to_remission = remission_time - time) |>
    rename(episode_time = time) |>
    mutate(ed = x) |> 
    ungroup() |>
    select(7, 1, 5, 4, 6)
}

#' Create a sequence of numbers. Internal function
#' @param x a tibble containing a column named ID and one named time
create_sequence <- function(x, to = 36) {
  months <- NULL
  ids <- NULL
  
  for(i in 1:nrow(x)){
    tmp_time <- seq(from = x$time[i], to = to)
    tmp_id <- rep(x$id[i], length(tmp_time))
    months <- c(months, tmp_time) 
    ids <- c(ids, tmp_id)
  }
  df <- tibble(id = ids, time = months)
  df
}

#' Calculate recurrence
#' @param x the name of eating disorder
calc_recurrence <- function(x) {
  ## id participants that remitted and extract
  ## id and time of remission
  cr <- calc_remission(x)
  tmp <- cr |>
    drop_na() |>
    select(id, remission_time) |>
    rename(time = remission_time)

  ## create vector of time of remit to 36 months
  df <- create_sequence(tmp)

  ## identify time of recurrence
  parts_that_recurred <- bp |> 
    filter(id %in% tmp$id) |>
    select(id, starts_with(x)) |>
    select(id, paste0(x, ".00"):paste0(x, ".36")) |>
    pivot_longer(cols = -id)  |>
    mutate(time = str_split_fixed(name, "[.]", 2)[,2] |> as.numeric()) |>
    inner_join(df, by = join_by(id, time)) |>
    filter(value == 1) |>
    group_by(id) |>
    slice(1) |>
    rename(recurred = value,
      recurred_time = time)

  ## recode 0 if they didn't recur but did remit
  cr <- cr |>
    full_join(parts_that_recurred 
      |> select(-name), by = join_by(id)) |>
    mutate(recurred = ifelse(is.na(recurred) & !is.na(time_to_remission), 0, recurred))
  
  ## return recurred, time of recurrence, and time-to-recurrence
  cr |>
    mutate(ed = x,
      time_to_recurred = recurred_time - remission_time ) |>
    select(ed, id, recurred, recurred_time, time_to_recurred)
}

#' Calculate diagnostic progression
calc_diag_progress <- function(x) {

  ## full
  fx <- paste0("f", x)
  parts_w_fx <- bp |> 
    select(id, starts_with(fx)) |>
    select(id, paste0(fx, ".00"):paste0(fx, ".36")) |>
    pivot_longer(cols = -id) |>
    mutate(time = rep(0:36, n_distinct(id))) |>
    filter(value == 1) |>
    group_by(id) |>
    slice(1) |>
    rename(full_time = time,
      full_value = value) |>
    select(-name)
  
  ## partial 
  px <- paste0("p", x)
  parts_w_px <- bp |> 
    select(id, starts_with(px)) |>
    select(id, paste0(px, ".00"):paste0(px, ".36")) |>
    pivot_longer(cols = -id) |>
    mutate(time = rep(0:36, n_distinct(id))) |>
    filter(value == 1) |>
    group_by(id) |>
    slice(1) |>
    rename(part_time = time,
      part_value = value) |>
    select(-name)
  
  parts_w_px |>
    left_join(parts_w_fx) |>
    arrange(desc(full_value)) |>
    mutate(
      ed = x,
      progress = ifelse(part_time < full_time, 1, 0),
      progress = ifelse(is.na(progress), 0, progress))
}

#' Determine what the episodic wave is
#' @param x a variable containing months since the start of the study
find_wave <- function(x) {
  c(0, (diff(x) != 1) |> cumsum())
}

#' Calculate length of episode
#' https://stackoverflow.com/questions/76632703/conditional-cumulative-sum-in-r-with-na-and-keeping-track-of-episodic-wave/76632818?noredirect=1#comment135117892_76632818
calc_len_episodes <- function(z) {
  r <- rle(replace(z, is.na(z), Inf))
  for (na in which(is.infinite(r$values))) {
    if (na == 1L) {
      r$values[na] <- 0
    } else if (na == length(r$values)) {
      r$values[na] <- 0
    } else {
      r$values[na] <- ifelse(r$values[na - 1] == r$values[na + 1], r$values[na - 1], 0.5)
    }
  }
  z2 <- inverse.rle(r)
  ave(z2, cumsum(z2 == 0), FUN = cumsum)
}

#' Determine if a series of numbers is sequential and number them. Used internally by calc_epi_dur()
#' @param x numeric variable
zz_count_seq <- function(x) {
  cnt <- NULL
  wave <- NULL
  epi <- 1
  i <- 1
  for(j in 1:length(x)) {
    if(j == 1) {
      cnt[j] <- 1
      wave[j] <- 1
    } else {
      cnt[j] <- ifelse(x[j] - x[j - 1] == 1, i + 1, 1)
      if(x[j] - x[j - 1] == 1) {
        i <- i + 1
        wave[j] <- epi
      } else {
        i <- 1
        epi <- epi + 1
        wave[j] <- epi
      }
    }
  }
  data.frame(cnt, wave)
}

