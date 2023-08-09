## examine age ---- 
## obtain frequency table of age 
age_freq <- table(bp$w1age) |> as.vector()

## find distinct ages and add 3 to represent 36-month
age <- names(table(bp$w1age)) |>
  as.numeric() |>
  (\(x) (x + 3))()

## combine into a tibble
age_tmp <- tibble(age,
  age_freq)
  
## 
y <- NULL
for (i in 1:nrow(age_tmp)) {
  y <- c(y, seq(1, age_tmp$age_freq[i]))
}
x <- rep(age_tmp$age, age_tmp$age_freq) |> as.numeric()
age_dist <- data.frame(x = x, y = y)

age_dist |>
  ggplot(aes(x = x, y = y)) +
  geom_point(alpha = 1/2, pch = 21) + 
  geom_label(aes(x = 60, y = 40, label = paste("Mean:", mean(x) |> round(1)))) +
  geom_label(aes(x = 60, y = 60, label = paste("Median:", median(x) |> round(1)))) +
  theme_bw() +
  xlab("Age at 36 Months") +
  ylab("Count") +
  ggtitle("Distribution of BP1234 Participants' Ages at 36 Months") 
## QUESTION: 

## Attrition ----


## TABLE 3 ----
### lifetime prevalence
disorders <- c(fan = "fan", fbn = "fbn", fbe = "fbe", fpu = "fpu", aan = "aan", pan = "pan", pbn = "pbn", pbe = "pbe")
disorders |>
  map_dbl(calc_prevalence, sum = TRUE)
disorders |> 
  map_chr(calc_prevalence_ci)

disorders |>
  map_df(calc_prevalence) |>
  filter(dev_ed == 1) |>
  pull(id) |>
  n_distinct() 

disorders |>
  map_df(calc_prevalence) |>
  filter(dev_ed == 1) |>
  pull(id) |>
  n_distinct() / nrow(bp)

disorders |>
  map_df(calc_prevalence) |>
  filter(dev_ed == 1 & ed %in% c("fan", "fbn", "fbe")) |>
  pull(id) |>
  n_distinct() / nrow(bp)

### cumulative incidence
disorders |> 
  map_dbl(\(x) calc_prevalence(x, cumulative = TRUE, sum = TRUE))
disorders |> 
  map_chr(\(x) calc_prevalence_ci(x, cumulative = TRUE))

### incidence per 100,000 person years
disorders |>
  map_dbl(calc_incidence_person_yrs)

### annual prevalence ----
disorders |>
  map_df(annual_prevalence)

## episode duration ----
## an episode was defined as length of the period in which a participant met all diagnostic criteria in sequential months.
disorders |> 
  map_df(\(x) calc_epi_dur(x, summarize = TRUE))

epi_dur <- disorders |> 
  map_df(calc_epi_dur)

## episode remission ----
## remission was defined as not meeting all criteria for an eating disorder for at least a 1-month period
epi_remi <- disorders |> 
  map_df(calc_remission)

epi_remi |>
  mutate(remit_yr1 = time_to_remission <= 12,
    remit_yr2 = time_to_remission <= 24,
    remit_yr3 = time_to_remission <= 36,
    remit_yr1 = ifelse(is.na(remit_yr1), FALSE, remit_yr1),
    remit_yr2 = ifelse(is.na(remit_yr2), FALSE, remit_yr2),
    remit_yr3 = ifelse(is.na(remit_yr3), FALSE, remit_yr3)) |>
  group_by(ed) |>
  summarize(remit_rate_yr1 = mean(remit_yr1) |> round(2),
    remit_rate_yr2 = mean(remit_yr2) |> round(2),
    remit_rate_yr3 = mean(remit_yr3) |> round(2),
    totl_remit_yr1 = sum(remit_yr1),
    totl_remit_yr2 = sum(remit_yr2),
    totl_remit_yr3 = sum(remit_yr3),
    totl_remit = sum(!is.na(time_to_remission)),
    n = n())

## episode recurrence ----
## defined as meeting criteria for another episode of the same eating disorder after exhibiting remission for at least 1 month
epi_dur |>
  select(ed, id, wave) |>
  group_by(ed, id) |>
  summarize(max_wave = max(wave)) |>
  select(ed, max_wave) |>
  table() |>
  as_tibble() |>
  mutate(ed = factor(ed, levels = disorders)) |>
  arrange(ed) |>
  filter(n != 0) |>
  group_by(ed) |>
  mutate(prop = n / sum(n),
    totl = sum(n),
    recur_rate = 1 - prop[1] |> round(2)) |>
  select(-prop) |>
  pivot_wider(names_from = c(max_wave),
    values_from = n)

epi_rec <-  disorders |> 
  map_df(calc_recurrence)
epi_rec |>
  group_by(ed) |>
  summarize(remitted_parts = sum(!is.na(recurred)),
    recurred_parts = sum(recurred, na.rm = TRUE),
    recur_rate = mean(recurred, na.rm = TRUE),
    ave_length_to_recur = mean(time_to_recurred, na.rm = TRUE))

## diagnostic progression ----
subthres <- c("an", "bn", "be")
diag_prog <- subthres |> 
  map_df(calc_diag_progress)

diag_prog |>
  group_by(ed) |>
  summarize(prop = mean(progress) |> round(2),
    totl = sum(progress))

## diagnostic crossover ----
ed_months <- paste(rep(disorders, each = 37), c(paste0("0", 0:9), 10:36), sep = "."); months

all_eds <- bp |> 
  select(id, all_of(ed_months)) |>
  pivot_longer(cols = -id) |>
  mutate(time = str_split_fixed(name, "[.]", 2))  |>
  mutate(ed = time[, 1],
    months = time[, 2] |> as.numeric()) |>
  select(-time) |>
  filter(value == 1) |>
  group_by(ed, id) |>
  arrange(months) |>
  slice(1) |> 
  select(-name, -value) |>
  ungroup()

all_eds_wide <- all_eds |>
  pivot_wider(names_from = ed, values_from = months)

all_eds_wide$frst_ed <- names(all_eds_wide)[-1][apply(all_eds_wide |> select(-id), 1, which.min)]

all_eds_wide |>
  pivot_longer(cols = 2:9) |>
  filter(!is.na(value)) |>
  select(frst_ed, name) |>
  table() 
## rows do not need to sum up to ED totals because they are not mutually exclusive!!! In other words, a participant could occur in multiple columns.
