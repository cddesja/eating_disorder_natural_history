## table 2 missing data
my <- tibble::tibble(month = c(paste0("0",0:9), 
  paste0("1", 0:9),
  paste0("2", 0:9),
  paste0("3", 0:6)),
  year = rep(c(0, 1, 2, 3), c(1, 12, 12, 12)))

miss_across_study <- bp |>
  select(id,
    fan.01:fan.36, 
    pan.01:pan.36,
    fbn.01:fbn.36,
    pbn.01:pbn.36,
    fbe.01:fbe.36, 
    pbe.01:pbe.36,
    fpu.01:fpu.36,
    aan.01:aan.36) |>
  pivot_longer(cols = -1) |>
  mutate(
    ed = str_split_fixed(name, "[.]", n = 2)[, 1],
    month = str_split_fixed(name, "[.]", n = 2)[, 2]) |>
  left_join(my) |>
  group_by(id, year, ed) |>
  summarize(missing = sum(is.na(value))) |>
  ungroup() |>
  group_by(id, year) |>
  summarize(months_missed_across_ed = sum(missing)) |>
  ungroup()
  
  miss_across_study |>  
  group_by(year) |>
  summarize(prop_miss = mean(months_missed_across_ed == 96))

bp |>
  select(id,
    fan.00, 
    pan.00,
    fbn.00,
    pbn.00,
    fbe.00, 
    pbe.00,
    fpu.00,
    aan.00) |>
  pivot_longer(cols = -1) |>
  mutate(
    ed = str_split_fixed(name, "[.]", n = 2)[, 1],
    month = str_split_fixed(name, "[.]", n = 2)[, 2]) |>
  left_join(my) |>
  group_by(id) |>
  summarize(months_missed_across_ed = sum(is.na(value))) |>
  ungroup() |>
  summarize(prop_miss = mean(months_missed_across_ed <= 4))

bpae <- haven::read_spss("/Users/cdesjardins/Library/CloudStorage/Dropbox/chris/ori/data/bp1234/BodyProjectsAllEpisodes.sav")
names(bpae) <- tolower(names(bpae))

## --- mental health utilization --- ##
#- mental health care is missing for BP-IV
table(bpae$sample, is.na(bpae$mentx1))
table(bpae$sample, is.na(bpae$mentx2))
table(bpae$sample, is.na(bpae$mentx3))

## --- emotional distress --- ##
#- not sure what this variable is called

## --- psychosocial functioning --- ##
table(bp$sample, is.na(bp$w1socf))
table(bp$sample, is.na(bp$w2socf)) ### missing for BP 1, 2, 4
table(bp$sample, is.na(bp$w3socf)) ### missing for BP 1 and 3

## --- suicidality --- ##
#- not sure what this variable is called

## --- body mass index --- ##
table(bp$sample, is.na(bp$w1intbmi))
table(bp$sample, is.na(bp$w2intbmi)) 
table(bp$sample, is.na(bp$w3intbmi)) ### missing for bp 3


