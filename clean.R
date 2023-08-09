## recode all 999 as NA
bp[bp == 999] <- NA

## convert variable names to lower case
names(bp) <- tolower(names(bp))
