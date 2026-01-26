
# 1. load required packages
# install.packages("ggdist")
library(plyr)
library(tidyverse)
library(brms)
library(tidybayes)
library(rstantools)
library(ggnewscale)  # used for Figure 4
library(patchwork)   # used for Figure 4
library(bayesplot)   # used for caterpillar plots
library(rnaturalearth)        # used for Figure 1
library(rnaturalearthhires)   # used for Figure 1
library(ggrepel)              # used for Figure 1
library(ggh4x)                # used for Figure 2
library(ggpubr)               # used for Figure 3
library(ggdist)               # used for Figure 3

# 2. set options
options(scipen = 999,                     # digits before scientific notation
        dplyr.summarise.inform = FALSE)   # silence grouping messages

# 3. set ggplot themes & color palettes
# 3. a) default theme for ggplot
theme_set(theme_bw())
theme_update(panel.grid   = element_blank(),
             legend.title = element_text(size = 9,  color = "black"),
             legend.text  = element_text(size = 9,  color = "black"),
             axis.title   = element_text(size = 10, color = "black"),
             axis.text    = element_text(size = 9,  color = "black"),
             strip.text   = element_text(size = 9,  color = "black"),
             plot.tag.position = c(0, 1))
# 3. b) custom colour palette for plots
pal11 <- scales::viridis_pal(option = "plasma")(11)
pal3  <- scales::viridis_pal(option = "plasma", end = 0.8)(3)


# 4. import raw dataset and select needed columns
fecundityRaw <- readRDS(
  "C:/Users/IMLAYT/OneDrive - DFO-MPO/Atlantic salmon body size project/MS fecundity/Analysis/fecundity.rds") %>% 
  select(DU, RIVER_NO, RIVER, YEAR, LATITUDE, LONGITUDE, EGG_COUNT_COR, 
         EGG_COUNT_RETAIN, EGG_DIAM_COR, STAGE_EGG_COUNT, FORK_LENGTH, SIZE, 
         RAGE, LAST_SPAWN, FULTONS_K, COLLECT_DOY, SAMPLE_DOY, STRIP_DOY,
         SAMPLE_DATE_RANGE, CITATION) %>% 
  droplevels()
saveRDS(fecundityRaw, "data - fecundity.rds")

# 6. set priors for models
priorList <- c(prior(normal(0, 10), class = b, coef = "Intercept"), # explicitly set prior for non-centered intercept
               prior(normal(0, 1),  class = b),
               prior(normal(0, 1),  class = sd),
               prior(normal(0, 1),  class = sigma),
               prior(lkj(2),        class = cor))

# 7. allow for parallel processing
options(mc.cores = parallel::detectCores())




# DATASET SUMMARY INFORMATION ---------------------------------------------

# 1. sample sizes
# 1. a) number of salmon
dim(fecundityRaw)
# 1. b) number of rivers
summary(fecundityRaw$RIVER) %>% length()
# 1. c) number of DUs
summary(fecundityRaw$DU) %>% length()

# 2. range of collection and stripping dates
# 2. a) NL rivers
filter(fecundityRaw, if_any(c(COLLECT_DOY, STRIP_DOY), ~!is.na(.)),
       DU %in% c("DU04", "DU05", "DU06", "DU07", "DU08", "DU09")) %>% #View()
  group_by(STAGE_EGG_COUNT) %>% 
  summarise(minCD = min(COLLECT_DOY, na.rm = TRUE), maxCD = max(COLLECT_DOY, na.rm = TRUE),
            minSD = min(STRIP_DOY, na.rm = TRUE), maxSD = max(STRIP_DOY, na.rm = TRUE)) %>% 
  mutate(across(c(minCD, maxCD, minSD, maxSD), 
                ~case_when(abs(.) == Inf ~ NA_character_,
                           .<=(31+28+31+30+31+30)               ~ paste(.-(31+28+31+30+31), "June", sep = " "),
                           .<=(31+28+31+30+31+30+31*2+30)       ~ paste(.-(31+28+31+30+31+30+31*2), "September", sep = " "),
                           .<=(31+28+31+30+31+30+31*2+30+31)    ~ paste(.-(31+28+31+30+31+30+31*2+30), "October", sep = " "),
                           .<=(31+28+31+30+31+30+31*2+30+31+30) ~ paste(.-(31+28+31+30+31+30+31*2+30+31), "November", sep = " ")),
                .names = "{.col}ay"))
# 2. b) QU rivers 
filter(fecundityRaw, !is.na(COLLECT_DOY), DU %in% c("DU12", "DU15")) %>% 
  summarise(minCD = min(COLLECT_DOY), maxCD = max(COLLECT_DOY)) %>% 
  mutate(across(c(minCD, maxCD), 
                ~case_when(.<=(31+28+31+30+31+30)               ~ paste(.-(31+28+31+30+31), "June", sep = " "),
                           .<=(31+28+31+30+31+30+31*2+30)       ~ paste(.-(31+28+31+30+31+30+31*2), "September", sep = " ")),
                .names = "{.col}ay"))
# 2. c) NB rivers
filter(fecundityRaw, DU == "DU16", RIVER != "Margaree River",
       if_any(c(COLLECT_DOY, STRIP_DOY), ~!is.na(.))) %>%
  group_by(str_sub(YEAR, 1, 3), RIVER) %>% 
  summarise(minCD = min(COLLECT_DOY, na.rm = TRUE), maxCD = max(COLLECT_DOY, na.rm = TRUE),
            minSD = min(STRIP_DOY, na.rm = TRUE), maxSD = max(STRIP_DOY, na.rm = TRUE)) %>% 
  mutate(across(c(minCD, maxCD, minSD, maxSD), 
                ~case_when(abs(.) == Inf ~ NA_character_,
                           .<=(31+28+31+30+31)                  ~ paste(.-(31+28+31+30), "May", sep = " "),
                           .<=(31+28+31+30+31+30)               ~ paste(.-(31+28+31+30+31), "June", sep = " "),
                           .<=(31+28+31+30+31+30+31*2+30)       ~ paste(.-(31+28+31+30+31+30+31*2), "September", sep = " "),
                           .<=(31+28+31+30+31+30+31*2+30+31)    ~ paste(.-(31+28+31+30+31+30+31*2+30), "October", sep = " "),
                           .<=(31+28+31+30+31+30+31*2+30+31+30) ~ paste(.-(31+28+31+30+31+30+31*2+30+31), "November", sep = " ")),
                .names = "{.col}ay"))
# 2. d) Margaree
filter(fecundityRaw, RIVER == "Margaree River") %>%
  summarise(minSD = min(STRIP_DOY, na.rm = TRUE), maxSD = max(STRIP_DOY, na.rm = TRUE)) %>% 
  mutate(across(c(minSD, maxSD), 
                ~case_when(.<=(31+28+31+30+31+30+31*2+30+31)    ~ paste(.-(31+28+31+30+31+30+31*2+30), "October", sep = " "),
                           TRUE ~ paste(.-(31+28+31+30+31+30+31*2+30+31+30), "December", sep = " ")),
                .names = "{.col}ay"))
# 2. e) rivers with estimated sample dates
filter(fecundityRaw, !is.na(SAMPLE_DATE_RANGE)) %>% 
  summarise(round(sum(SAMPLE_DATE_RANGE == 1)/5331*100, 1))
filter(fecundityRaw, SAMPLE_DATE_RANGE > 1) %>% 
  group_by(RIVER_NO, RIVER) %>% 
  summarise(nMid = sum(SAMPLE_DATE_RANGE > 1), 
            minSD = min(SAMPLE_DATE_RANGE), 
            maxSD = max(SAMPLE_DATE_RANGE)) #%>% 
  filter(maxSD > 1)
names(fecundityRaw)

# 3. NL rivers compare fecundity relationships & egg retention rates
# 3. a) compare relationships for Rocky River & Flat Bay Brook
summary(lm(EGG_COUNT ~ FORK_LENGTH, data = filter(fecundityRaw, grepl("Rocky", RIVER)))) # matches Bourgeois et al 1997
summary(lm(EGG_COUNT ~ FORK_LENGTH, data = filter(fecundityRaw, grepl("Flat", RIVER))))
# 3. b) calculate number of females with egg retention number
filter(fecundityRaw, STAGE_EGG_COUNT == "mature + retained") %>% 
  group_by(RIVER) %>% summarise(n = n(), nRet = sum(!is.na(EGG_COUNT_RETAIN))) %>% 
  summarise(sum(nRet), sum(n), sum(nRet)/sum(n))
# 3. c) calculate proportion of retained eggs in fecundity value
filter(fecundityRaw, !is.na(EGG_COUNT_RETAIN), STAGE_EGG_COUNT != "immature") %>%   ### DOUBLE CHECK THAT CONNE IMMATURE EGGS ARE NOT INCLUDED!!!!
  mutate(propRet = EGG_COUNT_RETAIN/EGG_COUNT_COR * 100) %>% 
  group_by(RIVER_NO, RIVER, LAST_SPAWN) %>% summarise(n= n(), 
                                          m = round(mean(propRet), 1), 
                                          sd = round(sd(propRet), 1), 
                                          min = round(min(propRet), 1), 
                                          max = round(max(propRet), 1))


filter(fecundityRaw, grepl("Conn|Expl|Flat|Indian|Pipers|Rocky|Terra", RIVER)) %>% 
  mutate(eggs = if_else(!is.na(EGG_COUNT_RETAIN), EGG_COUNT_COR-EGG_COUNT_RETAIN, EGG_COUNT_COR),
         method = if_else(!is.na(EGG_COUNT_RETAIN), "mature", STAGE_EGG_COUNT)) %>% 
  count(RIVER, method)

filter(fecundityRaw, !is.na(EGG_COUNT_RETAIN), STAGE_EGG_COUNT != "immature") %>%   ### DOUBLE CHECK THAT CONNE IMMATURE EGGS ARE NOT INCLUDED!!!!
  mutate(propRet = EGG_COUNT_RETAIN/EGG_COUNT_COR * 100) %>% 
  ggplot(aes(x = SAMPLE_DOY, y = propRet)) +
  stat_summary(aes(yintercept = after_stat(y), x = 300, col = RIVER), fun = mean, geom = "hline", linewidth = 0.2) +
  geom_jitter(aes(col = RIVER), shape = 1) +
  # geom_smooth() +
  scale_y_continuous(breaks = seq(0, 75, 5)) +
  labs(y = "prop. retained (%)", x = "sample DOY")


# 4. missing age data strategies spawning strategy summaries
# 4. a) proportion of salmon missing RAGE or LAST_SPAWN
filter(fecundityRaw, if_any(c(RAGE, LAST_SPAWN), ~is.na(.))) %>% 
  summarise(round(n()/5331*100, 1))
# 4. b) spawning strategies in size classes
filter(fecundityRaw, !is.na(LAST_SPAWN)) %>% count(RIVER, SIZE, LAST_SPAWN) %>%
  group_by(RIVER, SIZE) %>% mutate(tot = sum(n)) %>% spread(LAST_SPAWN, n) %>%
  mutate(sm1 = round(`1`/tot*100, 1),
         lg2 = round(`2`/tot*100, 1),
         lgR  = round((A + C)/tot*100, 1)) %>% group_by(SIZE) %>%
  summarise(meanSDsm1 = paste(round(mean(sm1, na.rm = TRUE), 1),  
                              round(sd(sm1, na.rm = TRUE), 1), sep = " ± "),
            sm1range  = paste(min(sm1, na.rm = TRUE), max(sm1, na.rm = TRUE), sep = "-"),
            
            meanSDlg2 = paste(round(mean(lg2, na.rm = TRUE), 1),  
                              round(sd(lg2, na.rm = TRUE), 1), sep = " ± "),
            lg2range  = paste(min(lg2, na.rm = TRUE), max(lg2, na.rm = TRUE), sep = "-"),
            
            meanSDlfR = paste(round(mean(lgR, na.rm = TRUE), 1),  
                              round(sd(lgR, na.rm = TRUE), 1), sep = " ± "),
            lgRrange  = paste(min(lgR, na.rm = TRUE), max(lgR, na.rm = TRUE), sep = "-"))



# TABLE 1 -----------------------------------------------------------------

table1 <- fecundityRaw %>% 
  mutate(CITATION = case_when(RIVER == "LaHave River"       ~ "Cutting et al. 1987; Cutting & Gray 1984", 
                              RIVER == "NW Miramichi River" ~ "Reid & Chaput 2012; Roth 2024",
                              RIVER == "SW Miramichi River" ~ "Reid & Chaput 2012",
                              TRUE ~ CITATION)) %>% 
  group_by(DU, RIVER_NO, RIVER, STAGE_EGG_COUNT, CITATION) %>% 
  summarise(minY  = min(YEAR), 
            maxY  = max(YEAR),
            nY    = n_distinct(YEAR),
            Small = sum(SIZE == "small", na.rm = TRUE),
            Large = sum(SIZE == "large", na.rm = TRUE),
            SD    = sum(!is.na(SAMPLE_DOY)) - sum(SAMPLE_DATE_RANGE > 7, na.rm = TRUE),
            FK    = sum(!is.na(FULTONS_K)),
            R1    = sum(RAGE == 1, na.rm = TRUE),
            R2    = sum(RAGE == 2, na.rm = TRUE),
            R3    = sum(RAGE == 3, na.rm = TRUE),
            R4    = sum(RAGE == 4, na.rm = TRUE),
            R5    = sum(RAGE == 5, na.rm = TRUE),
            R6    = sum(RAGE == 6, na.rm = TRUE),
            `1SW` = sum(LAST_SPAWN == 1, na.rm = TRUE),
            `2SW` = sum(LAST_SPAWN == 2, na.rm = TRUE),
            `3SW` = sum(LAST_SPAWN == 3, na.rm = TRUE),
            C     = sum(LAST_SPAWN == "C", na.rm = TRUE),
            A     = sum(LAST_SPAWN == "A", na.rm = TRUE),
            diam  = sum(!is.na(EGG_DIAM_COR))) %>% 
  arrange(RIVER_NO, minY) %>% 
  mutate(Years = case_when(minY == 1981 & RIVER %in% c("LaHave River", "Medway River") ~ paste("1980", maxY, sep = "-"),
                           minY == maxY ~ as.character(minY), 
                           TRUE         ~ paste(minY, maxY, sep = "-")),
         nY    = if_else(RIVER %in% c("LaHave River", "Medway River"), nY + 1, nY),
         across(c(Small:A), ~if_else(. == 0, "", as.character(.)))) %>% 
  select(DU, `No.` = RIVER_NO, River = RIVER, Years, N = nY,
         Method = STAGE_EGG_COUNT, Small, Large, contains(c("R", "SW")), C, A,  
         `sample date` = SD, `Fulton's K` = FK, `Egg diameter` = diam)



# FECUNDITY ~ FORK LENGTH: MODELS -----------------------------------------

# 1. prepare dataset for modelling
fecundity <- count(fecundityRaw, DU, RIVER, STAGE_EGG_COUNT) %>% 
  filter(n >= 5) %>% left_join(fecundityRaw) %>% droplevels() %>%   # min. sample size = 5
  mutate(RIVER_CAT = if_else(RIVER_NO %in% c(2:15), "1SW", "MSW")) %>% 
  mutate(lEGG_COUNT_COR = log(EGG_COUNT_COR),                       # log fecundity 
         lFORK_LENGTH   = log(FORK_LENGTH),                         # log fork length
         clFORK_LENGTH  = scale(lFORK_LENGTH),                      # center & scale log fork length
         SIZE           = fct_relevel(SIZE, "large", after = 1),    # re-order size 
         RIVER          = fct_reorder(RIVER, LATITUDE))             # re-order rivers


# 2. hierarchical candidate models
# 2. a) M0: vary intercept with maturity, no term for river
M0 <- brm(lEGG_COUNT_COR ~ 0 + Intercept + clFORK_LENGTH + STAGE_EGG_COUNT,
          prior  = filter(priorList, class %in% c("b", "sigma")),
          family = gaussian(link = "identity"),
          data   = fecundity, 
          chains = 4, 
          iter   = 2000,  # half discarded as warm-up
          seed   = 42)
M0
plot(M0)
# 2. b) M1. co-vary intercept with DU (and maturity)
M1 <- brm(lEGG_COUNT_COR ~ 0 + Intercept + clFORK_LENGTH + STAGE_EGG_COUNT + (1 | DU),
          prior  = filter(priorList, class %in% c("b", "sd", "sigma")),
          family = gaussian(link = "identity"), 
          data   = fecundity, 
          chains = 4, 
          iter   = 2000,  # half discarded as warm-up
          seed   = 42)
M1
plot(M1)
# 2. c) M2. co-vary intercept and slope with DU (and maturity)
M2 <- brm(lEGG_COUNT_COR ~ 0 + Intercept + clFORK_LENGTH + STAGE_EGG_COUNT + (clFORK_LENGTH | DU),
          prior  = priorList,
          family = gaussian(link = "identity"), 
          data   = fecundity, 
          chains = 4, 
          iter   = 2000,  # half discarded as warm-up
          seed   = 42,
          control = list(adapt_delta = 0.98)) #, max_treedepth = 12))
M2
plot(M2)
# 2. d) M1. co-vary intercept with river (and maturity)
M3 <- brm(lEGG_COUNT_COR ~ 0 + Intercept + clFORK_LENGTH + STAGE_EGG_COUNT + (1 | RIVER),
          prior  = filter(priorList, class %in% c("Intercept", "b", "sd", "sigma")),
          family = gaussian(link = "identity"), 
          data   = fecundity, 
          chains = 4, 
          iter   = 2000,  # half discarded as warm-up
          seed   = 42)
M3
plot(M3)
prior_summary(M3)
# 2. e) co-vary intercept and slope with river (and maturity)
M4 <- brm(lEGG_COUNT_COR ~ 0 + Intercept + clFORK_LENGTH + STAGE_EGG_COUNT + (clFORK_LENGTH | RIVER),
          prior  = priorList,
          family = gaussian(link = "identity"), 
          data   = fecundity, 
          chains = 4, 
          iter   = 2000,  # half discarded as warm-up
          seed   = 42)
M4
plot(M4)


# 3. leave-one-out cross-validation for model comparison
looM0 <- loo(M0)
looM1 <- loo(M1)
looM2 <- loo(M2)
looM3 <- loo(M3)
looM4 <- loo(M4)
loo_compare(looM0, looM1, looM2, looM3, looM4)

looM4 <- loo(topModel)
looM4t <- loo(M4test)
loo_compare(looM4, looM4t)


# 4. select top model for subsequent summaries and plots
# 5. a) identify top model
topModel <- M4
# 5. b) save models as RDS file
saveRDS(M0, "./models/M0.rds")
saveRDS(M1, "./models/M1.rds")
saveRDS(M2, "./models/M2.rds")
saveRDS(M3, "./models/M3.rds")
saveRDS(M4, "./models/M4.rds")
# 5. c) read models from RDS file
M0 <- readRDS("./models/M0.rds")
M1 <- readRDS("./models/M1.rds")
M2 <- readRDS("./models/M2.rds")
M3 <- readRDS("./models/M3.rds")
M4 <- readRDS("./models/M4.rds")



# FECUNDITY ~ FORK LENGTH: INTERPRETATION ---------------------------------

# 1. obtain posterior distribution from top model for median-sized salmon
# 1. a) bring in top model
topModel <- readRDS("./Analysis/M4.rds")
# 1. b) centering and scaling parameters used prior to modelling
FLcenter <- attr(fecundity$clFORK_LENGTH, "scaled:center")
FLscale  <- attr(fecundity$clFORK_LENGTH, "scaled:scale")
# 1. c) median fork length for size classes rounded to nearest cm & centered/scaled
FLmed <- group_by(fecundity, SIZE) %>% 
  summarize(FORK_LENGTH = median(round(FORK_LENGTH, 0))) %>%        # median FL rounded to the nearest cm
  mutate(clFORK_LENGTH = (log(FORK_LENGTH) - FLcenter) / FLscale)   # median FL centered & scaled
# 1. d) get predicted values from posterior distribution
pred <- distinct(fecundity, RIVER_NO, RIVER, SIZE, STAGE_EGG_COUNT) %>%   # single value/river/size/maturity
  left_join(FLmed) %>%                                                    # add median FL values
  # add_predicted_draws(topModel) %>%                                       # draw from posterior distribution
  # mutate(EGG_COUNT_pred = exp(.prediction))                               # transform predicted fecundity
  add_epred_draws(topModel) %>%                                       # draw from posterior distribution
  mutate(EGG_COUNT_pred = exp(.epred))                                # transform predicted fecundity

# 2. determine differences in fecundity for same size/STAGE_EGG_COUNT
group_by(pred, SIZE, STAGE_EGG_COUNT, RIVER) %>% 
  summarise(medEGGS = round(median(EGG_COUNT_pred), 0)) %>%   # determine median fecundity
  group_by(SIZE, STAGE_EGG_COUNT) %>% slice(c(1, n())) %>%    # select min/max median fecundity
  summarise(diff = max(medEGGS) - min(medEGGS),               # difference in median fecundity
            text = paste(min(medEGGS), " (", RIVER[which.min(medEGGS)], ") to ", 
                         max(medEGGS), " (", RIVER[which.max(medEGGS)], ")", sep = ""))


# 3. STAGE_EGG_COUNT coefficients with 95% credible intervals
# 3. a) check model summary
posterior_summary(topModel)
# posterior_summary(topModel) %>% data.frame() %>%                        # obtain all coefficients
#   rownames_to_column("STAGE_EGG_COUNT") %>% filter(grepl("STAGE", STAGE_EGG_COUNT)) %>%   # filter to STAGE_EGG_COUNTs
#   mutate(across(where(is.numeric), ~round(., 2)),                       # round values to 2 decimal places
#          STAGE_EGG_COUNT = str_sub(STAGE_EGG_COUNT, 18, -1),                              # shorten default nanmes
#          CIs    = paste("(", Q2.5, "-", Q97.5, ")", sep = "")) %>%      # create string for Results text
#   select(STAGE_EGG_COUNT, Estimate, CIs)
# 3. b) conditional means (95% CIs) 
as_draws_df(topModel) %>% dplyr::select(contains(c("b_Intercept", "b_STAGE_EGG_COUNT"))) %>%
  mutate(across(c(2:3), ~. + b_Intercept)) %>%
  pivot_longer(cols = everything()) %>% 
  group_by(name) %>% 
  summarise(mv = round(exp(mean(value)), 0),
            lq = round(exp(quantile(value, probs = 0.025)), 0),
            hq = round(exp(quantile(value, probs = 0.975)), 0)) %>% 
  mutate(text = paste0(mv, " (", lq, "-", hq, ")"),
         perc = round((max(mv) - mv) / max(mv), 2) * 100)

# 4. posterior distribution density plot (within 95% quantiles) & median values as lines
png("./Figures/Fig2.png", width = 7.5, height = 6, units = "in", res = 600)
group_by(pred, RIVER, SIZE, STAGE_EGG_COUNT) %>%
  filter(EGG_COUNT_pred >= quantile(EGG_COUNT_pred, probs = 0.025),       # filter <2.5% quantile
         EGG_COUNT_pred <= quantile(EGG_COUNT_pred, probs = 0.975)) %>%   # filter >97.5% quantile
  mutate(EGG_COUNT_pred = EGG_COUNT_pred/1000,
         RIVER          = fct_reorder(RIVER, desc(RIVER_NO))) %>% 
  ggplot(aes(x = EGG_COUNT_pred, 
             y = fct_reorder(RIVER, desc(RIVER_NO)),
             group = interaction(STAGE_EGG_COUNT, RIVER),
             fill  = STAGE_EGG_COUNT,
             col   = STAGE_EGG_COUNT)) +
  geom_boxplot(position = position_dodge(preserve = "single"),    # maintain same width across all boxes
               coef = NULL, alpha = 0.5) +                        # extent whiskers to min/max values (already within quantiles)
  # scale_fill_viridis_d(option = "plasma") +
  # scale_color_viridis_d(option = "plasma") +
  scale_fill_manual(values = pal3) + scale_color_manual(values = pal3) +
  scale_x_continuous(breaks = seq(1, 15, 1)) +
  facet_grid(~SIZE, scales = "free_x", space = "free_x",
           labeller = as_labeller(c(`large` = "large (76 cm)",
                                    `small` = "small (52 cm)"))) +
  labs(x = "predicted fecundity (x 1,000 eggs)", y = NULL, fill = NULL, col = NULL) +
  theme(legend.position  = "top",
        legend.direction = "horizontal")
dev.off()

# 5. conditional fecundity at size
# 5. a) coefficients from draws
coefM4 <- as_draws_df(topModel) %>% select(starts_with("b_"))
# 5. b) combine coefficients and fork length data
fl <- distinct(fecundity, STAGE_EGG_COUNT, FORK_LENGTH = round(FORK_LENGTH, 0)) %>% 
  group_by(STAGE_EGG_COUNT) %>%
  reframe(FORK_LENGTH   = seq(min(FORK_LENGTH), max(FORK_LENGTH), 1)) %>% 
  mutate(clFORK_LENGTH = (log(FORK_LENGTH) - FLcenter) / FLscale)
coefM4fl <- data.frame()
for(i in 1:nrow(coefM4)) {
  coefM4fl <- bind_rows(coefM4fl, bind_cols(fl, coefM4[i,]))
}
nrow(fl) * nrow(coefM4) == nrow(coefM4fl)  # check for sufficient rows
# 5. c) calculate conditional fecundity at fork lengths
names(coefM4fl)
predF_FL <- coefM4fl %>% 
  mutate(lEGG_COUNT_pred = b_Intercept + clFORK_LENGTH * b_clFORK_LENGTH +  
           case_when(STAGE_EGG_COUNT == "immature" ~ 0,
                     STAGE_EGG_COUNT == "mature"   ~ b_STAGE_EGG_COUNTmature,
                     TRUE                          ~ b_STAGE_EGG_COUNTmaturePretained)) %>% 
  group_by(STAGE_EGG_COUNT, FORK_LENGTH) %>% 
  # filter(lEGG_COUNT_pred > quantile(lEGG_COUNT_pred, probs = 0.025))
  summarise(EGG_COUNT_pred = round(exp(mean(lEGG_COUNT_pred)), 0),
            EGG_COUNT_lo = round(exp(quantile(lEGG_COUNT_pred, probs = 0.025)), 0),
            EGG_COUNT_hi = round(exp(quantile(lEGG_COUNT_pred, probs = 0.975)), 0)) %>% 
  mutate(across(contains("EGG"), ~./1000))

# 6. Figure S1
png("./Figures/FigS1.png", width = 7.5/2, height = 6/2, units = "in", res = 600)
ggplot(data = predF_FL, aes(x = FORK_LENGTH, y = EGG_COUNT_pred, 
                            col = STAGE_EGG_COUNT, fill = STAGE_EGG_COUNT)) +
  geom_ribbon(aes(ymin = EGG_COUNT_lo, ymax = EGG_COUNT_hi), 
              col = "transparent", alpha = 0.3) +
  geom_line() +
  scale_color_manual(values = pal3) + scale_fill_manual(values = pal3) +
  scale_x_continuous(breaks = seq(40, 110, 10)) +
  scale_y_continuous(breaks = seq(0, 25, 2)) +
  labs(x = "fork length (cm)", 
       y = "conditional fecundity ± 95% CIs\n(x 1,000 eggs)",
       color = NULL, fill = NULL) +
  theme(legend.position = c(0.25, 0.8))
dev.off()




# COVARIATES: MODELS ------------------------------------------------------

# 1. prepare dataset for covariate modes
# 1. a) filter missing values and centre/scale variables
fecundityCov <- fecundity %>% 
  filter(if_all(c(RAGE, LAST_SPAWN, FULTONS_K, SAMPLE_DOY), ~!is.na(.)),    # remove salmon missing covariates
         SAMPLE_DATE_RANGE <= 5) %>% droplevels() %>%                       # remove sample dates with range > 5 days
  # group_by(STAGE_EGG_COUNT) %>% 
  # mutate(    = scale(SAMPLE_DOY, center = TRUE, scale = TRUE)) %>%  # centre & scale sample date for immature vs mature
  # ungroup() %>% 
  mutate(cFULTONS_K     = scale(FULTONS_K, center = TRUE, scale = TRUE),    # centre & scale sample date
         cSAMPLE_DOY    = if_else(STAGE_EGG_COUNT == "immature",
                                  (SAMPLE_DOY - mean(SAMPLE_DOY[STAGE_EGG_COUNT == "immature"])) / 
                                     sd(SAMPLE_DOY[STAGE_EGG_COUNT == "immature"]),
                                  (SAMPLE_DOY - mean(SAMPLE_DOY[STAGE_EGG_COUNT != "immature"])) / 
                                    sd(SAMPLE_DOY[STAGE_EGG_COUNT != "immature"])),
         fRAGE          = as.factor(RAGE),                                  # set smolt ages as factor
         LAST_SPAWN     = fct_relevel(LAST_SPAWN, c(1, 2, 3, "C", "A")),    # re-order spawning strategies
         RIVER          = fct_reorder(RIVER, RIVER_NO))                     # re-order rivers
str(fecundityCov)
# 1. b) check sample sizes
count(fecundityCov, RIVER)                         # 20 rivers with min. 17 salmon 
count(fecundityCov, STAGE_EGG_COUNT, RAGE)         # n = 1 for one combination
count(fecundityCov, STAGE_EGG_COUNT, LAST_SPAWN)   # n = 2 for one combination
count(fecundityCov, RAGE, LAST_SPAWN)              # n < 3 for six combination
ggplot(fecundityCov, aes(x = cSAMPLE_DOY, fill = STAGE_EGG_COUNT)) +
  geom_histogram()
# 


# 2. candidate models
# 2. a) base model from previous analysis
M4S <- brm(lEGG_COUNT_COR ~ 0 + Intercept + clFORK_LENGTH + STAGE_EGG_COUNT + (clFORK_LENGTH | RIVER),
           prior  = priorList,
           family = gaussian(link = "identity"), 
           data   = fecundityCov, 
           chains = 4, 
           iter   = 4000,  # half discarded as warm-up
           seed   = 42)
M4S
plot(M4S)
# 2. b) add intercept for river age
M4Sa <- brm(lEGG_COUNT_COR ~ 0 + Intercept + clFORK_LENGTH + STAGE_EGG_COUNT + (clFORK_LENGTH | RIVER) + fRAGE,
            prior  = priorList,
            family = gaussian(link = "identity"), 
            data   = fecundityCov, 
            chains = 4, 
            iter   = 4000,  # half discarded as warm-up
            seed   = 42)
M4Sa
plot(M4Sa)
# 2. c) add intercept for last spawn
M4Sb <- brm(lEGG_COUNT_COR ~ 0 + Intercept + clFORK_LENGTH + STAGE_EGG_COUNT + (clFORK_LENGTH | RIVER) + LAST_SPAWN,
            prior  = priorList,
            family = gaussian(link = "identity"), 
            data   = fecundityCov, 
            chains = 4, 
            iter   = 4000,  # half discarded as warm-up
            seed   = 42)
M4Sb
plot(M4Sb)
# 2. d) add slope for body condition
M4Sc <- brm(lEGG_COUNT_COR ~ 0 + Intercept + clFORK_LENGTH + STAGE_EGG_COUNT + (clFORK_LENGTH | RIVER) + cFULTONS_K,
            prior  = priorList,
            family = gaussian(link = "identity"), 
            data   = fecundityCov, 
            chains = 4, 
            iter   = 4000,  # half discarded as warm-up
            seed   = 42)
M4Sc
plot(M4Sc)
# 2. e) add slope for sample date with STAGE_EGG_COUNT 
M4Sd <- brm(lEGG_COUNT_COR ~ 0 + Intercept + clFORK_LENGTH + STAGE_EGG_COUNT + (clFORK_LENGTH | RIVER) + 
              cSAMPLE_DOY + cSAMPLE_DOY : STAGE_EGG_COUNT,
            prior  = priorList,
            family = gaussian(link = "identity"), 
            data   = fecundityCov, 
            chains = 4, 
            iter   = 4000,  # half discarded as warm-up
            seed   = 42)
M4Sd
plot(M4Sd)


# 3. leave-one-out cross-validation for model comparison
# 3. a) re-run all models 
looM4S  <- loo(M4S)
looM4Sa <- loo(M4Sa)
looM4Sb <- loo(M4Sb)
looM4Sc <- loo(M4Sc)
looM4Sd <- loo(M4Sd)
# 3. b) initial comparison to determine if covariates improve the model
loo_compare(looM4S, looM4Sa, looM4Sb, looM4Sc, looM4Sd)


# 4. model with all covariates
M4Se <- brm(lEGG_COUNT_COR ~ 0 + Intercept + clFORK_LENGTH + STAGE_EGG_COUNT + (clFORK_LENGTH | RIVER) + 
              fRAGE + LAST_SPAWN + cFULTONS_K + cSAMPLE_DOY +  + cSAMPLE_DOY : STAGE_EGG_COUNT,
            prior  = priorList,
            family = gaussian(link = "identity"), 
            data   = fecundityCov,
            chains = 4, 
            iter   = 4000,  # half discarded as warm-up
            seed   = 42)
M4Se
plot(M4Se)


# 5. final leave-one-out cross-validation for model comparison
# 5. a) re-run latest model
looM4Se <- loo(M4Se)
# 5. b) additional comparison model
loo_compare(looM4S, looM4Sa, looM4Sb, looM4Sc, looM4Sd, looM4Se)


# 6. select top model for subsequent summaries and plots
# 6. a) identify top model
topCovModel <- M4Se
# 6. b) save models as RDS file
saveRDS(M4S,   "./models/M4S.rds")
saveRDS(M4Sa,  "./models/M4Sa.rds")
saveRDS(M4Sb,  "./models/M4Sb.rds")
saveRDS(M4Sc,  "./models/M4Sc.rds")
saveRDS(M4Sd,  "./models/M4Sd.rds")
saveRDS(M4Se,  "./models/M4Se.rds")
# 6. c) read models from RDS file
M4S  <- readRDS("./models/M4S.rds")
M4Sa <- readRDS("./models/M4Sa.rds")
M4Sb <- readRDS("./models/M4Sb.rds")
M4Sc <- readRDS("./models/M4Sc.rds")
M4Sd <- readRDS("./models/M4Sd.rds")
M4Se <- readRDS("./models/M4Se.rds")



# COVARIATES: INTERPRETATION ----------------------------------------------

# 1. set fixed values for plotting
# 1. a) centering and scaling parameters used prior to modelling
Kcenter <- round(attr(fecundityCov$cFULTONS_K, "scaled:center"), 2)
Kscale  <- round(attr(fecundityCov$cFULTONS_K, "scaled:scale"), 2)
DOYcenterI <- round(mean(fecundityCov$SAMPLE_DOY[fecundityCov$STAGE_EGG_COUNT == "immature"]), 1)
DOYscaleI  <- round(sd(fecundityCov$SAMPLE_DOY[fecundityCov$STAGE_EGG_COUNT == "immature"]), 1)
DOYcenterM <- round(mean(fecundityCov$SAMPLE_DOY[fecundityCov$STAGE_EGG_COUNT != "immature"]), 1)
DOYscaleM  <- round(sd(fecundityCov$SAMPLE_DOY[fecundityCov$STAGE_EGG_COUNT != "immature"]), 1)
# 1. b) median condition for size classes
fixedK <- bind_cols(q     = c("FULTONS_K", "FULTONS_K70"), 
                    value = quantile(fecundityCov$FULTONS_K, probs = c(0.5, 0.7))) %>%    # median & 70% quantiles      
  pivot_wider(names_from = q, values_from = "value") %>% 
  mutate(across(everything(), ~(. - Kcenter) / Kscale, .names = "c{.col}"))               # median FL centered & scaled
# 1. c) median sample dates
fixedDOY <- fecundityCov %>% 
  mutate(STAGE_SIMPLE = if_else(STAGE_EGG_COUNT == "immature", "immature", "mature")) %>%   # simplify STAGE_EGG_COUNT
  group_by(STAGE_SIMPLE) %>% 
  mutate(SAMPLE_DOY = round(median(SAMPLE_DOY), 0)) %>%                                     # median DOY rounded to the nearest date
  ungroup() %>% distinct(STAGE_EGG_COUNT, SAMPLE_DOY) %>% 
  mutate(SAMPLE_DOY2wk = SAMPLE_DOY + 14,
         across(starts_with("SAMPLE"),                         # median DOY centered & scaled
                ~if_else(STAGE_EGG_COUNT == "immature",
                         (. - DOYcenterI) / DOYscaleI,
                         (. - DOYcenterM) / DOYscaleM),
                .names = "c{.col}"))
# 2. model coefficients
sumCovModel <- posterior_summary(topCovModel) %>% data.frame() %>% 
  rownames_to_column("B") %>% filter(grepl("b_", B))


# 3. differences in fecundity by river age
# 3. a) 1SW predicted fecundity at each smolt age
fRAGE1SW2 <- sumCovModel$Estimate[1] +
  sumCovModel$Estimate[sumCovModel$B == "b_clFORK_LENGTH"] * FLmed$clFORK_LENGTH[FLmed$FORK_LENGTH == 52] +
  sumCovModel$Estimate[sumCovModel$B == "b_cFULTONS_K"]    * fixedK$cFULTONS_K +
  sumCovModel$Estimate[sumCovModel$B == "b_cSAMPLE_DOY"]   * fixedDOY$cSAMPLE_DOY[fixedDOY$STAGE_EGG_COUNT == "immature"]
fRAGE1SW4 <- fRAGE1SW2 + sumCovModel$Estimate[sumCovModel$B == "b_fRAGE4"]
fRAGE1SW6 <- fRAGE1SW2 + sumCovModel$Estimate[sumCovModel$B == "b_fRAGE6"]
# 3. b) 1SW differences in fecundity between smolt ages
round(exp(fRAGE1SW4) - exp(fRAGE1SW2), 0)
round((exp(fRAGE1SW4) - exp(fRAGE1SW2)) / exp(fRAGE1SW2) * 100, 1)
round(exp(fRAGE1SW6) - exp(fRAGE1SW4), 0)
round((exp(fRAGE1SW6) - exp(fRAGE1SW4)) / exp(fRAGE1SW4) * 100, 1)
# # 3. c) 2SW predicted fecundity at each smolt age
# fRAGE2SW2 <- sumCovModel$Estimate[1] +
#   sumCovModel$Estimate[sumCovModel$B == "b_LAST_SPAWN2"] +
#   sumCovModel$Estimate[sumCovModel$B == "b_clFORK_LENGTH"] * FLmed$clFORK_LENGTH[FLmed$FORK_LENGTH == 76] +
#   sumCovModel$Estimate[sumCovModel$B == "b_cFULTONS_K"]    * fixedK$cFULTONS_K +
#   sumCovModel$Estimate[sumCovModel$B == "b_cSAMPLE_DOY"]   * fixedDOY$cSAMPLE_DOY[fixedDOY$STAGE_EGG_COUNT == "immature"]
# fRAGE2SW4 <- fRAGE2 + sumCovModel$Estimate[sumCovModel$B == "b_fRAGE4"]
# fRAGE2SW6 <- fRAGE2 + sumCovModel$Estimate[sumCovModel$B == "b_fRAGE6"]
# # 3. d) 2SW differences in fecundity between smolt ages
# round(exp(fRAGE4) - exp(fRAGE2), 0)
# round(exp(fRAGE6) - exp(fRAGE4), 0)
# round((exp(fRAGE4) - exp(fRAGE2)) / exp(fRAGE2) * 100, 1)
# round((exp(fRAGE6) - exp(fRAGE4)) / exp(fRAGE4) * 100, 1)


# 4. differences in fecundity by spawning strategy
# 4. a) overlap in body size across spawning strategies
group_by(fecundityCov, LAST_SPAWN) %>% 
  reframe(paste(round(quantile(FORK_LENGTH, probs = 0.005), 0),
                round(quantile(FORK_LENGTH, probs = 0.995), 0), sep = "-"))
# 4. b) predicted fecundity for 1SW and C
SS1 <- sumCovModel$Estimate[1] +
  sumCovModel$Estimate[sumCovModel$B == "b_fRAGE3"] +
  sumCovModel$Estimate[sumCovModel$B == "b_clFORK_LENGTH"] * FLmed$clFORK_LENGTH[FLmed$FORK_LENGTH == 52] +
  sumCovModel$Estimate[sumCovModel$B == "b_cFULTONS_K"]    * fixedK$cFULTONS_K +
  sumCovModel$Estimate[sumCovModel$B == "b_cSAMPLE_DOY"]   * fixedDOY$cSAMPLE_DOY[fixedDOY$STAGE_EGG_COUNT == "immature"]
SS1C <- SS1 + sumCovModel$Estimate[sumCovModel$B == "b_LAST_SPAWNC"]
# 4. c) differences in fecundity between 1SW and C
round(exp(SS1C) - exp(SS1), 0)
round((exp(SS1C) - exp(SS1)) / exp(SS1) * 100, 1)
# 4. d) predicted fecundity for C and A
SS2C <- sumCovModel$Estimate[1] +
  sumCovModel$Estimate[sumCovModel$B == "b_LAST_SPAWNC"] +
  sumCovModel$Estimate[sumCovModel$B == "b_fRAGE3"] +
  sumCovModel$Estimate[sumCovModel$B == "b_clFORK_LENGTH"] * FLmed$clFORK_LENGTH[FLmed$FORK_LENGTH == 76] +
  sumCovModel$Estimate[sumCovModel$B == "b_cFULTONS_K"]    * fixedK$cFULTONS_K +
  sumCovModel$Estimate[sumCovModel$B == "b_cSAMPLE_DOY"]   * fixedDOY$cSAMPLE_DOY[fixedDOY$STAGE_EGG_COUNT == "immature"]
SS2A <- sumCovModel$Estimate[1] +
  sumCovModel$Estimate[sumCovModel$B == "b_LAST_SPAWNA"] +
  sumCovModel$Estimate[sumCovModel$B == "b_fRAGE3"] +
  sumCovModel$Estimate[sumCovModel$B == "b_clFORK_LENGTH"] * FLmed$clFORK_LENGTH[FLmed$FORK_LENGTH == 76] +
  sumCovModel$Estimate[sumCovModel$B == "b_cFULTONS_K"]    * fixedK$cFULTONS_K +
  sumCovModel$Estimate[sumCovModel$B == "b_cSAMPLE_DOY"]   * fixedDOY$cSAMPLE_DOY[fixedDOY$STAGE_EGG_COUNT == "immature"]
# 4. e) differences in fecundity between C and A
round(exp(SS2A) - exp(SS2C), 0)
round((exp(SS2A) - exp(SS2C)) / exp(SS2C) * 100, 1)


# 5. differences in fecundity by body condition
# 5. a)  determine % increase between mean and 70th quantile condition
round((fixedK$FULTONS_K70 - fixedK$FULTONS_K) / fixedK$FULTONS_K * 100, 1)
# 5. b) 1SW calculate fecundity for different body conditions
fK1SW50 <- sumCovModel$Estimate[sumCovModel$B == "b_Intercept"] +
  sumCovModel$Estimate[sumCovModel$B == "b_fRAGE3"] +
  sumCovModel$Estimate[sumCovModel$B == "b_clFORK_LENGTH"] * FLmed$clFORK_LENGTH[FLmed$FORK_LENGTH == 52] +
  sumCovModel$Estimate[sumCovModel$B == "b_cFULTONS_K"]    * fixedK$cFULTONS_K +
  sumCovModel$Estimate[sumCovModel$B == "b_cSAMPLE_DOY"]   * fixedDOY$cSAMPLE_DOY[fixedDOY$STAGE_EGG_COUNT == "immature"]
fK1SW70 <- sumCovModel$Estimate[sumCovModel$B == "b_Intercept"] +
  sumCovModel$Estimate[sumCovModel$B == "b_fRAGE3"] +
  sumCovModel$Estimate[sumCovModel$B == "b_clFORK_LENGTH"] * FLmed$clFORK_LENGTH[FLmed$FORK_LENGTH == 52] +
  sumCovModel$Estimate[sumCovModel$B == "b_cFULTONS_K"]    * fixedK$cFULTONS_K70 +
  sumCovModel$Estimate[sumCovModel$B == "b_cSAMPLE_DOY"]   * fixedDOY$cSAMPLE_DOY[fixedDOY$STAGE_EGG_COUNT == "immature"]
# 5. c) 1SW calculate difference in condition
round(exp(fK1SW70) - exp(fK1SW50), 0)
round((exp(fK1SW70) - exp(fK1SW50)) / exp(fK1SW50) * 100, 1)
# 5. d) 2SW calculate fecundity for different body conditions
fK2SW50 <- sumCovModel$Estimate[sumCovModel$B == "b_Intercept"] +
  sumCovModel$Estimate[sumCovModel$B == "b_LAST_SPAWN2"] +
  sumCovModel$Estimate[sumCovModel$B == "b_fRAGE3"] +
  sumCovModel$Estimate[sumCovModel$B == "b_clFORK_LENGTH"] * FLmed$clFORK_LENGTH[FLmed$FORK_LENGTH == 76] +
  sumCovModel$Estimate[sumCovModel$B == "b_cFULTONS_K"]    * fixedK$cFULTONS_K +
  sumCovModel$Estimate[sumCovModel$B == "b_cSAMPLE_DOY"]   * fixedDOY$cSAMPLE_DOY[fixedDOY$STAGE_EGG_COUNT == "immature"]
fK2SW70 <- sumCovModel$Estimate[sumCovModel$B == "b_Intercept"] +
  sumCovModel$Estimate[sumCovModel$B == "b_LAST_SPAWN2"] +
  sumCovModel$Estimate[sumCovModel$B == "b_fRAGE3"] +
  sumCovModel$Estimate[sumCovModel$B == "b_clFORK_LENGTH"] * FLmed$clFORK_LENGTH[FLmed$FORK_LENGTH == 76] +
  sumCovModel$Estimate[sumCovModel$B == "b_cFULTONS_K"]    * fixedK$cFULTONS_K70 +
  sumCovModel$Estimate[sumCovModel$B == "b_cSAMPLE_DOY"]   * fixedDOY$cSAMPLE_DOY[fixedDOY$STAGE_EGG_COUNT == "immature"]
# 5. d) 2SW calculate difference in condition
round(exp(fK2SW70) - exp(fK2SW50), 0)
round((exp(fK2SW70) - exp(fK2SW50)) / exp(fK2SW50) * 100, 1)


# 6. differences in fecundity by sample date - IMMATURE
# 6. a) calculate fecundity for different sample dates
fSD1SW <- sumCovModel$Estimate[sumCovModel$B == "b_Intercept"] +
  sumCovModel$Estimate[sumCovModel$B == "b_fRAGE3"] +
  sumCovModel$Estimate[sumCovModel$B == "b_clFORK_LENGTH"] * FLmed$clFORK_LENGTH[FLmed$FORK_LENGTH == 52] +
  sumCovModel$Estimate[sumCovModel$B == "b_cFULTONS_K"]    * fixedK$cFULTONS_K +
  sumCovModel$Estimate[sumCovModel$B == "b_cSAMPLE_DOY"]   * fixedDOY$cSAMPLE_DOY[fixedDOY$STAGE_EGG_COUNT == "immature"]
fSD1SW14 <- sumCovModel$Estimate[sumCovModel$B == "b_Intercept"] +
  sumCovModel$Estimate[sumCovModel$B == "b_fRAGE3"] +
  sumCovModel$Estimate[sumCovModel$B == "b_clFORK_LENGTH"] * FLmed$clFORK_LENGTH[FLmed$FORK_LENGTH == 52] +
  sumCovModel$Estimate[sumCovModel$B == "b_cFULTONS_K"]    * fixedK$cFULTONS_K +
  sumCovModel$Estimate[sumCovModel$B == "b_cSAMPLE_DOY"]   * fixedDOY$cSAMPLE_DOY2wk[fixedDOY$STAGE_EGG_COUNT == "immature"]
# 6. b) calculate difference in condition
round(exp(fSD1SW14) - exp(fSD1SW), 0)
round((exp(fSD1SW14) - exp(fSD1SW)) / exp(fSD1SW) * 100, 1)
# 6. c) calculate fecundity for different sample dates
fSD2SW <- sumCovModel$Estimate[sumCovModel$B == "b_Intercept"] +
  sumCovModel$Estimate[sumCovModel$B == "b_LAST_SPAWN2"] +
  sumCovModel$Estimate[sumCovModel$B == "b_fRAGE3"] +
  sumCovModel$Estimate[sumCovModel$B == "b_clFORK_LENGTH"] * FLmed$clFORK_LENGTH[FLmed$FORK_LENGTH == 76] +
  sumCovModel$Estimate[sumCovModel$B == "b_cFULTONS_K"]    * fixedK$cFULTONS_K +
  sumCovModel$Estimate[sumCovModel$B == "b_cSAMPLE_DOY"]   * fixedDOY$cSAMPLE_DOY[fixedDOY$STAGE_EGG_COUNT == "immature"]
fSD2SW14 <- sumCovModel$Estimate[sumCovModel$B == "b_Intercept"] +
  sumCovModel$Estimate[sumCovModel$B == "b_LAST_SPAWN2"] +
  sumCovModel$Estimate[sumCovModel$B == "b_fRAGE3"] +
  sumCovModel$Estimate[sumCovModel$B == "b_clFORK_LENGTH"] * FLmed$clFORK_LENGTH[FLmed$FORK_LENGTH == 76] +
  sumCovModel$Estimate[sumCovModel$B == "b_cFULTONS_K"]    * fixedK$cFULTONS_K +
  sumCovModel$Estimate[sumCovModel$B == "b_cSAMPLE_DOY"]   * fixedDOY$cSAMPLE_DOY2wk[fixedDOY$STAGE_EGG_COUNT == "immature"]
# 6. c) calculate difference in condition
round(exp(fSD2SW14) - exp(fSD2SW), 0)
round((exp(fSD2SW14) - exp(fSD2SW)) / exp(fSD2SW) * 100, 1)


# 7. differences in fecundity by sample date - MATURE
# 7. a) calculate fecundity for different sample dates
fSD1SWm <- sumCovModel$Estimate[sumCovModel$B == "b_Intercept"] +
  sumCovModel$Estimate[sumCovModel$B == "b_STAGE_EGG_COUNTmature"] +
  sumCovModel$Estimate[sumCovModel$B == "b_fRAGE3"] +
  sumCovModel$Estimate[sumCovModel$B == "b_clFORK_LENGTH"] * FLmed$clFORK_LENGTH[FLmed$FORK_LENGTH == 52] +
  sumCovModel$Estimate[sumCovModel$B == "b_cFULTONS_K"]    * fixedK$cFULTONS_K +
  sumCovModel$Estimate[sumCovModel$B == "b_cSAMPLE_DOY"]   * fixedDOY$cSAMPLE_DOY[fixedDOY$STAGE_EGG_COUNT == "mature"] +
  sumCovModel$Estimate[sumCovModel$B == "b_STAGE_EGG_COUNTmature:cSAMPLE_DOY"] * 
  fixedDOY$cSAMPLE_DOY[fixedDOY$STAGE_EGG_COUNT == "mature"]
fSD1SW14m <- sumCovModel$Estimate[sumCovModel$B == "b_Intercept"] +
  sumCovModel$Estimate[sumCovModel$B == "b_STAGE_EGG_COUNTmature"] +
  sumCovModel$Estimate[sumCovModel$B == "b_fRAGE3"] +
  sumCovModel$Estimate[sumCovModel$B == "b_clFORK_LENGTH"] * FLmed$clFORK_LENGTH[FLmed$FORK_LENGTH == 52] +
  sumCovModel$Estimate[sumCovModel$B == "b_cFULTONS_K"]    * fixedK$cFULTONS_K +
  sumCovModel$Estimate[sumCovModel$B == "b_cSAMPLE_DOY"]   * fixedDOY$cSAMPLE_DOY[fixedDOY$STAGE_EGG_COUNT == "mature"] +
  sumCovModel$Estimate[sumCovModel$B == "b_STAGE_EGG_COUNTmature:cSAMPLE_DOY"] * 
  fixedDOY$cSAMPLE_DOY2wk[fixedDOY$STAGE_EGG_COUNT == "mature"]
# 7. b) calculate difference in condition
round(exp(fSD1SW14m) - exp(fSD1SWm), 0)
round((exp(fSD1SW14m) - exp(fSD1SWm)) / exp(fSD1SWm) * 100, 1)
# 7. c) calculate fecundity for different sample dates
fSD2SWm <- sumCovModel$Estimate[sumCovModel$B == "b_Intercept"] +
  sumCovModel$Estimate[sumCovModel$B == "b_STAGE_EGG_COUNTmature"] +
  sumCovModel$Estimate[sumCovModel$B == "b_LAST_SPAWN2"] +
  sumCovModel$Estimate[sumCovModel$B == "b_fRAGE3"] +
  sumCovModel$Estimate[sumCovModel$B == "b_clFORK_LENGTH"] * FLmed$clFORK_LENGTH[FLmed$FORK_LENGTH == 76] +
  sumCovModel$Estimate[sumCovModel$B == "b_cFULTONS_K"]    * fixedK$cFULTONS_K +
  sumCovModel$Estimate[sumCovModel$B == "b_cSAMPLE_DOY"]   * fixedDOY$cSAMPLE_DOY[fixedDOY$STAGE_EGG_COUNT == "mature"] +
  sumCovModel$Estimate[sumCovModel$B == "b_STAGE_EGG_COUNTmature:cSAMPLE_DOY"] * 
  fixedDOY$cSAMPLE_DOY[fixedDOY$STAGE_EGG_COUNT == "mature"]
fSD2SW14m <- sumCovModel$Estimate[sumCovModel$B == "b_Intercept"] +
  sumCovModel$Estimate[sumCovModel$B == "b_STAGE_EGG_COUNTmature"] +
  sumCovModel$Estimate[sumCovModel$B == "b_LAST_SPAWN2"] +
  sumCovModel$Estimate[sumCovModel$B == "b_fRAGE3"] +
  sumCovModel$Estimate[sumCovModel$B == "b_clFORK_LENGTH"] * FLmed$clFORK_LENGTH[FLmed$FORK_LENGTH == 76] +
  sumCovModel$Estimate[sumCovModel$B == "b_cFULTONS_K"]    * fixedK$cFULTONS_K +
  sumCovModel$Estimate[sumCovModel$B == "b_cSAMPLE_DOY"]   * fixedDOY$cSAMPLE_DOY[fixedDOY$STAGE_EGG_COUNT == "mature"] +
  sumCovModel$Estimate[sumCovModel$B == "b_STAGE_EGG_COUNTmature:cSAMPLE_DOY"] * 
  fixedDOY$cSAMPLE_DOY2wk[fixedDOY$STAGE_EGG_COUNT == "mature"]
# 7. c) calculate difference in condition
round(exp(fSD2SW14m) - exp(fSD2SWm), 0)
round((exp(fSD2SW14m) - exp(fSD2SWm)) / exp(fSD2SWm) * 100, 1)



# 7. FIGURE 3 - conditional means
# conditional_effects(topCovModel, prob = 0.95, method = "posterior_epred", effects = "fRAGE")
# conditional_effects(topCovModel, prob = 0.95, method = "posterior_epred", effects = "LAST_SPAWN")
# conditional_effects(topCovModel, prob = 0.95, method = "posterior_epred", effects = "cFULTONS_K")#, spaghetti = TRUE)
# conditional_effects(topCovModel, prob = 0.95, method = "posterior_epred", effects = "cSAMPLE_DOY")
# 7. a) sample sizes for plots
nRAGE       <- count(fecundityCov, RAGE)
nLAST_SPAWN <- count(fecundityCov, LAST_SPAWN)
# 7. c) smolt age panel
pSmolt <- as_draws_df(topCovModel) %>% select(b_Intercept, b_fRAGE3:b_fRAGE6) %>% 
  mutate(across(c(b_fRAGE3:b_fRAGE6), ~. + b_Intercept )) %>%
  pivot_longer(cols = b_Intercept:b_fRAGE6, names_to = "age", values_to = "pred") %>% 
  mutate(age = if_else(age == "b_Intercept", "2", str_sub(age, -1, -1))) %>% 
  ggplot(aes(y = age, x = exp(pred))) + 
  stat_pointinterval(aes(x = age, y = exp(pred)),
                     point_interval = "median_qi", .width = c(0.025, 0.975), linewidth = 1) +
  geom_text(data = nRAGE, aes(y = 4600, x = RAGE-1, label = n), size = 3) +
  scale_y_continuous(breaks = seq(2500, 6000, 500)) +
  labs(y = "conditional median fecundity\n(95% CI)", x = "smolt age", tag = "A")
# 7. d) spawning strategy panel
pSpawn <- as_draws_df(topCovModel) %>% select(b_Intercept, b_LAST_SPAWN2:b_LAST_SPAWNA) %>% 
  mutate(across(c(b_LAST_SPAWN2:b_LAST_SPAWNA), ~. + b_Intercept )) %>%
  pivot_longer(cols = b_Intercept:b_LAST_SPAWNA, names_to = "SS", values_to = "pred") %>% 
  mutate(SS = case_when(SS == "b_Intercept"                  ~ "1SW",
                        str_sub(SS, -1, -1) %in% c("C", "A") ~ str_sub(SS, -1, -1),
                        TRUE ~ paste(str_sub(SS, -1, -1), "SW", sep = "")),
         SS = fct_relevel(SS, c("1SW", "2SW", "3SW", "C", "A"))) %>%
  ggplot(aes(x = SS, y = exp(pred))) + 
  stat_pointinterval(.width = 0.95, linewidth = 1) +
  geom_text(data = nLAST_SPAWN, 
            aes(y = 5500, x = if_else(LAST_SPAWN %in% c("C", "A"), LAST_SPAWN, paste0(LAST_SPAWN, "SW")), label = n),
            size = 3) +
  scale_y_continuous(breaks = seq(2500, 6000, 500)) +
  labs(y = "conditional median fecundity\n(95% CI)", 
       x = "spawning strategy", tag = "B")
# 7. e) condition panel
pK <- as_draws_df(topCovModel) %>% select(b_Intercept, b_cFULTONS_K) %>% 
  reframe(FULTONS_K = seq(min(fecundityCov$FULTONS_K), max(fecundityCov$FULTONS_K), 0.1),
          .by = c("b_Intercept", "b_cFULTONS_K")) %>% 
  mutate(cFULTONS_K = (FULTONS_K - Kcenter) / Kscale,
         EGG_COUNT_PRED = exp(b_Intercept + b_cFULTONS_K * cFULTONS_K)) %>%
  group_by(FULTONS_K, cFULTONS_K) %>%
  summarise(mEGG_COUNT_PRED = mean(EGG_COUNT_PRED),
            loEGG_COUNT_PRED = quantile(EGG_COUNT_PRED, prob = 0.025),
            hiEGG_COUNT_PRED = quantile(EGG_COUNT_PRED, prob = 0.975)) %>%
  ggplot(aes(x = FULTONS_K)) +
  geom_ribbon(aes(x = FULTONS_K, ymin = loEGG_COUNT_PRED, ymax = hiEGG_COUNT_PRED), fill = "grey80") +
  geom_line(aes(x = FULTONS_K, y = mEGG_COUNT_PRED)) +
  geom_rug(data = fecundityCov) +
  scale_x_continuous(breaks = seq(0, 2.6, 0.5)) +
  scale_y_continuous(breaks = seq(2500, 7000, 500)) +
  labs(y = "conditional median fecundity\n(95% CI)", 
       x = expression(paste("body condition (kg/cm"^3, ")", sep = "")), tag = "C")
# 7. f) sample date panel
# doyI <- 
pDOY <- as_draws_df(topCovModel) %>% 
  select(b_Intercept, b_cSAMPLE_DOY, `b_STAGE_EGG_COUNTmature:cSAMPLE_DOY`,
         `b_STAGE_EGG_COUNTmaturePretained:cSAMPLE_DOY`) %>% 
  reframe(SAMPLE_DOY = seq(min(fecundityCov$SAMPLE_DOY), max(fecundityCov$SAMPLE_DOY), 1),
        .by = c("b_Intercept", "b_cSAMPLE_DOY", "b_STAGE_EGG_COUNTmature:cSAMPLE_DOY",
                "b_STAGE_EGG_COUNTmaturePretained:cSAMPLE_DOY")) %>% 
  mutate(
    # immature
    cDOY_I   = if_else(SAMPLE_DOY %in% c(130:262), 
                                 (SAMPLE_DOY - DOYcenterI) / DOYscaleI, NA),
    immature = exp(b_Intercept + b_cSAMPLE_DOY * cDOY_I),
    # mature
    cDOY_M   = if_else(SAMPLE_DOY %in% c(283:313), 
                       (SAMPLE_DOY - DOYcenterM) / DOYscaleM, NA),
    mature   = exp(b_Intercept + b_cSAMPLE_DOY * cDOY_M + 
                     `b_STAGE_EGG_COUNTmature:cSAMPLE_DOY` * cDOY_M),
    # mature + retained
    cDOY_MR  = if_else(SAMPLE_DOY %in% c(292:317), 
                       (SAMPLE_DOY - DOYcenterM) / DOYscaleM, NA),
    `mature + retained` = exp(b_Intercept + b_cSAMPLE_DOY * cDOY_MR + 
                            `b_STAGE_EGG_COUNTmaturePretained:cSAMPLE_DOY` * cDOY_MR)) %>% 
  select(SAMPLE_DOY, immature, mature, `mature + retained`) %>% 
  pivot_longer(cols = 2:4, names_to = "STAGE_EGG_COUNT", values_to = "EGG_COUNT_PRED") %>% 
  filter(!is.na(EGG_COUNT_PRED)) %>% 
  group_by(STAGE_EGG_COUNT, SAMPLE_DOY) %>% 
  summarise(mEGG_COUNT_PRED  = mean(EGG_COUNT_PRED),
            loEGG_COUNT_PRED = quantile(EGG_COUNT_PRED, prob = 0.025),
            hiEGG_COUNT_PRED = quantile(EGG_COUNT_PRED, prob = 0.975)) %>% 
  ggplot(aes(x = SAMPLE_DOY, col = STAGE_EGG_COUNT, fill = STAGE_EGG_COUNT)) +
  geom_ribbon(aes(ymin = loEGG_COUNT_PRED, ymax = hiEGG_COUNT_PRED), 
              col = "transparent", alpha = 0.5) +
  geom_line(aes(y = mEGG_COUNT_PRED)) +
  scale_fill_manual(values = pal3) + scale_color_manual(values = pal3) +
  geom_rug(data = fecundityCov, col = "black") +
  scale_x_continuous(breaks = seq(130, 320, 30)) +
  scale_y_continuous(breaks = seq(2500, 6000, 500)) +
  labs(y = "conditional median fecundity\n(95% CI)", x = "sample date", tag = "D",
       col = NULL, fill = NULL) +
  theme(legend.position = c(0.6, 0.8))

# 7. b) multi-panelled plot
png("./figures/Fig3.png", width = 7.5, height = 6, units = "in", res = 600)
ggarrange(pSmolt, pSpawn, pK, pDOY)
dev.off()



# EGG DIAMETER: MODELS ----------------------------------------------------

# 1. prepare dataset
fecundityED <-  fecundity %>%
  filter(!is.na(EGG_DIAM_COR),                                          # remove missing egg diameters
         STAGE_EGG_COUNT != "immature",                                 # remove immature eggs
         !RIVER %in% c("LaHave River", "Little Salmonier River")) %>%   # remove river with < 5 obs
  group_by(RIVER) %>% 
  mutate(cEGG_DIAM_COR = scale(EGG_DIAM_COR)) %>% ungroup()   # center & scale egg diameter within rivers

count(fecundityED, RIVER)

# 2. egg diameter models
# 2. a) base model
M4ED <- brm(lEGG_COUNT_COR ~ 0 + Intercept + clFORK_LENGTH + STAGE_EGG_COUNT + (clFORK_LENGTH | RIVER),
            prior  = priorList,
            family = gaussian(link = "identity"), 
            data   = fecundityED,
            chains = 4, 
            iter   = 4000,  # half discarded as warm-up
            seed   = 42,
            control = list(adapt_delta = 0.99))
plot(M4ED)
summary(M4ED)
# 2. b) add slope for egg diameter
M4EDa <- brm(lEGG_COUNT_COR ~ 0 + Intercept + clFORK_LENGTH + STAGE_EGG_COUNT + (clFORK_LENGTH | RIVER) + 
               cEGG_DIAM_COR,
             prior  = priorList,
             family = gaussian(link = "identity"), 
             data   = fecundityED,
             chains = 4, 
             iter   = 4000,  # half discarded as warm-up
             seed   = 42,
             control = list(adapt_delta = 0.99))
plot(M4EDa)
prior_summary(M4EDa)


# 3. leave-one-out cross-validation for model comparison
# 3. a) re-run all models 
looED  <- loo(M4ED)
looEDa <- loo(M4EDa)
looEDb <- loo(M4EDb, moment_match = TRUE, reloo = TRUE)
# 3. b) initial comparison to determine if covariates improve the model
loo_compare(looED, looEDa)


# 4. select top model for subsequent summaries and plots
# 4. a) identify top model
topEDModel <- M4EDa
# 4. b) save models as RDS files
saveRDS(M4ED,  "./models/M4ED.rds")
saveRDS(M4EDa, "./models/M4EDa.rds")
# 4. c) read models from RDS files
M4ED  <- readRDS("./models/M4ED.rds")
M4EDa <- readRDS("./models/M4EDa.rds")



# EGG DIAMETER: INTERPRETATION ------------------------------------------

# 1. model coefficients
sumEDModel <- posterior_summary(topEDModel) %>% data.frame() %>% 
  rownames_to_column("B") %>% filter(grepl("b_", B))


# 2. conditional fecundity plot
png("./Figures/Fig4.png", width = 7.5/2, height = 5/2, units = "in", res = 600)
as_draws_df(topEDModel) %>% select(b_Intercept, b_cEGG_DIAM_COR) %>% 
  reframe(EGG_DIAM_COR = seq(min(fecundityED$EGG_DIAM_COR), max(fecundityED$EGG_DIAM_COR), 
                             length.out = 30),
          .by = c("b_Intercept", "b_cEGG_DIAM_COR")) %>% 
  mutate(cEGG_DIAM_COR  = (EGG_DIAM_COR - mean(fecundityED$EGG_DIAM_COR)) /
                             sd(fecundityED$EGG_DIAM_COR),
         EGG_COUNT_PRED = exp(b_Intercept + b_cEGG_DIAM_COR * cEGG_DIAM_COR)) %>% 
  group_by(EGG_DIAM_COR) %>% 
  summarise(mEGG_COUNT_PRED = mean(EGG_COUNT_PRED),
            loEGG_COUNT_PRED = quantile(EGG_COUNT_PRED, prob = 0.025),
            hiEGG_COUNT_PRED = quantile(EGG_COUNT_PRED, prob = 0.975)) %>%
  ggplot(aes(x = EGG_DIAM_COR)) +
  geom_ribbon(aes(ymin = loEGG_COUNT_PRED, ymax = hiEGG_COUNT_PRED), fill = "grey80") +
  geom_line(aes(y = mEGG_COUNT_PRED)) +
  geom_rug(data = fecundityED) +
  scale_x_continuous(breaks = seq(0, 11, 1)) +
  scale_y_continuous(breaks = seq(2000, 5000, 500)) +
  labs(y = "conditional mean fecundity\n(95% CI)", 
       x = "egg diameter (mm)")
dev.off()


# 3. differences in fecundity by egg diameter for small salmon
# 3. a) determine range of egg diameters at different quantiles
iqrED  <- round(quantile(fecundityED$EGG_DIAM_COR, probs = c(0.5, 0.7)), 1) 
# 1. a) centering and scaling parameters used prior to modelling
EDcenter <- round(mean(fecundityED$EGG_DIAM_COR), 2)
EDscale  <- round(sd(fecundityED$EGG_DIAM_COR), 2)
# 1. b) median condition for size classes
fixedK <- bind_cols(q     = c("FULTONS_K", "FULTONS_K70"), 
                    value = quantile(fecundityCov$FULTONS_K, probs = c(0.5, 0.7))) %>%    # median & 70% quantiles      
  pivot_wider(names_from = q, values_from = "value") %>% 
  mutate(across(everything(), ~(. - Kcenter) / Kscale, .names = "c{.col}"))               # median FL centered & scaled

# 3. b) predicted fecundity at varying egg diameters
smED50 <- sumEDModel$Estimate[sumEDModel$B == "b_Intercept"] +
  sumEDModel$Estimate[sumEDModel$B == "b_clFORK_LENGTH"] * FLmed$clFORK_LENGTH[FLmed$FORK_LENGTH == 52] +
  sumEDModel$Estimate[sumEDModel$B == "b_cEGG_DIAM_COR"] * ((iqrED[1] - EDcenter) / EDscale)
smED70 <- sumEDModel$Estimate[sumEDModel$B == "b_Intercept"] +
  sumEDModel$Estimate[sumEDModel$B == "b_clFORK_LENGTH"] * FLmed$clFORK_LENGTH[FLmed$FORK_LENGTH == 52] +
  sumEDModel$Estimate[sumEDModel$B == "b_cEGG_DIAM_COR"] * ((iqrED[2] - EDcenter) / EDscale)
# 3. c) differences in fecundity between smolt ages
round(exp(smED70) - exp(smED50), 0)
round((exp(smED70) - exp(smED50)) / exp(smED50) * 100, 1)
# 3. d) predicted fecundity at varying egg diameters
smED50 <- sumEDModel$Estimate[sumEDModel$B == "b_Intercept"] +
  sumEDModel$Estimate[sumEDModel$B == "b_clFORK_LENGTH"] * FLmed$clFORK_LENGTH[FLmed$FORK_LENGTH == 76] +
  sumEDModel$Estimate[sumEDModel$B == "b_cEGG_DIAM_COR"] * ((iqrED[1] - EDcenter) / EDscale)
smED70 <- sumEDModel$Estimate[sumEDModel$B == "b_Intercept"] +
  sumEDModel$Estimate[sumEDModel$B == "b_clFORK_LENGTH"] * FLmed$clFORK_LENGTH[FLmed$FORK_LENGTH == 76] +
  sumEDModel$Estimate[sumEDModel$B == "b_cEGG_DIAM_COR"] * ((iqrED[2] - EDcenter) / EDscale)
# 3. c) differences in fecundity between smolt ages
round(exp(smED70) - exp(smED50), 0)
round((exp(smED70) - exp(smED50)) / exp(smED50) * 100, 1)




# RESIDUAL FECUNDITY: MODELS ----------------------------------------------


# 1. determine number of rives with fecundity data over at least two decades
# 1. a) calculate residual fecundity
summary(fecundityRes)
fecundityRes <- fecundity %>% #filter(fecundity, SAMPLE_DATE_RANGE <= 5) %>%   # remove data without sample date
  mutate(DECADE = round_any(YEAR, 10, f = floor)) %>%           # calculate decade
  dplyr::count(RIVER, DECADE) %>%                               # determine rivers with min. 15 obs/decade 
  filter(n >= 15) %>% 
  distinct(RIVER, DECADE) %>% 
  dplyr::count(RIVER, name = "nD") %>% filter(nD > 1) %>%       # determine rivers with data in 2+ decades
  left_join(fecundity, by = "RIVER") %>% droplevels() %>%       # bring back original data
  add_epred_draws(M4) %>%                                       # add draws from posterior distribution
  group_by(RIVER_NO, RIVER, YEAR, EGG_COUNT_COR, lEGG_COUNT_COR, SIZE,
           FORK_LENGTH, clFORK_LENGTH, SAMPLE_DOY) %>% 
  dplyr::summarise(lEGGS_PRED = median(.epred),
                   EGGS_PRED  = exp(lEGGS_PRED)) %>% ungroup() %>%    # calculate median predicted fecundity
  mutate(EGG_COUNT_RES  = EGGS_PRED - EGG_COUNT_COR,                  # calculate residual fecundity
         fYEAR          = factor(YEAR, levels = seq(min(YEAR), max(YEAR), 1)),
         DECADE         = factor(plyr::round_any(YEAR, 10, f = floor),
                                 levels = c("1980", "1990", "2000", "2010", "2020"))) %>% 
  droplevels()
# 1. b) min/max years included
summarise(fecundityRes, min(YEAR), max(YEAR))
# 1. c) number of rivers
dplyr::count(fecundityRes, RIVER, DECADE) %>% pivot_wider(values_from = n, names_from = DECADE)


# 2. residual fecundity model
# 2. a) base model
R0 <- brm(EGG_COUNT_RES ~ 0 + Intercept + (1 | RIVER),
          prior  = filter(priorList, class != "cor"),
          family = gaussian(link = "identity"), 
          data   = fecundityRes,
          chains = 4, 
          iter   = 2000,  # half discarded as warm-up
          seed   = 41,
          save_pars = save_pars(all = TRUE))
R0
plot(R0)
# 2. b) annual variation
R1 <- brm(EGG_COUNT_RES ~ 0 + Intercept + fYEAR + (1 | RIVER),
          prior  = filter(priorList, class != "cor"),
          family = gaussian(link = "identity"), 
          data   = fecundityRes,
          chains = 4, 
          iter   = 2000,  # half discarded as warm-up
          seed   = 41,
          control = list(adapt_delta = 0.85),
          save_pars = save_pars(all = TRUE))
R1
plot(R1)
# 2. b) decadal and annual variation
R2 <- brm(EGG_COUNT_RES ~ 0 + Intercept + DECADE + (1 | fYEAR) + (1 | RIVER),
          prior  = filter(priorList, class != "cor"),
          family = gaussian(link = "identity"), 
          data   = fecundityRes,
          chains = 4, 
          iter   = 2000,  # half discarded as warm-up
          seed   = 41,
          control = list(adapt_delta = 0.85),
          save_pars = save_pars(all = TRUE))
R2
plot(R2)


# 3. leave-one-out cross-validation for model comparison
# 3. a) re-run all models 
looR0 <- loo(R0, moment_match = TRUE)
looR1 <- loo(R1, moment_match = TRUE), reloo = TRUE)
looR2 <- loo(R2, moment_match = TRUE)
# looR3 <- loo(R3)
# 3. b) initial comparison to determine if covariates improve the model
loo_compare(looR0, looR1, looR2) #, looR3


# 4. select top model for subsequent summaries and plots
# 4. a) identify top model
topRFModel <- R2
# 4. b) save models as RDS files
saveRDS(R0, "./Analysis/R0.rds")
saveRDS(R1, "./Analysis/R1.rds")
saveRDS(R2, "./Analysis/R2.rds")
# 4. c) read models from RDS files
R0 <- readRDS("./Analysis/R0.rds")
R1 <- readRDS("./Analysis/R1.rds")
R2 <- readRDS("./Analysis/R2.rds")



# RESIDUAL FECUNDITY: INTERPRETATION --------------------------------------

# predicted values
predRes <- distinct(fecundityRes) %>% 
  add_epred_draws(topRFModel) %>%                         # draw from posterior distribution
  group_by(.row) %>% 
  filter(.epred >= quantile(.epred, probs = 0.025),       # filter <2.5% quantile
         .epred <= quantile(.epred, probs = 0.975)) %>%   # filter >97.5% quantile
  ungroup()


# Figure 5
png("./Figures/FigS2.png", width = 7.5, height = 5, units = "in", res = 600)
ggplot(data = predRes, aes(x = YEAR, y = .epred, col = DECADE)) +
  geom_hline(yintercept = 0, col = "grey90") +
  geom_boxplot(aes(group = YEAR)) +
  scale_color_viridis_d(option = "plasma") +
  facet_wrap(~fct_reorder(RIVER, RIVER_NO), ncol = 4) +
  scale_x_continuous(breaks = seq(1980, 2025, 5)) +
  scale_y_continuous(breaks = seq(-200, 300, 100)) +
  labs(y = "predicted residual fecudity", x = NULL, col = "decade") +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position  = c(0.75, 0.1),
        legend.direction = "horizontal") +
  guides(colour = guide_legend(nrow = 2))
dev.off()



# FIGURE 1 ----------------------------------------------------------------

# 1. obtain basemap
basemap <- rnaturalearth::ne_states("canada")

  
# 2. plot locations of rivers and DUs
png("./Figures/Fig1.png", width = 4, height = 3.5, units = "in", res = 600)
distinct(fecundityRaw, DU, RIVER_NO, RIVER, LATITUDE, LONGITUDE) %>% 
  arrange(RIVER_NO) %>%
  ggplot() + 
  geom_sf(data = basemap, col = "grey75", fill = "grey95") +
  geom_point(aes(x = LONGITUDE, y = LATITUDE, fill = DU, shape = DU), 
             size = 1.5, alpha = 0.8) +
  geom_text_repel(aes(x = LONGITUDE, y = LATITUDE, label = RIVER_NO), 
                  direction = "both", max.iter = 20000, size = 2, 
                  min.segment.length = 0.2) + # max.overlaps = 20
  scale_fill_manual(values = c(pal11[1], pal11[3], pal11[5], pal11[7], pal11[9], 
                               pal11[11], pal11[2], pal11[4], pal11[6], pal11[8], 
                               pal11[10])) +
  scale_shape_manual(values = c(21, rep(22, 5), rep(21, 5))) +
  coord_sf(xlim = c(-70, -50), ylim = c(43, 54)) +
  labs(y = NULL, x = NULL, fill = NULL) + guides(fill = "none", shape = "none")
dev.off()


# 3. citation for rnaturalearth
citation("rnaturalearth")



# TABLE 2 -----------------------------------------------------------------

# 4. create model comparison table (Table 2)
# 4. a) first set of models
set1 <- data.frame(loo_compare(looM0, looM1, looM2, looM3, looM4)) %>% 
  rownames_to_column("model") %>% arrange(model) %>% 
  mutate(across(where(is.numeric), ~round(., 2))) %>% 
  add_row(.before = 1, 
          model = "Hierarchical fecundity-fork length relationships (25 rivers)")
# 4. b) second set of models (covariates)
set2 <- data.frame(loo_compare(looM4S, looM4Sa, looM4Sb, looM4Sc, looM4Sd, looM4Se)) %>% 
  rownames_to_column("model") %>% arrange(model) %>% 
  mutate(across(where(is.numeric), ~round(., 2))) %>% 
  add_row(.before = 1, 
          model = "Variation in fecundity due to life history, body condition, and sample date (20 rivers)")
# 4. c) third set of models (egg diameter)
set3 <- data.frame(loo_compare(looED, looEDa)) %>% 
  rownames_to_column("model") %>% arrange(model) %>% 
  mutate(across(where(is.numeric), ~round(., 2))) %>% 
    add_row(.before = 1, 
            model = "Variation in fecundity due to egg diameter (6 rivers)")
# 4. d) fourth set of models (residual fecundity)
set4 <- data.frame(loo_compare(looR0, looR1, looR2)) %>% 
  rownames_to_column("model") %>% arrange(model) %>% 
  mutate(across(where(is.numeric), ~round(., 2))) %>% 
  add_row(.before = 1, 
          model = "Temporal variation in residual fecundity (10 rivers)")
# 4. e) bind data across model sets and final formatting
# bind_rows(set1, set2, set3, set4) %>% 
set3 %>% 
  mutate(`Effective # of parameters (SE)` = paste(p_loo, " (", se_p_loo, ")", sep = ""),
         `LOO ELPD (SE)`                  = paste(elpd_loo, " (", se_elpd_loo, ")", sep = ""),
         `ΔELPD (SE)`                     = paste(elpd_diff, " (", se_diff, ")", sep = ""),
         `LOOIC (SE)`                     = paste(looic, " (", se_looic, ")", sep = ""),
         across(c(`Effective # of parameters (SE)`:`LOOIC (SE)`), ~if_else(. == "NA (NA)", "", .))) %>% 
  select(-c(elpd_diff:se_looic)) %>% View()



# TABLE S2 ----------------------------------------------------------------

# 1. get coefficients for two models
coefM4 <- distinct(fecundity, RIVER_NO, RIVER) %>% 
  left_join(coef(M4) %>% data.frame() %>% rownames_to_column("RIVER")) %>% 
  mutate(across(where(is.numeric), ~round(., 2)),
         β0 = paste(RIVER.Estimate.Intercept, " (", 
                    RIVER.Est.Error.Intercept, ")", sep = ""),
         βf = paste(RIVER.Estimate.clFORK_LENGTH, " (", 
                    RIVER.Est.Error.clFORK_LENGTH, ")", sep = ""))
coefM4se<- distinct(fecundity, RIVER) %>% 
  left_join(coef(M4Se) %>% data.frame() %>% rownames_to_column("RIVER")) %>% 
  mutate(across(where(is.numeric), ~round(., 2)),
         β0 = paste(RIVER.Estimate.Intercept, " (", 
                    RIVER.Est.Error.Intercept, ")", sep = ""),
         βf = paste(RIVER.Estimate.clFORK_LENGTH, " (", 
                    RIVER.Est.Error.clFORK_LENGTH, ")", sep = ""))

# 2. Table S2
left_join(select(coefM4, -contains(".")),
          select(coefM4se, -contains(".")),
          by = "RIVER") %>% arrange(RIVER_NO) %>% 
  mutate(across(c(β0.x:βf.y), ~if_else(. == "NA (NA)", "", .))) %>% View()
names(tableS2)

# 3. Table S3 footnotes
# 3. a) model M4 footnotes
mutate(coefM4, across(where(is.numeric), ~round(., 2)),
       βmM = paste0(RIVER.Estimate.STAGE_EGG_COUNTmature, " (",
                   RIVER.Est.Error.STAGE_EGG_COUNTmature, ")"),
       βmMR = paste0(RIVER.Estimate.STAGE_EGG_COUNTmaturePretained, " (",
                    RIVER.Est.Error.STAGE_EGG_COUNTmaturePretained, ")")) %>% 
  select(contains("βm")) %>% slice(1)

# 3. b) model M4se footnotes
posterior_summary(M4Se) %>% data.frame() %>% rownames_to_column(var = "coef") %>% 
  filter(grepl("STAGE_EGG_COUNT|fRAGE|LAST_SPAWN|FULTONS_K|SAMPLE_DOY", coef)) %>% 
  mutate(across(where(is.numeric), ~round(., 2)),
         coef = case_when(coef == "b_STAGE_EGG_COUNTmature"          ~ "βm = mature =",
                          coef == "b_STAGE_EGG_COUNTmaturePretained" ~ "βm = mature + retained =",
                          coef == "b_STAGE_EGG_COUNTmature"          ~ "βm = mature =",
                          coef == "b_cFULTONS_K"                     ~ "βk =",
                          coef == "b_cSAMPLE_DOY"                    ~ "βt =",
                          grepl("fRAGE", coef)      ~ paste0("βa = ", str_sub(coef, -1, -1), " = "),
                          grepl("LAST_SPAWN", coef) ~ paste0("βp = ", str_sub(coef, -1, -1), " = ")),
       β = paste0(Estimate, " (", Est.Error, ")")) %>%
  select(coef, β)
