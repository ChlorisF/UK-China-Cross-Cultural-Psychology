### ANALASES AND FIGURES
# Cultural generalisability of psychological measures
# 20-08-23 (updated 25-02-10)
# Boyin Feng

# Prepared with R Studio version 2021.09.0
# Prepared with R version 4.1.3

##############################################################################
#clean up the global environment
rm(list=ls())

#### ==== Package Installation ==== 

if (!require(tidyverse)) install.packages("tidyverse"); require(tidyverse)
if (!require(dplyr)) install.packages("dplyr"); require(dplyr)
if (!require(lme4)) install.packages("lme4"); require(lme4)
if (!require(lmerTest)) install.packages("lmerTest"); require(lmerTest)
if (!require(effsize)) install.packages("effsize"); require(effsize)
if (!require(metafor)) install.packages("metafor"); require(metafor)
if (!require(emmeans)) install.packages("emmeans"); require(emmeans)
if (!require(psych)) install.packages("psych"); require(psych)
if (!require(meta)) {install.packages("meta"); require(meta)}
if (!require(corrplot)) {install.packages("corrplot"); require(corrplot)}

#### ==== Data Setup ==== 
# set up working directory [CHANGE into own directory]
setwd("D:/analysis/data files")
dat <- read.csv("full_data.csv")

#### ==== Main Analyses ==== 
# check for data normality
## Function to create Q-Q plots for each measure
create_qq_plot <- function(var) {
  ggplot(dat, aes(sample = .data[[var]])) +  # Fixed: `data` → `dat`
    stat_qq() +
    stat_qq_line(color = "red") +
    ggtitle(paste("Q-Q Plot for", var)) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),  
          axis.title.y = element_blank()) 
}

## Generate and arrange Q-Q plots
# jpeg(file='QQplot.jpeg', height=1000, width=2000, res = 250) #Only use if want to export the plot
qq_plots <- lapply(var, create_qq_plot)

layout_matrix <- rbind(c(1, 2, 3, 4),  # First row: 4 plots
                       c(5, 6, 7, NA))  # Second row: 3 plots (last cell is empty)

grid.arrange(grobs = qq_plots, layout_matrix = layout_matrix) # Arrange plots using the layout matrix

# dev.off() #Use if want to export the plot
## Q-Q plots show skewness, refer to Supplementary Analyses for further confirmation using Mann-Whitney U-tests


# forest plot for all the measures between groups
dat_UK <- subset(dat, nationality == "UK")
dat_CN <- subset(dat, nationality == "China")

## multiple pairwise t-tests
dat_compare <- dat %>%
  select(!c(X,ID,sex,age, ethnicity, education,fixed_go_acc,random_go_acc,fixed_RT,random_RT)) %>%
  mutate(nationality = recode(nationality, 
                              "UK" = "British", 
                              "China" = "Chinese")) %>%
  rename(SART_fixed = fixed_nogo_acc, 
         SART_random = random_nogo_acc,
         Asian_congruent = EF_Chi_c,
         Asian_incongruent = EF_Chi_inc,
         West_congruent = EF_West_c,
         West_incongruent = EF_West_inc) #corrected hyphen error 2025-03-04


## Transform the data into long format
### Put all variables in the same column except `nationality`, the grouping variable
compareItem <- c("BDI","BIS","OCIR","STAI",
                 "Asian_Angry","Asian_Fearful","Asian_Happy",
                 "Western_Angry","Western_Fearful","Western_Happy",
                 "Asian_congruent","Asian_incongruent","West_congruent","West_incongruent",
                 "SART_fixed","SART_random")

dat.long <- dat_compare %>%
  pivot_longer(-nationality, names_to = "variables", values_to = "value") %>%
  filter((variables %in% compareItem))

stat.test <- dat.long %>%
  group_by(variables) %>%
  t_test(value ~ nationality, detailed = TRUE) %>% #"detailed" shows group means (estimate 1&2)
  adjust_pvalue(p.col = "p", output.col = "p.adj",method = "bonferroni") %>%
  add_significance()

## remove irrelevant columns and trailing zeros
stat.test <- stat.test %>%
  select(variables,estimate1,estimate2,estimate,statistic,df,p,p.adj,p.adj.signif) %>%
  rename(UK=estimate1,CN=estimate2,Diff=estimate,t = statistic) %>%
  mutate(df = as.numeric(as.character(df)),
         p = as.numeric(as.character(p)),
         p.adj = as.numeric(as.character(p.adj)))

## convert scientific notations to decimals
options(scipen = 999)

## calculate effect sizes
cohens.d <- dat.long %>%
  group_by(variables) %>%
  cohens_d(value ~ nationality,
           ci = TRUE)

## combine effect sizes with t-tests
cohens.d <- cohens.d %>%
  select(variables, effsize, conf.low, conf.high, magnitude) %>%
  mutate(conf.low = as.numeric(as.character(conf.low)),
         conf.high = as.numeric(as.character(conf.high)))


### define subgroups based on measures
ER_S <- c("Asian_Angry", "Asian_Fear", "Asian_Happy", "West_Angry", "West_Fear", "West_Happy")
ER_F <- c("Asian-congruent","Asian-incongruent","West-congruent","West-incongruent")
SART <- c("SART-fixed","SART-random")
MH <- c("BDI","BIS","OCIR","STAI")

forest <- stat.test %>%
  full_join(cohens.d, by = "variables") %>%
  rename(Item = variables,lower = conf.low, upper = conf.high) %>% #make column name readable for forest function
  mutate(type = case_when(Item %in% ER_S ~ "ER-S Accuracy, %",
                          Item %in% ER_F ~ "ER-F Accuracy, %",
                          Item %in% SART ~ "SART Accuracy, %",
                          Item %in% MH ~ "Mental Health Scores"))

## draw the forestplot
library(forestplot)
library(meta)
library(metafor)
forest <- read.csv("forest_all.csv")
forest <- forest %>%
  rename(Item = 锘縄tem) # skip this step if the dataset reads as Item already
forest <- tibble::tibble(mean  = forest$effsize,
                         lower = forest$lower,
                         upper = forest$upper,
                         Item = forest$Item,
                         type = forest$type,
                         UK = forest$UK,
                         CN = forest$CN,
                         p = forest$p,
                         adjusted_p = forest$p.adj)
### build random-effects model
m.gen <- metagen(TE = mean,
                 studlab = Item,
                 data = forest,
                 lower = forest$lower,
                 upper = forest$upper,
                 fixed = FALSE,
                 random = TRUE,
                 method.tau = "REML",
                 sm = "SMD",
                 hakn = TRUE,
                 title = "Effect Sizes for Items")

summary(m.gen)

sub <- update(m.gen, subgroup = forest$type,
              print.subgroup.name = FALSE)
update(m.gen, subgroup = forest$type,
       print.subgroup.name = FALSE)

forplot <- forest.meta(sub, 
                       prediction = TRUE, 
                       print.tau2 = FALSE,
                       leftcols = c("studlab", "CN"),
                       leftlabs = c("Item", "CN"),
                       rightcols = c("UK", "p","adjusted_p"),
                       rightlabs = c("UK", "P","adjusted P"),
                       test.subgroup = FALSE)
drapery(m.gen,
        labels = "studlab",
        type = "pvalue", # or zvalue
        legend = F)

# correlation matrices between variables
dat <- read.csv("full_data.csv")
dat_UK <- subset(dat, nationality == "UK")
dat_CN <- subset(dat, nationality == "China")

## select columns for correlations
UK_cor <- dat_UK %>% 
  select(BDI, BIS, OCIR, STAI, 
         ER_Happy, ER_Angry, ER_Fearful, EF_Happy_c, EF_Happy_inc, EF_Angry_c, EF_Angry_inc,
         fixed_nogo_acc, random_nogo_acc) %>% 
  rename(SART_fixed = fixed_nogo_acc, SART_random = random_nogo_acc,
         ES_Happy = ER_Happy, ES_Angry = ER_Angry, ES_Fear = ER_Fearful)
## compute correlation coefficients and p-values
cor_UK <- cor(UK_cor) 

### define cor.mtest function
### mat : is a matrix of data
### ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat) 
  p.mat
}
p.mat.UK <- cor.mtest(UK_cor, conf.level = 0.95)

## draw correlation matrix
jpeg(file='corrplot_UK_spearman.jpeg', height=2300, width=2500, res = 250)
corrplot(cor_UK, 
         type="lower", method="circle", 
         addCoef.col ='black', number.cex = 1.1, 
         tl.pos="lt", tl.col="black",  tl.offset=1, tl.srt=45,
         # col = COL2('RdBu'),
         # sig.level = c(0.001, 0.01, 0.05), pch.cex = 2.3,
         tl.cex = 1.25,cl.cex = 1.25,
         # is.corr = FALSE, 
         pch.col = 'grey20')$corrPos -> p1
text(p1$x, p1$y, round(p1$corr, 2), cex = 1.1)

corrplot(cor_UK, p.mat = p.mat.UK,add=T, 
         type="upper", method="color",
         diag=F, tl.pos="n", 
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 2,
         cl.cex = 1.25,
         insig = 'label_sig', pch.col = 'grey20')

n <- nrow(cor_UK)
dev.off()

# China plot
CN_cor <- dat_CN %>% 
  select(BDI, BIS, OCIR, STAI, 
         ER_Happy, ER_Angry, ER_Fearful, EF_Happy_c, EF_Happy_inc, EF_Angry_c, EF_Angry_inc,
         fixed_nogo_acc, random_nogo_acc) %>% 
  rename(SART_fixed = fixed_nogo_acc, SART_random = random_nogo_acc,
         ES_Happy = ER_Happy, ES_Angry = ER_Angry, ES_Fear = ER_Fearful)
## compute correlation coefficients and p-values
cor_CN <- cor(CN_cor, method = "spearman") # change method from default pearson to spearman
p.mat.CN <- cor.mtest(CN_cor, conf.level = 0.95)

## draw correlation matrix
jpeg(file='corrplot_CN_spearman.jpeg', height=2300, width=2500, res = 250)
corrplot(cor_CN, 
         type="lower", method="circle", 
         addCoef.col ='black', number.cex = 1.1, 
         tl.pos="lt", tl.col="black",  tl.offset=1, tl.srt=45,
         # col = COL2('RdBu'),
         # sig.level = c(0.001, 0.01, 0.05), pch.cex = 2.3,
         tl.cex = 1.25,cl.cex = 1.25,
         pch.col = 'grey20')$corrPos -> p1
text(p1$x, p1$y, round(p1$corr, 2), cex = 1.1)

corrplot(cor_CN, p.mat = p.mat.CN,add=T, 
         type="upper", method="color",
         diag=F, tl.pos="n", 
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 2,
         cl.cex = 1.25,
         insig = 'label_sig', pch.col = 'grey20')

n <- nrow(cor_CN)
dev.off()

# demographic comparisons
library(rcompanion)
library(effsize)

## vector type setting
dat$age <- as.numeric(dat$age)
dat$education <- as.factor(dat$education)
dat$sex <- as.factor(dat$sex)
orders <- c('Primary school level','Secondary school level','Undergraduate degree','Postgraduate degree')
unique(orders)
dat$education <- c("Primary school level"=6, 
                   "Secondary school level"=12,
                   "Undergraduate degree"=15,
                   "Postgraduate degree"=18)[dat$education]
dat$education <- as.numeric(dat$education)

## subset data
demo_UK <- subset(dat, nationality == "UK")
demo_CH <- subset(dat, nationality == "China")

## age
t.test(demo_UK$age, demo_CH$age)
cohen.d(demo_UK$age, demo_CH$age)

## sex
sex = table(dat$nationality,dat$sex)                
chisq.test(sex)

## education
t.test(demo_UK$education, demo_CH$education)
cohen.d(demo_UK$education, demo_CH$education)

#### ==== Supplementary Analyses ==== 
# Perform Mann-Whitney U-test for confirmation due to nonnormal distribution
stat.test <- dat.long %>%
  group_by(variables) %>%
  wilcox_test(value ~ nationality, detailed = TRUE) %>% # Perform the U-test
  adjust_pvalue(p.col = "p", output.col = "p.adj", method = "bonferroni") %>% # Adjust p-values
  add_significance() # Add significance labels

## Calculate group medians
group_medians <- dat.long %>%
  group_by(variables, nationality) %>%
  summarise(median = median(value, na.rm = TRUE)) %>% # Calculate medians for each group
  pivot_wider(names_from = nationality, values_from = median, names_prefix = "median_") # Reshape to wide format

## Merge medians with the U-test results
stat.test <- stat.test %>%
  left_join(group_medians, by = "variables") %>%  # Merge medians into the main results
  select(variables, median_British, median_Chinese, statistic, p, p.adj, p.adj.signif) %>% # Remove irrelevant columns and rename for clarity
  rename(UK = median_British, CN = median_Chinese, W = statistic) # Rename columns

## Convert p-values to numeric
stat.test <- stat.test %>%
  mutate(p = as.numeric(as.character(p)),
         p.adj = as.numeric(as.character(p.adj)))

## View the final results
print(stat.test) #The U-tests show consistent results with t-tests


# forest plot for OCIR question-wise comparison
for_OCIR <- read.csv("forest_OCIR.csv")

for_OCIR <- tibble::tibble(mean  = for_OCIR$effsize,
                         lower = for_OCIR$conf.low,
                         upper = for_OCIR$conf.high,
                         Question = for_OCIR$Question,
                         UK = for_OCIR$UK,
                         CN = for_OCIR$CN,
                         p = for_OCIR$p,
                         adjusted_p = for_OCIR$p.adj)


## build models
m.gen.ocir <- metagen(TE = mean,
                 studlab = Question,
                 data = for_OCIR,
                 lower = for_OCIR$lower,
                 upper = for_OCIR$upper,
                 fixed = FALSE,
                 random = TRUE,
                 method.tau = "REML",
                 sm = "SMD",
                 hakn = TRUE,
                 title = "Effect Sizes for Each Question")

summary(m.gen.ocir)

forplot_OCIR <- forest.meta(m.gen.ocir, 
                       prediction = TRUE, 
                       print.tau2 = TRUE,
                       xlim = c(-0.8, 0.2), at = c(-0.8,-0.5,-0.2,0,0.2),
                       leftcols = c("studlab", "CN"),
                       leftlabs = c("Question", "CN"),
                       rightcols = c("UK", "p","adjusted_p"),
                       rightlabs = c("UK", "P","adjusted P"),
                       test.subgroup = FALSE)


# ANCOVA models
library(lsr)
library(lmtest)

## BDI
m_BDI <- lm(BDI_Score ~ sex + age + education, data = dat)
summary(m_BDI)
m2_BDI <- lm(BDI_Score ~ sex + age + education + nationality, data = dat)
summary(m2_BDI)
### log likelihood ratio test
lrtest(m2_BDI, m_BDI)
### Compare model metrics
AIC(m_BDI); AIC(m2_BDI) 


## BIS
m_BIS <- lm(BIS_Score ~ sex + age + education, data = dat)
summary(m_BIS)
m2_BIS <- lm(BIS_Score ~ sex + age + education + nationality, data = dat)
summary(m2_BIS)
### log likelihood ratio test
lrtest(m2_BIS, m_BIS)
### Compare model metrics
AIC(m_BIS); AIC(m2_BIS) 

## OCIR
m_OCIR <- lm(OCIR_Score ~ sex + age + education, data = dat)
summary(m_OCIR)
m2_OCIR <- lm(OCIR_Score ~ sex + age + education + nationality, data = dat)
summary(m2_OCIR)
### Compare model metrics
AIC(m_OCIR); AIC(m2_OCIR) 
### log likelihood ratio test
lrtest(m2_OCIR, m_OCIR)

## STAI
m_STAI <- lm(STAI_Score ~ sex + age + education, data = dat)
summary(m_STAI)
m2_STAI <- lm(STAI_Score ~ sex + age + education + nationality, data = dat)
summary(m2_STAI) 
### log likelihood ratio test
lrtest(m2_STAI, m_STAI)
### Compare model metrics
AIC(m_STAI); AIC(m2_STAI) 

## ER-S
m1_ER <- aov(ER_meanAcc ~ ER_item * ER_img_Nat + sex + age + education, data = ER_ex)
m2_ER <- aov(ER_meanAcc ~ ER_item * ER_img_Nat * nationality + sex + age + education, data = ER_ex)
lrtest(m2_ER, m1_ER)

## ER-F
m1_EF <- aov(EF_meanAcc ~ EF_item * T_nat + Match + sex + age + education, data = EF_ex)
m2_EF <- aov(EF_meanAcc ~ EF_item * T_nat * nationality + Match + sex + age + education, data = EF_ex)
lrtest(m2_EF, m1_EF) 

## SART
m_sart_f1 <- aov(fixed_nogo_acc ~ sex + age + education, data = dat)
summary(m_sart_f1)
m_sart_f2 <- aov(fixed_nogo_acc ~ sex + age + education + nationality, data = dat)
summary(m_sart_f2)
lrtest(m_sart_f2, m_sart_f1)

m_sart_r1 <- aov(random_nogo_acc ~ sex + age + education, data = dat)
summary(m_sart_r1)
m_sart_r2 <- aov(random_nogo_acc ~ sex + age + education + nationality, data = dat)
summary(m_sart_r2)
lrtest(m_sart_r2, m_sart_r1)
