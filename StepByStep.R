###############################################################################X
#############         Alligator body condition method         #################X
###############################################################################X
## Citation: Balaguera-Reina et al. Redefining the Body Condition Approach:  ##X 
## American Alligators as a Case Study. Ecology. Submitted.                  ##X
###############################################################################X
#############        Last date modified: may 5th, 2025           ##############X          
###############################################################################X

rm(list = ls())
dev.off()

## Libraries----
library(tidyverse)
library(MASS)
library(emmeans)
library(patchwork)
library(ggpmisc)
library(ggpubr)

## 1. Call and format data for analysis----
Gators_All <- readxl::read_xlsx("American_Alligators_Database_20240605.xlsx", 
                          sheet = "Captures_Individuals")
Gators <- Gators_All %>%
  mutate(`Capture Date` = as.Date(`Capture Date`)) %>%
  subset(!`Total Length` < 125 & !Marsh_River_Canal %in% "Canal" & 
           !Route %in% c("ENP-NESSE", "Panther-Cochran Lake", "Panther-Tram", "ENP-EST")) %>%
  dplyr::select(Consecutive, `Capture Date`, Area, Route, Marsh_River_Canal, UTM_Easting, UTM_Northing, 
                `Water Depth`, `General Condition`, `SV length`, `Total Length`, Weight, Sex, 
                Deformities_Notes, Deformities_List, Notes, Exclude, `Was the animal caught for size estimates only?`) %>%
  mutate(lnWeight = log(Weight)) %>%
  mutate(lnSVL = log(`SV length`)) %>%
  mutate(SizeClass = case_when(`Total Length` < 175 ~ "Subadult",
                               `Total Length` >= 175 ~ "Adult")) %>%
  mutate(month = as.integer(format(`Capture Date`, "%m"))) %>%
  mutate(Season = case_when(month < 5 ~ "Spring",
                            month > 8 & month < 12 ~ "Fall")) %>%
  mutate(year = as.integer(format(`Capture Date`, "%Y"))) %>%
  mutate(WY = ifelse(month >= 6 & month <= 12, year + 1, year))

Gators <- Gators[complete.cases(Gators$Season),]
Gators <- Gators[complete.cases(Gators$Weight),]
Gators <- Gators[complete.cases(Gators$`SV length`),]

## 2. Summarise data (create and export table 1, Appendix 1)----
Gators %>%
  group_by(Sex) %>%
  summarise(
    n()
  )

Table1 <- Gators %>%
  group_by(Route) %>%
  summarise(
    `Number of alligators` = n(),
    Spring = length(which(Season %in% "Spring")),
    Fall = length(which(Season %in% "Fall")),
    Male = length(which(Sex %in% "Male")),
    Female = length(which(Sex %in% "Female")),
    Unknown = length(which(Sex %in% "U")),
    Subadult = length(which(SizeClass %in% "Subadult")),
    Adult = length(which(SizeClass %in% "Adult")),
    WY = paste(min(WY), max(WY), sep = "-")
  )

write.csv(Table1, "./Tables/Table1.csv")

sum(Table1$Subadult)
sum(Table1$Adult)

sum(Table1$Fall)
sum(Table1$Spring)

Appendix1 <- Gators %>%
  group_by(WY, Season, Route) %>%
  summarise(
    `Number of alligators` = n(),
    Male = length(which(Sex %in% "Male")),
    Female = length(which(Sex %in% "Female")),
    Unknown = length(which(Sex %in% "U")),
    Subadult = length(which(SizeClass %in% "Subadult")),
    Adult = length(which(SizeClass %in% "Adult"))
    )

Appendix1 %>%
  group_by(WY, Route) %>%
  summarise(
    n = sum(`Number of alligators`)
  ) %>%
  arrange(n)

Appendix1 %>%
  group_by(Route, Season) %>%
  summarise(
    n = sum(`Number of alligators`)
  ) %>%
  subset(Season %in% "Fall") %>%
  arrange(n)

Appendix1 %>%
  group_by(Route, Season) %>%
  summarise(
    n = sum(`Number of alligators`)
  ) %>%
  subset(Season %in% "Spring") %>%
  arrange(n)

range(Gators$`Total Length`, na.rm = T)
mean(Gators$`Total Length`, na.rm = T)
range(Gators$`SV length`)
mean(Gators$`SV length`)
range(Gators$Weight)
mean(Gators$Weight)

More10 <- Appendix1 %>% #Define sampling sites with more than 10 gators per season / WY
  subset(`Number of alligators` >= 10) %>%
  mutate(code = paste(WY, Season, Route, sep = "_"))

Less10 <- Appendix1 %>% #Define sampling sites with more than 10 gators per season / WY
  subset(`Number of alligators` < 10) %>%
  mutate(code = paste(WY, Season, Route, sep = "_"))

Less10 %>%
  group_by(Route) %>%
  summarise(
    n()
  )

## 3. Define allometric coefficients----
Gators$code <- paste(Gators$WY, Gators$Season, Gators$Route, sep = "_")
More10 <- data.frame(code = More10$code)
SitesMore10 <- left_join(More10, Gators, #Create dataframe with only sampling sites with 10 or more gators
               by = 'code')

a <- SitesMore10 %>% # Test it works
  group_by(WY, Season, Route) %>%
  summarise(
    `Number of alligators` = n(),
    Male = length(which(Sex %in% "Male")),
    Female = length(which(Sex %in% "Female")),
    Unknown = length(which(Sex %in% "U")),
    Subadult = length(which(SizeClass %in% "Subadult")),
    Adult = length(which(SizeClass %in% "Adult"))
  )

  ## 3.1 Test for outliers----
W_SVL <- lm(lnWeight ~ lnSVL, data = SitesMore10)
summary(W_SVL)
cor.test(x = W_SVL$residuals, y = W_SVL$fitted.values, # Test for heteroscedasticity 
         method = "spearman") 
stud_resid <- studres(W_SVL) # identify outliers 
SitesMore10 <- cbind(SitesMore10, stud_resid)
Outliers <- subset(SitesMore10, stud_resid < -3.0 | stud_resid > 3.0)
Outliers %>%
  group_by(Route) %>%
  summarise(
    n()
  )

p.outliers <- ggplot() +
  geom_point(SitesMore10, mapping = aes(y = stud_resid, x = lnSVL)) +
  geom_point(Outliers, mapping = aes(y = stud_resid, x = lnSVL, colour = Route)) +
  geom_hline(yintercept = c(0, -3, 3), col = c('red', 'blue', 'blue')) +
  labs( y = 'Studentized Residuals', x = 'Displacement') +
  scale_color_manual(name = 'Sampling site', labels = c("BICY", "ENP-FC", "ENP-SS", "LOX", "WCA2A",
                                                          "WCA3A-HD", "WCA3A-N41", "WCA3A-TW", "WCA3B"), 
                       values = c("#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F", "#00B9E3",
                                  "#619CFF", "#DB72FB", "#FF61C3")) +
  theme_bw()

Gators_NO <- subset(SitesMore10, !stud_resid < -3.0 & !stud_resid > 3.0) # 50 outliers (larger than 3 and lower than -3)
W_SVL <- lm(lnWeight ~ lnSVL, data = Gators_NO)
summary(W_SVL)
cor.test(x = W_SVL$residuals, y = W_SVL$fitted.values, 
         method = "spearman")
coefficients(W_SVL)
confint(W_SVL, level=0.95)

  ## 3.2. Assess effects data without outliers----
Appendix2 <- data.frame(contrast = NULL,
                     estimate = NULL,
                     SE = NULL,
                     df = NULL,
                     t.ratio = NULL,
                     p.value = NULL)
#By Sex
Gators_NOS <- subset(Gators_NO, !Sex %in% "U")
a <- lm(lnWeight ~ lnSVL*Sex, data = Gators_NOS)
summary(a)
anova(a)
a <- lstrends(a, "Sex", var = "lnSVL")
a <- data.frame(pairs(a)) # Moderate evidence by sex
Appendix2 <- rbind(Appendix2, a)

p.sexNO <- ggplot(Gators_NOS, aes(y = lnWeight, x = lnSVL, color = Sex)) +
  geom_point(alpha = 2/10) +
  geom_smooth(method = "lm", se = T, formula = y ~ x) +
  labs(y = "ln weight", x = "", colour = "Sex") +
  theme_bw()
  
#By Size Class
a <- lm(lnWeight ~ lnSVL*SizeClass, data = Gators_NO)
summary(a)
anova(a)
a <- lstrends(a, "SizeClass", var = "lnSVL")
a <- data.frame(pairs(a)) # Strong evidence by Size Class
Appendix2 <- rbind(Appendix2, a)

p.SCNO <- ggplot(Gators_NO, aes(y = lnWeight, x = lnSVL, color = SizeClass)) +
  geom_point(alpha = 2/10) +
  geom_smooth(method = "lm", se = T, formula = y ~ x) +
  labs(y = "", x = "", colour = "Size Class") +
  theme_bw()

#By season
a <- lm(lnWeight ~ lnSVL*Season, data = Gators_NO)
summary(a)
anova(a)
a <- lstrends(a, "Season", var = "lnSVL")
a <- data.frame(pairs(a)) # No evidence of an effect by season
Appendix2 <- rbind(Appendix2, a)

p.SNO <- ggplot(Gators_NO, aes(y = lnWeight, x = lnSVL, color = Season)) +
  geom_point(alpha = 2/10) +
  geom_smooth(method = "lm", se = T, formula = y ~ x) +
  labs(y = "ln weight", x = "", colour = "Season") +
  theme_bw()

#By Route
a <- lm(lnWeight ~ lnSVL*Route, data = Gators_NO)
summary(a)
anova(a)
confint(lstrends(a, "Route", var = "lnSVL"), adjust = "Bonferroni")
a <- test(pairs(lstrends(a, "Route", var = "lnSVL"), adjust = "Bonferroni"))
Appendix2 <- rbind(Appendix2, a) # Strong evidence of an effect by route

p.RNO <- ggplot(Gators_NO, aes(y = lnWeight, x = lnSVL, color = Route)) +
  geom_point(alpha = 2/10) +
  geom_smooth(method = "lm", se = F, formula = y ~ x) +
  stat_regline_equation(aes(y = lnWeight, x = lnSVL), 
                        label.x = c(4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.8, 4.8, 4.8), 
                        label.y = c(5, 4.75, 4.5, 4.25, 4, 3.75, 2.5, 2.25, 2)) +
  scale_color_manual(name = 'Sampling site', labels = c("BICY", "ENP-FC", "ENP-SS", "LOX", "WCA2A",
                                                        "WCA3A-HD", "WCA3A-N41", "WCA3A-TW", "WCA3B"), 
                     values = c("#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F", "#00B9E3",
                                "#619CFF", "#DB72FB", "#FF61C3")) +
  labs(y = "ln Weight", x = "ln Snout Vent Length") +
  theme_bw() +
  theme(legend.position = 'none')

#By WY
a <- lm(lnWeight ~ lnSVL*as.factor(WY), data = Gators_NO)
summary(a)
anova(a)
confint(lstrends(a, "WY", var = "lnSVL"), adjust = "Bonferroni")
a <- test(pairs(lstrends(a, "WY", var = "lnSVL"), adjust = "Bonferroni")) 
Appendix2 <- rbind(Appendix2, a) # Strong evidence of an effect by route

write.csv(Appendix2, "./Tables/Table2.csv")

p.WYNO <- ggplot(Gators_NO, aes(y = lnWeight, x = lnSVL, color = as.factor(WY))) +
  geom_point(alpha = 2/10) +
  geom_smooth(method = "lm", se = T, formula = y ~ x) +
  labs(y = "lm weight", x = "ln SVL", colour = "Water year") +
  theme_bw()

Appendix2 <- Appendix2 %>%
  mutate(across(where(is.numeric), ~ round(., digits = 3)))

  ## 3.3. Remove most different routes----
Gators_NO_Central <- subset(Gators_NO, !Route %in% c("Lox_Marsh", "BICY_Marsh", "ENP-SS"))

  ## 3.4. Assess effects data without most different routes (Central)----
Appendix2.1 <- data.frame(contrast = NULL,
                       estimate = NULL,
                       SE = NULL,
                       df = NULL,
                       t.ratio = NULL,
                       p.value = NULL)
#By Sex
Gators_1S <- subset(Gators_NO_Central, !Sex %in% "U")
a <- lm(lnWeight ~ lnSVL*Sex, data = Gators_1S)
summary(a)
anova(a)
a <- lstrends(a, "Sex", var = "lnSVL")
a <- data.frame(pairs(a)) # No difference by sex
Appendix2.1 <- rbind(Appendix2.1, a)

#By Size Class
a <- lm(lnWeight ~ lnSVL*SizeClass, data = Gators_NO_Central)
summary(a)
anova(a)
a <- lstrends(a, "SizeClass", var = "lnSVL")
a <- data.frame(pairs(a)) # Moderate evidence by Size Class

#By season
a <- lm(lnWeight ~ lnSVL*Season, data = Gators_NO_Central)
summary(a)
anova(a)
a <- lstrends(a, "Season", var = "lnSVL")
a <- data.frame(pairs(a)) # No evidence by season
Appendix2.1 <- rbind(Appendix2.1, a)

#By Route
a <- lm(lnWeight ~ lnSVL*Route, data = Gators_NO_Central)
summary(a)
anova(a)
confint(lstrends(a, "Route", var = "lnSVL"), adjust = "Bonferroni")
a <- test(pairs(lstrends(a, "Route", var = "lnSVL"), adjust = "Bonferroni")) #No evidence differences 
Appendix2.1 <- rbind(Appendix2.1, a)

#By WY
a <- lm(lnWeight ~ lnSVL*as.factor(WY), data = Gators_NO_Central)
summary(a)
anova(a)
b <- confint(lstrends(a, "WY", var = "lnSVL"), adjust = "Bonferroni")
b$Group <- "Central"
Slope_Groups <- b
a <- test(pairs(lstrends(a, "WY", var = "lnSVL"), adjust = "Bonferroni")) #Strong differences
Appendix2.1 <- rbind(Appendix2.1, a)

#By size class and route
a <- lm(lnWeight ~ lnSVL*Route*SizeClass, data = Gators_NO_Central)
summary(a)
anova(a)
confint(lstrends(a, c("Route", "SizeClass"), var = "lnSVL"), adjust = "Bonferroni")
a <- test(pairs(lstrends(a, c("Route", "SizeClass"), var = "lnSVL"), adjust = "Bonferroni")) #No evidence differences 

#No covariates
a <- lm(lnWeight ~ lnSVL, data = Gators_NO_Central)
summary(a)
coefficients(a)
confint(a)

Gators_NO_Central$Allometric_elevation <- (Gators_NO_Central$Weight / Gators_NO_Central$`SV length`^3.17)*10^5

Appendix2.1 <- Appendix2.1 %>%
  mutate(across(where(is.numeric), ~ round(., digits = 3)))

  ## 3.5. Group sampling areas based on evidence----
Gators_NO <- Gators_NO %>% #Group data based on slopes
  mutate(Groups = case_when(Route %in% c("ENP-FC", "WCA3B_Marsh", "WCA3A-HD_Marsh", "WCA3A-TW_Marsh", 
                                         "WCA3A-N41_Marsh", "WCA2A_Marsh") ~ "Central",
                            Route %in% c("Lox_Marsh") ~ "Northeast",
                            Route %in% c("ENP-SS") ~ "South",
                            Route %in% c("BICY_Marsh") ~ "West"))

writexl::write_xlsx(Gators_NO, "./Data_Used_Models.xlsx")
Groups <- split(Gators_NO, ~ Gators_NO$Groups)

  ## 3.6. Assess effects data grouped (South)----
Appendix2.2 <- data.frame(contrast = NULL,
                     estimate = NULL,
                     SE = NULL,
                     df = NULL,
                     t.ratio = NULL,
                     p.value = NULL)

#By Sex
Gators_S <- subset(Groups$South, !Sex %in% "U")
a <- lm(lnWeight ~ lnSVL*Sex, data = Gators_S)
summary(a)
anova(a)
a <- lstrends(a, "Sex", var = "lnSVL")
a <- data.frame(pairs(a)) # Weak difference by sex
Appendix2.2 <- rbind(Appendix2.2, a)

#By Size Class
a <- lm(lnWeight ~ lnSVL*SizeClass, data = Groups$South)
summary(a)
anova(a)
a <- lstrends(a, "SizeClass", var = "lnSVL")
a <- data.frame(pairs(a)) # No evidence by Size Class

#By season
a <- lm(lnWeight ~ lnSVL*Season, data = Groups$South)
summary(a)
anova(a)
a <- lstrends(a, "Season", var = "lnSVL")
a <- data.frame(pairs(a)) # No differences by season
Appendix2.2 <- rbind(Appendix2.2, a)

#By WY
a <- lm(lnWeight ~ lnSVL*as.factor(WY), data =Groups$South)
summary(a)
anova(a)
b <- confint(lstrends(a, "WY", var = "lnSVL"), adjust = "Bonferroni")
b$Group <- "South"
Slope_Groups <- rbind(Slope_Groups, b)
a <- test(pairs(lstrends(a, "WY", var = "lnSVL"), adjust = "Bonferroni")) 
Appendix2.2 <- rbind(Appendix2.2, a)

Appendix2.2 <- Appendix2.2 %>%
  mutate(across(where(is.numeric), ~ round(., digits = 3)))

#No covariates
a <- lm(lnWeight ~ lnSVL, data = Groups$South)
summary(a)
coefficients(a)
confint(a)

  ## 3.7. Assess effects data grouped (Northeast)----
Appendix2.3 <- data.frame(contrast = NULL,
                          estimate = NULL,
                          SE = NULL,
                          df = NULL,
                          t.ratio = NULL,
                          p.value = NULL)

#By Sex
Gators_S <- subset(Groups$Northeast, !Sex %in% "U")
a <- lm(lnWeight ~ lnSVL*Sex, data = Gators_S)
summary(a)
anova(a)
a <- lstrends(a, "Sex", var = "lnSVL")
a <- data.frame(pairs(a)) # No difference by sex
Appendix2.3 <- rbind(Appendix2.3, a)

#By Size Class
a <- lm(lnWeight ~ lnSVL*SizeClass, data = Groups$Northeast)
summary(a)
anova(a)
a <- lstrends(a, "SizeClass", var = "lnSVL")
a <- data.frame(pairs(a)) # No evidence by Size Class

#By season
a <- lm(lnWeight ~ lnSVL*Season, data = Groups$Northeast)
summary(a)
anova(a)
a <- lstrends(a, "Season", var = "lnSVL")
a <- data.frame(pairs(a)) # No differences by season
Appendix2.3 <- rbind(Appendix2.3, a)

#By WY
a <- lm(lnWeight ~ lnSVL*as.factor(WY), data =Groups$Northeast)
summary(a)
anova(a)
b <- confint(lstrends(a, "WY", var = "lnSVL"), adjust = "Bonferroni")
b$Group <- "Northeast"
Slope_Groups <- rbind(Slope_Groups, b)
a <- test(pairs(lstrends(a, "WY", var = "lnSVL"), adjust = "Bonferroni")) 
Appendix2.3 <- rbind(Appendix2.3, a)

Appendix2.3 <- Appendix2.3 %>%
  mutate(across(where(is.numeric), ~ round(., digits = 3)))

#No covariates
a <- lm(lnWeight ~ lnSVL, data = Groups$Northeast)
summary(a)
coefficients(a)
confint(a)

  ## 3.8. Assess effects data grouped (West)----
Appendix2.4 <- data.frame(contrast = NULL,
                          estimate = NULL,
                          SE = NULL,
                          df = NULL,
                          t.ratio = NULL,
                          p.value = NULL)

#By Sex
Gators_S <- subset(Groups$West, !Sex %in% "U")
a <- lm(lnWeight ~ lnSVL*Sex, data = Gators_S)
summary(a)
anova(a)
a <- lstrends(a, "Sex", var = "lnSVL")
a <- data.frame(pairs(a)) # No difference by sex
Appendix2.4 <- rbind(Appendix2.4, a)

#By Size Class
a <- lm(lnWeight ~ lnSVL*SizeClass, data = Groups$West)
summary(a)
anova(a)
a <- lstrends(a, "SizeClass", var = "lnSVL")
a <- data.frame(pairs(a)) # No evidence by Size Class

#By season
a <- lm(lnWeight ~ lnSVL*Season, data = Groups$West)
summary(a)
anova(a)
a <- lstrends(a, "Season", var = "lnSVL")
a <- data.frame(pairs(a)) # Weak differences by season
Appendix2.4 <- rbind(Appendix2.4, a)

#By WY
a <- lm(lnWeight ~ lnSVL*as.factor(WY), data =Groups$West)
summary(a)
anova(a)
b <- confint(lstrends(a, "WY", var = "lnSVL"), adjust = "Bonferroni")
b$Group <- "West"
Slope_Groups <- rbind(Slope_Groups, b)
a <- test(pairs(lstrends(a, "WY", var = "lnSVL"), adjust = "Bonferroni")) 
Appendix2.4 <- rbind(Appendix2.4, a)

#No covariates
a <- lm(lnWeight ~ lnSVL, data = Groups$West)
summary(a)
coefficients(a)
confint(a)

Appendix2.4 <- Appendix2.4 %>%
  mutate(across(where(is.numeric), ~ round(., digits = 3)))

## 4. Export the supplementary material table----
SupMat <- list(Appendix1, Appendix2, Appendix2.1, Appendix2.2, Appendix2.3, Appendix2.4)
names(SupMat) <- c("Summary Gators", "All data no outliers", "Central", "South", "Northeast", "West")
writexl::write_xlsx(SupMat, "./Supplementary material/SupMat.xlsx")

## 5. Residual analysis----
a <- lm(lnWeight ~ lnSVL, data = Gators_NO)
summary(a)
coef(a)
resid(a)

## 6. Figure 1 (No map) ----
(p.outliers + ggtitle("B") + theme(legend.position = "bottom")) / (p.RNO + ggtitle("C"))
ggsave("./Figures/Figure1.1.jpeg", plot = last_plot(), height = 3300, width = 2500, dpi = 400, 
       unit = "px")

## 7. Figure 2 ----
a <- ggplot(Gators_NO, aes(x = lnSVL, y = lnWeight, color = Groups)) +
  geom_abline(slope = 3, intercept = -10.8, col = "black", linewidth = 1) +
  geom_abline(slope = 3.12, intercept = -11.32, col = "yellow3", linewidth = 1) +
  geom_smooth(formula = y ~ x, method = "lm") +
  geom_point(alpha = 0.05) +
  stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = T, coef.digits = 3, f.digits = 3, p.digits = 3) +
  annotate(geom = "text", x = 4.9, y = 2.05, label = "CI = 3.14 - 3.19", col = "#F8766D") +
  annotate(geom = "text", x = 4.9, y = 1.84, label = "CI = 2.86 - 2.94", col = "#7CAE00") +
  annotate(geom = "text", x = 4.9, y = 1.63, label = "CI = 2.99 - 3.12", col = "#00BFC4") +
  annotate(geom = "text", x = 4.9, y = 1.42, label = "CI = 3.19 - 3.34", col = "#C77CFF") +
  labs(title = "A", y = "ln Weight", x = "ln Snout-Vent Length") +
  theme_bw() +
  theme(legend.position = "none")

b <- ggplot(Slope_Groups, aes(y = lnSVL.trend, x = as.factor(WY), group = Group, 
                                          colour = Group)) +
  geom_point() +
  geom_line(linewidth = 0.7) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = .3, alpha = 4/10) +
  labs(title = "B", y = "Marginal Mean Slope", x = "Water Year") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")

a / (b + theme(legend.position = "top"))

ggsave("./Figures/Figure2.jpeg", plot = last_plot(), height = 2000, width = 1500, dpi = 250, units = "px")

## 8. Figure 3----
a <- Gators_NO_Central %>%
  group_by(WY, Route) %>%
  summarise(
    Allometric_elevation = mean(Allometric_elevation),
  ) %>%
  
  ggplot() +
  geom_point(mapping = aes(y = Allometric_elevation, x = as.factor(WY), 
                           group = Route, colour = Route)) +
  geom_line(mapping = aes(y = Allometric_elevation, x = as.factor(WY), 
                          group = Route, colour = Route), linewidth = 0.7) +
  labs(title = "A", y = "Mean Allometric Elevation", x = "Water Year", color = "Sampling site") +
  scale_color_discrete(labels = c("ENP-FC", "WCA2A", "WCA3A-HD", "WCA3A-N41", "WCA3A-TW", "WCA3B")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

b <- Gators_NO_Central %>%
  group_by(WY, Season) %>%
  summarise(
    Allometric_elevation = mean(Allometric_elevation),
  ) %>%
  ggplot() +
  geom_point(mapping = aes(y = Allometric_elevation, x = as.factor(WY), 
                           group = Season, colour = Season)) +
  geom_line(mapping = aes(y = Allometric_elevation, x = as.factor(WY), 
                          group = Season, colour = Season), linewidth = 0.7) +
  labs(title = "B", y = "Mean Allometric Elevation", x = "Water Year") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

a / b

ggsave("./Figures/Figure3.jpeg", plot = last_plot(), height = 1500, width = 2000, dpi = 250, 
       units = "px")

## 9. Figure 4----
Gators_NO_Central %>%
  subset(!Sex %in% "U") %>%
  ggplot(aes(x = Allometric_elevation, color = as.factor(Sex))) +
  geom_histogram(fill = "white", binwidth = 0.05) +
  labs(y = "American alligator count", x = "Body Condition") +
  geom_vline(aes(xintercept = mean(Allometric_elevation)), color = "blue", linetype = "dashed", size = 1) +
  guides(color = guide_legend(title = "Sex")) +
  theme_bw()

Gators_NO_Central %>%
  subset(!Sex %in% "U") %>%
  summarise(
    mean(Allometric_elevation),
    min(Allometric_elevation),
    max(Allometric_elevation)
  )

ggsave("./Figures/Figure4.jpeg", plot = last_plot(), height = 1500, width = 2500, units = "px",
       dpi = 250)