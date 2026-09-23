# ======================================================================================
#                               DATA ANALYSIS SCRIPT
#         Brain State Dynamics in Ketamine-Induced Dissociation Resemble Those in PTSD
# ======================================================================================

# -----------------------------------------
# 1. Load Necessary Libraries
# -----------------------------------------

# Data Manipulation and Preparation Libraries
library(readxl)    # For reading Excel files
library(tidyr)     # For tidying data
library(dplyr)     # For data manipulation

# Statistical Analysis and Tests Libraries
library(afex)      # Methods for analyzing factorial experiments
library(emmeans)   # Tools for estimated marginal means (LSMs)
library(Hmisc)     # Miscellaneous functions for analysis and presentation of data

# Visualization Libraries
library(ggplot2)       # For creating graphics based on the Grammar of Graphics
library(pheatmap)      # For drawing pretty heatmaps
library(corrplot)      # For visualizing correlation matrices
library(ggstatsplot)   # For creating graphics with statistical details
library(gridExtra)     # For arranging multiple grid-based plots

# String Operations Library
library(stringr)       # For string manipulation

# -----------------------------------------
# 2. Define Root Directory and Subject IDs
# -----------------------------------------

# Define the root directory where your data files are located
root_dir <- "/Users/noamgoldway/Library/CloudStorage/Box-Box/Goldway, Noam/ketamine brain states/manuscript/data"

# Ketamine group subject IDs
k_subji_with_brain_data <- c(
  "0435", "0614", "0642", "0839", "0863", "1008", "1137", "1175", "1405", "1842",
  "2494", "2711", "3440", "3571", "3911", "5005", "5071", "5407", "7826", "7856",
  "7932", "8295", "8298", "8446", "8754", "9018", "9235", "9364", "9501", "9616"
)

# PTSD group subject IDs
ptsd_subji_with_brain_data <- c(
  "001", "002", "007", "008", "010", "012", "013", "016", "018", "023",
  "031", "034", "036", "038", "039", "040", "043", "045", "047", "049",
  "054", "056", "059", "060", "061", "063", "064", "584", "586", "589",
  "590", "591", "594", "595", "596", "599", "601", "604", "606", "607",
  "608", "610", "612", "616", "618", "619", "621", "622", "625", "630",
  "633", "834", "836", "837", "838", "950", "534", "535", "537", "538",
  "540", "541", "542", "543", "544", "545", "546", "548", "549", "554",
  "555", "556", "557", "558", "564", "565", "568", "569"
)

# -----------------------------------------
# 3. Load and Process Demographic Data
# -----------------------------------------

## 3.1 Ketamine Group Demographics

# Load ketamine demographic data
demog_k <- read.csv(file.path(root_dir, 'demog_k.csv'))

# Ensure subji values have leading zeros if necessary
demog_k$subji <- ifelse(nchar(demog_k$subji) == 3, paste0("0", demog_k$subji), demog_k$subji)

# Convert to character type and trim whitespace
demog_k$subji <- trimws(as.character(demog_k$subji))
k_subji_with_brain_data <- trimws(as.character(k_subji_with_brain_data))

# Filter to keep only subjects with brain data
demog_k <- demog_k[demog_k$subji %in% k_subji_with_brain_data, ]

# Calculate and print demographic statistics
mean_age_k <- mean(demog_k$age, na.rm = TRUE)
sd_age_k <- sd(demog_k$age, na.rm = TRUE)
gender_count_k <- table(demog_k$gender)

cat("Ketamine Group:\n")
cat("Mean age:", mean_age_k, "\n")
cat("Standard deviation of age:", sd_age_k, "\n")
cat("Gender count:\n")
print(gender_count_k)

## 3.2 PTSD Group Demographics

# Load PTSD demographic data
demog_ptsd <- read.csv(file.path(root_dir, 'demog_ptsd.csv'))

# Ensure subji values have exactly 3 digits
demog_ptsd$subji <- sprintf("%03s", demog_ptsd$subji)

# Convert to character type and trim whitespace
demog_ptsd$subji <- trimws(as.character(demog_ptsd$subji))
ptsd_subji_with_brain_data <- trimws(as.character(ptsd_subji_with_brain_data))

# Filter to keep only subjects with brain data
demog_ptsd <- demog_ptsd[demog_ptsd$subji %in% ptsd_subji_with_brain_data, ]

# Calculate and print demographic statistics
mean_age_ptsd <- mean(demog_ptsd$age, na.rm = TRUE)
sd_age_ptsd <- sd(demog_ptsd$age, na.rm = TRUE)
gender_count_ptsd <- table(demog_ptsd$gender)

cat("\nPTSD Group:\n")
cat("Mean age:", mean_age_ptsd, "\n")
cat("Standard deviation of age:", sd_age_ptsd, "\n")
cat("Gender count:\n")
print(gender_count_ptsd)

# -----------------------------------------
# 4. Load and Process CADSS Data (Ketamine Group)
# -----------------------------------------

# Load CADSS data
CADSS_wide <- read.csv(file.path(root_dir, 'CADSS.csv'))

# Ensure subji values have leading zeros if necessary
CADSS_wide$subji <- ifelse(nchar(CADSS_wide$subji) == 3, paste0("0", CADSS_wide$subji), CADSS_wide$subji)

# Filter to keep only subjects with brain data
CADSS_wide <- CADSS_wide[CADSS_wide$subji %in% k_subji_with_brain_data, ]

# Reshape CADSS data to long format
CADSS <- CADSS_wide %>%
  pivot_longer(
    cols = -subji,
    names_to = c("scale", "sub_scale", "session", "time_point"),
    names_sep = "_",
    values_to = "rating"
  )

# Data cleaning and factor conversion
CADSS <- CADSS %>%
  mutate(
    sub_scale = as.factor(sub_scale),
    session = ifelse(session == '2', 'Placebo', 'Ketamine'),
    time_point = case_when(
      time_point == 'baseline' ~ 'Pre-infusion',
      time_point == 'post' ~ 'Post-bolus',
      time_point == 'end' ~ 'End of infusion',
      TRUE ~ NA_character_
    )
  ) %>%
  drop_na()

CADSS$time_point <- factor(CADSS$time_point, levels = c("Pre-infusion", "Post-bolus", "End of infusion"))
CADSS$session <- factor(CADSS$session, levels = c("Ketamine", "Placebo"))

# Calculate total CADSS score per subject per time point
CADSS_total <- CADSS %>%
  group_by(subji, session, time_point) %>%
  summarise(total_score = sum(rating, na.rm = TRUE), .groups = 'drop')

# -----------------------------------------
# 5. Plot CADSS Total Scores
# -----------------------------------------

# Create CADSS total score plot
CADSS_total_plot <- ggplot(CADSS_total, aes(x = time_point, y = total_score, fill = session)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(0.9), color = "black", size = 0.7) +
  geom_point(position = position_jitterdodge(0.25, 0.9), size = 3, alpha = 0.7) +
  scale_fill_brewer(palette = "Pastel1", name = "Session") +
  labs(
    title = "CADSS Total Scores",
    x = "Time Point",
    y = "Total Score"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.title = element_text(face = "bold", size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 19),
    strip.text = element_text(face = "bold", size = 26)
  )

# Print the plot
print(CADSS_total_plot)

# -----------------------------------------
# 6. Statistical Analysis: CADSS Data
# -----------------------------------------

# Mixed-effects model for CADSS total scores
CADSS_m <- mixed(
  total_score ~ time_point * session + (1 | subji),
  data = CADSS_total,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e9))
)
nice(CADSS_m, ddf = "Kenward-Roger")

# Post-hoc pairwise comparisons
CADSS_effects <- emmeans(CADSS_m, ~ session | time_point, adjust = "none")
contrast(CADSS_effects, method = "pairwise")

# -----------------------------------------
# 7. Load and Process Clinical Measures (PTSD Group)
# -----------------------------------------

# Load clinical data
file_path <- file.path(root_dir, "fo_results+des.xlsx")
fo_clinic <- read_excel(file_path, sheet = "ptsd")

# Clean subject IDs
PTSD_subid <- gsub("[^0-9]", "", fo_clinic$`subji TP1`)

# Pivot data to long format
fo_clinic_long <- fo_clinic %>%
  pivot_longer(
    cols = -c("subji TP1", "subji TP2"),
    names_to = c("tp", ".value"),
    names_pattern = "^(TP\\d)(.*)"
  ) %>%
  mutate(
    tp = recode(tp, "TP1" = "1", "TP2" = "2"),
    subj_id = if_else(tp == "1", `subji TP1`, `subji TP2`)
  ) %>%
  select(-c("subji TP1", "subji TP2"))

# Trim whitespace from column names
names(fo_clinic_long) <- trimws(names(fo_clinic_long))

# Select relevant clinical measures
ptsd_clinic_long <- fo_clinic_long %>%
  select(tp, caps5, Diss, subj_id) %>%
  mutate(
    tp = recode(tp, `1` = "Pre-treatment", `2` = "Post-treatment"),
    tp = factor(tp, levels = c("Pre-treatment", "Post-treatment"))
  ) %>%
  arrange(tp)

# -----------------------------------------
# 8. Plot Clinical Measures
# -----------------------------------------

# Reshape data for plotting
fo_clinic_long_melted <- ptsd_clinic_long %>%
  pivot_longer(cols = c(caps5, Diss), names_to = "Measure", values_to = "Value")

# Define labels for measures
measure_labels <- c(caps5 = "CAPS-5 Total Score", Diss = "CAPS-5 Dissociation")

# Create plot for clinical measures
ptsd_clinic_p <- ggplot(fo_clinic_long_melted, aes(x = tp, y = Value, fill = tp)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(0.8), color = "black", size = 0.7) +
  geom_point(position = position_jitterdodge(0.25, 0.8), size = 3, alpha = 0.7) +
  facet_wrap(~ Measure, scales = "free_y", labeller = labeller(Measure = measure_labels)) +
  scale_fill_brewer(palette = "Pastel1", name = "Time Point") +
  labs(
    x = "Time Point",
    y = "Score"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold", size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 26)
  )

# Print the plot
print(ptsd_clinic_p)

# -----------------------------------------
# 9. Statistical Analysis: Clinical Measures
# -----------------------------------------

# Mixed-effects model for CAPS-5 total score
caps_m <- mixed(
  caps5 ~ tp + (1 | subj_id),
  data = ptsd_clinic_long,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e9))
)
nice(caps_m, ddf = "Kenward-Roger")

# Wilcoxon signed-rank test for Dissociation scores
diss_tp1 <- filter(ptsd_clinic_long, tp == 'Pre-treatment')
diss_tp2 <- filter(ptsd_clinic_long, tp == 'Post-treatment')

# Merge data on subj_id
diss_data <- inner_join(
  select(diss_tp1, subj_id, Diss),
  select(diss_tp2, subj_id, Diss),
  by = "subj_id",
  suffix = c("_pre", "_post")
) %>%
  drop_na()

# Perform the Wilcoxon signed-rank test
wilcox_test <- wilcox.test(diss_data$Diss_pre, diss_data$Diss_post, paired = TRUE, exact = FALSE, correct = TRUE)
print(wilcox_test)

# -----------------------------------------
# 10. Load and Process FO Data (Ketamine Group)
# -----------------------------------------

# Load FO data for ketamine group
fo_ket <- read_excel(file_path, sheet = "ket")

# Pivot data to long format
fo_ket_long <- fo_ket %>%
  pivot_longer(
    cols = -subji,
    names_to = c("session", "meta_state", "sign"),
    names_pattern = "([^ ]+) ([A-Z]+)([-+])",
    values_to = "value"
  )

# Calculate mean FO values per meta-state per session
fo_ket_mean <- fo_ket_long %>%
  group_by(subji, session, meta_state) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = 'drop') %>%
  mutate(session = factor(session, levels = c("ketamine", "placebo")))

# -----------------------------------------
# 11. Plot FO Data: Ketamine vs Placebo
# -----------------------------------------

# Create FO plot for ketamine vs placebo
fo_ket_plot <- ggplot(fo_ket_mean, aes(x = session, y = mean_value, fill = session)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(0.8), color = "black", size = 0.7) +
  geom_point(position = position_jitterdodge(0.25, 0.8), size = 3, alpha = 0.7) +
  facet_wrap(~ meta_state, scales = "free_y", labeller = as_labeller(c(SOM = "SOM", VIS = "VIS", DMN = "DMN"))) +
  scale_fill_brewer(palette = "Pastel1", name = "Session") +
  labs(
    title = "Ketamine vs Placebo",
    x = "",
    y = "Fractional Occupancy"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold", size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 26)
  )

# Print the plot
print(fo_ket_plot)

# -----------------------------------------
# 12. Statistical Analysis: FO Data (Ketamine)
# -----------------------------------------

# Mixed-effects model for FO data
ket_fo_m <- mixed(
  mean_value ~ session * meta_state + (1 | subji),
  data = fo_ket_mean,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e9))
)
nice(ket_fo_m, ddf = "Kenward-Roger")

# Post-hoc pairwise comparisons
ket_fo_post_hoc <- emmeans(ket_fo_m, ~ session | meta_state, adjust = "none")
contrast(ket_fo_post_hoc, method = "pairwise")

# -----------------------------------------
# 13. Process FO Data (PTSD Group)
# -----------------------------------------

# Select relevant FO data for PTSD group
ptsd_FO <- fo_clinic_long %>%
  select(tp, `VIS-`, `SOM+`, `SOM-`, `DMN-`, `DMN+`, `VIS+`, subj_id) %>%
  mutate(tp = recode(tp, `1` = "Pre-treatment", `2` = "Post-treatment"))

# Calculate mean FO values per meta-state per time point
ptsd_FO_means <- ptsd_FO %>%
  group_by(subj_id, tp) %>%
  mutate(
    SOM = (`SOM+` + `SOM-`) / 2,
    VIS = (`VIS+` + `VIS-`) / 2,
    DMN = (`DMN+` + `DMN-`) / 2
  ) %>%
  select(subj_id, tp, SOM, VIS, DMN) %>%
  ungroup()

# Convert data to long format
ptsd_FO_long <- ptsd_FO_means %>%
  pivot_longer(cols = SOM:DMN, names_to = "Meta_state", values_to = "Value") %>%
  mutate(tp = factor(tp, levels = c("Pre-treatment", "Post-treatment")))

# -----------------------------------------
# 14. Plot FO Data: PTSD Pre vs Post Treatment
# -----------------------------------------

# Define custom labels
measure_labels <- c(SOM = "SOM", VIS = "VIS", DMN = "DMN")

# Create FO plot for PTSD group
ptsd_FO_plot <- ggplot(ptsd_FO_long, aes(x = tp, y = Value, fill = tp)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(0.8), color = "black", size = 0.7) +
  geom_point(position = position_jitterdodge(0.25, 0.8), size = 3, alpha = 0.7) +
  facet_wrap(~ Meta_state, scales = "free_y", labeller = labeller(Meta_state = measure_labels)) +
  scale_fill_brewer(palette = "Pastel1", name = "Time Point") +
  labs(
    title = "PTSD: Pre vs Post Treatment",
    x = "",
    y = "Fractional Occupancy"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold", size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 26)
  )

# Print the plot
print(ptsd_FO_plot)

# -----------------------------------------
# 15. Statistical Analysis: FO Data (PTSD)
# -----------------------------------------

# Mixed-effects model for FO data
ptsd_fo_m <- mixed(
  Value ~ tp * Meta_state + (1 | subj_id),
  data = ptsd_FO_long,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e9))
)
nice(ptsd_fo_m, ddf = "Kenward-Roger")

# Post-hoc pairwise comparisons
ptsd_fo_post_hoc <- emmeans(ptsd_fo_m, ~ tp | Meta_state, adjust = "none")
contrast(ptsd_fo_post_hoc, method = "pairwise")

# -----------------------------------------
# 16. Combined Mixed Model Analysis
# -----------------------------------------

# Prepare PTSD data for combined analysis
ptsd_FO_long_with_source <- ptsd_FO_long %>%
  mutate(
    session = NA,
    source = "ptsd"
  )

# Prepare Ketamine data for combined analysis
fo_ket_mean_with_source <- fo_ket_mean %>%
  rename(subj_id = subji, Value = mean_value) %>%
  mutate(
    tp = if_else(session == "ketamine", "Pre-treatment", "Post-treatment"),
    source = "ket",
    Meta_state = meta_state
  ) %>%
  select(subj_id, tp, Meta_state, Value, session, source)

# Merge datasets
merged_data <- bind_rows(ptsd_FO_long_with_source, fo_ket_mean_with_source) %>%
  filter(!is.na(Meta_state) & !is.na(Value)) %>%
  mutate(
    tp = factor(tp, levels = c("Pre-treatment", "Post-treatment")),
    Meta_state = factor(Meta_state, levels = c("SOM", "VIS", "DMN")),
    source = factor(source, levels = c("ptsd", "ket")),
    subj_id = factor(subj_id)
  ) %>%
  droplevels()

# Mixed-effects model for combined data
ptsd_ket_fo_m <- mixed(
  Value ~ tp * Meta_state * source + (1 | subj_id),
  data = merged_data,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e9))
)
nice(ptsd_ket_fo_m, ddf = "Kenward-Roger")

# Post-hoc pairwise comparisons
ptsd_ket_fo_post_hoc <- emmeans(ptsd_ket_fo_m, ~ tp | Meta_state, adjust = "none")
contrast(ptsd_ket_fo_post_hoc, method = "pairwise")

# -----------------------------------------
# 17. Correlation Analysis between Clinical Measures and FO (PTSD Group)
# -----------------------------------------

# Clean column names
ptsd_clinic_long <- ptsd_clinic_long %>% janitor::clean_names()

# Calculate differences in Dissociation scores
ptsd_clinic_long <- ptsd_clinic_long %>%
  group_by(subj_id) %>%
  mutate(
    diss_pre = first(diss[tp == "Pre-treatment"]),
    diss_post = first(diss[tp == "Post-treatment"]),
    diss_diff = if_else(!is.na(diss_pre) & !is.na(diss_post), diss_pre - diss_post, NA_real_),
    diss_bl = if_else(any(diss[tp == "Pre-treatment"] > 0), 1, 0)
  ) %>%
  ungroup()

# Print updated data
print(select(ptsd_clinic_long, subj_id, tp, diss, diss_pre, diss_post, diss_diff, diss_bl))

# Calculate FO differences
ptsd_FO_long_diff <- ptsd_FO_long %>%
  group_by(subj_id, Meta_state) %>%
  summarise(
    Value_diff = Value[tp == "Pre-treatment"] - Value[tp == "Post-treatment"],
    .groups = 'drop'
  ) %>%
  filter(Meta_state != "VIS") %>%
  pivot_wider(names_from = Meta_state, values_from = Value_diff)

# Merge clinical and FO data
ptsd_fo_clinic_diff <- left_join(ptsd_clinic_long, ptsd_FO_long_diff, by = "subj_id") %>%
  filter(tp == "Pre-treatment" & diss_bl == 1) %>%
  select(diss_diff, DMN, SOM) %>%
  drop_na()

# Perform Spearman correlation
cor_diss_diff_DMN <- cor.test(ptsd_fo_clinic_diff$diss_diff, ptsd_fo_clinic_diff$DMN, method = "spearman")
cor_diss_diff_SOM <- cor.test(ptsd_fo_clinic_diff$diss_diff, ptsd_fo_clinic_diff$SOM, method = "spearman")

# Print correlation results
cat("Spearman correlation between CAPS-5 Dissociation Change and DMN FO Change:\n")
cat("Correlation coefficient (rho):", cor_diss_diff_DMN$estimate, "\n")
cat("p-value:", cor_diss_diff_DMN$p.value, "\n\n")

cat("Spearman correlation between CAPS-5 Dissociation Change and SOM FO Change:\n")
cat("Correlation coefficient (rho):", cor_diss_diff_SOM$estimate, "\n")
cat("p-value:", cor_diss_diff_SOM$p.value, "\n")

# Create correlation plot
cor_diss_diff_DMN_plot <- ggplot(ptsd_fo_clinic_diff, aes(x = diss_diff, y = DMN)) +
  geom_point(position = position_jitter(width = 0.25), size = 3, alpha = 0.3, color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(
    x = "CAPS-5 Dissociation Difference (Pre - Post)\nHigher Values = Reduction in Symptoms",
    y = "DMN FO Difference (Pre - Post)\nHigher Values = Reduction in DMN FO"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold", size = 30),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 24),
    axis.text.y = element_text(size = 24),
    legend.position = "none"
  )

# Print the plot
print(cor_diss_diff_DMN_plot)

# -----------------------------------------
# 18. Analyze Brain State Transition Energy Data
# -----------------------------------------

# Load energy data
energy_data <- read.csv(file.path(root_dir, 'energy.csv'), check.names = FALSE)

# Reshape data to long format
energy_data_long <- energy_data %>%
  pivot_longer(
    cols = -c(subji, tp, group, From),
    names_to = "to",
    values_to = "value"
  ) %>%
  mutate(
    subji = trimws(subji),
    From = as.character(From),
    to = as.character(to)
  ) %>%
  arrange(subji, From, to, group, tp)

# Remove subject with abnormal imaging finding
energy_data_long <- filter(energy_data_long, subji != "K9653")

# Filter out 'stay' transitions
energy_data_long_no_stay <- filter(energy_data_long, From != to)

# Convert columns to factors
energy_data_long_no_stay <- energy_data_long_no_stay %>%
  mutate(
    From = as.factor(From),
    to = as.factor(to),
    group = as.factor(group),
    subji = as.factor(subji)
  )

# Split data into ketamine and PTSD groups
energy_data_ketamine <- filter(energy_data_long_no_stay, group == "ketamine") %>%
  mutate(value_z = scale(value))

energy_data_ptsd <- filter(energy_data_long_no_stay, group == "ptsd") %>%
  mutate(value_z = scale(value))

## Ketamine Group Energy Analysis

# Calculate mean energy per start state
mean_value_ketamine <- energy_data_ketamine %>%
  mutate(From = str_remove(From, "[-+]")) %>%
  group_by(subji, tp, From, group) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    mean_value_z = mean(value_z, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  rename(start_state = From, session = tp)

# Calculate group means and standard errors
mean_ketamine <- mean_value_ketamine %>%
  group_by(session, start_state) %>%
  summarise(
    mean_z = mean(mean_value_z, na.rm = TRUE),
    se = sd(mean_value_z, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  mutate(session = factor(session, levels = c(1, 2), labels = c("placebo", "ketamine"))) %>%
  drop_na()


# Mixed-effects model for energy data (ketamine)
ketamine_e_m <- mixed(
  mean_value_z ~ start_state * session + (1 | subji),
  data = mean_value_ketamine,
  control = list(optimizer = "bobyqa", optCtrl = list(maxfun = 1e9))
)
nice(ketamine_e_m)

# Post-hoc pairwise comparisons
ketamine_e_post_hoc <- emmeans(ketamine_e_m, ~ session | start_state, adjust = "none")
contrast(ketamine_e_post_hoc, method = "pairwise")

## PTSD Group Energy Analysis

# Calculate mean energy per start state
mean_value_ptsd <- energy_data_ptsd %>%
  mutate(From = str_remove(From, "[-+]")) %>%
  group_by(subji, tp, From, group) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    mean_value_z = mean(value_z, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  rename(start_state = From)

# Calculate group means and standard errors
mean_ptsd <- mean_value_ptsd %>%
  group_by(tp, start_state) %>%
  summarise(
    mean_z = mean(mean_value_z, na.rm = TRUE),
    se = sd(mean_value_z, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

mean_ptsd$tp <- as.factor(mean_ptsd$tp)

# Mixed-effects model for energy data (PTSD)
ptsd_e_m <- mixed(
  mean_value_z ~ start_state * tp + (1 | subji),
  data = mean_value_ptsd,
  control = list(optimizer = "bobyqa", optCtrl = list(maxfun = 1e9))
)
nice(ptsd_e_m)

# Post-hoc pairwise comparisons
ptsd_e_post_hoc <- emmeans(ptsd_e_m, ~ tp | start_state, adjust = "none")
contrast(ptsd_e_post_hoc, method = "pairwise")

# -----------------------------------------
# 19. Plot Combined Energy Differences
# -----------------------------------------

# Load necessary libraries
library(gridExtra)

# Define y-axis limits
y_limits <- c(-0.15, 0.3)

# PTSD Difference Plot
mean_ptsd_wide <- mean_ptsd %>%
  pivot_wider(names_from = tp, values_from = c(mean_z, se), names_prefix = "TP") %>%
  mutate(
    mean_diff = mean_z_TP2 - mean_z_TP1,
    se_diff = sqrt(se_TP1^2 + se_TP2^2)
  )

ptsd_diff_plot <- ggplot(mean_ptsd_wide, aes(x = start_state, y = mean_diff)) +
  geom_bar(stat = "identity", fill = "lightgray", color = "black", size = 0.7) +
  geom_errorbar(aes(ymin = mean_diff - se_diff, ymax = mean_diff + se_diff), width = 0.25, color = "black") +
  labs(
    title = "PTSD Post vs Pre Treatment",
    x = "Start State",
    y = "Difference (Post - Pre Treatment)"
  ) +
  scale_y_continuous(limits = y_limits) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
    axis.title = element_text(face = "bold", size = 24),
    axis.text.x = element_text(size = 24),
    axis.text.y = element_text(size = 24)
  )

# Ketamine Difference Plot
mean_ketamine_wide <- mean_ketamine %>%
  pivot_wider(names_from = session, values_from = c(mean_z, se)) %>%
  mutate(
    mean_diff = mean_z_placebo - mean_z_ketamine,
    se_diff = sqrt(se_placebo^2 + se_ketamine^2)
  )

ketamine_diff_plot <- ggplot(mean_ketamine_wide, aes(x = start_state, y = mean_diff)) +
  geom_bar(stat = "identity", fill = "white", color = "black", size = 0.7) +
  geom_errorbar(aes(ymin = mean_diff - se_diff, ymax = mean_diff + se_diff), width = 0.25, color = "black") +
  labs(
    title = "Ketamine vs Placebo",
    x = "Start State",
    y = "Difference (Ketamine - Placebo)"
  ) +
  scale_y_continuous(limits = y_limits) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
    axis.title = element_text(face = "bold", size = 24),
    axis.text.x = element_text(size = 24),
    axis.text.y = element_text(size = 24)
  )

# Print both plots side by side with the same y-axis limits
grid.arrange(ketamine_diff_plot, ptsd_diff_plot, ncol = 2)
