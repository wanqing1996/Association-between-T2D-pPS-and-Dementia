#1.Compute Overall T2D PRS (Including APOE)
auto_PRS_calc_v3 \
    -g <path_to_genotype_data>/merged_imputed_info0.3_maf0.01_2638ref \
    -s <path_to_summary_statistics>/T2D_GWAS_All_Metal_LDSC-CORR_Neff.v2.hg38.txt \
    -o T2D_overall \
    -r 0.2 \
    --colname_chr Chromsome \
    --colname_pos pos_hg38 \
    --colname_ea EffectAllele \
    --colname_nea NonEffectAllele \
    --colname_beta Beta \
    --colname_p Pval

#2.Compute T2D PRS (Excluding APOE rs429358)
awk '$2 != "rs429358"' <path_to_summary_statistics>/T2D_GWAS_All_Metal_LDSC-CORR_Neff.v2.hg38.txt > T2D_GWAS_noAPOE.txt
auto_PRS_calc_v3 \
    -g <path_to_genotype_data>/merged_imputed_info0.3_maf0.01_2638ref \
    -s T2D_GWAS_noAPOE.txt \
    -o T2D_noAPOE \
    -r 0.2 \
    --colname_chr Chromsome \
    --colname_pos pos_hg38 \
    --colname_ea EffectAllele \
    --colname_nea NonEffectAllele \
    --colname_beta Beta \
    --colname_p Pval

3.Specific PRS Computation Commands for 12 Genetic Clusters
#3. 1. ALP Negative pPS
auto_PRS_calc_v3 -g <path_to_genotype_data>/merged_imputed_info0.3_maf0.01_2638ref \
                 -s <path_to_summary_statistics>/ALP_Negative_pPS.txt \
                 -o ALP_Negative \
                 -r 0.2 -p 0.01 \
                 --colname_chr=chr --colname_pos=pos \
                 --colname_ea=EffectAllele --colname_nea=NonEffectAllele \
                 --colname_beta=Beta --colname_p=Pval
# 3.2. Beta Cell 1 pPS
auto_PRS_calc_v3 -g <path_to_genotype_data>/merged_imputed_info0.3_maf0.01_2638ref \
                 -s <path_to_summary_statistics>/Beta_Cell_1_pPS.txt \
                 -o Beta_Cell_1 \
                 -r 0.2 -p 0.01 \
                 --colname_chr=chr --colname_pos=pos \
                 --colname_ea=EffectAllele --colname_nea=NonEffectAllele \
                 --colname_beta=Beta --colname_p=Pval
# 3.3. Beta Cell 2 pPS
auto_PRS_calc_v3 -g <path_to_genotype_data>/merged_imputed_info0.3_maf0.01_2638ref \
                 -s <path_to_summary_statistics>/Beta_Cell_2_pPS.txt \
                 -o Beta_Cell_2 \
                 -r 0.2 -p 0.01 \
                 --colname_chr=chr --colname_pos=pos \
                 --colname_ea=EffectAllele --colname_nea=NonEffectAllele \
                 --colname_beta=Beta --colname_p=Pval
#3. 4. Bilirubin pPS
auto_PRS_calc_v3 -g <path_to_genotype_data>/merged_imputed_info0.3_maf0.01_2638ref \
                 -s <path_to_summary_statistics>/Bilirubin_pPS.txt \
                 -o Bilirubin \
                 -r 0.2 -p 0.01 \
                 --colname_chr=chr --colname_pos=pos \
                 --colname_ea=EffectAllele --colname_nea=NonEffectAllele \
                 --colname_beta=Beta --colname_p=Pval
#3. 5. Cholesterol pPS
auto_PRS_calc_v3 -g <path_to_genotype_data>/merged_imputed_info0.3_maf0.01_2638ref \
                 -s <path_to_summary_statistics>/Cholesterol_pPS.txt \
                 -o Cholesterol \
                 -r 0.2 -p 0.01 \
                 --colname_chr=chr --colname_pos=pos \
                 --colname_ea=EffectAllele --colname_nea=NonEffectAllele \
                 --colname_beta=Beta --colname_p=Pval
# 3.6. Hyper Insulin pPS
auto_PRS_calc_v3 -g <path_to_genotype_data>/merged_imputed_info0.3_maf0.01_2638ref \
                 -s <path_to_summary_statistics>/Hyper_Insulin_pPS.txt \
                 -o Hyper_Insulin \
                 -r 0.2 -p 0.01 \
                 --colname_chr=chr --colname_pos=pos \
                 --colname_ea=EffectAllele --colname_nea=NonEffectAllele \
                 --colname_beta=Beta --colname_p=Pval
# 3.7. Lipodystrophy 1 pPS
auto_PRS_calc_v3 -g <path_to_genotype_data>/merged_imputed_info0.3_maf0.01_2638ref \
                 -s <path_to_summary_statistics>/Lipodystrophy_1_pPS.txt \
                 -o Lipodystrophy_1 \
                 -r 0.2 -p 0.01 \
                 --colname_chr=chr --colname_pos=pos \
                 --colname_ea=EffectAllele --colname_nea=NonEffectAllele \
                 --colname_beta=Beta --colname_p=Pval
# 3.8. Lipodystrophy 2 pPS
auto_PRS_calc_v3 -g <path_to_genotype_data>/merged_imputed_info0.3_maf0.01_2638ref \
                 -s <path_to_summary_statistics>/Lipodystrophy_2_pPS.txt \
                 -o Lipodystrophy_2 \
                 -r 0.2 -p 0.01 \
                 --colname_chr=chr --colname_pos=pos \
                 --colname_ea=EffectAllele --colname_nea=NonEffectAllele \
                 --colname_beta=Beta --colname_p=Pval
# 3.9. Liver-Lipid pPS
auto_PRS_calc_v3 -g <path_to_genotype_data>/merged_imputed_info0.3_maf0.01_2638ref \
                 -s <path_to_summary_statistics>/Liver_Lipid_pPS.txt \
                 -o Liver_Lipid \
                 -r 0.2 -p 0.01 \
                 --colname_chr=chr --colname_pos=pos \
                 --colname_ea=EffectAllele --colname_nea=NonEffectAllele \
                 --colname_beta=Beta --colname_p=Pval
# 3.10. Obesity pPS
auto_PRS_calc_v3 -g <path_to_genotype_data>/merged_imputed_info0.3_maf0.01_2638ref \
                 -s <path_to_summary_statistics>/Obesity_pPS.txt \
                 -o Obesity \
                 -r 0.2 -p 0.01 \
                 --colname_chr=chr --colname_pos=pos \
                 --colname_ea=EffectAllele --colname_nea=NonEffectAllele \
                 --colname_beta=Beta --colname_p=Pval
# 3.11. Proinsulin pPS
auto_PRS_calc_v3 -g <path_to_genotype_data>/merged_imputed_info0.3_maf0.01_2638ref \
                 -s <path_to_summary_statistics>/Proinsulin_pPS.txt \
                 -o Proinsulin \
                 -r 0.2 -p 0.01 \
                 --colname_chr=chr --colname_pos=pos \
                 --colname_ea=EffectAllele --colname_nea=NonEffectAllele \
                 --colname_beta=Beta --colname_p=Pval
# 3.12. SHBG-LpA pPS
auto_PRS_calc_v3 -g <path_to_genotype_data>/merged_imputed_info0.3_maf0.01_2638ref \
                 -s <path_to_summary_statistics>/SHBGLpA_pPS.txt \
                 -o SHBGLpA \
                 -r 0.2 -p 0.01 \
                 --colname_chr=chr --colname_pos=pos \
                 --colname_ea=EffectAllele --colname_nea=NonEffectAllele \
                 --colname_beta=Beta --colname_p=Pval

#4.Table 1 Association of Polygenic Risk Scores for Type 2 Diabetes With all-cause dementia, Alzheimer's disease, Vascular dementia and Other & unspecified dementia Risk
library(survival)
data <- read.csv("data.csv")
run_cox_model <- function(event_col, prs_col, covariates) {  
    surv_obj <- Surv(time = data$follow_up_years, event = data[[event_col]])  
    cox_model <- coxph(surv_obj ~ ., data = data[, c(event_col, prs_col, covariates)])  
    summary(cox_model)  
}  
adjust_age_sex <- c("AGE", "SEX")  
adjust_age_sex_T2D <- c("AGE", "SEX", "T2D_status")  
prs_col <- "T2D_PRS"
cox_age_sex <- run_cox_model("dementia", prs_col, adjust_age_sex)  
cox_age_sex_T2D <- run_cox_model("dementia", prs_col, adjust_age_sex_T2D)  
cox_age_sex <- run_cox_model("VD", prs_col, adjust_age_sex)  
cox_age_sex_T2D <- run_cox_model("VD", prs_col, adjust_age_sex_T2D)  
cox_age_sex <- run_cox_model("Alzheimers", prs_col, adjust_age_sex)  
cox_age_sex_T2D <- run_cox_model("Alzheimers", prs_col, adjust_age_sex_T2D)  
cox_age_sex <- run_cox_model("Other_Dementia", prs_col, adjust_age_sex)  
cox_age_sex_T2D <- run_cox_model("Other_Dementia", prs_col, adjust_age_sex_T2D)  
prs_col_non_apoe <- "T2D_PRS_non_apoe"
cox_age_sex <- run_cox_model("dementia", prs_col_non_apoe, adjust_age_sex)  
cox_age_sex_T2D <- run_cox_model("dementia", prs_col_non_apoe, adjust_age_sex_T2D)  
cox_age_sex <- run_cox_model("VD", prs_col_non_apoe, adjust_age_sex)  
cox_age_sex_T2D <- run_cox_model("VD", prs_col_non_apoe, adjust_age_sex_T2D)  
cox_age_sex <- run_cox_model("Alzheimers", prs_col_non_apoe, adjust_age_sex)  
cox_age_sex_T2D <- run_cox_model("Alzheimers", prs_col_non_apoe, adjust_age_sex_T2D)  
cox_age_sex <- run_cox_model("Other_Dementia", prs_col_non_apoe, adjust_age_sex)  
cox_age_sex_T2D <- run_cox_model("Other_Dementia", prs_col_non_apoe, adjust_age_sex_T2D)  


#5.Table 2 Association between T2D genetic clusters  and Risks of All-cause Dementia, Alzheimer's Disease, Vascular Dementia, and Other/Unspecified Dementia
library(survival)
data <- read.csv("data.csv")
run_cox_model <- function(event_col, prs_col, covariates) {  
    surv_obj <- Surv(time = data$follow_up_years, event = data[[event_col]])  
    cox_model <- coxph(surv_obj ~ ., data = data[, c(event_col, prs_col, covariates)])  
    summary(cox_model)  
}  
adjust_covariates <- c("AGE", "SEX", "APOE_e4")
genetic_clusters <- c("ALP_Negative_pPS", "Beta_Cell_1_pPS", "Beta_Cell_2_pPS", 
                      "Bilirubin_pPS", "Cholesterol_pPS", "Hyper_Insulin_pPS", 
                      "Lipodystrophy_1_pPS", "Lipodystrophy_2_pPS", "Liver_Lipid_pPS",
                      "Obesity_pPS", "Proinsulin_pPS", "SHBG_LpA_pPS")
dementia_outcomes <- c("dementia", "Alzheimers", "VD", "Other_Dementia")
for (prs in genetic_clusters) {
    for (outcome in dementia_outcomes) {
        result <- run_cox_model(outcome, prs, adjust_covariates)
        print(result)
    }
}

#6.Table 3 Associations Between clusters Tertiles for Hyperinsulinemia and SHBGLpA and Risk of Vascular Dementia Across Adjusted Models
library(survival)
data <- read.csv("data.csv")
data$HyperInsulin_tertile <- cut(data$HyperInsulin_pPS, 
                                 breaks = quantile(data$HyperInsulin_pPS, probs = c(0, 1/3, 2/3, 1)), 
                                 labels = c("Low", "Medium", "High"), 
                                 include.lowest = TRUE)
data$SHBGLpA_tertile <- cut(data$SHBGLpA_pPS, 
                            breaks = quantile(data$SHBGLpA_pPS, probs = c(0, 1/3, 2/3, 1)), 
                            labels = c("Low", "Medium", "High"), 
                            include.lowest = TRUE)
model1_cov <- c("AGE", "SEX", "EDUCATION", "MARITAL_STATUS", "RESIDENTIAL_LOCATION")
model2_cov <- c(model1_cov, "SMOKING", "ALCOHOL")
model3_cov <- c(model2_cov, "LDL_C", "HDL_C", "TOTAL_CHOLESTEROL", "TRIGLYCERIDES", "BLOOD_PRESSURE", "BMI", "FBG")
model4_cov <- c(model3_cov, "CORONARY_ARTERY_DISEASE", "CEREBROVASCULAR_DISEASE", "SLEEP_DISORDERS", "DEPRESSION")
run_cox_model <- function(event_col, tertile_col, covariates) {  
    surv_obj <- Surv(time = data$follow_up_years, event = data[[event_col]])  
    cox_model <- coxph(surv_obj ~ relevel(as.factor(data[[tertile_col]]), ref = "Low") + ., 
                       data = data[, c(event_col, tertile_col, covariates)])  
    summary(cox_model)  
}  
cox_hyperinsulin_m1 <- run_cox_model("VD", "HyperInsulin_tertile", model1_cov)
cox_hyperinsulin_m2 <- run_cox_model("VD", "HyperInsulin_tertile", model2_cov)
cox_hyperinsulin_m3 <- run_cox_model("VD", "HyperInsulin_tertile", model3_cov)
cox_hyperinsulin_m4 <- run_cox_model("VD", "HyperInsulin_tertile", model4_cov)
cox_shbglpa_m1 <- run_cox_model("VD", "SHBGLpA_tertile", model1_cov)
cox_shbglpa_m2 <- run_cox_model("VD", "SHBGLpA_tertile", model2_cov)
cox_shbglpa_m3 <- run_cox_model("VD", "SHBGLpA_tertile", model3_cov)
cox_shbglpa_m4 <- run_cox_model("VD", "SHBGLpA_tertile", model4_cov)


#7.Figure 1： Effect Sizes of Hyper Insulin on Vascular Dementia
library(ggplot2)
library(dplyr)
data <- read.csv("data.csv")
p <- ggplot(data, aes(x = BETA_CHIDAI, y = P_CHIDAI)) +
  geom_point(aes(color = -log10(P_CHIDAI)), size = 2) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) + # Threshold line at 0.05
  scale_color_gradient(low = "blue", high = "red", name = "-log10(P_CHIDAI)") +
  theme_minimal() +
  labs(
    title = "Effect Sizes of HyperInsulin pRS on Vascular Dementia",
    x = "Beta",
    y = "P value"
  ) +
  theme(
    text = element_text(size = 12),
    legend.position = "right"
  )
print(p)

#8. Table 4 Associations Between rs11129735 and rs329118 Genotypes and Risk of Vascular Dementia
library(survival)
data <- read.csv("data.csv")
model1_cov <- c("AGE", "SEX")  
model2_cov <- c(model1_cov, "T2D_status")  
model4_cov <- c(model2_cov, "EDUCATION", "MARITAL_STATUS", "RESIDENTIAL_LOCATION", 
                 "SMOKING", "ALCOHOL", "LDL_C", "HDL_C", "TOTAL_CHOLESTEROL", 
                 "TRIGLYCERIDES", "BLOOD_PRESSURE", "BMI", "FBG",
                 "CORONARY_ARTERY_DISEASE", "CEREBROVASCULAR_DISEASE", "SLEEP_DISORDERS", "DEPRESSION")  
run_cox_model <- function(event_col, genotype_col, covariates) {  
    surv_obj <- Surv(time = data$follow_up_years, event = data[[event_col]])  
    cox_model <- coxph(surv_obj ~ relevel(as.factor(data[[genotype_col]]), ref = "reference") + ., 
                       data = data[, c(event_col, genotype_col, covariates)])  
    summary(cox_model)  
}  
cox_trank1_m1 <- run_cox_model("VD", "rs11129735", model1_cov)  
cox_trank1_m2 <- run_cox_model("VD", "rs11129735", model2_cov)  
cox_trank1_m4 <- run_cox_model("VD", "rs11129735", model4_cov)  
cox_jade2_m1 <- run_cox_model("VD", "rs329118", model1_cov)  
cox_jade2_m2 <- run_cox_model("VD", "rs329118", model2_cov)  
cox_jade2_m4 <- run_cox_model("VD", "rs329118", model4_cov)  
cox_shbg_m1 <- run_cox_model("VD", "rs858519", model1_cov)  
cox_shbg_m2 <- run_cox_model("VD", "rs858519", model2_cov)  
cox_shbg_m4 <- run_cox_model("VD", "rs858519", model4_cov)  
cox_combined_m1 <- run_cox_model("VD", "Combined_Risk_Alleles", model1_cov)  
cox_combined_m2 <- run_cox_model("VD", "Combined_Risk_Alleles", model2_cov)  
cox_combined_m4 <- run_cox_model("VD", "Combined_Risk_Alleles", model4_cov)  

#9.eFigure 2：Distribution of T2D_PRS by Diabetes Status
library(ggplot2)
data <- read.csv("data.csv")
p1 <- ggplot(data, aes(x = T2D_PRS, fill = Diabetes_Status)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Distribution of T2D_PRS by Diabetes Status",
    x = "T2D_PRS",
    y = "Density",
    fill = "Diabetes Status"
  ) +
  scale_fill_manual(values = c("red", "blue"), labels = c("No diabetes", "Type 2 diabetes")) +
  theme_minimal()
p2 <- ggplot(data, aes(x = T2D_PRS_no_APOE, fill = Diabetes_Status)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Distribution of T2D_PRS (Excluding APOE) by Diabetes Status",
    x = "T2D_PRS_no_APOE",
    y = "Density",
    fill = "Diabetes Status"
  ) +
  scale_fill_manual(values = c("red", "blue"), labels = c("No diabetes", "Type 2 diabetes")) +
  theme_minimal()
library(gridExtra)
grid.arrange(p1, p2, ncol = 1)

#10.eTable 4   Associations Between Type 2 Diabetes Polygenic Risk Score (PRS) and Diabetes Risk
library(dplyr)
library(broom)
data <- read.csv("data.csv")
data$PRS_Category <- factor(data$PRS_Category, levels = c("Low", "Medium", "High"))
model_t2d_prs_continuous <- glm(T2D_Status ~ T2D_PRS, data = data, family = binomial)
model_t2d_prs_categorical <- glm(T2D_Status ~ PRS_Category, data = data, family = binomial)
model_t2d_prs_no_apoe_continuous <- glm(T2D_Status ~ T2D_PRS_non_apoe, data = data, family = binomial)
model_t2d_prs_no_apoe_categorical <- glm(T2D_Status ~ PRS_Category, data = data, family = binomial)
summary_continuous <- function(model) {
  tidy(model, exponentiate = TRUE) %>% 
    select(term, estimate, conf.low, conf.high, p.value)
}
summary_categorical <- function(model) {
  tidy(model, exponentiate = TRUE) %>% 
    filter(term != "(Intercept)") %>%
    select(term, estimate, conf.low, conf.high, p.value)
}
result_t2d_prs_continuous <- summary_continuous(model_t2d_prs_continuous)
result_t2d_prs_categorical <- summary_categorical(model_t2d_prs_categorical)
result_t2d_prs_no_apoe_continuous <- summary_continuous(model_t2d_prs_no_apoe_continuous)
result_t2d_prs_no_apoe_categorical <- summary_categorical(model_t2d_prs_no_apoe_categorical)
print(result_t2d_prs_continuous)
print(result_t2d_prs_categorical)
print(result_t2d_prs_no_apoe_continuous)
print(result_t2d_prs_no_apoe_categorical)

#11.eTable 5 Baseline Characteristics of Participants by T2D PRS
library(tableone)
df <- read.csv("data.csv")
df$T2D_PRS_Cat <- factor(df$T2D_PRS_Cat, 
                         levels = c("Low", "Medium", "High", "Missing"))
vars <- c(
  "Age", 
  "Male_sex", 
  "Female_sex",
  "Age_at_dementia_diagnosis",
  "Follow_up_years",
  "Dementia_rate_1000py", 
  "Alzheimers_rate_1000py",
  "Vascular_dementia_rate_1000py",
  "Other_dementia_rate_1000py",
  "Total_T2D_events",
  "Education_level", 
  "Smoking_status", 
  "Alcohol_consumption",
  "Marital_status",
  "Residence",
  "BMI",
  "FBG", 
  "LDL_C", 
  "HDL_C", 
  "TC", 
  "TG", 
  "SBP", 
  "DBP",
  "Coronary_heart_disease",
  "Stroke",
  "Heart_failure",
  "Sleep_disorders",
  "Depression"
)
factorVars <- c(
  "Male_sex", 
  "Female_sex", 
  "Education_level", 
  "Smoking_status", 
  "Alcohol_consumption",
  "Marital_status",
  "Residence",
  "Coronary_heart_disease",
  "Stroke",
  "Heart_failure",
  "Sleep_disorders",
  "Depression"
)
table1 <- CreateTableOne(
  vars = vars,
  strata = "T2D_PRS_Cat",
  data = df,
  factorVars = factorVars
)
print(table1, showAllLevels = TRUE)

#12.eTable 6 Association between T2D genetic clusters  and Risks of All-cause Dementia, Alzheimer's Disease, Vascular Dementia, and Other/Unspecified Dementia
library(survival)
library(broom)
library(dplyr)
df <- read.csv("data.csv")
run_cox_model <- function(data, time_col, event_col, predictor, covariates) {
   surv_obj <- Surv(time = data[[time_col]], event = data[[event_col]])
   formula_str <- as.formula(paste("surv_obj ~", predictor, "+", paste(covariates, collapse = " + ")))
  cox_model <- coxph(formula_str, data = data)
  tidy(cox_model, exponentiate = TRUE) %>% 
    filter(term == predictor) %>% 
    select(term, estimate, conf.low, conf.high, p.value)
}
time_col <- "follow_up_years"  
dementia_outcomes <- c("all_cause_dementia", "Alzheimers", "Vascular_dementia", "Other_dementia") 
genetic_clusters <- c(
  "ALP_Negative_pPS",
  "Beta_Cell_1_pPS",
  "Beta_Cell_2_pPS",
  "Bilirubin_pPS",
  "Cholesterol_pPS",
  "Hyper_Insulin_pPS",
  "Lipodystrophy_1_pPS",
  "Lipodystrophy_2_pPS",
  "Liver_Lipid_pPS",
  "Obesity_pPS",
  "Proinsulin_pPS",
  "SHBG_LpA_pPS"
)
covariates <- c("SEX", "AGE", "T2D_status", "APOE_e4")
results <- list()
for (cluster in genetic_clusters) {
  for (outcome in dementia_outcomes) {
    model_result <- run_cox_model(
      data = df,
      time_col = time_col,
      event_col = outcome,
      predictor = cluster,
      covariates = covariates
    )
        results <- append(results, list(data.frame(
      Genetic_Cluster = cluster,
      Dementia_Type = outcome,
      HR = model_result$estimate,
      CI = paste0("(", round(model_result$conf.low, 2), "-", round(model_result$conf.high, 2), ")"),
      P_value = model_result$p.value
    )))
  }
}
final_results <- bind_rows(results)
print(final_results)

#13.eTable 7 Associations Between clusters Tertiles for Hyperinsulinemia and SHBGLpA and Risk of All-Cause Dementia, Alzheimer’s Disease, Vascular Dementia, and Other/Unspecified Dementia.
library(survival)
library(dplyr)
library(broom)
df <- read.csv("data.csv")
df$HyperInsulin_tertile <- factor(df$HyperInsulin_tertile, levels = c("Low", "Medium", "High"))
df$SHBGLpA_tertile <- factor(df$SHBGLpA_tertile, levels = c("Low", "Medium", "High"))
run_cox_model <- function(data, time_col, event_col, predictor, covariates) {
   surv_obj <- Surv(time = data[[time_col]], event = data[[event_col]])
    formula_str <- as.formula(paste("surv_obj ~", predictor, "+", paste(covariates, collapse = " + ")))
  cox_model <- coxph(formula_str, data = data)
  tidy(cox_model, exponentiate = TRUE) %>% 
    filter(term != "(Intercept)") %>% 
    select(term, estimate, conf.low, conf.high, p.value)
}
time_col <- "follow_up_years"  dementia_outcomes <- c("all_cause_dementia", "Alzheimers", "Vascular_dementia", "Other_dementia") 
genetic_clusters <- c("HyperInsulin_tertile", "SHBGLpA_tertile")
covariates <- c("SEX", "AGE")
results <- list()
for (cluster in genetic_clusters) {
  for (outcome in dementia_outcomes) {
    model_result <- run_cox_model(
      data = df,
      time_col = time_col,
      event_col = outcome,
      predictor = cluster,
      covariates = covariates
    )
    results <- append(results, list(data.frame(
      Genetic_Cluster = cluster,
      Dementia_Type = outcome,
      HR = model_result$estimate,
      CI = paste0("(", round(model_result$conf.low, 2), "-", round(model_result$conf.high, 2), ")"),
      P_value = model_result$p.value
    )))
  }
}
final_results <- bind_rows(results)
print(final_results)

#14.eTable 8 Associations Between clusters Tertiles for Hyperinsulinemia and SHBGLpA and Risk of All-Cause Dementia, Alzheimer’s Disease, Vascular Dementia, and Other/Unspecified Dementia.
library(survival)
library(dplyr)
library(broom)
df <- read.csv("data.csv")
df$HyperInsulin_tertile <- factor(df$HyperInsulin_tertile, levels = c("Low", "Medium", "High"))
df$SHBGLpA_tertile <- factor(df$SHBGLpA_tertile, levels = c("Low", "Medium", "High"))
run_cox_model <- function(data, time_col, event_col, predictor, covariates) {
  surv_obj <- Surv(time = data[[time_col]], event = data[[event_col]])
    formula_str <- as.formula(paste("surv_obj ~", predictor, "+", paste(covariates, collapse = " + ")))
  cox_model <- coxph(formula_str, data = data)
  tidy(cox_model, exponentiate = TRUE) %>% 
    filter(term != "(Intercept)") %>% 
    select(term, estimate, conf.low, conf.high, p.value)
}
time_col <- "follow_up_years"  
dementia_outcomes <- c("all_cause_dementia", "Alzheimers", "Vascular_dementia", "Other_dementia") 
genetic_clusters <- c("HyperInsulin_tertile", "SHBGLpA_tertile")
covariates <- c("SEX", "AGE", "T2D_status")
results <- list()
for (cluster in genetic_clusters) {
  for (outcome in dementia_outcomes) {
    model_result <- run_cox_model(
      data = df,
      time_col = time_col,
      event_col = outcome,
      predictor = cluster,
      covariates = covariates
    )
    results <- append(results, list(data.frame(
      Genetic_Cluster = cluster,
      Dementia_Type = outcome,
      HR = model_result$estimate,
      CI = paste0("(", round(model_result$conf.low, 2), "-", round(model_result$conf.high, 2), ")"),
      P_value = model_result$p.value
    )))
  }
}
final_results <- bind_rows(results)
print(final_results)

#15.eTable 9: Association of Hyper Insulin Cluster-Specific Genetic Variants with Vascular Dementia Risk
auto_PRS_calc_v3 \
    -g <path_to_genotype_data>/merged_imputed_info0.3_maf0.01_2638ref \
    -s <path_to_summary_statistics>/HyperInsulin_VD_GWAS.txt \
    -o HyperInsulin_VD \
    -r 0.2 \
    -p 0.01 \
    --colname_chr=chr \
    --colname_pos=pos \
    --colname_ea=EffectAllele \
    --colname_nea=NonEffectAllele \
    --colname_beta=Beta \
--colname_p=Pval

#16.eTable 10: Association of SHBG-LpA Cluster-Specific Genetic Variants with Vascular Dementia Risk
auto_PRS_calc_v3 \
    -g <path_to_genotype_data>/merged_imputed_info0.3_maf0.01_2638ref \
    -s <path_to_summary_statistics>/SHBGLpA_VD_GWAS.txt \
    -o SHBGLpA_VD \
    -r 0.2 \
    -p 0.01 \
    --colname_chr=chr \
    --colname_pos=pos \
    --colname_ea=EffectAllele \
    --colname_nea=NonEffectAllele \
    --colname_beta=Beta \
    --colname_p=Pval
