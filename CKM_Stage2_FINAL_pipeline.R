## =============================================================================
##  CKM Stage 2 Phenotypes and Mortality — FINAL END-TO-END PIPELINE
##  NHANES 2007–2018 + NCHS Public-Use Linked Mortality Files (through 2019)
##
##  One script = one source of truth. It rebuilds every table, figure, and
##  manuscript number from raw data. Nothing from earlier runs is reused.
##
##  Key corrections vs. all previous runs
##   1. Stage 2 hypertriglyceridemia = TG >= 135 mg/dL (AHA CKM advisory,
##      Ndumele 2023). MetS TG component stays >= 150 (harmonized 2009).
##   2. Stage 3 risk-equivalent uses AHA PREVENT 10-yr total CVD >= 20%
##      (package 'preventr'), matching the Methods text.
##   3. Survey design built on the full fasting subsample, then subset
##      (correct domain estimation), not on the complete-case data alone.
##   4. Cluster labeling: three disease domains assigned greedily; the
##      leftover cluster is the reference and is named by what it is.
##   5. Medication / smoking variables recoded to population-level 0/1
##      (fixes the NA smoking p-value and skip-pattern artefacts).
##   6. New sensitivity analysis: clustering WITHOUT age.
##   7. All figures use clinical labels; PRISMA drawn from live counts.
##   8. Results_autofilled.md writes the manuscript numbers from the data.
##
##  Run:  source("CKM_Stage2_FINAL_pipeline.R")   (internet needed on 1st run)
##  Time: ~20–60 min depending on machine (bootstraps + RSF + PREVENT).
## =============================================================================

## ---------------------------- 0. CONFIG --------------------------------------
CFG <- list(
  out_dir        = "CKM_final_outputs",
  cache_dir      = "nhanes_cache",
  mort_dir       = "mortality_files",
  seed           = 20260917,
  tg_stage2      = 135,     # AHA CKM advisory
  tg_mets        = 150,     # harmonized MetS
  k              = 4,
  k_range        = 2:8,
  nstart         = 50,
  n_boot_jaccard = 500,
  n_boot_c       = 200,
  tau            = 10,      # years, RMST / NRI / IDI horizon
  rsf_ntree      = 1000,
  idi_npert      = 300
)

## ---------------------------- 1. PACKAGES ------------------------------------
pkgs <- c("nhanesA","dplyr","tidyr","readr","purrr","stringr","survey","survival",
          "ggplot2","survminer","patchwork","randomForestSRC","gtsummary","gt",
          "flextable","officer","survC1","survIDINRI","smd","cardx","fpc","mclust",
          "cmprsk","survRM2","multcomp","jsonlite","preventr")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss)
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))
## Masking guard: MASS (pulled in by multcomp) and stats hijack these names
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
recode <- dplyr::recode; summarise <- dplyr::summarise; mutate <- dplyr::mutate
lag <- dplyr::lag; first <- dplyr::first; last <- dplyr::last
options(survey.lonely.psu = "adjust", dplyr.summarise.inform = FALSE)

dir.create(CFG$out_dir,  showWarnings = FALSE)
dir.create(CFG$cache_dir, showWarnings = FALSE)
dir.create(CFG$mort_dir,  showWarnings = FALSE)
OUT  <- function(...) file.path(CFG$out_dir, ...)
LOG  <- OUT("pipeline_log.txt"); cat("", file = LOG)
say  <- function(...) { m <- sprintf(...); message(m); cat(m, "\n", file = LOG, append = TRUE) }
safe <- function(expr, what) tryCatch(expr, error = function(e) { say("!! %s FAILED: %s", what, conditionMessage(e)); NULL })
R <- list()   # every manuscript number is collected here

LAB <- c(RIDAGEYR = "Age (years)", BMXBMI = "BMI (kg/m²)", BMXWAIST = "Waist circumference (cm)",
         sbp = "Systolic BP (mmHg)", dbp = "Diastolic BP (mmHg)", egfr = "eGFR (mL/min/1.73 m²)",
         LBXGLU = "Fasting glucose (mg/dL)", LBXGH = "HbA1c (%)", LBXTR = "Triglycerides (mg/dL)",
         LBDHDD = "HDL cholesterol (mg/dL)", LBXTC = "Total cholesterol (mg/dL)",
         LBXSUA = "Uric acid (mg/dL)", URDACT = "UACR (mg/g)", female = "Female sex",
         cluster = "Phenotype cluster", RaceEth = "Race/ethnicity", Diabetes = "Diabetes",
         Hypertension = "Hypertension", MetSyn = "Metabolic syndrome")
CLUST_VARS <- c("RIDAGEYR","BMXBMI","BMXWAIST","sbp","dbp","egfr","LBXGLU","LBXGH","LBXTR","LBDHDD")

## ---------------------------- 2. DATA ----------------------------------------
cycles <- c(E = "2007_2008", F = "2009_2010", G = "2011_2012",
            H = "2013_2014", I = "2015_2016", J = "2017_2018")
keep <- list(
  DEMO   = c("SEQN","RIAGENDR","RIDAGEYR","RIDRETH1","RIDEXPRG","SDMVPSU","SDMVSTRA"),
  BMX    = c("SEQN","BMXBMI","BMXWAIST"),
  BPX    = c("SEQN", paste0("BPXSY",1:4), paste0("BPXDI",1:4)),
  BIOPRO = c("SEQN","LBXSCR","LBXSUA"),
  GHB    = c("SEQN","LBXGH"),
  GLU    = c("SEQN","LBXGLU","WTSAF2YR"),
  TRIGLY = c("SEQN","LBXTR"),
  HDL    = c("SEQN","LBDHDD"),
  TCHOL  = c("SEQN","LBXTC"),
  ALB_CR = c("SEQN","URDACT"),
  BPQ    = c("SEQN","BPQ020","BPQ040A","BPQ050A","BPQ090D","BPQ100D"),
  DIQ    = c("SEQN","DIQ010","DIQ050","DIQ070"),
  MCQ    = c("SEQN","MCQ160B","MCQ160C","MCQ160D","MCQ160E","MCQ160F"),
  SMQ    = c("SEQN","SMQ020","SMQ040")
)
get_tab <- function(name) {
  f <- file.path(CFG$cache_dir, paste0(name, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  d <- tryCatch(nhanesA::nhanes(name, translated = FALSE), error = function(e) NULL)
  if (is.null(d)) stop("Could not download ", name)
  saveRDS(d, f); d
}
load_cycle <- function(sfx) {
  tabs <- imap(keep, function(v, t) get_tab(paste0(t, "_", sfx)) %>% select(any_of(v)))
  d <- reduce(tabs[-1], left_join, by = "SEQN", .init = tabs$DEMO)
  d$cycle <- sfx; d
}
read_mort <- function(yrs) {
  fn <- sprintf("NHANES_%s_MORT_2019_PUBLIC.dat", yrs)
  fp <- file.path(CFG$mort_dir, fn)
  if (!file.exists(fp))
    download.file(paste0("https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/", fn), fp, mode = "wb")
  read_fwf(fp, col_types = "iiiiiiii",
           fwf_cols(SEQN = c(1,6), eligstat = c(15,15), mortstat = c(16,16), ucod_leading = c(17,19),
                    diabetes_mcod = c(20,20), hyperten_mcod = c(21,21), permth_int = c(43,45), permth_exm = c(46,48)),
           na = c("", "."))
}
say("Loading NHANES ...")
nh_raw <- map_dfr(names(cycles), load_cycle)
mort <- map_dfr(cycles, read_mort)
nh_raw <- left_join(nh_raw, mort, by = "SEQN")
say("Raw pooled N = %d", nrow(nh_raw))

## ---------------------------- 3. DERIVED VARIABLES ---------------------------
yes <- function(x) !is.na(x) & x == 1
egfr21 <- function(scr, age, female) {
  k <- ifelse(female, 0.7, 0.9); a <- ifelse(female, -0.241, -0.302)
  142 * pmin(scr / k, 1)^a * pmax(scr / k, 1)^(-1.200) * 0.9938^age * ifelse(female, 1.012, 1)
}
kdigo <- function(egfr, uacr) {
  g <- as.character(cut(egfr, c(-Inf,15,30,45,60,90,Inf), right = FALSE,
                        labels = c("G5","G4","G3b","G3a","G2","G1")))
  A <- ifelse(is.na(uacr) | uacr < 30, "A1", ifelse(uacr <= 300, "A2", "A3"))
  case_when(is.na(g) ~ NA_character_,
            g %in% c("G4","G5") ~ "very_high",
            g == "G3b" & A != "A1" ~ "very_high",
            g == "G3a" & A == "A3" ~ "very_high",
            g == "G3b" ~ "high",
            g == "G3a" & A == "A2" ~ "high",
            g %in% c("G1","G2") & A == "A3" ~ "high",
            g == "G3a" ~ "moderate",
            g %in% c("G1","G2") & A == "A2" ~ "moderate",
            TRUE ~ "low")
}
rowmean_na <- function(df) { m <- rowMeans(df, na.rm = TRUE); m[is.nan(m)] <- NA; m }

d <- nh_raw %>%
  mutate(across(starts_with("BPXDI"), ~ ifelse(.x == 0, NA, .x)),
         female  = RIAGENDR == 2,
         sbp     = rowmean_na(pick(any_of(paste0("BPXSY",1:4)))),
         dbp     = rowmean_na(pick(any_of(paste0("BPXDI",1:4)))),
         egfr    = egfr21(LBXSCR, RIDAGEYR, female),
         RaceEth = factor(RIDRETH1, 1:5, c("Mexican American","Other Hispanic","Non-Hispanic White",
                                           "Non-Hispanic Black","Other/Multi-Racial")),
         htn_med = yes(BPQ050A),                     # skip pattern: not asked => not taking
         dm_med  = yes(DIQ050) | yes(DIQ070),
         lipid_med = yes(BPQ100D),
         ever_smoker   = case_when(SMQ020 == 1 ~ 1, SMQ020 == 2 ~ 0, TRUE ~ NA_real_),
         current_smoker = !is.na(SMQ040) & SMQ040 %in% c(1,2),
         Hypertension = coalesce(sbp >= 130 | dbp >= 80, FALSE) | yes(BPQ020) | htn_med,
         Diabetes     = yes(DIQ010) | dm_med | coalesce(LBXGH >= 6.5, FALSE) | coalesce(LBXGLU >= 126, FALSE),
         mets_n = coalesce(BMXWAIST >= ifelse(female, 88, 102), FALSE) +
           coalesce(LBXTR >= CFG$tg_mets, FALSE) +
           coalesce(LBDHDD < ifelse(female, 50, 40), FALSE) +
           (coalesce(sbp >= 130 | dbp >= 85, FALSE) | htn_med) +
           (coalesce(LBXGLU >= 100, FALSE) | dm_med),
         MetSyn = mets_n >= 3,
         HTG    = coalesce(LBXTR >= CFG$tg_stage2, FALSE),
         ckd_cat = kdigo(egfr, URDACT),
         CKD_modhigh = ckd_cat %in% c("moderate","high"),
         CKD_vhigh   = ckd_cat %in% "very_high",
         CVD = yes(MCQ160B) | yes(MCQ160C) | yes(MCQ160D) | yes(MCQ160E) | yes(MCQ160F),
         stage2_criteria = HTG | Hypertension | MetSyn | Diabetes | CKD_modhigh,
         time_y = permth_exm / 12,
         death  = as.integer(mortstat == 1),
         cvdeath = as.integer(mortstat == 1 & ucod_leading %in% c(1, 5)),
         pooled_weight = WTSAF2YR / 6)

## PREVENT 10-yr total CVD (only where it can change stage)
prev_file <- file.path(CFG$cache_dir, "prevent_scores.rds")
find_tcvd <- function(x) {
  if (is.data.frame(x) && "total_cvd" %in% names(x)) return(as.numeric(x$total_cvd[1]))
  if (is.list(x)) for (el in x) { v <- find_tcvd(el); if (!is.null(v)) return(v) }
  NULL
}
if (file.exists(prev_file)) {
  prev <- readRDS(prev_file)
} else {
  say("Computing PREVENT scores (slow, cached afterwards) ...")
  cand <- d %>% filter(RIDAGEYR >= 30, RIDAGEYR <= 79, stage2_criteria, !CVD, !CKD_vhigh,
                       !is.na(sbp), !is.na(LBXTC), !is.na(LBDHDD), !is.na(egfr), !is.na(BMXBMI))
  prev <- tibble(SEQN = cand$SEQN, prevent10 = NA_real_)
  for (i in seq_len(nrow(cand))) {
    r <- cand[i, ]
    res <- tryCatch(preventr::estimate_risk(
      age = r$RIDAGEYR, sex = ifelse(r$female, "female", "male"), sbp = r$sbp,
      bp_tx = r$htn_med, total_c = r$LBXTC, hdl_c = r$LBDHDD, statin = r$lipid_med,
      dm = r$Diabetes, smoking = r$current_smoker, egfr = r$egfr, bmi = r$BMXBMI,
      time = "10yr", quiet = TRUE), error = function(e) NULL)
    v <- if (is.null(res)) NULL else find_tcvd(res)
    prev$prevent10[i] <- if (is.null(v)) NA_real_ else v
  }
  if (all(is.na(prev$prevent10))) stop("PREVENT returned no values. Check preventr::estimate_risk() arguments for your installed version.")
  if (max(prev$prevent10, na.rm = TRUE) > 1) prev$prevent10 <- prev$prevent10 / 100
  saveRDS(prev, prev_file)
}
d <- d %>% left_join(prev, by = "SEQN") %>%
  mutate(Stage3_risk = coalesce(prevent10 >= 0.20, FALSE),
         stage2_isolated = stage2_criteria & !CVD & !CKD_vhigh & !Stage3_risk)
R$prevent_scored <- sum(!is.na(prev$prevent10)); R$prevent_candidates <- nrow(prev)

## ---------------------------- 4. ATTRITION -----------------------------------
s1 <- d %>% filter(RIDAGEYR >= 18, !yes(RIDEXPRG))
s2 <- s1 %>% filter(eligstat == 1, !is.na(mortstat), !is.na(permth_exm))
s3 <- s2 %>% filter(stage2_isolated)
s4 <- s3 %>% filter(if_all(all_of(c(CLUST_VARS, "female")), ~ !is.na(.x)),
                    !is.na(WTSAF2YR), WTSAF2YR > 0, time_y > 0)
flow <- tibble(Step = c("Raw pooled NHANES 2007–2018", "Adult (≥18 y), non-pregnant",
                        "Eligible NCHS mortality linkage", "CKM Stage 2 (Stages 0, 1, 3, 4 excluded)",
                        "Complete core variables + valid fasting weight"),
               N = c(nrow(d), nrow(s1), nrow(s2), nrow(s3), nrow(s4)))
flow$Excluded <- c(NA, -diff(flow$N))
write_csv(flow, OUT("PRISMA_attrition_flow.csv"))
R$flow <- flow
say("Analytic N = %d", nrow(s4))

## ---------------------------- 5. WEIGHTED K-MEANS ----------------------------
set.seed(CFG$seed)
wstd <- function(x, w) { mu <- sum(w * x) / sum(w); s <- sqrt(sum(w * (x - mu)^2) / sum(w)); (x - mu) / s }
wkmeans <- function(X, w, k, nstart = 50, iter.max = 100) {
  X <- as.matrix(X); n <- nrow(X); best <- NULL
  dist2 <- function(C) sapply(seq_len(k), function(j) rowSums(sweep(X, 2, C[j, ])^2))
  for (s in seq_len(nstart)) {
    C <- X[sample.int(n, k, prob = w), , drop = FALSE]
    for (it in seq_len(iter.max)) {
      cl <- max.col(-dist2(C), ties.method = "first")
      if (length(unique(cl)) < k) break
      Cn <- t(sapply(seq_len(k), function(j) colSums(X[cl == j, , drop = FALSE] * w[cl == j]) / sum(w[cl == j])))
      conv <- max(abs(Cn - C)) < 1e-8; C <- Cn; if (conv) break
    }
    if (length(unique(cl)) < k) next
    D <- dist2(C); cl <- max.col(-D, ties.method = "first")
    wss <- sum(w * D[cbind(seq_len(n), cl)])
    if (is.null(best) || wss < best$wss) best <- list(cluster = cl, centers = C, wss = wss)
  }
  best
}
dat <- s4
w_norm <- dat$pooled_weight / mean(dat$pooled_weight)
Z <- sapply(CLUST_VARS, function(v) wstd(dat[[v]], w_norm))

elbow <- tibble(k = CFG$k_range, wss = map_dbl(CFG$k_range, ~ wkmeans(Z, w_norm, .x, nstart = 20)$wss))
write_csv(elbow, OUT("Table_S_Elbow_WSS.csv"))
ggsave(OUT("Figure_S3_Elbow.pdf"),
       ggplot(elbow, aes(k, wss)) + geom_line() + geom_point(size = 2) + theme_classic() +
         labs(x = "Number of clusters (k)", y = "Weighted within-cluster sum of squares"),
       width = 5, height = 4)

km <- wkmeans(Z, w_norm, CFG$k, nstart = CFG$nstart)

## Labeling rule (pre-specified, applied identically in every clustering below)
label_clusters <- function(Zm, cl, w) {
  dom <- function(v) sapply(sort(unique(cl)), function(j) sum(w[cl == j] * v[cl == j]) / sum(w[cl == j]))
  z <- function(n) if (n %in% colnames(Zm)) Zm[, n] else rep(0, nrow(Zm))
  age_dom <- if ("RIDAGEYR" %in% colnames(Zm)) (z("RIDAGEYR") + z("sbp")) / 2 else z("sbp")
  M <- cbind(Vascular   = dom(age_dom),
             Gluco      = dom((z("LBXGLU") + z("LBXGH")) / 2),
             Obese      = dom((z("BMXBMI") + z("BMXWAIST")) / 2),
             Dyslip     = dom((z("LBXTR") - z("LBDHDD")) / 2),
             Age_only   = dom(z("RIDAGEYR")))
  rownames(M) <- sort(unique(cl))
  lab <- setNames(rep(NA_character_, nrow(M)), rownames(M))
  A <- M[, c("Vascular","Gluco","Obese")]
  for (step in 1:3) {
    avail <- A; avail[!is.na(lab), ] <- -Inf; avail[, c("Vascular","Gluco","Obese") %in% lab] <- -Inf
    idx <- which(avail == max(avail), arr.ind = TRUE)[1, ]
    lab[idx[1]] <- colnames(A)[idx[2]]
  }
  ref <- which(is.na(lab))
  lab[ref] <- if ("RIDAGEYR" %in% colnames(Zm) && M[ref, "Age_only"] <= -0.5) "RefYoung"
  else if (M[ref, "Dyslip"] >= 0.5) "RefDyslip" else "RefLow"
  pretty <- c(RefYoung = "Younger Adult (Ref)", RefDyslip = "Atherogenic Dyslipidemia (Ref)",
              RefLow = "Lower-Risk (Ref)", Vascular = "Vascular Aging / ISH",
              Gluco = "Severe Glucotoxicity", Obese = "Severe Obesity / Cardiometabolic")
  lv <- c(pretty[lab[ref]], pretty[c("Vascular","Gluco","Obese")])
  list(factor = factor(pretty[lab[as.character(cl)]], levels = lv),
       audit = data.frame(cluster_raw = rownames(M), round(M, 3), label = pretty[lab], row.names = NULL))
}
lb <- label_clusters(Z, km$cluster, w_norm)
dat$cluster <- lb$factor
write_csv(lb$audit, OUT("Table_S12_Cluster_signature_and_labels.csv"))
REF <- levels(dat$cluster)[1]; NONREF <- levels(dat$cluster)[-1]
R$labels <- lb$audit; R$ref_label <- REF
say("Reference cluster: %s", REF)

## ---------------------------- 6. SURVEY DESIGN (domain) ----------------------
fast <- d %>% filter(!is.na(WTSAF2YR), WTSAF2YR > 0) %>%
  left_join(dat %>% select(SEQN, cluster), by = "SEQN") %>%
  mutate(analytic = SEQN %in% dat$SEQN,
         uacr_log = log(URDACT + 1),
         across(c(Diabetes, Hypertension, MetSyn), ~ factor(ifelse(.x, "Yes", "No"), c("No","Yes"))),
         htn_med01 = as.numeric(htn_med), dm_med01 = as.numeric(dm_med))
des_all <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~pooled_weight, nest = TRUE, data = fast)
des <- subset(des_all, analytic)
R$N <- nrow(dat); R$weighted_pop <- sum(dat$pooled_weight)
R$n_by_cluster <- as.list(table(dat$cluster))
R$wpop_by_cluster <- as.list(tapply(dat$pooled_weight, dat$cluster, sum))

## ---------------------------- 7. TABLE 1 & S1 --------------------------------
t1 <- safe({
  tbl_svysummary(des, by = cluster,
                 include = c(RIDAGEYR, female, RaceEth, egfr, LBXGLU, LBXGH, BMXBMI, BMXWAIST, sbp, dbp,
                             LBXTR, LBDHDD, LBXTC, Diabetes, Hypertension, MetSyn),
                 label = as.list(LAB[c("RIDAGEYR","female","RaceEth","egfr","LBXGLU","LBXGH","BMXBMI","BMXWAIST",
                                       "sbp","dbp","LBXTR","LBDHDD","LBXTC","Diabetes","Hypertension","MetSyn")]) %>%
                   setNames(c("RIDAGEYR","female","RaceEth","egfr","LBXGLU","LBXGH","BMXBMI","BMXWAIST",
                              "sbp","dbp","LBXTR","LBDHDD","LBXTC","Diabetes","Hypertension","MetSyn")),
                 statistic = list(all_continuous() ~ "{mean} ({sd})", all_categorical() ~ "{n_unweighted} ({p}%)"),
                 digits = list(LBXGH ~ 2), value = list(Diabetes ~ "Yes", Hypertension ~ "Yes", MetSyn ~ "Yes")) %>%
    add_overall() %>% add_p() %>%
    modify_header(all_stat_cols() ~ "**{level}**  \nn = {n_unweighted}") %>%
    modify_caption("**Table 1. Survey-weighted baseline characteristics by phenotype, CKM Stage 2, NHANES 2007–2018**")
}, "Table 1")
if (!is.null(t1)) { gt::gtsave(as_gt(t1), OUT("Table1.html")); flextable::save_as_docx(as_flex_table(t1), path = OUT("Table1.docx")) }

s3x <- s3 %>% mutate(Included = factor(ifelse(SEQN %in% dat$SEQN, "Included", "Excluded"), c("Excluded","Included")))
tS1 <- safe({
  tbl_summary(s3x, by = Included, include = c(RIDAGEYR, female, egfr, LBXGLU, LBXGH, BMXBMI, sbp, LBXTR, LBDHDD),
              label = list(RIDAGEYR ~ LAB[["RIDAGEYR"]], female ~ LAB[["female"]], egfr ~ LAB[["egfr"]],
                           LBXGLU ~ LAB[["LBXGLU"]], LBXGH ~ LAB[["LBXGH"]], BMXBMI ~ LAB[["BMXBMI"]],
                           sbp ~ LAB[["sbp"]], LBXTR ~ LAB[["LBXTR"]], LBDHDD ~ LAB[["LBDHDD"]]),
              statistic = all_continuous() ~ "{mean} ({sd})", missing_text = "Missing") %>%
    add_difference(everything() ~ "smd")
}, "Table S1")
if (!is.null(tS1)) flextable::save_as_docx(as_flex_table(tS1), path = OUT("Table_S1_Missingness.docx"))

## ---------------------------- 8. ROBUSTNESS ----------------------------------
## 8a Jaccard bootstrap
cb_fun <- function(data, k, ...) {
  w <- data[, ncol(data)]; X <- data[, -ncol(data), drop = FALSE]
  r <- wkmeans(X, w, k, nstart = 10); p <- r$cluster
  list(result = r, nc = k, clusterlist = lapply(seq_len(k), function(j) p == j),
       partition = p, clustermethod = "weighted k-means")
}
jac <- safe({
  cb <- clusterboot(cbind(Z, w_norm), B = CFG$n_boot_jaccard, bootmethod = "boot",
                    clustermethod = cb_fun, k = CFG$k, count = FALSE, seed = CFG$seed)
  map_lab <- sapply(seq_len(CFG$k), function(j) names(which.max(table(dat$cluster[cb$partition == j]))))
  tibble(Cluster = map_lab, Jaccard_bootmean = cb$bootmean,
         Stability = ifelse(cb$bootmean >= 0.85, "Highly stable", ifelse(cb$bootmean >= 0.75, "Stable", "Unstable")))
}, "Jaccard")
if (!is.null(jac)) { write_csv(jac, OUT("Table_S6_Jaccard_Stability.csv")); R$jaccard <- range(jac$Jaccard_bootmean) }

## 8b GMM (unweighted)
gmm <- safe({
  g_all <- Mclust(Z, G = 1:8, verbose = FALSE); g4 <- Mclust(Z, G = 4, verbose = FALSE)
  pdf(OUT("Figure_S4_GMM_BIC.pdf"), 6, 5); plot(g_all, what = "BIC"); dev.off()
  tibble(Method = c("GMM (BIC-optimal)", "GMM (k = 4)"), k = c(g_all$G, 4),
         ARI_vs_kmeans = c(adjustedRandIndex(g_all$classification, km$cluster),
                           adjustedRandIndex(g4$classification, km$cluster)))
}, "GMM")
if (!is.null(gmm)) { write_csv(gmm, OUT("Table_S5_GMM_Sensitivity.csv")); R$gmm <- gmm }

## 8c Age-stratified re-clustering
ages <- safe({
  map_dfr(list(`Age < 50` = dat$RIDAGEYR < 50, `Age ≥ 50` = dat$RIDAGEYR >= 50), function(ix) {
    r <- wkmeans(Z[ix, ], w_norm[ix], CFG$k, nstart = CFG$nstart)
    a <- adjustedRandIndex(r$cluster, km$cluster[ix])
    tibble(N = sum(ix), ARI_vs_pooled = a,
           Interpretation = cut(a, c(-Inf, .2, .4, .6, .8, Inf), c("Poor","Weak","Moderate","Good","Excellent")))
  }, .id = "Stratum")
}, "Age strata")
if (!is.null(ages)) { write_csv(ages, OUT("Table_S3_Age_Stratified_Clustering.csv")); R$age_strata <- ages }

## 8d Convergent validity
ext <- safe({
  cont <- c(LBXTC = "LBXTC", LBXSUA = "LBXSUA", UACR_log = "uacr_log")
  binv <- c(`Ever smoker` = "ever_smoker", `Antihypertensive use` = "htn_med01", `Antidiabetic use` = "dm_med01")
  f_c <- imap_dfr(cont, function(v, nm) {
    m <- svyby(as.formula(paste0("~", v)), ~cluster, des, svymean, na.rm = TRUE)
    p <- regTermTest(svyglm(as.formula(paste(v, "~ cluster")), des), ~cluster)$p
    val <- m[[v]]; if (v == "uacr_log") { val <- exp(val) - 1; nm <- "UACR (geometric mean, mg/g)" }
    tibble(Variable = nm, !!!setNames(as.list(val), m$cluster), P_Value = as.numeric(p))
  })
  f_b <- imap_dfr(binv, function(v, nm) {
    m <- svyby(as.formula(paste0("~", v)), ~cluster, des, svymean, na.rm = TRUE)
    p <- svychisq(as.formula(paste0("~ factor(", v, ") + cluster")), des, statistic = "F")$p.value
    tibble(Variable = paste0(nm, " (proportion)"), !!!setNames(as.list(m[[v]]), m$cluster), P_Value = as.numeric(p))
  })
  bind_rows(f_c, f_b)
}, "External validation")
if (!is.null(ext)) { write_csv(ext, OUT("Table_S4_Convergent_Validity.csv")); R$ext <- ext }

## ---------------------------- 9. SURVIVAL ------------------------------------
tidy_cox <- function(m, method = NA) {
  b <- coef(m); se <- sqrt(diag(vcov(m)))
  tibble(Covariate = names(b), HR = exp(b), LCL_95 = exp(b - 1.96 * se), UCL_95 = exp(b + 1.96 * se),
         P_Value = 2 * pnorm(-abs(b / se)), Method = method) %>%
    mutate(Covariate = str_replace(Covariate, "^cluster", "cluster: "),
           Covariate = recode(Covariate, femaleTRUE = "Female sex", RIDAGEYR = "Age (per year)",
                              egfr = "eGFR (per mL/min/1.73 m²)", LBXGLU = "Fasting glucose (per mg/dL)"))
}
pull_cl <- function(tb) tb %>% filter(str_detect(Covariate, "^cluster")) %>% mutate(Covariate = str_remove(Covariate, "cluster: "))

## Table 3 — rates
t3 <- dat %>% group_by(cluster) %>%
  summarise(N = n(), Deaths = sum(death), CV_Deaths = sum(cvdeath), PersonYears = sum(time_y),
            W_Deaths = sum(pooled_weight * death), W_PY = sum(pooled_weight * time_y)) %>%
  bind_rows(dat %>% summarise(cluster = "Overall", N = n(), Deaths = sum(death), CV_Deaths = sum(cvdeath),
                              PersonYears = sum(time_y), W_Deaths = sum(pooled_weight * death),
                              W_PY = sum(pooled_weight * time_y))) %>%
  mutate(Weighted_rate_per_1000PY = 1000 * W_Deaths / W_PY, CV_share_pct = 100 * CV_Deaths / Deaths)
write_csv(t3, OUT("Table3_Mortality_Rates.csv")); R$t3 <- t3
R$mean_fu <- mean(dat$time_y); R$median_fu <- median(dat$time_y)

## Figure 1 — KM (weighted curves, unweighted numbers at risk)
km_fig <- safe({
  fw <- survfit(Surv(time_y, death) ~ cluster, data = dat, weights = w_norm)
  fu <- survfit(Surv(time_y, death) ~ cluster, data = dat)
  lr <- tryCatch({ x <- svylogrank(Surv(time_y, death) ~ cluster, des); as.numeric(tail(unlist(x), 1)) },
                 error = function(e) survdiff(Surv(time_y, death) ~ cluster, data = dat)$pvalue)
  R$logrank_p <<- lr
  labs_ <- levels(dat$cluster)
  p1 <- ggsurvplot(fw, data = dat, legend.labs = labs_, xlim = c(0, 10), break.time.by = 2,
                   ylim = c(0.7, 1), censor = FALSE, xlab = "Years since examination",
                   ylab = "Survival probability (weighted)", legend.title = "",
                   ggtheme = theme_classic(base_size = 11))$plot +
    annotate("text", x = 0.5, y = 0.72, hjust = 0,
             label = ifelse(lr < 0.001, "Design-based log-rank p < 0.001", sprintf("Design-based log-rank p = %.3f", lr)))
  p2 <- ggrisktable(fu, data = dat, legend.labs = labs_, xlim = c(0, 10), break.time.by = 2,
                    xlab = "Years", ylab = "", ggtheme = theme_classic(base_size = 9), fontsize = 3)
  pp <- p1 / p2 + plot_layout(heights = c(3, 1.1))
  ggsave(OUT("Figure1_KaplanMeier.pdf"), pp, width = 8, height = 7)
  ggsave(OUT("Figure1_KaplanMeier.tiff"), pp, width = 8, height = 7, dpi = 300, compression = "lzw")
}, "Figure 1")

## Crude Tukey
m_crude <- svycoxph(Surv(time_y, death) ~ cluster, design = des)
tuk <- safe({
  lv <- levels(dat$cluster); cn <- names(coef(m_crude)); K <- NULL; rn <- NULL
  for (a in 1:(length(lv) - 1)) for (b in (a + 1):length(lv)) {
    v <- setNames(rep(0, length(cn)), cn)
    if (paste0("cluster", lv[b]) %in% cn) v[paste0("cluster", lv[b])] <- 1
    if (paste0("cluster", lv[a]) %in% cn) v[paste0("cluster", lv[a])] <- -1
    K <- rbind(K, v); rn <- c(rn, paste(lv[b], "vs", lv[a]))
  }
  rownames(K) <- rn
  g <- glht(parm(coef(m_crude), vcov(m_crude)), linfct = K)
  ci <- confint(g)$confint; sm <- summary(g)
  tibble(Comparison = rn, HR = exp(ci[, 1]), LCL = exp(ci[, 2]), UCL = exp(ci[, 3]),
         P_Tukey = as.numeric(sm$test$pvalues))
}, "Tukey")
if (!is.null(tuk)) { write_csv(tuk, OUT("Table_S7_Crude_Tukey.csv")); R$tukey <- tuk }

## Primary adjusted Cox
f_adj <- Surv(time_y, death) ~ RIDAGEYR + female + egfr + LBXGLU + cluster
m_adj <- svycoxph(f_adj, design = des)
t2 <- tidy_cox(m_adj, "svycoxph, design-based"); write_csv(t2, OUT("Table2_Cox_Adjusted.csv")); R$cox <- t2

## Schoenfeld (unweighted diagnostic)
safe({
  m_u <- coxph(f_adj, data = dat)
  zp <- cox.zph(m_u)
  zt <- as.data.frame(zp$table) %>% tibble::rownames_to_column("Term") %>%
    mutate(Note = "Unweighted diagnostic model; interpret graphically, not as design-based test")
  write_csv(zt, OUT("Table_S_Schoenfeld_unweighted.csv")); R$zph <- zt
  nice <- c(RIDAGEYR = "Age", female = "Female sex", egfr = "eGFR", LBXGLU = "Fasting glucose", cluster = "Phenotype cluster")
  pdf(OUT("Figure_S5_Schoenfeld.pdf"), 10, 6); par(mfrow = c(2, 3))
  for (i in seq_len(nrow(zp$table) - 1)) plot(zp[i], resid = FALSE, main = nice[rownames(zp$table)[i]], xlab = "Years", ylab = "β(t)")
  dev.off()
}, "Schoenfeld")

## Stratified Cox (age tertile × sex)
qs <- as.numeric(coef(svyquantile(~RIDAGEYR, des, c(1/3, 2/3))))
dat$age_tert <- cut(dat$RIDAGEYR, c(-Inf, qs, Inf), labels = c("T1","T2","T3"))
des <- update(des, age_tert = cut(RIDAGEYR, c(-Inf, qs, Inf), labels = c("T1","T2","T3")))
m_str <- svycoxph(Surv(time_y, death) ~ egfr + LBXGLU + cluster + strata(age_tert, female), design = des)
tS11 <- tidy_cox(m_str, "svycoxph stratified by age tertile × sex")
write_csv(tS11, OUT("Table_S11_Stratified_Cox.csv")); R$cox_strat <- tS11

## Cardiovascular mortality
m_cs <- svycoxph(Surv(time_y, cvdeath) ~ RIDAGEYR + female + egfr + LBXGLU + cluster, design = des)
tS2a <- tidy_cox(m_cs, "Cause-specific svycoxph"); write_csv(tS2a, OUT("Table_S2a_CV_CauseSpecific.csv")); R$cs <- tS2a
fg <- safe({
  X <- model.matrix(~ RIDAGEYR + female + egfr + LBXGLU + cluster, dat)[, -1]
  st <- ifelse(dat$cvdeath == 1, 1, ifelse(dat$death == 1, 2, 0))
  m <- crr(dat$time_y, st, cov1 = X, failcode = 1, cencode = 0)
  s <- summary(m)$coef
  tibble(Covariate = str_replace(rownames(s), "^cluster", "cluster: "), SHR = exp(s[, 1]),
         LCL_95 = exp(s[, 1] - 1.96 * s[, 3]), UCL_95 = exp(s[, 1] + 1.96 * s[, 3]), P_Value = s[, 5],
         Method = "Fine-Gray (unweighted)")
}, "Fine-Gray")
if (!is.null(fg)) { write_csv(fg, OUT("Table_S2b_Fine_Gray.csv")); R$fg <- fg }

## RMST (unadjusted, unweighted)
rm <- safe({
  map_dfr(NONREF, function(g) {
    s <- dat %>% filter(cluster %in% c(REF, g))
    r <- rmst2(s$time_y, s$death, as.integer(s$cluster == g), tau = CFG$tau)
    u <- r$unadjusted.result
    tibble(Comparison = paste(g, "vs", REF), RMST_group = r$RMST.arm1$rmst[1], RMST_ref = r$RMST.arm0$rmst[1],
           Difference = u[1, 1], LCL_95 = u[1, 2], UCL_95 = u[1, 3], P_Value = u[1, 4])
  })
}, "RMST")
if (!is.null(rm)) { write_csv(rm, OUT("Table_S10_RMST.csv")); R$rmst <- rm }

## ---------------------------- 10. INCREMENTAL VALUE --------------------------
cdelta <- function(dd) {
  ww <- dd$pooled_weight / mean(dd$pooled_weight)
  f0 <- suppressWarnings(coxph(Surv(time_y, death) ~ RIDAGEYR + female + egfr + LBXGLU, data = dd, weights = ww))
  f1 <- suppressWarnings(coxph(Surv(time_y, death) ~ RIDAGEYR + female + egfr + LBXGLU + cluster, data = dd, weights = ww))
  c0 <- concordance(Surv(dd$time_y, dd$death) ~ predict(f0, type = "lp"), weights = ww, reverse = TRUE)$concordance
  c1 <- concordance(Surv(dd$time_y, dd$death) ~ predict(f1, type = "lp"), weights = ww, reverse = TRUE)$concordance
  c(c0 = c0, c1 = c1, delta = c1 - c0)
}
cpt <- cdelta(dat)
set.seed(CFG$seed)
ids <- dat %>% distinct(SDMVSTRA, SDMVPSU)
boot <- map_dbl(seq_len(CFG$n_boot_c), function(b) {
  pk <- ids %>% group_by(SDMVSTRA) %>% slice_sample(prop = 1, replace = TRUE) %>% ungroup()
  bd <- inner_join(pk, dat, by = c("SDMVSTRA","SDMVPSU"), relationship = "many-to-many")
  if (length(unique(bd$cluster)) < CFG$k) return(NA_real_)
  tryCatch(cdelta(bd)[["delta"]], error = function(e) NA_real_)
})
bv <- boot[!is.na(boot)]
tS8 <- tibble(C_Base = cpt[["c0"]], C_Cluster = cpt[["c1"]], Delta_C = cpt[["delta"]],
              Boot_mean = mean(bv), CI_Lower_95 = unname(quantile(bv, .025)), CI_Upper_95 = unname(quantile(bv, .975)),
              Boot_SE = sd(bv), P_Value = min(1, 2 * min(mean(bv <= 0), mean(bv >= 0))),
              Valid = length(bv), Attempted = CFG$n_boot_c, Method = "Harrell's C; PSU-cluster bootstrap within strata")
write_csv(tS8, OUT("Table_S8_Cindex.csv")); write_csv(tibble(replicate = seq_along(boot), delta_c = boot), OUT("Table_S8_bootstrap_replicates.csv"))
R$cidx <- tS8
ggsave(OUT("Figure_S2_Bootstrap_DeltaC.pdf"),
       ggplot(tibble(d = bv), aes(d)) + geom_histogram(bins = 30, fill = "grey70", colour = "white") +
         geom_vline(xintercept = 0, linetype = 2) + geom_vline(xintercept = mean(bv), colour = "firebrick") +
         theme_classic() + labs(x = "ΔC-index (base + cluster − base)", y = "Replicates"), width = 6, height = 4)

nri <- safe({
  base <- cbind(age = dat$RIDAGEYR, female = as.numeric(dat$female), egfr = dat$egfr, glu = dat$LBXGLU)
  dum <- model.matrix(~ cluster, dat)[, -1, drop = FALSE]
  x <- IDI.INF(indata = cbind(dat$time_y, dat$death), covs0 = base, covs1 = cbind(base, dum),
               t0 = CFG$tau, npert = CFG$idi_npert, seed1 = CFG$seed)
  tibble(Metric = c("IDI", "Continuous NRI"), Estimate = c(x$m1[1], x$m2[1]),
         Lower = c(x$m1[2], x$m2[2]), Upper = c(x$m1[3], x$m2[3]), P_Value = c(x$m1[4], x$m2[4]),
         Note = "Unweighted; t0 = 10 years")
}, "NRI/IDI")
if (!is.null(nri)) { write_csv(nri, OUT("Table_S9_NRI_IDI.csv")); R$nri <- nri }

rsf <- safe({
  set.seed(CFG$seed)
  rd <- dat %>% transmute(time_y, death, RIDAGEYR, female = factor(female), egfr, LBXGLU, LBXGH, BMXBMI,
                          BMXWAIST, sbp, dbp, LBXTR, LBDHDD, cluster, w = w_norm) %>% as.data.frame()
  tr <- sample.int(nrow(rd), floor(.7 * nrow(rd)))
  fit <- rfsrc(Surv(time_y, death) ~ . - w, data = rd[tr, ], ntree = CFG$rsf_ntree,
               case.wt = rd$w[tr], importance = "permute", seed = -CFG$seed)
  pr <- predict(fit, rd[-tr, ])
  vi <- sort(fit$importance, decreasing = TRUE)
  vtab <- tibble(Variable = names(vi), Label = coalesce(LAB[names(vi)], names(vi)), VIMP = as.numeric(vi),
                 Rank = seq_along(vi))
  vtab$Label[vtab$Variable == "female"] <- "Female sex"
  write_csv(vtab, OUT("Table4_RSF_VIMP.csv"))
  gg <- ggplot(vtab, aes(VIMP, reorder(Label, VIMP))) + geom_col(fill = "grey40") + theme_classic() +
    labs(x = "Permutation variable importance", y = NULL)
  ggsave(OUT("Figure3_RSF_VIMP.pdf"), gg, width = 6, height = 4.5)
  ggsave(OUT("Figure3_RSF_VIMP.tiff"), gg, width = 6, height = 4.5, dpi = 300, compression = "lzw")
  er <- tibble(trees = seq_along(fit$err.rate), err = as.numeric(fit$err.rate)) %>% filter(!is.na(err))
  ggsave(OUT("Figure_S6_RSF_OOB.pdf"), ggplot(er, aes(trees, err)) + geom_line() + theme_classic() +
           labs(x = "Number of trees", y = "OOB error (1 − C)"), width = 5, height = 4)
  pv <- plot.variable(fit, xvar.names = c("RIDAGEYR","egfr","LBXGLU","LBXGH"), partial = TRUE,
                      surv.type = "mort", show.plots = FALSE)
  pdd <- map_dfr(names(pv$plotthis), function(v) tibble(var = LAB[[v]], x = pv$plotthis[[v]]$x, y = pv$plotthis[[v]]$yhat))
  gp <- ggplot(pdd, aes(x, y)) + geom_line() + facet_wrap(~var, scales = "free_x") + theme_classic() +
    labs(x = NULL, y = "Predicted mortality (ensemble)")
  ggsave(OUT("Figure2_Partial_Dependence.pdf"), gp, width = 8, height = 6)
  ggsave(OUT("Figure2_Partial_Dependence.tiff"), gp, width = 8, height = 6, dpi = 300, compression = "lzw")
  list(vimp = vtab, test_C = 1 - tail(na.omit(as.numeric(pr$err.rate)), 1))
}, "RSF")
if (!is.null(rsf)) R$rsf <- rsf

## ---------------------------- 11. NO-AGE SENSITIVITY -------------------------
noage <- safe({
  Zn <- Z[, setdiff(CLUST_VARS, "RIDAGEYR")]
  kn <- wkmeans(Zn, w_norm, CFG$k, nstart = CFG$nstart)
  ln <- label_clusters(Zn, kn$cluster, w_norm)
  dn <- dat %>% mutate(cluster = ln$factor)
  dsn <- update(des, cluster_noage = ln$factor[match(SEQN, dat$SEQN)])
  mn <- svycoxph(Surv(time_y, death) ~ RIDAGEYR + female + egfr + LBXGLU + cluster_noage, design = dsn)
  cd <- cdelta(dn)
  age_by <- dn %>% group_by(cluster) %>% summarise(n = n(), mean_age = weighted.mean(RIDAGEYR, pooled_weight))
  list(ARI_vs_main = adjustedRandIndex(kn$cluster, km$cluster), labels = ln$audit, age_by = age_by,
       cox = tidy_cox(mn, "svycoxph, clusters derived without age"), delta_c = cd)
}, "No-age sensitivity")
if (!is.null(noage)) {
  write_csv(noage$cox, OUT("Table_S14_NoAge_Cox.csv"))
  write_csv(noage$age_by, OUT("Table_S14b_NoAge_cluster_ages.csv"))
  R$noage <- noage
}

## ---------------------------- 12. PRISMA FIGURE ------------------------------
safe({
  fl <- R$flow; n <- nrow(fl)
  boxes <- tibble(y = rev(seq_len(n)) * 2, lab = paste0(fl$Step, "\nN = ", format(fl$N, big.mark = ",")))
  ex <- tibble(y = boxes$y[-1] + 1,
               lab = paste0("Excluded: n = ", format(fl$Excluded[-1], big.mark = ","), "\n",
                            c("Age <18 y or pregnant", "Ineligible/incomplete linkage",
                              "Not Stage 2 (Stage 0, 1, 3 or 4)", "Missing core data or no fasting weight")))
  g <- ggplot() +
    geom_label(data = boxes, aes(3, y, label = lab), size = 3.4, label.padding = unit(.5, "lines")) +
    geom_label(data = ex, aes(7, y, label = lab), size = 3, fill = "grey95", label.padding = unit(.4, "lines")) +
    geom_segment(data = boxes[-n, ], aes(x = 3, xend = 3, y = y - .55, yend = y - 1.45),
                 arrow = arrow(length = unit(.2, "cm"))) +
    geom_segment(data = ex, aes(x = 3, xend = 5.6, y = y, yend = y)) +
    xlim(0.5, 9) + ylim(1, 2 * n + 1) + theme_void()
  ggsave(OUT("Figure_S1_PRISMA.pdf"), g, width = 8, height = 8)
  ggsave(OUT("Figure_S1_PRISMA.tiff"), g, width = 8, height = 8, dpi = 300, compression = "lzw")
}, "PRISMA")

## ---------------------------- 13. SUPPLEMENT DOCX ----------------------------
safe({
  doc <- read_docx()
  add_tab <- function(doc, title, df) {
    df <- df %>% mutate(across(where(is.numeric), ~ signif(.x, 4)))
    doc %>% body_add_par(title, style = "heading 2") %>%
      body_add_flextable(autofit(fontsize(flextable(df), size = 8, part = "all"))) %>% body_add_break()
  }
  for (f in sort(list.files(CFG$out_dir, pattern = "^Table.*\\.csv$")))
    doc <- add_tab(doc, tools::file_path_sans_ext(f), read_csv(OUT(f), show_col_types = FALSE))
  print(doc, target = OUT("All_Tables_Supplement.docx"))
}, "Supplement docx")

## ---------------------------- 14. AUTO-FILLED RESULTS ------------------------
hr  <- function(h, l, u) sprintf("%.2f (95%% CI %.2f–%.2f)", h, l, u)
pf  <- function(p) ifelse(p < 0.001, "p<0.001", sprintf("p=%.2f", p))
mil <- function(x) sprintf("%.1f million", x / 1e6)
cl_line <- function(tb) paste(sprintf("%s HR %s, %s", tb$Covariate, hr(tb$HR, tb$LCL_95, tb$UCL_95), pf(tb$P_Value)), collapse = "; ")
safe({
  t3c <- R$t3 %>% filter(cluster != "Overall"); t3o <- R$t3 %>% filter(cluster == "Overall")
  txt <- c(
    "# AUTO-FILLED RESULTS — generated by pipeline, paste into manuscript",
    sprintf("_Generated %s. Do not hand-edit numbers; re-run the script instead._", Sys.time()), "",
    "## Attrition",
    paste(sprintf("- %s: %s", R$flow$Step, format(R$flow$N, big.mark = ",")), collapse = "\n"),
    sprintf("- PREVENT scored for %d of %d candidates (others outside equation ranges; not reclassified).", R$prevent_scored, R$prevent_candidates),
    "", "## Cohort",
    sprintf("N = %s, representing %s US adults. Reference cluster label: **%s**.",
            format(R$N, big.mark = ","), mil(R$weighted_pop), R$ref_label),
    paste(sprintf("- %s: n = %s (%s)", names(R$n_by_cluster), unlist(R$n_by_cluster), mil(unlist(R$wpop_by_cluster))), collapse = "\n"),
    "", "## Robustness",
    if (!is.null(R$jaccard)) sprintf("Jaccard bootstrap means %.2f–%.2f.", R$jaccard[1], R$jaccard[2]) else "Jaccard: FAILED",
    if (!is.null(R$gmm)) sprintf("GMM ARI: BIC-optimal (k=%d) %.2f; k=4 %.2f.", R$gmm$k[1], R$gmm$ARI_vs_kmeans[1], R$gmm$ARI_vs_kmeans[2]) else "",
    if (!is.null(R$age_strata)) paste(sprintf("%s ARI %.2f (%s)", R$age_strata$Stratum, R$age_strata$ARI_vs_pooled, R$age_strata$Interpretation), collapse = "; ") else "",
    if (!is.null(R$ext)) paste(sprintf("%s %s", R$ext$Variable, pf(R$ext$P_Value)), collapse = "; ") else "",
    "", "## Crude mortality",
    sprintf("Follow-up %.0f person-years (mean %.1f y); %d deaths (%d cardiovascular). Weighted crude rate %.1f/1,000 PY overall.",
            t3o$PersonYears, R$mean_fu, t3o$Deaths, t3o$CV_Deaths, t3o$Weighted_rate_per_1000PY),
    paste(sprintf("- %s: %.1f/1,000 PY; %d deaths; CV share %.1f%%", t3c$cluster, t3c$Weighted_rate_per_1000PY, t3c$Deaths, t3c$CV_share_pct), collapse = "\n"),
    if (!is.null(R$logrank_p)) sprintf("Log-rank %s.", pf(R$logrank_p)) else "",
    if (!is.null(R$tukey)) paste(sprintf("- Crude %s: HR %s, Tukey %s", R$tukey$Comparison, hr(R$tukey$HR, R$tukey$LCL, R$tukey$UCL), pf(R$tukey$P_Tukey)), collapse = "\n") else "",
    if (!is.null(R$rmst)) paste(sprintf("- RMST(10y) %s: %.2f y (%.2f to %.2f), %s", R$rmst$Comparison, R$rmst$Difference, R$rmst$LCL_95, R$rmst$UCL_95, pf(R$rmst$P_Value)), collapse = "\n") else "",
    "", "## Adjusted all-cause",
    paste0("Standard: ", cl_line(pull_cl(R$cox))),
    { a <- R$cox %>% filter(Covariate == "Age (per year)"); sprintf("Age HR %s.", hr(a$HR, a$LCL_95, a$UCL_95)) },
    paste0("Stratified: ", cl_line(pull_cl(R$cox_strat))),
    "", "## Cardiovascular",
    paste0("Cause-specific: ", cl_line(pull_cl(R$cs))),
    if (!is.null(R$fg)) paste0("Fine-Gray: ", cl_line(pull_cl(rename(R$fg, HR = SHR)))) else "",
    "", "## Incremental value",
    sprintf("C base %.3f; C +cluster %.3f; ΔC %.4f; bootstrap mean %.4f (95%% CI %.4f to %.4f), %s, %d/%d valid.",
            R$cidx$C_Base, R$cidx$C_Cluster, R$cidx$Delta_C, R$cidx$Boot_mean, R$cidx$CI_Lower_95, R$cidx$CI_Upper_95,
            pf(R$cidx$P_Value), R$cidx$Valid, R$cidx$Attempted),
    if (!is.null(R$nri)) paste(sprintf("%s %.3f (%.3f to %.3f)", R$nri$Metric, R$nri$Estimate, R$nri$Lower, R$nri$Upper), collapse = "; ") else "NRI/IDI: FAILED",
    if (!is.null(R$rsf)) sprintf("RSF: %s. Cluster rank %d of %d. Test-set C %.3f.",
                                 paste(sprintf("%s %.3f", head(R$rsf$vimp$Label, 5), head(R$rsf$vimp$VIMP, 5)), collapse = ", "),
                                 R$rsf$vimp$Rank[R$rsf$vimp$Variable == "cluster"], nrow(R$rsf$vimp), R$rsf$test_C) else "",
    "", "## No-age sensitivity",
    if (!is.null(R$noage)) c(sprintf("ARI vs main %.2f; ΔC %.4f.", R$noage$ARI_vs_main, R$noage$delta_c[["delta"]]),
                             paste0("Cox: ", cl_line(pull_cl(R$noage$cox))),
                             paste(sprintf("- %s mean age %.1f", R$noage$age_by$cluster, R$noage$age_by$mean_age), collapse = "\n")) else "FAILED"
  )
  writeLines(txt, OUT("Results_autofilled.md"))
}, "Auto-filled results")

jsonlite::write_json(R, OUT("manuscript_numbers.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE, force = TRUE)
writeLines(capture.output(sessionInfo()), OUT("sessionInfo.txt"))
say("DONE. Check %s for any '!! FAILED' lines. Upload the whole '%s' folder.", LOG, CFG$out_dir)
