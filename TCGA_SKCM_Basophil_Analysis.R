# ============================================================
# TCGA-SKCM 嗜碱性粒细胞预后与转录组分析
# ============================================================
# 作者：乔雨慧
# 日期：2026年4月
# 描述：本研究基于TCGA-SKCM队列，系统性评估外周血嗜碱性粒细胞
#       在黑色素瘤中的预后价值，并通过全转录组相关性分析和
#       功能富集探究其潜在免疫机制。
# ============================================================

# 0. 环境准备-----------------------------------------------------
# 首次运行请取消注释安装所需包
# install.packages(c("tidyverse", "glmnet", "survival", "survminer", 
#                    "ggplot2", "ggrepel", "rms", "Hmisc", "survivalROC"))
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))

library(tidyverse)
library(glmnet)
library(survival)
library(survminer)
library(ggplot2)
library(ggrepel)
library(rms)
library(Hmisc)
library(survivalROC)
library(clusterProfiler)
library(org.Hs.eg.db)

# 设置工作目录（请修改为你的实际路径）
# setwd("D:/datatjjm")

# ============================================================
# 第一部分：嗜碱性粒细胞浸润度计算
# ============================================================

# 1.1 读取并清洗转录组数据 -----------------------------------------------
cat("\n===== 第一部分：计算嗜碱性粒细胞浸润度 =====\n")
exp_raw <- read.delim("TCGA-SKCM.star_counts.tsv.gz", 
                      stringsAsFactors = FALSE, 
                      row.names = 1, 
                      check.names = FALSE)

# 过滤低表达基因：在至少10%样本中表达>0
exp_clean1 <- exp_raw[rowSums(exp_raw > 0) >= ncol(exp_raw) * 0.1, ]
# 去除重复基因名
exp_clean2 <- exp_clean1[!duplicated(rownames(exp_clean1)), ]
# 去除含NA的行
exp_clean <- exp_clean2[complete.cases(exp_clean2), ]

cat("原始基因数：", nrow(exp_raw), "\n")
cat("清洗后基因数：", nrow(exp_clean), "\n")

# 处理Ensembl ID（去掉版本号，确保行名唯一）
temp_rownames <- sapply(strsplit(rownames(exp_clean), "\\."), `[`, 1)
dup_rows <- which(duplicated(temp_rownames))
for (i in dup_rows) {
  temp_rownames[i] <- paste0(temp_rownames[i], "_", sum(temp_rownames[1:i] == temp_rownames[i]))
}
rownames(exp_clean) <- temp_rownames
cat("行名是否唯一：", !any(duplicated(rownames(exp_clean))), "\n")

# 1.2 定义嗜碱性粒细胞特征基因 -------------------------------------------
baso_genes <- c("CCL2", "CCL7", "CXCL8", "IL4", "IL5", "IL13",
                "FCER1A", "FCER1G", "CPA3", "MS4A2", "TPSAB1", "IL3RA")

# Symbol转Ensembl ID
ensg_map <- mapIds(org.Hs.eg.db, keys = baso_genes, 
                   keytype = "SYMBOL", column = "ENSEMBL", multiVals = "first")
ensg_map <- na.omit(ensg_map)
ensg_ids <- as.character(ensg_map)
ensg_ids <- sapply(strsplit(ensg_ids, "\\."), `[`, 1)

# 1.3 计算浸润度（特征基因平均表达）--------------------------------------
exp_clean_rownames <- sapply(strsplit(rownames(exp_clean), "\\."), `[`, 1)
use_ensg <- intersect(ensg_ids, exp_clean_rownames)
cat("实际使用的特征基因数：", length(use_ensg), "\n")

baso_expr <- exp_clean[exp_clean_rownames %in% use_ensg, ]
baso_data <- data.frame(
  submitter_id = colnames(baso_expr),
  basophils = colMeans(baso_expr, na.rm = TRUE)
) %>%
  mutate(submitter_id = substr(gsub("\\.", "-", submitter_id), 1, 12))

head(baso_data)

# ============================================================
# 第二部分：临床数据处理与生存分析
# ============================================================

cat("\n===== 第二部分：临床生存分析 =====\n")

# 2.1 读取并清洗临床数据 -------------------------------------------------
clin_raw <- read.delim("TCGA-SKCM.clinical.tsv.gz", 
                       stringsAsFactors = FALSE, 
                       check.names = FALSE)

df_clin_final <- clin_raw %>%
  dplyr::select(sample, case_id, disease_type,
                vital_status.demographic,
                age_at_index.demographic,
                days_to_death.demographic,
                days_to_last_follow_up.diagnoses) %>%
  dplyr::rename(sample_id = sample,
                vital_status = vital_status.demographic,
                days_to_death = days_to_death.demographic,
                days_to_last_follow_up = days_to_last_follow_up.diagnoses,
                age = age_at_index.demographic) %>%
  dplyr::filter(!is.na(vital_status),
                vital_status %in% c("Dead", "Alive"),
                !is.na(age), age > 0) %>%
  dplyr::mutate(
    OS.time = ifelse(vital_status == "Dead", days_to_death, days_to_last_follow_up),
    OS = as.integer(vital_status == "Dead")
  ) %>%
  dplyr::filter(!is.na(OS.time), OS.time > 0, !is.na(OS))

# 2.2 合并临床与嗜碱性粒细胞数据 -----------------------------------------
# 统一ID格式
df_clin_final$sample_id <- substr(gsub("\\.", "-", df_clin_final$sample_id), 1, 12)
colnames(baso_data)[colnames(baso_data) == "submitter_id"] <- "sample_id"
baso_data$sample_id <- substr(gsub("\\.", "-", baso_data$sample_id), 1, 12)

# 去重并取交集
baso_data <- baso_data[!duplicated(baso_data$sample_id), ]
df_clin_final <- df_clin_final[!duplicated(df_clin_final$sample_id), ]
common_ids <- intersect(baso_data$sample_id, df_clin_final$sample_id)
cat("共同样本数：", length(common_ids), "\n")

df_analysis <- merge(
  baso_data[baso_data$sample_id %in% common_ids, ],
  df_clin_final[df_clin_final$sample_id %in% common_ids, ],
  by = "sample_id"
)

# 2.3 LASSO-Cox回归构建风险评分 -----------------------------------------
set.seed(123)
x <- as.matrix(scale(df_analysis[, c("basophils", "age")]))
y <- Surv(time = df_analysis$OS.time, event = df_analysis$OS)

cv_fit <- cv.glmnet(x = x, y = y, family = "cox", nfolds = 10, alpha = 1)
plot(cv_fit)
title("LASSO-Cox Cross-Validation")

final_lasso <- glmnet(x = x, y = y, family = "cox", alpha = 1, lambda = cv_fit$lambda.min)
cat("LASSO系数：\n")
print(coef(final_lasso))

# 计算风险评分
df_analysis$risk_score <- predict(final_lasso, newx = x, type = "link")[, 1]

# 2.4 C-index及Bootstrap置信区间 ----------------------------------------
c_index_temp <- rcorr.cens(df_analysis$risk_score, Surv(df_analysis$OS.time, df_analysis$OS))[1]
if (c_index_temp < 0.5) {
  df_analysis$risk_score <- -1 * df_analysis$risk_score
  cat("模型方向已自动纠正\n")
}

boot_cindex <- function(i) {
  idx <- sample(nrow(df_analysis), replace = TRUE)
  rcorr.cens(df_analysis$risk_score[idx], Surv(df_analysis$OS.time, df_analysis$OS)[idx])["C Index"]
}
set.seed(123)
boot_results <- replicate(1000, boot_cindex())
ci_lower <- quantile(boot_results, 0.025)
ci_upper <- quantile(boot_results, 0.975)
cat("C-index =", round(mean(boot_results), 3), "\n")
cat("95%CI = [", round(ci_lower, 3), ",", round(ci_upper, 3), "]\n")

# 2.5 Kaplan-Meier生存曲线 ----------------------------------------------
df_analysis$risk_group <- ifelse(df_analysis$risk_score > median(df_analysis$risk_score),
                                 "Low Risk", "High Risk")
fit_km <- survfit(Surv(OS.time, OS) ~ risk_group, data = df_analysis)

ggsurvplot(fit_km, data = df_analysis,
           pval = TRUE, risk.table = TRUE,
           palette = c("#2E9FDF", "#E7B800"),
           legend.labs = c("Low Risk", "High Risk"),
           title = "Kaplan-Meier Survival Curve")

# 2.6 时间依赖性ROC曲线 --------------------------------------------------
roc_1year <- survivalROC(Stime = df_analysis$OS.time, status = df_analysis$OS,
                         marker = df_analysis$risk_score, predict.time = 365,
                         method = "NNE", span = 0.25)
roc_3year <- survivalROC(Stime = df_analysis$OS.time, status = df_analysis$OS,
                         marker = df_analysis$risk_score, predict.time = 365 * 3,
                         method = "NNE", span = 0.25)
roc_5year <- survivalROC(Stime = df_analysis$OS.time, status = df_analysis$OS,
                         marker = df_analysis$risk_score, predict.time = 365 * 5,
                         method = "NNE", span = 0.25)

plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "1-Specificity", ylab = "Sensitivity",
     main = "Time-dependent ROC Curves", las = 1)
abline(a = 0, b = 1, lty = 2, col = "gray60")
lines(roc_1year$FP, roc_1year$TP, col = "#E63946", lwd = 2, type = "s")
lines(roc_3year$FP, roc_3year$TP, col = "#1D3557", lwd = 2, type = "s")
lines(roc_5year$FP, roc_5year$TP, col = "#F0E442", lwd = 2, type = "s")
legend("bottomright", legend = c("1-Year", "3-Year", "5-Year"),
       col = c("#E63946", "#1D3557", "#F0E442"), lwd = 2, box.col = "black", bg = "white")

# 2.7 C-index森林图 -----------------------------------------------------
plot_data <- data.frame(model = "Prognostic Model",
                        c_index = mean(boot_results),
                        lower = ci_lower, upper = ci_upper)
ggplot(plot_data, aes(x = c_index, y = model)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_point(size = 6, color = "#E63946", shape = 15) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray60") +
  scale_x_continuous(limits = c(0.5, 0.7)) +
  labs(x = "C-index (95% Confidence Interval)", y = "",
       title = "Concordance Index of Prognostic Model") +
  theme_bw()

# ============================================================
# 第三部分：转录组相关性分析
# ============================================================

cat("\n===== 第三部分：全转录组相关性分析 =====\n")

# 3.1 读取TPM数据并匹配样本 ----------------------------------------------
tpm <- read.table("TCGA-SKCM.star_tpm.tsv.gz", header = TRUE, row.names = 1)

colnames(tpm) <- gsub("\\.", "-", colnames(tpm))
colnames(tpm) <- substr(colnames(tpm), 1, 12)
tpm <- tpm[, !duplicated(colnames(tpm))]

common_samples <- intersect(colnames(tpm), baso_data$sample_id)
cat("匹配样本数:", length(common_samples), "\n")

tpm_filtered <- tpm[, common_samples]
baso_filtered <- baso_data[match(common_samples, baso_data$sample_id), ]

# 3.2 过滤低表达基因 ----------------------------------------------------
keep_genes <- rowSums(tpm_filtered > 1) >= (0.3 * ncol(tpm_filtered))
tpm_filtered <- tpm_filtered[keep_genes, ]
cat("过滤后剩余基因数:", nrow(tpm_filtered), "\n")

# 3.3 批量Spearman相关性计算 --------------------------------------------
cor_results <- data.frame(gene = rownames(tpm_filtered), cor = NA, pval = NA)
baso_vals <- baso_filtered$basophils

cat("开始计算相关性...\n")
pb <- txtProgressBar(min = 0, max = nrow(tpm_filtered), style = 3)
for (i in 1:nrow(tpm_filtered)) {
  gene_exp <- as.numeric(tpm_filtered[i, ])
  if (sum(!is.na(gene_exp)) >= 10 && length(unique(gene_exp)) > 2) {
    test <- tryCatch({
      cor.test(gene_exp, baso_vals, method = "spearman", exact = FALSE)
    }, error = function(e) NULL)
    if (!is.null(test)) {
      cor_results$cor[i] <- test$estimate
      cor_results$pval[i] <- test$p.value
    }
  }
  setTxtProgressBar(pb, i)
}
close(pb)

cor_results <- cor_results[!is.na(cor_results$pval), ]
cat("\n完成！共计算", nrow(cor_results), "个基因的相关性\n")

# 3.4 筛选显著基因 ------------------------------------------------------
sig_pos <- cor_results[cor_results$pval < 0.05 & cor_results$cor > 0.3, ]
sig_pos <- sig_pos[order(sig_pos$cor, decreasing = TRUE), ]
sig_neg <- cor_results[cor_results$pval < 0.05 & cor_results$cor < -0.3, ]
sig_neg <- sig_neg[order(sig_neg$cor), ]

cat("显著正相关基因:", nrow(sig_pos), "\n")
cat("显著负相关基因:", nrow(sig_neg), "\n")

# ============================================================
# 第四部分：GO功能富集分析
# ============================================================

cat("\n===== 第四部分：GO功能富集分析 =====\n")

sig_pos$ensembl <- gsub("\\..*", "", sig_pos$gene)

entrez_ids <- bitr(sig_pos$ensembl[1:500], 
                   fromType = "ENSEMBL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

go_pos <- enrichGO(gene = entrez_ids$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05)

cat("Top10 GO富集通路：\n")
print(head(go_pos@result[, c("Description", "p.adjust", "Count")], 10))

# GO气泡图
go_df <- as.data.frame(go_pos@result)
go_df <- go_df[order(go_df$p.adjust), ][1:12, ]
go_df$GeneRatio <- sapply(go_df$GeneRatio, function(x) eval(parse(text = x)))

ggplot(go_df, aes(x = GeneRatio, y = reorder(Description, -p.adjust))) +
  geom_point(aes(size = Count, color = p.adjust), alpha = 0.8) +
  scale_color_gradient(low = "#E74C3C", high = "#3498DB", trans = "log10") +
  scale_size_continuous(range = c(4, 10)) +
  labs(x = "Gene Ratio", y = "", title = "GO Enrichment: Basophil-Correlated Genes") +
  theme_bw()

# ============================================================
# 第五部分：可视化
# ============================================================

cat("\n===== 第五部分：火山图与Top10基因 =====\n")

# 5.1 火山图 ------------------------------------------------------------
volcano_data <- cor_results
volcano_data$logP <- -log10(volcano_data$pval)
volcano_data$sig <- "Not Significant"
volcano_data$sig[volcano_data$pval < 0.05 & volcano_data$cor > 0.3] <- "Positive"
volcano_data$sig[volcano_data$pval < 0.05 & volcano_data$cor < -0.3] <- "Negative"

top10_genes <- head(sig_pos$gene, 10)
top10_ensembl <- gsub("\\..*", "", top10_genes)
top10_symbols <- AnnotationDbi::select(org.Hs.eg.db, keys = top10_ensembl,
                                        keytype = "ENSEMBL", columns = "SYMBOL")

key_genes <- c("CPA3", "CCL2", "TPSAB1", "MMP19", "PDPN", "CSF2RB")

volcano_data$symbol <- ""
for (i in 1:nrow(volcano_data)) {
  ens <- gsub("\\..*", "", volcano_data$gene[i])
  if (ens %in% top10_symbols$ENSEMBL) {
    sym <- top10_symbols$SYMBOL[match(ens, top10_symbols$ENSEMBL)]
    if (sym %in% key_genes) {
      volcano_data$symbol[i] <- sym
    }
  }
}

ggplot(volcano_data, aes(x = cor, y = logP, color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Positive" = "#E74C3C", 
                                "Negative" = "#3498DB",
                                "Not Significant" = "gray70")) +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  geom_text_repel(data = subset(volcano_data, symbol != ""),
                  aes(label = symbol), size = 5, max.overlaps = 20,
                  box.padding = 0.8, point.padding = 0.5, force = 2) +
  labs(x = "Spearman Correlation (r)", y = "-log10(P-value)",
       title = "Transcriptome-Wide Correlation with Basophil Levels") +
  theme_bw() +
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = "white", color = "gray80"),
        legend.title = element_blank())

# 5.2 Top10基因表格 -----------------------------------------------------
top10_info <- AnnotationDbi::select(org.Hs.eg.db, keys = top10_ensembl,
                                     keytype = "ENSEMBL",
                                     columns = c("SYMBOL", "GENENAME"))
top10_result <- data.frame(ensembl = top10_ensembl,
                           cor = head(sig_pos$cor, 10),
                           pval = head(sig_pos$pval, 10))
top10_result <- merge(top10_result, top10_info, by.x = "ensembl", by.y = "ENSEMBL")
top10_result <- top10_result[order(top10_result$cor, decreasing = TRUE), ]
cat("\nTop10正相关基因：\n")
print(top10_result[, c("SYMBOL", "GENENAME", "cor", "pval")])

cat("\n===== 全部分析完成 =====\n")
