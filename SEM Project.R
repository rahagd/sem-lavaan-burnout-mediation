# ============================================================
# مشروع تعريفي: النمذجة بالمعادلات البنائية في R
# الضغط الوظيفي -> الدعم -> الاحتراق الوظيفي -> الأداء
# تنبيه: البيانات محاكاة لأغراض تعليمية فقط
# ============================================================

library(lavaan)
library(semPlot)
library(psych)

# 1) توليد بيانات محاكاة --------------------------------------
set.seed(123)
n <- 300

stress  <- rnorm(n)
support <- -0.40*stress + rnorm(n, sd = 0.9)
burnout <-  0.45*stress - 0.35*support + rnorm(n, sd = 0.8)
perf    <- -0.40*burnout + 0.25*support - 0.10*stress + rnorm(n, sd = 0.8)

make_items <- function(f, prefix, loads = c(0.80, 0.75, 0.70)) {
  items <- sapply(loads, function(l) {
    x <- l*f + rnorm(length(f), sd = sqrt(1 - l^2))
    as.integer(cut(x, breaks = c(-Inf, -1.2, -0.4, 0.4, 1.2, Inf), labels = 1:5))
  })
  colnames(items) <- paste0(prefix, 1:3)
  items
}

data <- data.frame(
  make_items(stress,  "ST"),
  make_items(support, "SU"),
  make_items(burnout, "BO"),
  make_items(perf,    "PE")
)

# 2) نموذج القياس (CFA) ---------------------------------------
model_cfa <- '
  ST =~ ST1 + ST2 + ST3
  SU =~ SU1 + SU2 + SU3
  BO =~ BO1 + BO2 + BO3
  PE =~ PE1 + PE2 + PE3
'
fit_cfa <- cfa(model_cfa, data = data, estimator = "MLR")
summary(fit_cfa, fit.measures = TRUE, standardized = TRUE)

# 3) الثبات والصدق --------------------------------------------
alpha(data[, c("ST1","ST2","ST3")])$total$raw_alpha
alpha(data[, c("SU1","SU2","SU3")])$total$raw_alpha
alpha(data[, c("BO1","BO2","BO3")])$total$raw_alpha
alpha(data[, c("PE1","PE2","PE3")])$total$raw_alpha

est  <- standardizedSolution(fit_cfa)
load <- est[est$op == "=~", c("lhs", "est.std")]
res  <- sapply(split(load$est.std, load$lhs), function(l) {
  c(CR = sum(l)^2 / (sum(l)^2 + sum(1 - l^2)), AVE = mean(l^2))
})
round(t(res), 3)
round(lavInspect(fit_cfa, "cor.lv"), 2)
round(sqrt(res["AVE", ]), 3)

# 4) النموذج البنائي مع الوساطة -------------------------------
model_sem <- '
  ST =~ ST1 + ST2 + ST3
  SU =~ SU1 + SU2 + SU3
  BO =~ BO1 + BO2 + BO3
  PE =~ PE1 + PE2 + PE3

  SU ~ a1*ST
  BO ~ a2*ST + b1*SU
  PE ~ c*ST + d1*BO + d2*SU

  ST_SU_PE    := a1*d2
  ST_BO_PE    := a2*d1
  ST_SU_BO_PE := a1*b1*d1
  indirect    := a1*d2 + a2*d1 + a1*b1*d1
  total       := c + indirect
'
fit_sem <- sem(model_sem, data = data, estimator = "MLR")
summary(fit_sem, fit.measures = TRUE, standardized = TRUE)
fitMeasures(fit_sem, c("chisq","df","pvalue","cfi","tli","rmsea","srmr"))

# 5) Bootstrap للآثار غير المباشرة -----------------------------
fit_boot <- sem(model_sem, data = data, se = "bootstrap", bootstrap = 1000)
pe <- parameterEstimates(fit_boot, boot.ci.type = "perc", standardized = TRUE)
pe[pe$op == ":=", c("label","est","ci.lower","ci.upper","pvalue")]

# 6) المسارات المباشرة و R2 ------------------------------------
pe2 <- parameterEstimates(fit_sem, standardized = TRUE)
pe2[pe2$op == "~", c("lhs","rhs","est","se","pvalue","std.all")]
lavInspect(fit_sem, "rsquare")[c("SU","BO","PE")]

# 7) رسم النموذج ------------------------------------------------
semPaths(fit_sem, what = "std", layout = "tree2", edge.label.cex = 1.1,
         sizeMan = 7, sizeLat = 10, residuals = FALSE)
