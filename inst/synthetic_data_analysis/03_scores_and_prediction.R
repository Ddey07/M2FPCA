# =============================================================================
# 03_scores_and_prediction.R
#
# (a) Compute subject-level multivariate FPC scores (digital biomarkers) with
#     mpfpca_scores(), which predicts each domain's latent trajectory
#     (getLatentPreds) and projects onto the shared eigenfunctions.
# (b) Reconstruct the full multivariate covariance from the partially separable
#     fit and compare to the analytic truth used to generate the data.
# (c) Multivariate PREDICTION: use the leading component scores as features in a
#     multinomial classifier of diagnostic group (mirrors the paper's
#     multinomial-LASSO digital-phenotyping step), with a train/test split.
#
# Resumable: latent predictions are the slow part and are cached per domain.
# =============================================================================

suppressMessages(library(M2FPCA))
outdir <- Sys.getenv("SDA_OUT", unset = "output")
st  <- readRDS(file.path(outdir, "synthetic_data.rds"))
mpf <- readRDS(file.path(outdir, "mpf.rds"))
dat_list <- st$dat_list; type <- st$type; argvals <- st$argvals
m <- st$m; p <- st$p; n <- st$n; group <- st$group

t0 <- Sys.time(); budget <- as.numeric(Sys.getenv("SDA_BUDGET", "600"))
save_atomic <- function(o, path){tmp<-paste0(path,".tmp");saveRDS(o,tmp);file.rename(tmp,path)}

## ---- (a) latent predictions (cached per domain) --------------------------
gl <- getFromNamespace(".call_getLatentPreds", "M2FPCA")
latf <- file.path(outdir, "latent.rds")
latent <- if (file.exists(latf)) readRDS(latf) else vector("list", p)
for (v in 1:p) {
  if (!is.null(latent[[v]])) next
  if (as.numeric(Sys.time() - t0, units = "secs") > budget - 20) {
    cat("PAUSE before latent domain", v, "\n"); save_atomic(latent, latf); quit(status = 3)
  }
  ti <- Sys.time()
  # impute.missing = TRUE so the latent trajectory is complete on the grid;
  # required for irregular/sparse EMA data before projecting onto phi_l.
  z <- gl(X = dat_list[[v]], type = rep(type[v], m),
          lat.cov.est = mpf$cov_marginal[[v]], ridge = 0.05, impute.missing = TRUE)
  latent[[v]] <- pmin(pmax(z, -8), 8)
  save_atomic(latent, latf)
  cat(sprintf("  latent %d (%s) in %.1fs\n", v, type[v], as.numeric(Sys.time()-ti,units="secs")))
}

## ---- project onto shared eigenfunctions: theta_{i,v,l} -------------------
trapz <- getFromNamespace("trapz", "M2FPCA")
npc <- min(4, mpf$L)
mu_list <- lapply(latent, colMeans)
scores <- vector("list", npc)
for (l in 1:npc) {
  Sl <- matrix(NA, n, p)
  for (v in 1:p) {
    phi_l <- mpf$pfpca_results[[v]]$phi[, l]
    Sl[, v] <- apply(latent[[v]], 1, function(row) trapz(argvals, (row - mu_list[[v]]) * phi_l))
  }
  colnames(Sl) <- paste0(c("Mood","Anx","Energy","TLAC"), "_PC", l)
  scores[[l]] <- Sl
}
scores_wide <- do.call(cbind, scores)

cat("\n== Variance explained by shared components (FVE, %) ==\n")
fve <- mpf$pfpca_results[[1]]$cumFVE[1:npc]
print(round(data.frame(PC = 1:npc, cumFVE = fve,
                       lambda = mpf$pfpca_results[[1]]$lambda[1:npc]), 3))

cat("\n== Cross-domain score correlation, component 1 (should reflect coupling) ==\n")
print(round(cor(scores[[1]], use = "pairwise.complete.obs"), 2))

## ---- (b) reconstructed vs true multivariate covariance -------------------
Chat <- cov.from.mpfpca.dir(mpf)
# analytic truth: full latent covariance = Rbase (x) K, K = sum_l lambda_l phi phi'
phi <- st$phi; lambda <- st$lambda; Rbase <- st$Rbase
K_true <- Reduce("+", lapply(1:ncol(phi), function(l) lambda[l] * outer(phi[,l], phi[,l])))
R_true <- kronecker(Rbase, cov2cor(K_true))
Chat_c <- cov2cor(Chat)
rel_err <- norm(Chat_c - R_true, "F") / norm(R_true, "F")
vcor    <- cor(as.numeric(Chat_c), as.numeric(R_true))
cat(sprintf("\n== Covariance recovery ==\n  rel. Frobenius error: %.3f | cor(vec(true), vec(est)): %.3f\n",
            rel_err, vcor))

## ---- (c) multivariate prediction: multinomial classification -------------
suppressMessages(library(nnet))
feats <- as.data.frame(scores_wide)
feats$group <- group
ok <- stats::complete.cases(scores_wide)
if (any(!ok)) cat(sprintf("  dropping %d subject-days with incomplete latent projection\n", sum(!ok)))
feats <- feats[ok, ]
set.seed(7)
tr <- unlist(lapply(split(seq_len(nrow(feats)), feats$group), function(idx) sample(idx, max(1, round(0.7*length(idx))))))
train <- feats[tr, ]; test <- feats[-tr, ]

fit_mn <- multinom(group ~ ., data = train, trace = FALSE, maxit = 500)
pred   <- predict(fit_mn, newdata = test)
acc    <- mean(pred == test$group)
cm     <- table(true = test$group, predicted = pred)
base   <- max(table(train$group)) / nrow(train)   # majority-class baseline

cat(sprintf("\n== Multivariate prediction of diagnostic group ==\n"))
cat(sprintf("  features: %d component scores (%d PCs x %d domains)\n", ncol(scores_wide), npc, p))
cat(sprintf("  test accuracy: %.2f   (majority baseline %.2f, chance %.2f)\n",
            acc, base, 1/length(st$groups)))
cat("  confusion matrix (rows=true):\n"); print(cm)

# variable-importance proxy: sum |coef| across logits per feature
imp <- colSums(abs(coef(fit_mn)))[-1]
imp <- sort(imp, decreasing = TRUE)
cat("\n  top predictive component scores (|coef| across logits):\n")
print(round(head(imp, 6), 2))

## ---- persist results -----------------------------------------------------
res <- list(scores = scores, scores_wide = scores_wide, fve = fve, npc = npc,
            score_cor_pc1 = cor(scores[[1]]), Chat = Chat, R_true = R_true,
            rel_err = rel_err, vcor = vcor, latent = latent,
            pred = list(acc = acc, baseline = base, chance = 1/length(st$groups),
                        cm = cm, importance = imp, test_group = test$group, pred = pred))
save_atomic(res, file.path(outdir, "results.rds"))
cat("\nSaved:", file.path(outdir, "results.rds"), "\nSCORES_DONE\n")
