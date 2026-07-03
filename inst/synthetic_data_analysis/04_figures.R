# =============================================================================
# 04_figures.R  --  render result figures (base graphics, no extra deps)
# =============================================================================
suppressMessages(library(M2FPCA))
outdir <- Sys.getenv("SDA_OUT", unset = "output")
figdir <- Sys.getenv("SDA_FIG", unset = file.path(dirname(outdir), "figures"))
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

st  <- readRDS(file.path(outdir, "synthetic_data.rds"))
mpf <- readRDS(file.path(outdir, "mpf.rds"))
res <- readRDS(file.path(outdir, "results.rds"))
hours <- st$hours; argvals <- st$argvals; m <- st$m; p <- st$p
domains <- st$domains; groups <- st$groups; group <- st$group

pal_div <- colorRampPalette(c("#2166AC","#4393C3","#F7F7F7","#D6604D","#B2182B"))(64)
pal_seq <- colorRampPalette(c("#FFFFCC","#A1DAB4","#41B6C4","#225EA8","#0C2C84"))(64)
gcol <- c(Control="#4D4D4D", MDD="#1B9E77", `Bipolar I`="#D95F02", `Bipolar II`="#7570B3")
im4  <- function(M, main, col=pal_div, zlim=NULL, ax=TRUE) {
  if (is.null(zlim)) { r <- max(abs(M), na.rm=TRUE); zlim <- c(-r, r) }
  image(1:nrow(M), 1:ncol(M), M, col=col, zlim=zlim, xlab="", ylab="",
        main=main, axes=FALSE); box()
  if (ax) for (k in 1:p) { at <- (k-0.5)*m; axis(1,at,domains[k],tick=FALSE,cex.axis=.7)
                           axis(2,at,domains[k],tick=FALSE,cex.axis=.7,las=1) }
  abline(v=(1:(p-1))*m+.5, h=(1:(p-1))*m+.5, col="grey40", lty=1)
}

## --- Fig 1: example synthetic diurnal trajectories, by domain x group ------
png(file.path(figdir,"fig1_example_trajectories.png"), 1100, 850, res=130)
op <- par(mfrow=c(2,2), mar=c(4,4,3,1))
set.seed(3)
for (v in 1:p) {
  X <- st$dat_list[[v]]
  yl <- if (v==4) "Total Log Activity" else "ordinal level"
  plot(NA, xlim=range(hours), ylim=range(X, na.rm=TRUE),
       xlab="hour of day", ylab=yl, main=domains[v])
  for (g in seq_along(groups)) {
    idx <- sample(which(group==groups[g]), 6)
    for (i in idx) {
      xi <- X[i,]; ok <- !is.na(xi)
      lines(hours[ok], xi[ok], col=adjustcolor(gcol[g],.5), lwd=1)
      points(hours[ok], xi[ok], col=gcol[g], pch=16, cex=.4)
    }
  }
  if (v==1) legend("topright", names(gcol), col=gcol, lwd=2, cex=.7, bty="n")
}
par(op); dev.off()

## --- Fig 2: estimated shared eigenfunctions + FVE --------------------------
png(file.path(figdir,"fig2_eigenfunctions.png"), 1100, 480, res=130)
op <- par(mfrow=c(1,2), mar=c(4,4,3,1))
phi_hat <- mpf$pfpca_results[[1]]$phi
K <- min(3, ncol(phi_hat)); ecol <- c("#1B9E77","#D95F02","#7570B3")
matplot(hours, phi_hat[,1:K], type="l", lty=1, lwd=2.5, col=ecol,
        xlab="hour of day", ylab=expression(phi[l](t)),
        main="Estimated shared eigenfunctions")
abline(h=0, col="grey70", lty=3)
legend("topright", c("PC1: diurnal level","PC2: morning-evening","PC3: circadian"),
       col=ecol, lwd=2.5, cex=.7, bty="n")
fve <- mpf$pfpca_results[[1]]$cumFVE[1:6]
barplot(fve, names.arg=paste0("PC",1:6), col="#41B6C4", border=NA,
        ylab="cumulative FVE (%)", main="Variance explained", ylim=c(0,100))
abline(h=90, lty=2, col="grey30")
par(op); dev.off()

## --- Fig 3: true vs estimated multivariate covariance ---------------------
png(file.path(figdir,"fig3_covariance.png"), 1200, 560, res=130)
op <- par(mfrow=c(1,2), mar=c(3,5,3,1))
im4(res$R_true, "True correlation (4 domains)")
im4(cov2cor(res$Chat), sprintf("M2FPCA estimate (cor=%.2f)", res$vcor))
par(op); dev.off()

## --- Fig 4: cross-domain score correlation (component 1) ------------------
png(file.path(figdir,"fig4_score_correlation.png"), 620, 560, res=130)
op <- par(mar=c(6,6,3,2))
Rs <- res$score_cor_pc1
image(1:p,1:p, Rs[,p:1], col=pal_div, zlim=c(-1,1), axes=FALSE,
      xlab="", ylab="", main="Cross-domain PC1 score correlation")
axis(1, 1:p, domains, las=2, cex.axis=.7, tick=FALSE)
axis(2, 1:p, rev(domains), las=1, cex.axis=.7, tick=FALSE)
for (a in 1:p) for (b in 1:p) text(a, p-b+1, sprintf("%.2f", Rs[a,b]), cex=.8)
box(); par(op); dev.off()

## --- Fig 5: digital biomarkers -- score distributions by diagnosis --------
png(file.path(figdir,"fig5_biomarkers_by_group.png"), 1100, 520, res=130)
op <- par(mfrow=c(1,3), mar=c(7,4,3,1))
sw <- as.data.frame(res$scores_wide); sw$group <- group
show <- c("TLAC_PC2","Energy_PC1","Energy_PC3")
labs <- c("Activity morning-evening contrast (PC2)",
          "Energy diurnal level (PC1)","Energy circadian (PC3)")
for (k in seq_along(show)) {
  boxplot(sw[[show[k]]] ~ sw$group, col=gcol[levels(group)], las=2,
          xlab="", ylab="score", main=labs[k], cex.axis=.8)
}
par(op); dev.off()

## --- Fig 6: prediction -- confusion matrix --------------------------------
png(file.path(figdir,"fig6_prediction.png"), 640, 600, res=130)
op <- par(mar=c(6,6,4,2))
cm <- res$pred$cm; cmn <- cm/rowSums(cm)
image(1:ncol(cm),1:nrow(cm), t(cmn)[,nrow(cm):1], col=pal_seq, zlim=c(0,1),
      axes=FALSE, xlab="predicted", ylab="true",
      main=sprintf("Diagnostic-group prediction\ntest acc = %.0f%%  (baseline %.0f%%)",
                    100*res$pred$acc, 100*res$pred$baseline))
axis(1, 1:ncol(cm), colnames(cm), las=2, cex.axis=.7, tick=FALSE)
axis(2, 1:nrow(cm), rev(rownames(cm)), las=1, cex.axis=.7, tick=FALSE)
for (a in 1:ncol(cm)) for (b in 1:nrow(cm)) text(a, nrow(cm)-b+1, cm[b,a], cex=1)
box(); par(op); dev.off()

cat("Figures written to", figdir, ":\n"); print(list.files(figdir, "\\.png$"))
