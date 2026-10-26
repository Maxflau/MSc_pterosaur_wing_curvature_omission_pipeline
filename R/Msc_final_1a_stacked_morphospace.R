#######################################################################################
# The lines are coding building the stacked stratigraphical biomechanical morphospace
#######################################################################################

##########################################################################################
# The block is lining distinguishing Pterodactyliform and non-pterodactyliform data point
##########################################################################################
if (!exists("PTERODACT")) {
PTERODACT <- c("Azhdarchoidea", "Euctenochasmatia", "Dsungaripteroidea","Ornithocheiromorpha", "Pteranodontia", "Basal pterodactyliform")
}
##############################################################################
# The block is setting the shared layout options for every stacked plot
# PANEL_ASPECT for the plane height / width in plot units,
# SHEAR for the horizontal shift from front to back edge of a plane,
# PLANE_GAP for the vertical spacing between planes
##############################################################################
GROUP_VAR    <- "clade"
PANEL_ASPECT <- 0.45   
SHEAR        <- 0.34 
PLANE_GAP    <- 0.62  
HULL_ALPHA   <- 0.38
PT_CEX       <- 0.50
grp_all <- sort(unique(as.character(temporal_data[[GROUP_VAR]])))
grp_all <- grp_all[!is.na(grp_all) & grp_all != ""]
PAL <- if (length(grp_all) <= 12) {
brewer.pal(max(3, min(12, length(grp_all))), "Set3")[seq_along(grp_all)]
} else colorRampPalette(brewer.pal(12, "Set3"))(length(grp_all))
names(PAL) <- grp_all
cat(sprintf("Groups: %d\n", length(grp_all)))

########################################################################################################
# The lines are drawing one stacked morphospace plot, with one plane correspond to one specific time bin
########################################################################################################
stacked_panels <- function(d, xv, yv, xlab, ylab, main_title) {
d <- d[is.finite(d[[xv]]) & is.finite(d[[yv]]) & !is.na(d$Time_Bin), ]
##################################################################################################
# The block is sorting out the order of the time bins from oldest to youngest
##################################################################################################
bin_age <- tapply(d$Midpoint, d$Time_Bin, mean, na.rm = TRUE)
bin_age <- bin_age[!is.na(bin_age)]
bins <- names(sort(bin_age, decreasing = TRUE))
nb <- length(bins)
if (nb < 2) { cat("SKIPPED", main_title, "- fewer than two bins\n"); return(invisible(NULL)) }
cat(sprintf("\n  %s\n  Stacking order (bottom to top):\n", main_title))
for (b in bins) cat(sprintf("    %-24s ~%.0f Ma  (n = %d)\n",b, bin_age[[b]], sum(d$Time_Bin == b)))
##################################################################################################
# The block is scaling  the whole dataset to make sure every plane shares the same axes
##################################################################################################
xr <- range(d[[xv]]); yr <- range(d[[yv]])
xpad <- diff(xr) * 0.10; ypad <- diff(yr) * 0.10
xr <- xr + c(-xpad, xpad); yr <- yr + c(-ypad, ypad)
ux <- (d[[xv]] - xr[1]) / diff(xr)                  
uy <- (d[[yv]] - yr[1]) / diff(yr) * PANEL_ASPECT  
total_h <- PANEL_ASPECT + PLANE_GAP * (nb - 1)
plot(NA, xlim = c(-0.16, 1 + SHEAR + 0.24), ylim = c(-0.14, total_h + 0.08), axes = FALSE, xlab = "", ylab = "", asp = 1)
title(main_title, line = 0.4, cex.main = 1.2, adj = 0)
for (k in seq_along(bins)) {
base <- PLANE_GAP * (k - 1)
tx <- function(x, y) x + (y / PANEL_ASPECT) * SHEAR
ty <- function(y) base + y
polygon(c(tx(0, 0), tx(1, 0), tx(1, PANEL_ASPECT), tx(0, PANEL_ASPECT)), c(ty(0), ty(0), ty(PANEL_ASPECT), ty(PANEL_ASPECT)),
border = "grey25", lwd = 1.1,col = adjustcolor("grey97", alpha.f = 0.85))
sel <- d$Time_Bin == bins[k]
px <- tx(ux[sel], uy[sel]); py <- ty(uy[sel])
gv <- as.character(d[[GROUP_VAR]][sel])
##################################################################################################
# The block is assigning a colour to a hull around each group when enough points are present
##################################################################################################
for (g in unique(gv)) {
if (is.na(g) || g == "") next
i <- which(gv == g)
if (length(i) < 3) next
h <- chull(px[i], py[i])
polygon(px[i][h], py[i][h], col = adjustcolor(PAL[[g]], alpha.f = HULL_ALPHA), border = PAL[[g]], lwd = 1.3)
}
cols <- PAL[gv]; cols[is.na(cols)] <- "grey70"
pts_pch <- ifelse(gv %in% PTERODACT, 22, 21)   
points(px, py, pch = pts_pch, bg = cols, col = "grey15", cex = PT_CEX, lwd = 0.45)
text(-0.035, ty(PANEL_ASPECT / 2),
labels = sprintf("%s (~%.0f Ma)", bins[k], bin_age[[bins[k]]]), srt = 90, cex = 0.68, font = 2, xpd = TRUE)
text(tx(1, PANEL_ASPECT) - 0.04, ty(PANEL_ASPECT) - 0.035, labels = sprintf("n=%d", sum(sel)), cex = 0.62, col = "grey40")
}
##################################################################################################
# The block is assigning to all planes  the same axes
##################################################################################################
at_x <- pretty(xr, 5); ax <- (at_x - xr[1]) / diff(xr)
ok <- ax >= 0 & ax <= 1
segments(ax[ok], 0, ax[ok], -0.022, xpd = TRUE)
text(ax[ok], -0.055, labels = at_x[ok], cex = 0.62, xpd = TRUE)
text(0.5, -0.105, labels = xlab, cex = 0.82, xpd = TRUE)
##################################################################################################
# The block is preparing the  right-hand axis for the y variable
##################################################################################################
tx0 <- function(x, y) x + (y / PANEL_ASPECT) * SHEAR
ty0 <- function(y) y
at_y <- pretty(yr, 4); ay <- (at_y - yr[1]) / diff(yr) * PANEL_ASPECT
ok <- ay >= 0 & ay <= PANEL_ASPECT
##################################################################################################
# The block is adding the y-axis tick marks and label on the slanted edge
##################################################################################################
ex <- tx0(1, ay[ok]); ey <- ty0(ay[ok])
segments(ex, ey, ex + 0.03, ey, xpd = TRUE)
text(ex + 0.06, ey, labels = at_y[ok], cex = 0.62, xpd = TRUE)
text(tx0(1, PANEL_ASPECT / 2) + 0.13, PANEL_ASPECT / 2,
labels = ylab, srt = 90, cex = 0.82, xpd = TRUE)
invisible(TRUE)
}
