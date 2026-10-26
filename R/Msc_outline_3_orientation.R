######################################################################
# The lines are normalising the orientation of the wing outlines
######################################################################

cat("Normalising outline orientation...\n")

######################################################################
# The block is rotating each outline so its longest axis lies flat
######################################################################

orient_principal_axis <- function(coords) {
cx <- mean(coords[,1]); cy <- mean(coords[,2])
centered <- cbind(coords[,1]-cx, coords[,2]-cy)
ev <- eigen(cov(centered))$vectors          # columns = eigenvectors
rotated <- centered %*% ev                   # project onto principal axes
######################################################################
# The lines are eigenning() returns axes with PC1 first; rotated[,1] is now the long axis (x)
######################################################################
rotated
}

######################################################################
# The block is making every wing point the same way
######################################################################
orient_tip_right <- function(coords) {
x <- coords[,1]
xr <- diff(range(x))
left_mask  <- x < (min(x) + 0.15*xr)
right_mask <- x > (max(x) - 0.15*xr)
width_left  <- diff(range(coords[left_mask,  2]))
width_right <- diff(range(coords[right_mask, 2]))
######################################################################################################################
# The block is flipping x when the right end is wider than the left, since that means the tip is on the left
######################################################################################################################
if (width_right > width_left) coords[,1] <- -coords[,1]
coords
}

######################################################################################
# The lines are flipping outlines so the same side always faces up 
# Convention check so wings don't appear upside-down relative to each other
######################################################################################
orient_dorsal_up <- function(coords) {
######################################################################################################
# The block is signedding area; if negative (clockwise), the y-convention is flipped -> mirror y
######################################################################################################
n <- nrow(coords)
a <- sum(coords[,1]*coords[c(2:n,1),2] - coords[c(2:n,1),1]*coords[,2]) / 2
if (a < 0) coords[,2] <- -coords[,2]
coords
}

######################################################################################################################################
# The block is smoothing the outline very gently Moving average over the closed contour, small window so shape is preserved
######################################################################################################################################
smooth_contour <- function(coords, window = 3) {
n <- nrow(coords)
if (window < 2) return(coords)
sm <- coords
half <- (window - 1) %/% 2
for (k in 1:n) {
idx <- ((k - half - 1):(k + half - 1)) %% n + 1
sm[k,1] <- mean(coords[idx,1]); sm[k,2] <- mean(coords[idx,2])
}
sm
}

######################################################################
# The block is running the full orientation steps for one outline
######################################################################
orient_outline <- function(coords) {
coords <- orient_principal_axis(coords)
coords <- orient_tip_right(coords)
coords <- orient_dorsal_up(coords)
coords <- smooth_contour(coords, window = 3)   # mild; set window=1 to disable
coords
}

######################################################################
# The block is applying the orientation steps to every saved outline
######################################################################
n_done <- 0
for (nm in names(outlines_list)) {
cm <- as.matrix(outlines_list[[nm]][, c("x","y")])
oriented <- orient_outline(cm)
outlines_list[[nm]] <- data.frame(x = oriented[,1], y = oriented[,2])
n_done <- n_done + 1
}

cat(paste("✓ Oriented + lightly smoothed", n_done, "outlines (tip right, dorsal up)\n\n"))
