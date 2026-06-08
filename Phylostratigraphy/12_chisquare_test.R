# =========================================================
# BRAIN CHI-SQUARE TEST
# =========================================================

brain_matrix <- matrix(
  c(
    9388, 823,
    1262, 121
  ),
  
  nrow = 2,
  byrow = TRUE
)

colnames(brain_matrix) <- c("Non_DEG", "DEG")
rownames(brain_matrix) <- c("Old", "New")

brain_matrix

chisq.test(brain_matrix)

# =========================================================
# FAT BODY CHI-SQUARE TEST
# =========================================================

fatbody_matrix <- matrix(
  c(
    8820, 1391,
    1224, 159
  ),
  
  nrow = 2,
  byrow = TRUE
)

colnames(fatbody_matrix) <- c("Non_DEG", "DEG")
rownames(fatbody_matrix) <- c("Old", "New")

fatbody_matrix

chisq.test(fatbody_matrix)

