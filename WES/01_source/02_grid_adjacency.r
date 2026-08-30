################################################################################
# @Project - WES Hair Follicle
# @Description - Pure helpers for grid distance and adjacency classification
################################################################################

calculate_distance <- function(x1, y1, x2, y2) pmax(abs(x1 - x2), abs(y1 - y2))

adjacency_class <- function(d, same_mouse) {
  if_else(same_mouse & d <= 1, "Adjacent", "Distal")
}
