# Function from the github of Imad Ali
# https://github.com/imadmali/bball-hmm/tree/master/data
# Associated with the vignette: https://mc-stan.org/users/documentation/case-studies/bball-hmm.html

plt_defense_example <- function(defense_data, ...) {
  list2env(defense_data, environment())
  plot(h[1], h[2], xlim = c(-2,2), ylim = c(-2,2), pch = 19,
       xlab = "x coordinate", ylab = "y coordinate", ...)
  points(b)
  points(o[1,,], pch = 3)
  points(o[2,,], pch = 3)
  points(o[3,,], pch = 3)
  points(o[4,,], pch = 3)
  points(o[5,,], pch = 3)
  text(d[1,1], d[1,2], "defender start", pos = 2, cex = 0.8)
  text(d[N,1], d[N,2], "defender end", pos = 4, cex = 0.8)
  text(h[1], h[2], "hoop", pos = 1, cex = 0.8)
  text(b[1,1], b[1,2], "ball", pos = 4, cex = 0.8)
  text(o[1,1,1], o[1,1,2], "o1", pos = 4, cex = 0.8)
  text(o[2,1,1], o[2,1,2], "o2", pos = 4, cex = 0.8)
  text(o[3,1,1], o[3,1,2], "o3", pos = 4, cex = 0.8)
  text(o[4,1,1], o[4,1,2], "o4", pos = 2, cex = 0.8)
  text(o[5,1,1], o[5,1,2], "o5", pos = 4, cex = 0.8) 
}