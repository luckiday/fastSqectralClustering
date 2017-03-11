library(mlbench)

set.seed(111)
obj <- mlbench.spirals(100,1,0.025)
my.data <-  4 * obj$x
plot(my.data)