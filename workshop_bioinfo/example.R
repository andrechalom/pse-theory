library(pse)
l <- LHS(model=function(x) {x[,1]+x[,2]*x[,3]+x[,4]*x[,5]/(1+x[,6])+rnorm(1,1,0.3)}, factors=7, N=200, nboot=50, res.names="Model results", repetitions=5)
plotprcc(l)
plotecdf(l)
plotscatter(l)
plotcv(l)

