source("pse.R")
library(sensitivity)
# Hipercubo gerado de tamanho 650
# Analises realizadas para modelo independente de densidade
load("Independent.Rdata")
# Analises realizadas para modelo dependente de densidade
load("Dependent.Rdata")

# Plots
par(mfrow=c(2,1), cex=0.8, mar=c(4,2,1.5,1))
mu = apply(imorr$ee, 2, function(x) mean(abs(x)))
plot(imorr, log='xy', xlim=c(min(mu), 1.2*max(mu)))
mtext('(a)', at=min(mu))
mu = apply(dmorr$ee, 2, function(x) mean(abs(x)))
plot(dmorr, log='xy', xlim=c(min(mu), 1.2*max(mu)))
mtext('(b)', at=min(mu))

par(mfrow=c(1,2))
plot(ecdf(ires), do.points=FALSE, main="ECDF, independent")
plot(ecdf(dres), do.points=FALSE, main="ECDF, dependent", log="x", xlim=c(1, 1e6))

corPlot(ivars, ires)
corPlot(ivars[,c(1,3,5,7)], ires)
corPlot(dvars, dres, log="y")

par(mfrow=c(2,1), cex=1, cex.axis=0.8)
plot(iprcc); abline(h=0, lty=3)
plot(dprcc); abline(h=0, lty=3)
# p77 e p11 pulam!

par(mfrow=c(2,1), cex.axis=0.8)
plot(ifast); title(main="eFAST")
plot(dfast); title(main="eFAST")
