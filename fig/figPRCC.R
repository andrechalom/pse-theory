PRCCg <- function(N, res, r, k, Time, x) {
resr <- rank(res); rr <- rank(r); kr <- rank(k); Timer <- rank(Time); xr <- rank(x)
PRCC<-matrix(0,1,4)
colnames(PRCC)<-c("r","k","Time","x")
#PRCC between res & r
PRCC[1]<-
cor(
    residuals(lm(rr~kr+Timer+xr)),
    residuals(lm(resr~kr+Timer+xr))
)
#PRCC between res & k
PRCC[2] <-
cor(
    residuals(lm(kr~rr+Timer+xr)),
    residuals(lm(resr~rr+Timer+xr))
)
#PRCC between res & Time
PRCC[3] <-
cor(
    residuals(lm(Timer~kr+rr+xr)),
    residuals(lm(resr~kr+rr+xr))
)
#PRCC between res & x
PRCC[4] <-
cor(
    residuals(lm(xr~kr+Timer+rr)),
    residuals(lm(resr~kr+Timer+rr))
)
barplot(PRCC, ylim=c(-1,1), main=paste("N =",N))
#Desenha as linhas a partir das quais o PRCC eh significativo
tval <- function (prcc, df) {
	return(prcc*(sqrt(df/(1-prcc*prcc))) - qt(0.975,df) ) 
}
psig <- uniroot(tval, c(0,1), N-5)$root ### df deve ser alterado NA MAO
abline(h=0);abline(h=psig,lty=2);abline(h=-psig,lty=2)
}
# OBS: when running from inside "./fig", use "../Rdata"
par(mfrow=c(2,2))
load('Rdata/LHS15.Rdata'); N<- 15;
PRCCg(N, res, r, k, Time, x)
load('Rdata/LHS45.Rdata'); N<- 45;
PRCCg(N, res, r, k, Time, x)
load('Rdata/LHS150.Rdata'); N<- 150;
PRCCg(N, res, r, k, Time, x)
load('Rdata/LHS450.Rdata'); N<- 450;
PRCCg(N, res, r, k, Time, x)
