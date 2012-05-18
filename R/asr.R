# Deve ser chamado apos o modelo estar constuido. Pode ser chamado mais de uma vez
# source("modelrun.R")
# EXTENDENDO o cubo
newr <- LHSextend(r, 0.25, 2)
newk <- LHSextend(k, 10, 50)
newx <- LHSextend(x, 1, 10)
newTime <- LHSextend(Time, 100, 40)
newN <- 2*N
newvars <- LHScorcorr (cbind(newr,newk, newTime, newx))
newk[1:newN] <- newvars[,2]
newTime[1:newN] <- newvars[,3]
newx[1:newN] <- newvars[,4]

# Vemos que as correlacoes continuam baixas:
print(cor(r,k)^2); print(cor(newr, newk)^2); print(cor(c(r,newr),c(k,newk))^2)
# Resultado do Teo.1
print(cor(c(r,newr),c(k,newk))^2 < (cor(r,k)^2 + cor(newr, newk)^2) )
# Propriedades de media e momentos se conservam:
mean(r); mean(newr); mean(c(newr,r))
sd(r); sd(newr); sd(c(newr,r))
newres <- mapply(modelRun, newx, newr, newk, newTime)
# Triplicamos o N original
N <- 3*N
res <- c(res,newres)
r <- c(r, newr); attributes(r) <- attributes(newr)
k <- c(k, newk); attributes(k) <- attributes(newk)
Time <- c(Time, newTime); attributes(Time) <- attributes(newTime)
x <- c(x, newx); attributes(x) <- attributes(newx)

