
# Calculo do numero de series para var. binarias
NumSeriesDisc <-function(y){
# y deve ser um vetor de 0s e 1s
n1=sum(y)
n2=length(y)-n1
series=1
for (i in 2:length(y)){
if(y[i]!=y[i-1]){series<-series+1}}
r=series
print(c("n1", "n2", "nº de séries"))
return(c(n1,n2,r))
}

# Calculo do numero de series para var. continuas (usando mediana)
NumSeriesCont <-function(y){
grupo=c()
for(i in 1:length(y)){
if (y[i]>median(y)) grupo[i]<-'+'  
if (y[i]<median(y)) grupo[i]<-'-' 
if (y[i]==median(y)) grupo[i]<-0
}
grupo=grupo[grupo!=0]
n1<-0
n2<-0
for(i in 1:length(grupo)){
if(grupo[i]=='+'){n1<-n1+1}
if(grupo[i]=='-'){n2<-n2+1}}
series<-1
for(i in 2:length(grupo)){
if(grupo[i]!=grupo[i-1]){series<-series+1}
}
r=series
print(grupo)
print(c("n1", "n2", "numero de series"))
return(c(n1,n2,r))
}

densRuns <-function(n1,n2,r){
n=n1+n2
# r: nº de séries ou corridas
if (r/2>as.integer(r/2)) fr=(choose(n1-1,(r-1)/2)*choose(n2-1,(r-3)/2)+choose(n1-1,(r-3)/2)*choose(n2-1,(r-1)/2))/choose(n,n1)
else fr=2*choose(n1-1,r/2-1)*choose(n2-1,r/2-1)/choose(n,n1)
return(fr)
}

cdfRuns <- function(n1,n2,qe){
#qe >=2
p=0
for (r in 2:qe) p = p+densRuns(n1,n2,r)
return(p)
}

qRuns <- function(n1,n2,p){
r=2
p1=densRuns(n1,n2,r)
if (p1>=p) qR=r else{
while (p1 < p) {
r=r+1
p1 = p1 + densRuns(n1,n2,r)
}
if (p<0.5) qR=r-1
if (p>=0.5)qR=r+1}
return(qR)
}


# Funcao distribuicao empirica e teorica 
Fn=ecdf(x)
plot(Fn)
xsort=sort(x) 
lines(xsort,pDIST(xsort,thetas))

# Funcao Quantil empirica e teorica
u=seq(0,1,by=1/n)         
Qn=quantile(x,u)
plot(u,qDIST(u,theta),type="l")
points(u,Qn,pch=19)
for (i in 1:n){  
lines(c(u[i],u[i+1]),c(Qn[i],Qn[i]))
}

# Para teste de aderencia var. discretas
# Exemplo
tab=matrix(0,nrow=2,ncol=8)
tab[1,]=c(29,19,18,25,17,10,15,11) # valores observados
tab[2,]=rep(18,8)  # valores esperados
estat<-sum((tab[1,]-tab[2,])^2/tab[2,])


# Teste de ponto de mudanca 
# var. binarias
TestPMDisc<-function(x){
m=sum(x)
N=length(x)
n=N-m
Sj=cbind(rep(0,N-1))
for (j in 1:(N-1)) Sj[j]=sum(x[1:j])
Sj=cbind(as.numeric(Sj))
jj=cbind(as.numeric(1:(N-1)))
Dt=N/(n*m)*abs(Sj-jj*m/N)
Dcal=max(Dt)
pm=which.max(Dt)
print(c("m","n","Estat. Dmn","ponto mudança"))
return(c(m,n,Dcal,pm))
}

Dmn_calc=TestPMDisc(x)
# Use Tabela G L_II do Apendice Siegel (m e n < 25)
# Se m ou n > 25, use Tabela L_III do Apendice

# var. continuas
TestePMCont <-function(x){
N=length(x)
posto=rank(x)
Wj=cbind(rep(0,N))
for (j in 1:N) Wj[j]=sum(posto[1:j])
cbind(posto,Wj)
Wj1=Wj[1:N-1]
jj=cbind(as.numeric(1:(N-1)))
Kt=abs(2*Wj1-jj*(N+1))
Kcal=max(Kt)
m=which.max(Kt)
Kmn=Kt[m]
n=N-m
print(c("m","n","Kmn"))
return(c(m,n,Kmn))
}


# Teste exato: Use Tabela J (Siegel)
# Para usar a tabela J: calcule a soma dos postos até a posição m e compare com a tabela
sum(rank(x[1:m]))


#### Amostra grande
TestePMContaprox <-function(x){
pmx=TestePMCont(x)
m=pmx[1]
N=length(x)
n=N-m
W=rank(x)
Wm=sum(W[1:m])
EW=m*(N+1)/2
VarW=m*n*(N+1)/12
Wz=(Wm-EW-1/2*sign(Wm-EW))/sqrt(VarW)
p_valor=2*pnorm(-abs(Wz))
print(c("pto mudança","Wz","p_valor"))
return(c(m,Wz,p_valor))
}


# Teste das Corridas para grande amostras (aproxima??o pela normal)
RunsTesteAprox <- function(y,n1,n2,r){
# hc: ccorre??o de continuidade
n=n1+n2
media=2*n1*n2/n+1
varR=2*n1*n2*(2*n1*n2-n)/(n^2*(n-1))
if (r>media) hc=-1/2 else hc=1/2
rz=(r+hc-media)/sqrt(varR)
p_valor=2*pnorm(-abs(rz))
print("p_valor")
return(p_valor)
}


dagostino.pearson.test <- function(x) {
    # from Zar (1999), implemented by Doug Scofield,  scofield at bio.indiana.edu
    DNAME <- deparse(substitute(x))
    n <- length(x)
    x2 <- x * x
    x3 <- x * x2
    x4 <- x * x3
    # compute Z_g1
    k3 <- ((n*sum(x3)) - (3*sum(x)*sum(x2)) + (2*(sum(x)^3)/n)) /
((n-1)*(n-2))
    g1 <- k3 / sqrt(var(x)^3)
    sqrtb1 <- ((n - 2)*g1) / sqrt(n*(n - 1))
    A <- sqrtb1 * sqrt(((n + 1)*(n + 3)) / (6*(n - 2)))
    B <- (3*(n*n + 27*n - 70)*(n+1)*(n+3)) / ((n-2)*(n+5)*(n+7)*(n+9))
    C <- sqrt(2*(B - 1)) - 1
    D <- sqrt(C)
    E <- 1 / sqrt(log(D))
    F <- A / sqrt(2/(C - 1))
    Zg1 <- E * log(F + sqrt(F*F + 1))
    # compute Z_g2
    G <- (24*n*(n-2)*(n-3)) / (((n+1)^2)*(n+3)*(n+5))
    k4 <- (((n*n*n + n*n)*sum(x4)) - (4*(n*n + n)*sum(x3)*sum(x)) -
(3*(n*n - n)*sum(x2)^2) + (12*n*sum(x2)*sum(x)^2) - (6*sum(x)^4)) /
(n*(n-1)*(n-2)*(n-3))
    g2 <- k4 / var(x)^2
    H <- ((n-2)*(n-3)*abs(g2)) / ((n+1)*(n-1)*sqrt(G))
    J <- ((6*(n*n - 5*n + 2)) / ((n+7)*(n+9))) * sqrt((6*(n+3)*(n+5)) /
(n*(n-2)*(n-3)))
    K <- 6 + (8/J)*(2/J + sqrt(1 + 4/(J*J)))
    L <- (1 - 2/K) / (1 + H*sqrt(2/(K-4)))
    Zg2 <- (1 - 2/(9*K) - (L^(1/3))) / (sqrt(2/(9*K)))
    K2 <- Zg1*Zg1 + Zg2*Zg2
    pk2 <- pchisq(K2, 2, lower.tail=FALSE)
    RVAL <- list(statistic = c(K2 = K2), p.value = pk2, method =
"D'Agostino-Pearson normality test\n\nK2 is distributed as Chi-squared
with df=2", alternative = "distribution is not normal", data.name =
DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
}


