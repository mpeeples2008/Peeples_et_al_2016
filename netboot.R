# create random table of nsim replicates
nsim <- 100
net.data <- AD1250cer

sim.mat <- function(x) {
  names <- row.names(x)
  x <- na.omit(x)
  x <- prop.table(as.matrix(x),1)*100
  rd <- dim(x)[1]
  results <- matrix(0,rd,rd)
  for (s1 in 1:rd) {
    for (s2 in 1:rd) {
      x1Temp <- as.numeric(x[s1, ])
      x2Temp <- as.numeric(x[s2, ])
      results[s1,s2] <- 200 - (sum(abs(x1Temp - x2Temp)))}}
  row.names(results) <- names
  colnames(results) <- names
  results <- results/200
  results <- round(results,3)
  return(results)}

x <- sim.mat(net.data)

net.x <- event2dichot(x,method='absolute',thresh=0.749)
cut <- 20

set.eg <- function(x) {
eg <- as.matrix(evcent(x,rescale=T))
eg <- as.matrix((eg^2))
eg <- as.matrix(eg * (1/sum(eg)))
eg <- eg*nrow(eg)
eg <- sqrt(eg)
return(eg)}
x.eg <- set.eg(x)



data.rowsums <- apply(net.data,1,sum)
for (i in 1:length(data.rowsums)) {
if (data.rowsums[i] > cut) {data.rowsums[i] <- cut}}
data.rowsums <- round(data.rowsums,0)

data.sim <- rmultinom(nsim,data.rowsums[1],prob=net.data[1,])

for (i in 2:nrow(net.data)) {
data.sim <- cbind(data.sim,rmultinom(nsim,data.rowsums[i],prob=net.data[i,]))}

data.sim <- data.sim
data.sim <- t(data.sim)

# resort table by original rows
data.sim2 <- matrix(rep(0,nsim*nrow(net.data)*ncol(net.data)),nrow=nsim*nrow(net.data))
for (k in 1:nsim) {
	for (i in 1:nrow(net.data)) {
		data.sim2[(k-1)*nrow(net.data)+i,] <- data.sim[k+(i-1)*nsim,]}}

# label rows and columns
colnames(data.sim2) <- colnames(net.data)
sites <- rownames(net.data)
site.sim <- sites
while (length(site.sim) < length(sites)*nsim) {
site.sim <- c(site.sim,sites)}
rownames(data.sim2) <- site.sim

rand.dg <- NULL
rand.eg <- NULL
graph.cor <- NULL
#plot(x,x,type='n',xlim=c(0,1),ylim=c(0,1),xlab='Subsample Eigenvector Centrality',ylab='Original Eigenvector Centrality')
for (i in 1:nsim) {
rand.table <- data.sim2[(nrow(net.data)*(i-1)+1):(nrow(net.data)*i),]
rand.sim <- 1-diss.mat(rand.table)
temp1 <- met.degree(rand.sim)
rand.dg[i] <- cor(temp1,met.degree(x))
temp2 <- set.eg(rand.sim)
#points(temp2[order(temp2)]/max(temp2),x.eg[order(x.eg)]/max(x.eg),type='l',col='#0000ff25')
rand.eg[i] <- cor(temp2,set.eg(x))
z <- as.vector(rand.sim)
z2 <- as.vector(x)
graph.cor[i] <- cor(z,z2)
}

mean(rand.dg)
mean(rand.eg)
mean(graph.cor)

