library(plot3D)

df <- read.table("5Y_douglas_Anew.txt", header = TRUE)
snum <- 20
vnum <- 4
rnum <- 10
bnum <- 5

bb <- bnum
i <- 10
j <- 2
k <- 7

df_s <- df[1:(snum+1), 4]
df_v <- df[seq(1,((snum+1)*(vnum+1)),(snum+1)), 3]
df_r <- df[seq(1,((snum+1)*(vnum+1)*(rnum+1)),((snum+1)*(vnum+1))),2]
df_b <- df[seq(1,((snum+1)*(vnum+1)*(rnum+1)*(bnum+1)),((snum+1)*(vnum+1)*(rnum+1))),1]

df_z <- c()
for (bb in 0:bnum) {
  df_z <- cbind(df_z, df[(bb*(rnum+1)*(vnum+1)*(snum+1) + i*(vnum+1)*(snum+1) + j*(snum+1) + 1):(bb*(rnum+1)*(vnum+1)*(snum+1) + i*(vnum+1)*(snum+1) + j*(snum+1) + snum + 1), 5])
}

persp3D(x=df_s[1:10], y=df_b, z=df_z[1:10,], theta=-50, phi=10, expand=0.75, ticktype="detailed", xlab="S", ylab="B", zlab="u value", main=paste("r = ", df_r[i+1], ", v = ", df_v[j+1], split = ""),axes=TRUE)
