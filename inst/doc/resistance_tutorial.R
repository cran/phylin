### R code from vignette source 'resistance_tutorial.Snw'
### Encoding: UTF-8

###################################################
### code chunk number 1: attach
###################################################
library(phylin)
data(simulations)


###################################################
### code chunk number 2: dgen
###################################################
size <- (simul.env$env.sur - min(simul.env$env.sur)) / 10
grid.full <- simul.env[,1:2]
plot(grid.full, cex = size)
points(simul.sample[,1:2], col=simul.sample$lineage, pch=16)


###################################################
### code chunk number 3: grid
###################################################
grid <- grid.full[abs(grid.full$x) < 1.05 & abs(grid.full$y) < 1.05,]


###################################################
### code chunk number 4: samples_tree
###################################################
plot(hclust(simul.gen.dist, method="average"), labels=FALSE)
abline(h = 1.25, col='red', lty=2)


###################################################
### code chunk number 5: geo.dist
###################################################
geo.dist <- dist(simul.sample[,1:2])


###################################################
### code chunk number 6: conductance
###################################################
conductance <- 1/(1+exp(-1*(simul.env$env.sur-5)))
plot(simul.env$env.sur, conductance, cex=0.25, 
     xlab="Environmental surface", ylab = "Conductance")


###################################################
### code chunk number 7: conductanceraster
###################################################
library(gdistance)
conductance <- rasterFromXYZ(data.frame(grid.full, conductance))
plot(conductance)


###################################################
### code chunk number 8: transition
###################################################
tr <- transition(conductance, mean, 8)
tr <- geoCorrection(tr, type="r")


###################################################
### code chunk number 9: resist
###################################################
res.dist <- commuteDistance(tr, as.matrix(simul.sample[,1:2]))


###################################################
### code chunk number 10: geovario
###################################################
gv.geo <- gen.variogram(geo.dist, simul.gen.dist, lag=0.01)
plot(gv.geo)


###################################################
### code chunk number 11: geomodel
###################################################
gv.geo <- gv.model(gv.geo, range=2, nugget=0.3)
plot(gv.geo)


###################################################
### code chunk number 12: resvario
###################################################
gv.res <- gen.variogram(res.dist, simul.gen.dist, lag=120)
plot(gv.res)


###################################################
### code chunk number 13: resmodel
###################################################
gv.res <- gv.model(gv.res, range=19000)
plot(gv.res)


###################################################
### code chunk number 14: geoFUN
###################################################
my.geo.dist <- function (from, to) {
    nf <- nrow(from)
    allcoords <- rbind(from, to)
    dist <- as.matrix(dist(allcoords))
    geo.dist <- dist[1:nf, (nf+1):ncol(dist)]
    return(geo.dist)
}


###################################################
### code chunk number 15: geoFUNtry
###################################################
# The 'x' and 'y' columns of the first 3 samples
sp <- simul.sample[1:3, 1:2]
# The first 6 locations in the grid
grd <- grid.full[1:6,]
# Calculate distances from samples to grid
my.geo.dist(sp, grd)


###################################################
### code chunk number 16: resFUN
###################################################
my.res.dist <- function (from, to, tr) {
    nf <- nrow(from)
    allcoords <- as.matrix(rbind(from, to))
    dist <- as.matrix(commuteDistance(tr, allcoords))
    my.dist <- dist[1:nf, (nf+1):ncol(dist)]
    return(my.dist)
}


###################################################
### code chunk number 17: resFUNtry
###################################################
# Calculate distances from samples to grid
my.res.dist(sp, grd, tr)


###################################################
### code chunk number 18: geoInterpol
###################################################
lin <- as.integer(simul.sample$lineage == 2)
intpl <- krig(lin, simul.sample[,1:2], grid, gv.geo, my.geo.dist, 
              neg.weights=FALSE, verbose=FALSE)
grid.image(intpl, grid)
points(simul.sample[,1:2], pch=lin+1)


###################################################
### code chunk number 19: resInterpol
###################################################
lin <- as.integer(simul.sample$lineage == 2)
intpl <- krig(lin, simul.sample[,1:2], grid, gv.res, my.res.dist, tr=tr, 
              neg.weights=FALSE, verbose=FALSE)
grid.image(intpl, grid)
points(simul.sample[,1:2], pch=lin+1)


###################################################
### code chunk number 20: loopInterpol
###################################################
geo <- matrix(NA, nrow(grid), 4) 
res <- matrix(NA, nrow(grid), 4) 
for (l in 1:4) {
    lin <- as.integer(simul.sample$lineage == l)
    geo[,l] <- krig(lin, simul.sample[,1:2], grid, gv.geo, my.geo.dist, 
                  neg.weights=FALSE, verbose=FALSE)$Z
    res[,l] <- krig(lin, simul.sample[,1:2], grid, gv.res, my.res.dist, tr=tr, 
                  neg.weights=FALSE, verbose=FALSE)$Z
}


###################################################
### code chunk number 21: rasterInterpol
###################################################
geo.raster <- rasterFromXYZ(data.frame(grid, geo))
res.raster <- rasterFromXYZ(data.frame(grid, res))


