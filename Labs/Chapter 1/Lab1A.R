##############################################=
# Spatial Econometrics
# Laboratory: Pollution in Metropolitan Region, Chile
# By: Mauricio Sarrias
##############################################=

## I: Clean directory and load packages ====
rm(list = ls(all = TRUE)) # Clean workspace

#install.packages("maptools", "spdep", "foreign", "RColorBrewer")

library("maptools") # Tools for reading and handling spatial objects
library("spdep")    # Tools for spatial analysis
library("foreign")  # Read stata files
library("RColorBrewer") # Creates nice looking palettes specially for thematic maps

# Change the following:
setwd("~/Dropbox/Mis Clases/Spatial Econometrics/Data/PollutionRM")

## II: Merge shape file with data.frame ====

# shapefile: popular geospatial vector data format that describes
# vector features: points, lines, and polygons.
# Mandatory files:
#     .shp = shape format; the feature geometry itself.
#     .shx = shape index format; a positional index of the feature geogmetry
#     .dbf = attribute format; columnar attributes for each shape in dBase IV format
# Other files:
#     .prj = projection format; the coordinates system and projection information. 


mr_shape <- readShapeSpatial("mr_chile")  # Load shape file
names(mr_shape)                           # Names of variables in dbf
mr_shape@data                             # Show data in dbf

poll_data <- read.dta("Pollution_RM.dta") # Load data is Stata (for example)
names(poll_data)                          # Names of variables

# Remark!: both shapefile and data have the same id variable as first variable

# Combine
poll_sh <- mr_shape                     # Create a new shape file
poll_sh@data <- merge(mr_shape@data,    # object of class Spatial
                      poll_data,        # object of class data.frame        
                      by.x  = "ID", 
                      by.y  = "ID", 
                      all.x = TRUE,     # It will have all rows of mr_shape@data
                      sort  = FALSE)    # sort = F is critical, otherwise the new data-table is re-ordered
summary(poll_sh)

## III: Plot Pollution Variables ====
poll_sh$lco2 <- log(poll_sh$co2)
poll_sh$lnox <- log(poll_sh$nox)
poll_sh$lcov <- log(poll_sh$cov)


spplot(poll_sh, 
       "lco2", 
       at = quantile(poll_sh$lco2, p = c(0, .25, .5, .75, 1), 
                     na.rm = TRUE),
       col.regions = brewer.pal(5, "Blues"), 
       main = expression(Log(C0[2])))
spplot(poll_sh, 
       "lnox", 
       at = quantile(poll_sh$lnox, p = c(0, .25, .5, .75, 1), 
                     na.rm = TRUE),
       col.regions = brewer.pal(5, "Blues"), 
       main = expression(Log(NOX)))
spplot(poll_sh, 
       "lcov", 
       at = quantile(poll_sh$lcov, p = c(0, .25, .5, .75, 1), 
                     na.rm = TRUE),
       col.regions = brewer.pal(5, "Blues"), 
       main = expression(Log(COV)))

# Equal breaks, quantile breaks, jenks breaks
library("classInt")
pal <- brewer.pal(7, "Greens")
brks.qt <- classIntervals(poll_sh$lco2, n = 7, style = "quantile")
brks.jk <- classIntervals(poll_sh$lco2, n = 7, style = "jenks")
brks.eq <- classIntervals(poll_sh$lco2, n = 7, style = "equal")


library("gridExtra") # To combine graphs
p1 <- spplot(poll_sh, "lco2", 
       at = brks.qt$brks,
       col.regions = pal, 
       col = "transparent",
       main = "Quantile breaks")
p2 <- spplot(poll_sh, "lco2", 
       at = brks.jk$brks,
       col.regions = pal, 
       col = "transparent",
       main = "Jenks breaks")
p3 <- spplot(poll_sh, "lco2", 
       at = brks.eq$brks,
       col.regions = pal, 
       col = "transparent",
       main = "Equal breaks")

grid.arrange(p1, p2, p3)

#Cumulative dist plot
par(mfrow = c(1, 3))
plot(brks.qt, pal = pal, main = "Quantile")
plot(brks.jk, pal = pal, main = "Jenks")
plot(brks.eq, pal = pal, main = "Equal")


## IV: W Matrices ====

# Queen
queen.w <- poly2nb(poll_sh, 
                   row.names = poll_sh$NAME, 
                   queen =  TRUE)
summary(queen.w)
attributes(queen.w)
is.symmetric.nb(queen.w)

listw2mat(nb2listw(queen.w, 
                   style = "W")) #See the matrix
rowSums(listw2mat(nb2listw(queen.w,
                           style = "W"))) #Row sums
rowSums(listw2mat(nb2listw(queen.w,
                           style = "B")))

# Rook
rook.w <- poly2nb(poll_sh, row.names = poll_sh$NAME, queen =  FALSE)
rowSums(listw2mat(nb2listw(rook.w)))

# K-neighboords
coords <- coordinates(poll_sh)

k1neigh <- knn2nb(knearneigh(coords, k = 1))
k2neigh <- knn2nb(knearneigh(coords, k = 2))
rowSums(listw2mat(nb2listw(k2neigh)))

k3neigh <- knn2nb(knearneigh(coords, k = 3))
k4neigh <- knn2nb(knearneigh(coords, k = 4))
k5neigh <- knn2nb(knearneigh(coords, k = 5))

# Distance
dist.mat <- as.matrix(dist(coords, method = "euclidean")) # Euclidean Distance
dist.mat.inv <- 1 / dist.mat
diag(dist.mat.inv) <- 0
dist.mat.inv[1:5, 1:5]

dist.mat.inve <- mat2listw(dist.mat.inv, 
                           style = "W", 
                           row.names = poll_sh$NAME) #Normalized
summary(dist.mat.inve)

rowSums(listw2mat(mat2listw(dist.mat.inv, style = "W", row.names = poll_sh$NAME)))
rowSums(listw2mat(mat2listw(dist.mat.inv, style = "B", row.names = poll_sh$NAME)))


# Plot W matrices
par(mfrow = c(2,2))
plot(poll_sh, border = "grey", main = "Queen")
plot(queen.w, coordinates(poll_sh), add =  TRUE, col = "red")
plot(poll_sh, border = "grey", main = "Rook")
plot(rook.w, coordinates(poll_sh), add =  TRUE, col = "red")
plot(poll_sh, border = "grey", main = "1-Neighboors")
plot(k1neigh, coordinates(poll_sh), add =  TRUE, col = "red")
plot(poll_sh, border = "grey", main = "2-Neighboors")
plot(k2neigh, coordinates(poll_sh), add =  TRUE, col = "red")
dev.off()

plot(poll_sh, border = "grey", main = "Inverse Diance")
plot(dist.mat.inve, coordinates(poll_sh), add =  TRUE, col = "red")

plot(poll_sh, border = "grey", main = "Queen")
plot(queen.w, coordinates(poll_sh), add =  TRUE, col = "red")
plot(rook.w, coordinates(poll_sh), add =  TRUE, col = "yellow")


## V: Moran Scatterplot ====
par(mfrow = c(2, 2))
moran.plot(poll_sh$lco2, 
           listw = nb2listw(queen.w),
           xlab = "Log(C02)",
           ylab = expression(WLog(C0[2])),
           labels = as.character(poll_sh$NAME), pch = 19,
           main = "Queen")

moran.plot(poll_sh$lco2, 
           listw = nb2listw(rook.w),
           xlab = "Log(C02)",
           ylab = expression(WLog(C0[2])),
           labels = as.character(poll_sh$NAME), pch = 19,
           main = "Rook")

moran.plot(poll_sh$lco2, 
           listw = nb2listw(k5neigh),
           xlab = "Log(C02)",
           ylab = expression(WLog(C0[2])),
           labels = as.character(poll_sh$NAME), pch = 19,
           main = "5-Neighboors")

moran.plot(poll_sh$lco2, 
           listw = dist.mat.inve,
           xlab = "Log(C02)",
           ylab = expression(WLog(C0[2])),
           labels = as.character(poll_sh$NAME), pch = 19,
           main = "Inverse Distance")
dev.off()


par(mfrow = c(2, 2))
moran.plot(poll_sh$lnox, 
           listw = nb2listw(queen.w),
           xlab = "Log(NOX)",
           ylab = expression(Wlog(NOX)),
           labels = as.character(poll_sh$NAME), pch = 19,
           main = "Queen")

moran.plot(poll_sh$lnox, 
           listw = nb2listw(rook.w),
           xlab = "Log(NOX)",
           ylab = expression(Wlog(NOX)),
           labels = as.character(poll_sh$NAME), pch = 19,
           main = "Rook")

moran.plot(poll_sh$lnox, 
           listw = nb2listw(k5neigh),
           xlab = "Log(NOX)",
           ylab = expression(Wlog(NOX)),
           labels = as.character(poll_sh$NAME), pch = 19,
           main = "5-Neighboors")

moran.plot(poll_sh$lnox, 
           listw = dist.mat.inve,
           xlab = "Log(NOX)",
           ylab = expression(Wlog(NOX)),
           labels = as.character(poll_sh$NAME), pch = 19,
           main = "Inverse Distance")
dev.off()


## VI: Moran Tests ====
moran.test(x = poll_sh$lco2, listw = nb2listw(queen.w), randomisation = FALSE)
moran.test(x = poll_sh$lco2, listw = nb2listw(rook.w), randomisation = FALSE)
moran.test(x = poll_sh$lco2, listw = nb2listw(k2neigh), randomisation = FALSE)
moran.test(x = poll_sh$lco2, listw = nb2listw(k3neigh), randomisation = FALSE)
moran.test(x = poll_sh$lco2, listw = nb2listw(k4neigh), randomisation = FALSE)
moran.test(x = poll_sh$lco2, listw = nb2listw(k5neigh), randomisation = FALSE)
moran.test(x = poll_sh$lco2, listw = dist.mat.inve, randomisation = FALSE)

# Spatial correllograms
plot(sp.correlogram(neighbours = queen.w,
                    var = poll_sh$lco2,
                    order = 4,
                    method = "I",
                    style = "W"))

# VII: LISA ====
# Local moran
local.mi <- localmoran(poll_sh$lco2,
                       listw = nb2listw(queen.w),
                       alternative = "two.sided",
                       p.adjust.method  = "fdr")

poll_sh$lmi   <- local.mi[, 1]
poll_sh$lmi.p <- local.mi[, 5]

poll_sh$lmi.p.sig <- as.factor(ifelse(local.mi[,5] < 0.001, "Sig p < 0.001", 
                              ifelse(local.mi[, 5] < 0.05, "Sig p< 0.05", "NS")))

spplot(poll_sh, "lmi", at = summary(poll_sh$lmi), 
       col.regions = brewer.pal(5, "RdBu"), 
       main = "Local Moran's I")

spplot(poll_sh, 
       "lmi.p.sig", 
       col.regions = c("white", "Red", "Yellow"))


