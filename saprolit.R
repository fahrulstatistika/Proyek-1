SAPROLIT
## 1.data_preparation ##

##plot data##
plot(saprolit$X,saprolit$Y, xlab = "Easting", ylab = "Northing", pch = 19)
text(saprolit$X, saprolit$Y, row.names(saprolit), cex=0.6, pos=4 , col="red")
      
##create the histogram of variable nilai##
summary(saprolit$Ni)
sd(saprolit$Ni)
var(saprolit$Ni)
quantile(saprolit$Ni)
hist(saprolit$Ni, xlab = "Ni", breaks = "Freedman-Diaconis", main = "Histogram Ni", col = "gray87")
qqPlot(saprolit$Ni, ylab = "kadar Ni")
library("car")
scatterplot(collar_vale$Z ~ collar_vale$X, data = collar_vale, xlab = "Easting", ylab = "Elevation")
scatterplot(collar_vale$Z ~ collar_vale$Y, data = collar_vale, xlab = "Northing", ylab = "Elevation")
model1 = lm(collar_vale$Z ~ collar_vale$X, data = collar_vale)
summary(model1)
model2 = lm(collar_vale$Z ~ collar_vale$Y, data = collar_vale)
summary(model2)
#regresion
plot(collar_vale$X, collar_vale$Z, xlab = "Easting", ylab = "Elevation")
abline(lm(collar_vale$Z ~ collar_vale$X, data = collar_vale), col = "red")
plot(collar_vale$Y, collar_vale$Z, xlab = "Northing", ylab = "Elevation")
abline(lm(collar_vale$Z ~ collar_vale$Y, data = collar_vale), col = "red")
library("scatterplot3d")
scatterplot3d(collar_vale$X, collar_vale$Y, collar_vale$Z, pch = 16, type="h", xlab = "Easting", 
              ylab = "Northing", zlab = "Elevation")

qqnorm(saprolit$Ni, pch = 1)
qqline(saprolit$Ni, col = "steelblue", lwd = 2)


persp(collar_vale$X, collar_vale$Y, collar_vale$Z, zlab = "Elevation", theta = 30, phi = 15,
      col = "springgreen", shade = 0.5)

Points_new <- data.frame(Points.sp$Long, Points.sp$Lat, Points.sp$Elevation)
Points_new <- as.matrix(Points_new)
p <- plot_ly(z = ~Points_new) %>% add_surface()
p

plot_ly(collar_vale$X, collar_vale$Z)
z <- runif(50,0,1) 
z
y <- runif(50,1,2)
x <- runif(50,3,6)
plot_ly(x = x, y = y, z = z, type = 'mesh3d')

#Creating Plot
histogram(~ Ni, data = saprolit,
          xlab = "Kadar Ni", 
          col = "slategray4",
          ylab = "Frequency",
          breaks = 40,
          type = "density", 
          panel = function(x, ...) {
          panel.histogram(x, ...)
      xn <- seq(min(x), max(x), length.out = 100)
      yn <- dnorm(xn, mean(x), sd(x))
          panel.lines(xn, yn, col = "black")
          })


a = density(saprolit$Ni)
hist(saprolit$Ni, breaks = 40, xlab = "Ni", freq = FALSE, ylab = "Frequency", main = "Histogram saprolit Ni", col = "gray47" )
lines(a)

#spatial interpolation#
coordinates(saprolit) <- ~X+Y
v.cloud = variogram(Ni ~ 1, data = saprolit, cloud= T)
head(v.cloud)
plot(v.cloud, xlab = "Separation distance (m)", 
     col = "black")
v = variogram(Ni ~ 1, data = saprolit, cloud = F)
v
plot(v, main = "Variogram - default", xlab = "Separation distance (m)")
v.map = variogram(Ni ~ 1, saprolit, map = TRUE, cutoff=200, width=50)
plot(v.map, col.regions = bpy.colors(64), main="Variogram Map", xlab="x", ylab="y")

#create variogram model#
vg_model = vgm(psill = 0.30, model = "Sph", range = 170, nugget = 0.18)
vg_model = fit.variogram(exp_var, model = vg_model)

## 2. Experimental Variogram ##
gs = gstat(id = "Ni", formula = Ni~1, locations = ~X+Y, data = saprolit)

##calculate the experimental variogram##
exp_var = variogram(gs)

##plot the experimental variogram##
plot(exp_var)

## 3. Variogram Modeling ##
#create variogram model#
vg_model4 = vgm(psill = 0.13, model = "Sph", range = 64, nugget = 0.18)
vg_model5 = vgm(psill = 0.13, model = "Exp", range = 64, nugget = 0.18)
vg_model6 = vgm(psill = 0.13, model = "Gau", range = 64, nugget = 0.18)

vg_model14 = fit.variogram(exp_var, model = vg_model4)
vg_model15 = fit.variogram(exp_var, model = vg_model5)
vg_model16 = fit.variogram(exp_var, model = vg_model6)

a = plot(exp_var, model = vg_model14, cloud=T, col = "black", pch = 19, type = "b")
b = plot(exp_var, model = vg_model15, cloud=T, col = "black", pch = 19, type = "b")
c = plot(exp_var, model = vg_model16, cloud=T, col = "black", pch = 19, type = "b")
d = c(a,b,c)
d

with(exp_var, plot(dist, gamma, pch=19, cex = 1.3))
lines(variogramLine(vg_model14, maxdist = 2), lwd = 2)
lines(variogramLine(vg_model15, maxdist = 2), col = 2, lwd = 2)
lines(variogramLine(vg_model16, maxdist = 2), col = 4, lwd = 2)
legend("bottomright", bty = "n", lty = 1, col = c(1, 2, 4), lwd = 2, 
       legend = c("Exponential", "Spherical", "Gaussian"))

#plot the variogram model#
plot(exp_var, model = vg_model, cloud=T, col = "blue")

#update the gstat object#
gs = gstat(id = "Ni", formula = Ni~1, locations = ~X+Y, 
           data = saprolit, nmax = 0.5, model = vg_model14)

## 4. kriging ##
#creating grid for kriging#
summary(saprolit)
x = seq(from = 8638, to = 9051)
y = seq(from = 12662, to = 13075)

#the variable X and Y are the same name of the gs object#
xy_grid = expand.grid(X=x, Y=y)

#do kriging#
hasil_kriging = predict(gs, newdata = xy_grid, debug.level = -1)

#plot the result#
aux = hasil_kriging$Ni.pred
#convert a vector to a matrix (for plotting purpose)#
dim(aux) = c(length(x),length(y))
image(x,y,aux,asp=1)
contour(x,y,aux, asp=1, drawlabels = TRUE, col = "black", add = TRUE) 
#plot the data location#
points(saprolit$X, saprolit$Y, pch = "+")

#other plotting methods#
#use lattice package#
levelplot(hasil_kriging$Ni.pred ~ hasil_kriging$X + hasil_kriging$Y, 
          col.regions = heat.colors(100), 
          main = "Ordinary Kriging Estimate",
          xlab = "X", ylab = "Y")

#plot the kriging variance result#
aux_var = hasil_kriging$Ni.var
#convert a vector to a matrix (for plotting purpose)
dim(aux_var) = c(length(x),length(y))
image(x,y,aux_var,asp=1)
contour(x,y,aux_var, asp=1, drawlabels = TRUE, col = "blue", add = TRUE) 
points(saprolit$X, saprolit$Y, pch = "+")

#other plotting methods#
#use lattice package#
levelplot(hasil_kriging$Ni.var ~ hasil_kriging$X + hasil_kriging$Y, 
          col.regions = heat.colors(100),
          main = "Ordinary Kriging Variances" ,
          xlab = "X", ylab = "Y")

#us sp package#
aux_var = hasil_kriging
coordinates(aux_var) <- c("X", "Y")
p <- spplot(aux_var["Ni.var"], asp = 1, colorkey = TRUE, xlim = c(8638,9051), ylim = c(12662,13075), 
            scales = list(draw = TRUE) )
p + layer(panel.points(X, Y, col = "black", pch = 1), data = saprolit)

#cropping the prediction#
aux = hasil_kriging$Ni.pred
aux[hasil_kriging$Ni.var>1.4]=NA
#convert a vector to a matrix (for plotting purpose)
dim(aux) = c(length(x),length(y))
image(x,y,aux,asp=1)
points(limonit$X, limonit$Y, pch = "+")
#other plotting methods#
#use lattice package#
levelplot(aux ~ hasil_kriging$X + hasil_kriging$Y, col.regions = heat.colors(100))
#use sp package#
aux = hasil_kriging$Ni.pred
aux[hasil_kriging$Ni.var>1.4]=NA
aux_sp = hasil_kriging
aux_sp["Ni.pred"]=aux
coordinates(aux_sp) <- c("X", "Y")
p <- spplot(aux_sp["Ni.pred"], asp = 1, colorkey = TRUE, xlim = c(8638,9051), ylim = c(12662,13075), 
            scales = list(draw = TRUE) )
p + layer(panel.points(X, Y, col = "black", pch = 1), data = limonit)

## 5. cross validation ##
cross_val = gstat.cv(gs, nfold = nrow(saprolit))
plot(observed ~ Ni.pred, data = cross_val, asp = 1)
abline(a = 0, b = 1, col = 2, lwd = 2)
#calculate RMSE (root mean squared error)
#sqrt (cross_val$residual^2)/length (cross_val$residual))
sqrt(mean(cross_val$residual^2))
#calculate mean error#
#should be close to 0#
mean(cross_val$residual)
#calculate MSDR (mean squared deviation ratio)#
#should be close to 1#
mean(cross_val$residual^2/cross_val$Ni.var)

## 6. simulation ##
hasil_simulasi = predict(gs, newdata = xy_grid, debug.level = -1, nsim = 10)
aux = hasil_simulasi$sim3
#convert a vector to a matrix (for plotting purpose)#
dim(aux) = c(length(x), length(y))
image(x,y,aux,asp = 1)
points(limonit$X, limonit$Y, pch = "+")
#another plotting method#
#use lattice package#
levelplot(aux ~ hasil_simulasi$X + hasil_simulasi$Y, col.regions = heat.colors(100))
#use sp package#
aux = hasil_simulasi
coordinates(aux) = c("X", "Y")
spplot(aux, colorkey = TRUE)J

## 7. variogram map ##
var_map = variogram(gs, cutoff = 200, width = 20, map = TRUE)
plot(var_map,col.regions = rainbow(100), main="Variogram Map", xlab="x", ylab="y") 

#anisotropic variograms#
exp_vg_aniso = variogram(gs, alpha = c(0,90,180,270))
plot(exp_vg_aniso)
vg_model_aniso = vgm(psill = 0.09, "Sph", range = 25, nugget = 0.15, anis = c(90,1))
plot(exp_vg_aniso, vg_model_aniso, as.table = TRUE, pch = 19, type = "b")

vg_model_anisofit = fit.variogram(exp_vg_aniso, model = vg_model_aniso)
plot(exp_vg_aniso, vg_model_anisofit, as.table = TRUE, pch = 19, type = "b")

#update the gstat object#
gs_aniso = gstat(id = "Ni", formula = Ni~1, locations = ~X+Y, data = limonit)

#kriging using anisotropic variograms#
hasil_kriging_aniso = predict(gs_aniso, newdata = xy_grid, debug.level = -1)

#plot use lattice package#
levelplot(hasil_kriging_aniso$Ni.pred ~ hasil_kriging_aniso$X + hasil_kriging_aniso$Y, col.regions = heat.colors(100))

scatr::scat(data = saprolit, x = X, y = Y, marg = "box")
