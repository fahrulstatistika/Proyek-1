# Geostatistika Studi Kasus pada Limonit
LIMONIT
## 1.data_preparation ##

##plot data##
plot(limonit$X,limonit$Y, xlab = "Easting", ylab = "Northing")


##create the histogram of variable nilai##
hist(limonit$Ni, xlab = "Ni", main = "Histogram Ni",col = "gray87")

a = density(limonit$Ni)
hist(limonit$Ni, breaks = 50, xlab = "Ni", freq = FALSE, 
     ylab = "Frequency", main = "Histogram limonit Ni", 
     col = "gray47" )
lines(a)

#Creating Plot
histogram(~ Ni, data = limonit,
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

#spatial interpolation#
coordinates(limonit) <- ~X+Y
v.cloud = variogram(Ni~ 1, data = limonit, cloud=T)
head(v.cloud)
plot(v.cloud, 
     xlab = "Separation distance (m)", 
     col = "black")
v = variogram(Ni ~ 1, data = limonit, cloud = F)
v
plot(v, main = "Variogram - default", xlab = "Separation distance (m)")
v.map = variogram(Ni ~ 1, limonit, map = TRUE, cutoff=200, width=50)
plot(v.map, col.regions = bpy.colors(64), main="Variogram Map", xlab="x", ylab="y")

#create variogram model#
vg_model = vgm(psill = 0.047, model = "Sph", range = 130, nugget = 0.001)
vg_model = fit.variogram(exp_var, model = vg_model)

## 2. Experimental Variogram ##
gs = gstat(id = "Ni", formula = Ni~1, locations = ~X+Y, data = limonit)

##calculate the experimental variogram##
exp_var = variogram(gs)

##plot the experimental variogram##
plot(exp_var)

## 3. Variogram Modeling ##
#create variogram model#
vg_model1 = vgm(psill = 0.044, model = "Sph", range = 130, nugget = 0.004)
vg_model2 = vgm(psill = 0.044, model = "Exp", range = 130, nugget = 0.004)
vg_model3 = vgm(psill = 0.044, model = "Gau", range = 130, nugget = 0.004)

vg_model11 = fit.variogram(exp_var, model = vg_model1)
vg_model12 = fit.variogram(exp_var, model = vg_model2)
vg_model13 = fit.variogram(exp_var, model = vg_model3)

#plot the variogram model#
par(mfrow = c(1,3))
par(mfrow=c(1,3), mar=c(3,3,.5,.5), mgp=c(1.5,.7,0))
a = plot(exp_var, model = vg_model11, cloud=T, col = "black", pch = 19, type = "b")
b = plot(exp_var, model = vg_model12, cloud=T, col = "black", pch = 19, type = "b")
c = plot(exp_var, model = vg_model13, cloud=T, col = "black", pch = 19, type = "b")
d = c(a,b,c)
d


with(exp_var, plot(dist, gamma, pch=19, cex = 1.3))
lines(variogramLine(vg_model11, maxdist = 2), lwd = 2)
lines(variogramLine(vg_model12, maxdist = 2), col = 2, lwd = 2)
lines(variogramLine(vg_model13, maxdist = 2), col = 4, lwd = 2)
legend("bottomright", bty = "n", lty = 1, col = c(1, 2, 4), lwd = 2, 
       legend = c("Exponential", "Spherical", "Gaussian"))

a = lines(variogramLine(vg_model12, maxdist = 2), col = 2, lwd = 2)
b = legend("bottomright", bty = "n", lty = 1, col = c(1, 2, 4), lwd = 2, 
       legend = c("Exponential"))
plot(b)


Co.v <- variogram(Co~1, Co.points.sp,cressie=TRUE)
Co.vm <- fit.variogram(Co.v, vgm(15, 'Exp', 1, nugget=1))
Co.vm

Co.vm.sph <- fit.variogram(Co.v,vgm(14, 'Sph', 1.5, 1))
Co.vm.mat <- fit.variogram(Co.v, vgm(16, 'Mat', 1, 1, k=1))

#plot the empirical binned variogram
with(Co.v, plot(dist, gamma, cex=1.3, pch=19))    
lines(variogramLine(Co.vm, maxdist=2), lwd=2)
lines(variogramLine(Co.vm.sph, maxdist=2), col=2, lwd=2)
lines(variogramLine(Co.vm.mat, maxdist=2), col=4, lwd=2)
legend('bottomright', bty='n', lty=1, col=c(1,2,4), lwd=2,
       legend=c("Exponential","Spherical","Gaussian"))

par(mfrow = c(1,3))
plot(a$x,a$tTF)
plot(a$y,a$tTF)

ols2 <- variofit(geodat.v1, ini=c(.8,30), nugget=.4, cov.model='exponential')
lines(ols2, col=2)
ols3 <- variofit(geodat.v1, ini=c(.8,30), nugget=.4, cov.model='gaussian')
lines(ols3, col=4)
ols4 <- variofit(geodat.v1, ini=c(.8,30), nugget=.4, cov.model='spherical')
lines(ols4, col=3)
ols5 <- variofit(geodat.v1, ini=c(.8,30), nugget=.4, cov.model='linear')
lines(ols5, col='gray70')
legend('bottomright', c('exponential', 'Gaussian', 'spherical', 'linear'), col=c(2:4, 'gray70'), lty=rep(1,4), cex=.8, bty='n')

ols2 = fit.variogram(exp_var,vgm(psill = 0.044, model = "Sph", range = 130, nugget = 0.004) )
lines(ols2, col=2)

#update the gstat object#
gs = gstat(id = "Ni", formula = Ni~1, locations = ~X+Y, 
           data = limonit, nmax = 0.05, 
           model = vg_model12)

## 4. kriging ##
#creating grid for kriging#
summary(limonit)
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
image(x,y, aux,asp=1)
contour(x,y, aux, asp=1, drawlabels = TRUE, col = "black", add = TRUE)
#plot the data location#
points(limonit$X, limonit$Y, pch = "+")

#other plotting methods#
#use lattice package#
levelplot(hasil_kriging$Ni.pred ~ hasil_kriging$X + hasil_kriging$Y, 
         col.regions = heat.colors(100), 
         main = "Ordinary Kriging Estimate",
          xlab = "X", ylab = "Y")

#us sp package#
aux = hasil_kriging
coordinates(aux) = c("X", "Y")
p = spplot(aux["Ni.pred"], asp = 1, colorkey = TRUE, xlim = c(8638,9051), ylim = c(12662,13075), 
           scales = list(draw = TRUE) )
p + layer(panel.points(X, Y, col = "black", pch = 1), data = limonit)

#plot the kriging variance result#
aux_var = hasil_kriging$Ni.var
#convert a vector to a matrix (for plotting purpose)
dim(aux_var) = c(length(x),length(y))
image(x,y,aux_var,asp=1)
points(limonit$X, limonit$Y, pch = "+")

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
p + layer(panel.points(X, Y, col = "black", pch = 1), data = limonit)

#cropping the prediction#
aux = hasil_kriging$Ni.pred
aux[hasil_kriging$Ni.var>1.4]=NA
#convert a vector to a matrix (for plotting purpose)
dim(aux) = c(length(x),length(y))
image(x,y,aux,asp=1)
points(limonit$X, limonit$Y, pch = "+")

#other plotting methods#
#use lattice package#
levelplot(aux ~ hasil_kriging$X + hasil_kriging$Y, col.regions = rainbow(100))

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
cross_val = gstat.cv(gs, nfold = nrow(limonit))
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
spplot(aux, colorkey = TRUE)

## 7. variogram map ##
var_map = variogram(gs, cutoff = 200, width = 100, map = TRUE)
plot(var_map)

#anisotropic variograms#
exp_vg_aniso = variogram(gs, alpha = c(0,45,225,135,315))
plot(exp_vg_aniso)
vg_model_aniso = vgm(psill = 0.02, "Exp", range = 50, nugget = 0, anis = c(45,0.89))
plot(exp_vg_aniso, vg_model_aniso, as.table = TRUE, pch = 19, type = "b")

vg_model_anisofit = fit.variogram(exp_vg_aniso, model = vg_model_aniso)
plot(exp_vg_aniso, vg_model_anisofit, as.table = TRUE, pch = 19, type = "b")

#update the gstat object#
gs_aniso = gstat(id = "Ni", formula = Ni~1, locations = ~X+Y, data = limonit)

#kriging using anisotropic variograms#
hasil_kriging_aniso = predict(gs_aniso, newdata = xy_grid, debug.level = -1)

#plot use lattice package#
levelplot(hasil_kriging_aniso$Ni.pred ~ hasil_kriging_aniso$X + hasil_kriging_aniso$Y, 
          col.regions = rainbow(100))

vg_model1 = vgm(psill = 0.044, model = "Sph", range = 130, nugget = 0.004)
vg_model11 = fit.variogram(exp_var, model = vg_model1)
plot(exp_var, model = vg_model11, cloud=T, col = "black", pch = 19, type = "b")
