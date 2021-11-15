#==================================================================================
# Medium Spatial Autocorrelation Data Generation
#==================================================================================

library(sp)
library(raster)
library(gstat)
library(lattice)
library(ncf)

set.seed(44)
rmvn <- function(n, mu = 0, V = matrix(1)){
  p <- length(mu)
  if (any(is.na(match(dim(V), p)))) 
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
}

# square lattice region
simgrid <- expand.grid(1:100, 1:100)
n <- nrow(simgrid)

# distance matrix
distance <- as.matrix(dist(simgrid))
# generation of random variable
phi <- 0.04
set.seed(44)
X1 <- rmvn(1, rep(0, n), exp(-phi * distance))
X2 <- rmvn(1, rep(10,n), exp(-phi * distance))

# visualizing the results of simulating a random variable 
Xraster1 <- raster::rasterFromXYZ(cbind(simgrid[, 1:2] - 0.5, X1))
par(mfrow = c(2, 2))
plot(1:200, exp(-phi * 1:200), type = "l", xlab = "Distance", ylab = "Correlation")
plot(Xraster1, main = 'X1') # variable 1 with high spatial autocorrelation

Xraster2 <- raster::rasterFromXYZ(cbind(simgrid[, 1:2] - 0.5, X2))
par(mfrow = c(1, 2))
plot(1:200, exp(-phi * 1:200), type = "l", xlab = "Distance", ylab = "Correlation")
plot(Xraster2, main = 'X2') # variable 2 with high spatial autocorrelation


set.seed(44)
Xvar3 <- raster(matrix(rnorm(n), 100, 100), xmn = 0, xmx = 50, ymn = 0, ymx = 50)
plot(Xvar3, main = "X3") # covariate 3 with no spatial autocorrelation, variable is random on plot


beta_0 <- 0.1
beta_1 <- 1
beta_2 <- -1

sigma <- 1 

set.seed(44)
eps <- rnorm(length(X1),0,sigma)
# linear function
Y <- (beta_0 + beta_1*values(Xvar3) + beta_2*values(Xvar3)^2 + 
        2*values(Xraster1) + 1*values((Xraster2))+ eps) # we need an error term attached to our y_variable
coords <- coordinates(Xvar3)
Y = cbind(Y,coords)
head(Y)

XrasterY <- raster::rasterFromXYZ(cbind(simgrid[, 1:2] - 0.5, Y[,1]))
par(mfrow = c(1, 4))
plot(1:200, exp(-phi * 1:200), type = "l", xlab = "Distance", ylab = "Correlation")
plot(XrasterY, main = 'Y')
mtext("y = (0.1 + (Xvar3) - (Xvar3)^2 + 2(Xvar1) + (Xvar2) + E)", side = 3, line = -2, outer = TRUE)


id <- sample(1:n, 10000)
coords <- coordinates(Y)[id, c("x","y")]
y_sampled = Y[id ,]

head(y_sampled)
data <- data.frame(coords, # this gives x,y coordinates of sample
                   Xvar3 = raster::extract(Xvar3, coords), # X1 variable
                   Xvar1 = raster::extract(Xraster1,coords), # X2 variable
                   Xvar2 = raster::extract(Xraster2,coords), # X3 variable
                   y_variable = y_sampled[,1]) # Y variable
head(data) 

library(sf)
# Convert data frame to sf object
my.sf.point <- sf::st_as_sf(x = data, 
                            coords = c("x", "y"),
                            crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

#plot(my.sf.point)

#========================================
# Moran's I Calculation 
#========================================

coordinates = st_coordinates(my.sf.point)
library(spdep)
dnbTresh1=dnearneigh(coordinates,0,2,longlat=FALSE) # change distance to account for all points (no islands)
summary(dnbTresh1)
plot(st_geometry(my.sf.point))
plot(dnbTresh1, coordinates, add=TRUE, col="blue")
dnbTresh1.listw=nb2listw(dnbTresh1,style="W",zero.policy=TRUE)
class(dnbTresh1.listw)
moran.test(my.sf.point$y_variable,dnbTresh1.listw)

#====================================
# RELATIONSHIP BETWEEN VARIABLES
#====================================
par(mfrow = c(1,3))
plot(data$Xvar1,(data$y_variable), ylab = "y-variable", xlab = "X-variable 1", col = "blue") 
plot(data$Xvar2,(data$y_variable),  ylab = "y-variable", xlab = "X-variable 2", col = "blue")
plot(data$Xvar3,(data$y_variable),  ylab = "y-variable", xlab = "X-variable 3", col = "blue")


corr_matrix <- 
  cor(data[,-c(1,2)]) %>% 
  round(., digits=2)
library(ggcorrplot)
ggcorrplot( corr = corr_matrix,               # correlation matrix to visualise
            ggtheme = theme_bw,     # simple theme for plot
            title = "Visualisation of Correlations Between Variables", 
            show.legend = TRUE,
            show.diag = TRUE,
            legend.title = "Correlation",
            hc.order = TRUE,         # order the coefficients       
            lab = TRUE,              # show correlation coefficients on the plot 
            lab_size = 3,            # size of displayed coefficients
            colors = c("tomato2", "white", "springgreen3")
)

#========================================
# Linear Model
#========================================

rp <- data %>%
  dplyr::select(Xvar1, Xvar2, Xvar3, y_variable)
head(rp, 2)
# co-ordinates dataframe, required for mlr 
co <- data %>%
  dplyr::select(x, y)
head(co, 2)

library(mlr)

task_lm = mlr::makeRegrTask(data = rp, target = "y_variable", coordinates = co)

lrns_lm = mlr::listLearners(task_lm, warn.missing.packages = FALSE)
dplyr:: select(lrns_lm, class, name, short.name, package)

lrn_lm = mlr:: makeLearner(cl = "regr.lm", predict.type = "response")

getLearnerPackages(lrn_lm)
helpLearner(lrn_lm)

getLearnerModel(train(lrn_lm, task_lm))

perf_level_lm = mlr::makeResampleDesc(method = "SpRepCV",
                                      folds = 5,
                                      reps = 50)

cv_sp_lm = mlr::resample(
  task = task_lm, 
  learner = lrn_lm,
  resampling = perf_level_lm,
  measures = mlr::rmse
) # 3.3221932

task_nsp_lm = mlr::makeRegrTask(data = rp, target = "y_variable")
perf_level_nsp_lm = mlr::makeResampleDesc(method = "RepCV", folds = 5, reps = 50)
cv_nsp_lm = mlr::resample(learner = lrn_lm,
                          task = task_nsp_lm,
                          resampling = perf_level_nsp_lm,
                          measures = mlr::rmse) # 2.8637226 

rmse_sp <- cv_sp_lm$measures.test$rmse
rmse_sp
rmse_nsp <- cv_nsp_lm$measures.test$rmse
rmse_nsp

library("RColorBrewer")
boxplot(rmse_sp, rmse_nsp,
        horizontal = TRUE,
        names= c("S-CV","NS-CV"),
        col=brewer.pal(n = 3, name = "RdBu"),
        xlab="RMSE")

(summary(cv_sp_lm$measures.test$rmse))
(summary(cv_nsp_lm$measures.test$rmse))

#========================================
# Random Forest
#========================================

library(dplyr)

rp <- data %>%
  dplyr::select(Xvar1, Xvar2, Xvar3, y_variable)
head(rp, 2)

co <- data %>%
  dplyr::select(x, y)
head(co, 2)

task_rf = mlr::makeRegrTask(data = rp, target = "y_variable", coordinates = co)

lrns_rf = mlr::listLearners(task_rf, warn.missing.packages = FALSE)
dplyr::select(lrns_rf, class, name, short.name, package)

lrn_rf = mlr::makeLearner(cl = "regr.randomForest", predict.type = "response") # random forest model
lrn_rf$par.set #OR
getParamSet("regr.randomForest")

library(ParamHelpers)
integer_ps = makeParamSet(
  makeIntegerParam("mtry", lower = 1, upper = 3),
  makeIntegerParam("sampsize", lower = 1, upper = 3000)
)

print(integer_ps)

tune_level_rf = mlr::makeResampleDesc(method = "SpCV", iters = 2)
ctrl = mlr::makeTuneControlRandom(maxit = 20) # random search with 50 iterations, to find the optimal hyperparameters 

wrapped_lrn_rf = 
  mlr::makeTuneWrapper(learner = lrn_rf,
                       #inner loop
                       resampling = tune_level_rf,
                       #hyperparameter search space
                       par.set = integer_ps,
                       #random search
                       control = ctrl,
                       show.info = TRUE,
                       #measuree of performance
                       measures = mlr::rmse)

perf_level_rf = mlr::makeResampleDesc(method = "SpRepCV",
                                      folds = 5,
                                      reps = 20) # change back to 100 


set.seed(44) # this might be slow, run it after doing time series 
cv_sp_rf = mlr::resample(learner = wrapped_lrn_rf,
                         task = task_rf,
                         resampling = perf_level_rf,
                         extract = getTuneResult,
                         measures = mlr::rmse)

(summary(cv_sp_rf$measures.test$rmse))

# normal (non-spatial) cross validation 

task_nsp = makeRegrTask(data = rp, target = "y_variable")
perf_level_nsp = makeResampleDesc(method = "RepCV", folds = 5, reps = 20) # change reps to 100

tune_level_nsp = mlr::makeResampleDesc(method = "CV", iters = 2)
ctrl_nsp = mlr::makeTuneControlRandom(maxit = 20) # random search with 50 iterations, to find the optimal hyperparameters 

#wrapping everything together now 
wrapped_lrn_rf_nsp = 
  mlr::makeTuneWrapper(learner = lrn_rf,
                       #inner loop
                       resampling = tune_level_nsp,
                       #hyperparameter search space
                       par.set = integer_ps,
                       #random search
                       control = ctrl_nsp,
                       show.info = TRUE,
                       #measuree of performance
                       measures = mlr::rmse)

set.seed(33)
cv_nsp_rf = mlr::resample(learner = wrapped_lrn_rf_nsp,
                          task = task_nsp,
                          resampling = perf_level_nsp,
                          extract = getTuneResult,
                          measures = mlr::rmse)

(summary(cv_nsp_rf$measures.test$rmse))

rmse_sp_rf <- cv_sp_rf$measures.test$rmse
# onject of 100 rmse values 
rmse_nsp_rf <- cv_nsp_rf$measures.test$rmse

library("RColorBrewer")
boxplot(rmse_sp_rf, rmse_nsp_rf,
        horizontal = TRUE,
        names= c("S-CV","NS-CV"),
        col=brewer.pal(n = 3, name = "RdBu"),
        xlab="RMSE")

#=============================================================================================
# Support vector machine
#=============================================================================================
library(mlr)
# SPATIAL CROSS-VALIDATION-----------
smv <- data %>%
  dplyr::select(Xvar1, Xvar2, Xvar3, y_variable)
head(smv, 2)
# co-ordinates dataframe, required for mlr 
co <- data %>%
  dplyr::select(x, y)
head(co, 2)

# Task
task = makeRegrTask(data = smv, target = "y_variable",
                    coordinates = co)
listLearners(task, warn.missing.packages = FALSE) %>%
  dplyr::select(class, name, short.name, package) %>%
  head(20)

#learner
lrn_svm <- makeLearner("regr.svm", predict.type = "response")
getParamSet("regr.svm")

# define the outer limits of the randomly selected hyperparameters

library(ParamHelpers)
kernels <- c("polynomial","radial","sigmoid")
?makeLearner
ps =  makeParamSet(
  makeDiscreteParam("kernel", values = kernels),
  makeIntegerParam("degree", lower =1, upper=3),
  makeNumericParam("cost",lower=5, upper =20),
  makeNumericParam("gamma", lower = 0, upper = 2)
)

# use 50 randomly selected hyperparameters
ctrl = makeTuneControlRandom(maxit = 20)

# performance estimation level
perf_level = makeResampleDesc(method = "SpRepCV", folds = 5, reps = 50)

?makeResampleDesc

# five spatially disjoint partitions
tune_level = makeResampleDesc("SpCV", iters = 20.)



dwrapped_lrn_svm = makeTuneWrapper(learner = lrn_svm, 
                                   resampling = tune_level,
                                  par.set = ps,
                                  control = ctrl, 
                                  show.info = TRUE,
                                  measures = mlr::rmse)



set.seed(44)
sp_cv_svm = mlr::resample(learner = wrapped_lrn_svm,
                          task = task,
                          resampling = perf_level,
                          extract = getTuneResult,
                          measures = mlr::rmse) #3.0832771


# 

# NORMAL CROSS-VALIDATION ------------

task_nsp = makeRegrTask(data = smv, target = "y_variable")
perf_level_nsp= makeResampleDesc(method = "RepCV", folds = 5, reps = 20) # change reps to 100

tune_level_nsp = mlr::makeResampleDesc(method = "CV", iters = 5)
ctrl_nsp = mlr::makeTuneControlRandom(maxit = 5) # random search with 20 iterations, to find the optimal hyperparameters 

#wrapping everything together now 
wrapped_lrn_svm_nsp = 
  mlr::makeTuneWrapper(learner = lrn_svm,
                       #inner loop
                       resampling = tune_level_nsp,
                       #hyperparameter search space
                       par.set = ps,
                       #random search
                       control = ctrl_nsp,
                       show.info = TRUE,
                       #measuree of performance
                       measures = mlr::rmse)

set.seed(44)
cv_nsp = mlr::resample(learner = wrapped_lrn_svm_nsp,
                       task = task_nsp,
                       resampling = perf_level_nsp,
                       extract = getTuneResult,
                       measures = mlr::rmse) #2.3781431

#===============================================================================
# Print results
#===============================================================================

rmse_sp <- sp_cv_svm$measures.test$rmse
rmse_sp
summary(rmse_sp)
rmse_nsp <- cv_nsp$measures.test$rmse
rmse_nsp
summary(rmse_nsp)
library("RColorBrewer")
boxplot(rmse_sp, rmse_nsp,
        horizontal = TRUE,
        names = c("S-CV","NS-CV"),
        col=brewer.pal(n=3,name="RdBu"),
        xlab="RMSE")

