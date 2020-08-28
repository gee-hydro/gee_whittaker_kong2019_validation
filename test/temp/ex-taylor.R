## in the examples below, most effort goes into making some artificial data
## the function itself can be run very simply
## Not run:
## dummy model data for 2003
dat <- selectByDate(mydata, year = 2003)
dat <- data.frame(date = mydata$date, obs = mydata$nox, mod = mydata$nox)

## now make mod worse by adding bias and noise according to the month
## do this for 3 different models
dat <- transform(dat, month = as.numeric(format(date, "%m")))
mod1 <- transform(dat, mod = mod + 10 * month + 10 * month * rnorm(nrow(dat)),
                  model = "model 1")
## lag the results for mod1 to make the correlation coefficient worse
## without affecting the sd
mod1 <- transform(mod1, mod = c(mod[5:length(mod)], mod[(length(mod) - 3) :
                                                            length(mod)]))

## model 2
mod2 <- transform(dat, mod = mod + 7 * month + 7 * month * rnorm(nrow(dat)),
                  model = "model 2")
## model 3
mod3 <- transform(dat, mod = mod + 3 * month + 3 * month * rnorm(nrow(dat)),
                  model = "model 3")

mod.dat <- rbind(mod1, mod2, mod3)

# date obs      mod month   model
# 1 1998-01-01 00:00:00 285 488.6313     1 model 1
# 2 1998-01-01 01:00:00  NA 268.8517     1 model 1
# 3 1998-01-01 02:00:00  NA 172.3369     1 model 1
# 4 1998-01-01 03:00:00 493 201.8693     1 model 1
# 5 1998-01-01 04:00:00 468 155.8600     1 model 1
# 6 1998-01-01 05:00:00 264 136.4204     1 model 1

## basic Taylor plot
TaylorDiagram(mod.dat[, c(2:3, 5)] %>% data.table(), obs = "obs", mod = "mod", group = "model")

# ## Taylor plot by season
TaylorDiagram(mod.dat, obs = "obs", mod = "mod", group = "model", type = "season")
#
# ## now show how to evaluate model improvement (or otherwise)
# mod1a <- transform(dat, mod = mod + 2 * month + 2 * month * rnorm(nrow(dat)),
#                    model = "model 1")
# mod2a <- transform(mod2, mod = mod * 1.3)
# mod3a <- transform(dat, mod = mod + 10 * month + 10 * month * rnorm(nrow(dat)),
#                    model = "model 3")
# mod.dat2 <- rbind(mod1a, mod2a, mod3a)
# mod.dat$mod2 <- mod.dat2$mod
#
# ## now we have a data frame with 3 models, 1 set of observations
# ## and TWO sets of model predictions (mod and mod2)
#
# ## do for all models
# TaylorDiagram(mod.dat, obs = "obs", mod = c("mod", "mod2"), group = "model")
#
# ## End(Not run)
# ## Not run:
# ## all models, by season
# TaylorDiagram(mod.dat, obs = "obs", mod = c("mod", "mod2"), group = "model",
#               type = "season")
#
# ## consider two groups (model/month). In this case all months are shown by model
# ## but are only differentiated by model.
#
# TaylorDiagram(mod.dat, obs = "obs", mod = "mod", group = c("model", "month"))
#
# ## End(Not run)
