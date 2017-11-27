library(data.table)
library(lubridate)
library(ggplot2)
library(glmnet)
library(xgboost)
library(caret)
library(pls)

datadir <- 'data/'
resdir <- 'data/results/'
datafiles <- list.files(datadir)
currdir <- getwd()


getSeason <- function(date) {
  month <- as.numeric(substr(date, 6, 7))
  if (month < 3) return (4) # Winter
  if (month < 6) return (1) # Spring
  if (month < 9) return (2) # Summer
  if (month < 12) return (3) # Fall
  else return (4) # Winter
}

saveGGPlot <- function(name, plot) {
  pdf(paste("data/plots/",name,".pdf",sep=""))
  print(plot)
  dev.off()
}

# Load data and create categoricals
setwd("C:/Users/boltz/mlhomework/f2017-cs584/project/data")
data <- fread("ALL_MWRD.csv")
names(data)[8] <- "NH3N"
names(data)[9] <- "PTOT"
data <- data[, names(data)[-1] := lapply(.SD, as.numeric), .SDcols=names(data)[-1]]
data <- data[, "DATE" := lapply(.SD, function(x) as.IDate(x, format='%d-%b-%y')), .SDcols="DATE"]
# Relatively few NAs for non BOD5 variables, so just remove
data <- data[ !is.na(data$pH) & !is.na(data$SS) & !is.na(data$TS) & !is.na(data$TKN) & !is.na(data$NH3N) & !is.na(data$PTOT)& data$FLOW != 0]
data$SEASON <- as.factor(sapply(data$DATE, getSeason))
data$YEAR <- as.factor(year(data$DATE))
data$MONTH <- as.factor(month(data$DATE))
data$LOCATION <- as.factor(data$LOCATION)

# Impute BOD5
bod5.train <- data[ !is.na(BOD5) ][,c(-1,-10,-11,-13)]
bod5.imp.lm <- lm(BOD5~., data = bod5.train)
data$BOD5[is.na(data$BOD5)] <- predict(bod5.imp.lm, data[is.na(BOD5)][,c(-1,-4,-10,-11,-13)])

stopifnot(!anyNA(data))

# Training is prior to 2016, testing is 2016 
train <- data[YEAR != 2016][,c(-1,-12)]
test <- data[YEAR == 2016][,c(-1,-12)]
stopifnot(NROW(train) + NROW(test) == NROW(data))

# Preliminary plots

plot.bod5 <- ggplot(data = data, mapping = aes(y=PTOT, x=BOD5)) + 
  geom_point(alpha=0.2)  # Huge outliers 
plot.bod5.nout <- ggplot(data = data[BOD5 < 1000], mapping = aes(y=PTOT, x=BOD5)) +
  geom_point(alpha=0.2)  

saveGGPlot('bod5_ptot', plot.bod5)
saveGGPlot('bod5_ptot_lessoutliers', plot.bod5.nout)

plot.flow <- ggplot(data=data, mapping=aes(y=PTOT,x=FLOW, color=LOCATION)) + 
  geom_point(alpha=0.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"), 
                     name="Location",
                     breaks=levels(data$LOCATION),
                     labels=c("Calumet", "Egan", "Hanover Park", "Kirie", "Lamont", "O'Brien")) +
  xlab("Flow (millions of gallons per day") + ylab("Total Phosphorous (mg / L)") + title("Total Phosphorous by Flow and Location (2011 - 2016)")
saveGGPlot('flow_ptot', plot.flow)

plot.date <- ggplot(data=data, mapping=aes(y=PTOT,x=DATE,color=LOCATION)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"), 
                     name="Location",
                     breaks=levels(data$LOCATION),
                     labels=c("Calumet", "Egan", "Hanover Park", "Kirie", "Lamont", "O'Brien")) +
  xlab("Date") + ylab("Total Phosphorous") + title("Total Phosphorous by Location (2011 - 2016)")
saveGGPlot('date_ptot', plot.date)

plot.tkn <- ggplot(data=data, mapping=aes(y=PTOT, x=TKN, color=LOCATION)) +
  geom_point(alpha=0.5) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"), 
                                           name="Location",
                                           breaks=levels(data$LOCATION),
                                           labels=c("Calumet", "Egan", "Hanover Park", "Kirie", "Lamont", "O'Brien")) +
  xlab("Total Kjeldahl Nitrogen") + ylab("Total Phosphorous") + title("Total Phosphorous by Total Kjeldahl Nitrogen (2011 - 2016)")
saveGGPlot('tkn_ptot', plot.tkn)

### Baseline

base.preds <- rep(mean(data$PTOT), NROW(test))
base.mse <- mean((base.preds - test$PTOT)^2)

### Linear model

all.lm <- lm(PTOT~., data=train)
lm.preds <- predict(all.lm, test)
lm.mse <- mean((lm.preds - test$PTOT)^2)
cat(title="Full Linear Regression Model",
    capture.output(summary(all.lm)),
    file="data/results/FullLM_Output.txt",
    sep="\n")
sstot <- sum((test$PTOT - mean(test$PTOT))^2)
ssres <- sum((test$PTOT - lm.preds)^2)
all.rsq <- 1 - (ssres/sstot)


## Single models
lm.rsq <- data.frame(testmse=rep(NA,NCOL(train)),rsq=rep(NA, NCOL(train)),names=names(train))
i <- 1
for (p in train) {
  if (names(train)[i] != "PTOT") {
    formula <- paste("PTOT~",names(train)[i],sep="")
    curr.lm <- lm(formula, data=train)
    print(summary(curr.lm))
    lm.rsq$rsq[i] <- summary(curr.lm)$adj.r.squared
    lm.rsq$testmse[i] <- mean((predict(curr.lm,test[,-8])-test$PTOT)^2)
    
    outname <- paste(resdir,names(train)[i],"output.txt",sep="")
    cat(title=paste(names(train)[i], "Linear Regression Output\n"),
        capture.output(summary(curr.lm)), 
        file=outname,
        sep="\n")
  }
  i <- i+1
}
lm.rsq <- lm.rsq[order(lm.rsq$rsq, decreasing = TRUE),]  # TKN and BOD5 are highest R-Squared single predictors

cat(title="Single Regression Summary",
    capture.output(lm.rsq),
    file="data/results/LMSummary.txt",
    sep="\n")

# Regression with highest R-Squared 

lm.opt <- lm("PTOT~TKN+BOD5+NH3N+LOCATION",data=train)
opt.mse <- mean((predict(lm.opt,test[,-8])-test$PTOT)^2)

# Log transform of response variable
all.resplot <- ggplot(data=train, mapping = aes(y=all.lm$residuals, x=PTOT)) + geom_point(alpha=0.2,color="blue") +
  ylab("Residuals of PTOT~.") + ylab("Total Phosphorous (mg/L)")
saveGGPlot('res_ptot', all.resplot)

lm.resplog <- lm("log(PTOT)~.", data = train)
mse.resplog <- mean((predict(lm.int, test[,-8])-log(test$PTOT))^2)

resplog.resplot <- ggplot(data=train, mapping = aes(y=lm.resplog$residuals, x=log(PTOT))) + geom_point(alpha=0.2,color="blue") +
  ylab("Residuals") + ylab("Total Phosphorous (mg/L)") 
saveGGPlot('resplog_ptot', resplog.resplot)

# Log transforms of independent variables
log.mean <- mean(log(data$PTOT)) 
log.baseline <- mean((log(data$PTOT) - log.mean)^2)

lm.log1 <- lm("log(PTOT)~TKN+log(BOD5)+NH3N+log(SS)+LOCATION", data = train)
lm.log2 <- lm("log(PTOT)~log(TKN)+log(BOD5)+NH3N+log(SS)+LOCATION+SEASON", data = train)
log1.mse <- mean((predict(lm.log1, test[,-8])-log(test$PTOT))^2)
log2.mse <- mean((predict(lm.log2, test[,-8])-log(test$PTOT))^2) ## WINNER

log2.resplot <- ggplot(data=train, mapping = aes(y=lm.log2$residuals, x=log(PTOT))) + geom_point(alpha=0.2,color="blue") +
  ylab("Residuals") + xlab("log(Total Phosphorous) in mg/L") + ggtitle("Residual Plot with Log Transformations\n of Response and Predictors")
saveGGPlot('log2_ptot', log2.resplot) 

# New log dataset to use for the rest of the analysis
log.data <- data.frame(flow=data$FLOW, ph=data$pH, logbod5=log(data$BOD5), logts=log(data$TS),
                       logss=log(data$SS), logtkn=log(data$TKN), nh3n=data$NH3N, logptot=log(data$PTOT), 
                       location=data$LOCATION, season=data$SEASON, year=data$YEAR)
log.train <- log.data[log.data$year != 2016,][,-11]
log.test <- log.data[log.data$year == 2016,][,-11]

# Ridge regression 
train.num <- model.matrix(~.-1,log.train)
test.num <- model.matrix(~.-1,log.test)

set.seed(53902)
cv.ridge <- cv.glmnet(x=train.num[,-8],y=train.num[,8])
ridge.mse <- mean((predict(cv.ridge,newx=test.num[,-8],s=cv.ridge$lambda.min) - log.test$logptot)^2)

#PCR
pcr.fit <- pcr(formula=logptot~.,
             data=log.train,
             scale=TRUE,
             validation="CV")
validationplot(pcr.fit,val.type="MSEP")
pcr.pred <- predict(pcr.fit, log.test[,-8])
pcr.mse <- mean((pcr.pred-log.test$logptot)^2)

# xgboost
xgb.train <- xgb.DMatrix(data=train.num[,-8],label=train.num[,8])
xgb.test <- xgb.DMatrix(data=test.num[,-8], label=test.num[,8])
wtch <- list(train=xgb.train, test=xgb.test)

xgb.fit <- xgb.train(data=xgb.train, max_depth=5, eta=0.01,nthread=4,nrounds=1000,watchlist=wtch, objective="reg:linear")
xgb.min.mse <- xgb.fit$evaluation_log$test_rmse[which(xgb.fit$evaluation_log$test_rmse == min(xgb.fit$evaluation_log$test_rmse))]^2
xgb.imp <- xgb.importance(model=xgb.fit)
xgb.plot.importance(importance_matrix = xgb.imp)

xgb.plot <- ggplot(data=xgb.imp, mapping = aes(x=Feature, y=Gain)) + geom_col(fill="#FF6666") + 
  scale_x_discrete(name="Feature",
                   labels=c("Flow","pH","locKirie","locLamont","locObrien","Summer", "Fall","Winter", "log(BOD5)", 
                            "log(TS)", "log(SS)", "log(TKN)", "NH3N", "locCalumet", "locEgan","locHanPark")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggtitle("Feature Importance by Gain in XGBoost")
saveGGPlot('xgb_imp', xgb.plot)
