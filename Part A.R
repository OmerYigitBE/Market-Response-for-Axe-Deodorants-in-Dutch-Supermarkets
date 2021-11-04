library(ggplot2)
library(dplyr)



#Importing data

deo <- read_xls("DataDeodorant.xls")

#Changing f+d into comb
colnames(deo)[43] = "FAComb"; colnames(deo)[44] = "NIVEAComb"; colnames(deo)[45] = "REXONAComb"; colnames(deo)[46] = "SANEXComb"; colnames(deo)[47] = "VOGUEComb"; colnames(deo)[48] = "8X4Comb"; colnames(deo)[49] = "DOVEComb"; colnames(deo)[50] = "AXEComb"

#Exploratory analysis
str(deo)
summary(deo)



#Detecting top 3 brands
totalSales <- data.frame(Brand = c("Dove","Fa","Nivea","Rexona","Sanex","Vogue","8x4","Axe"),
                         Sales = c(sum(deo$DOVESales), #15979.72
                                   sum(deo$FASales), #19193.77
                                   sum(deo$NIVEASales), #13052.5
                                   sum(deo$REXONASales), #37782.92
                                   sum(deo$SANEXSales), #19188.61
                                   sum(deo$VOGUESales), #11883.94
                                   sum(deo$`8X4Sales`), #18309.83
                                   sum(deo$AXESales))) #35887.45

### AXE ###

# Market share
deo$AXEMarket <- deo$AXESales/(deo$DOVESales+deo$FASales+deo$NIVEASales+deo$REXONASales+deo$SANEXSales+deo$VOGUESales+deo$`8X4Sales`+deo$AXESales)
#

share <- deo[,2:10] %>% group_by(WEEK) %>% summarise_each(sum)
share$Total <- rowSums(share[,2:9])
share$AXEShare <- share$AXESales/share$Total
share$Year <- 2000 + as.integer(substr(share$WEEK, 2,3))
share$Week <- as.integer(substr(share$WEEK, 5,6))

axe_ts <- ts(share$AXEShare, frequency = 52)
plot.ts(axe_ts)
plot(decompose(axe_ts))

# Price
price <- deo[, c(2,11:18)] %>% group_by(WEEK) %>% summarise_each(mean)
cor(price[,-1])

# Price Difference

price_diff <- data.frame("8X4" = deo$`8X4RPrice` - deo$`8X4Price`,
                         "AXE" = deo$AXERPrice - deo$AXEPrice,
                         "DOVE" = deo$DOVERPrice - deo$DOVEPrice,
                         "FA" = deo$FARPrice - deo$FAPrice,
                         "NIVEA" = deo$NIVEARPrice - deo$NIVEAPrice,
                         "REXONA" = deo$REXONARPrice - deo$REXONAPrice,
                         "SANEX" = deo$SANEXRPrice - deo$SANEXPrice,
                         "VOGUE" = deo$VOGUERPrice - deo$VOGUEPrice)
deo <- cbind(deo, price_diff)

promotion <- data.frame("8X4" = c(cor(deo$`8X4DISP`,deo$`8X4Sales`), cor(deo$`8X4FEAT`,deo$`8X4Sales`), cor(deo$`8X4Comb`,deo$`8X4Sales`)),
                        "AXE" = c(cor(deo$AXEDISP,deo$AXESales), cor(deo$AXEFEAT,deo$AXESales), cor(deo$AXEComb,deo$AXESales)),
                        "DOVE" = c(cor(deo$DOVEDISP,deo$DOVESales), cor(deo$DOVEFEAT,deo$DOVESales), cor(deo$DOVEComb,deo$DOVESales)),
                        "FA" = c(cor(deo$FADISP,deo$FASales), cor(deo$FAFEAT,deo$FASales), cor(deo$FAComb,deo$FASales)),
                        "NIVEA" = c(cor(deo$NIVEADISP,deo$NIVEASales), cor(deo$NIVEAFEAT,deo$NIVEASales), cor(deo$NIVEAComb,deo$NIVEASales)),
                        "REXONA" = c(cor(deo$REXONADISP,deo$REXONASales), cor(deo$REXONAFEAT,deo$REXONASales), cor(deo$REXONAComb,deo$REXONASales)),
                        "SANEX" = c(cor(deo$SANEXDISP,deo$SANEXSales), cor(deo$SANEXFEAT,deo$SANEXSales), cor(deo$SANEXComb,deo$SANEXSales)),
                        "VOGUE" = c(cor(deo$VOGUEDISP,deo$VOGUESales), cor(deo$VOGUEFEAT,deo$VOGUESales), cor(deo$VOGUEComb,deo$VOGUESales)))

promotion <- round(promotion, 4)

promotion_axe <- data.frame(Chain = deo$Chain, Sales = deo$AXESales, Display = deo$AXEDISP, Feature = deo$AXEFEAT, Combination = deo$AXEComb)
promotion_axe_retail <- promotion_axe %>% group_by(Chain) %>% summarize(corr_disp = cor(Sales,Display), corr_disp_test = cor.test(Sales,Display)$p.value,
                                                                        corr_feat = cor(Sales,Feature), corr_feat_test = cor.test(Sales,Feature)$p.value)
                                                                        #corr_comb = cor(Sales,Combination), corr_comb_test = cor.test(Sales,Combination)$p.value)

alb <- promotion_axe[promotion_axe$Chain == "ALBERT HEIJN",]
c1000 <- promotion_axe[promotion_axe$Chain == "C-1000",]
edah <- promotion_axe[promotion_axe$Chain == "EDAH",]
jumbo <- promotion_axe[promotion_axe$Chain == "JUMBO",]
boer <- promotion_axe[promotion_axe$Chain == "SUPER DE BOER",]


cor.test(deo$`8X4DISP`,deo$`8X4Sales`)$p.value
cor.test(deo$AXEFEAT,deo$AXESales)$p.value
cor.test(deo$DOVEComb, deo$DOVESales)$p.value
cor.test(deo$FAFEAT, deo$FASales)$p.value
cor.test(deo$NIVEAFEAT, deo$NIVEASales)$p.value
cor.test(deo$REXONADISP, deo$REXONASales)$p.value
cor.test(deo$SANEXFEAT, deo$SANEXSales)$p.value
cor.test(deo$VOGUEDISP, deo$VOGUESales)$p.value

#8X4 Promotions retail level
eight4<- data.frame(Chain = deo$Chain, Sales = deo$`8X4Sales`, Promotion = deo$X8X4)
retail_8x4 <- eight4 %>% group_by(Chain) %>% summarize(correlation = cor(Sales, Promotion), test = cor.test(Sales, Promotion)$p.value)

#AXE Promotions retail level
axe <- data.frame(Chain = deo$Chain, Sales = deo$AXESales, Promotion = deo$AXE)
retail_axe <- axe %>% group_by(Chain) %>% summarize(correlation = cor(Sales, Promotion), test = cor.test(Sales, Promotion)$p.value)

#DOVE Promotions retail level
dove <- data.frame(Chain = deo$Chain, Sales = deo$DOVESales, Promotion = deo$DOVE)
retail_dove <- dove %>% group_by(Chain) %>% summarize(correlation = cor(Sales, Promotion), test = cor.test(Sales, Promotion)$p.value)

#FA Promotions retail level
fa <- data.frame(Chain = deo$Chain, Sales = deo$FASales, Promotion = deo$FA)
retail_fa <- fa %>% group_by(Chain) %>% summarize(correlation = cor(Sales, Promotion), test = cor.test(Sales, Promotion)$p.value)

#NIVEA Promotions retail level
nivea <- data.frame(Chain = deo$Chain, Sales = deo$NIVEASales, Promotion = deo$NIVEA)
retail_nivea <- nivea %>% group_by(Chain) %>% summarize(correlation = cor(Sales, Promotion), test = cor.test(Sales, Promotion)$p.value)

#REXONA Promotions retail level
rexona <- data.frame(Chain = deo$Chain, Sales = deo$REXONASales, Promotion = deo$REXONA)
retail_rexona <- rexona %>% group_by(Chain) %>% summarize(correlation = cor(Sales, Promotion), test = cor.test(Sales, Promotion)$p.value)

#SANEX Promotions retail level
sanex <- data.frame(Chain = deo$Chain, Sales = deo$SANEXSales, Promotion = deo$SANEX)
retail_sanex <- sanex %>% group_by(Chain) %>% summarize(correlation = cor(Sales, Promotion), test = cor.test(Sales, Promotion)$p.value)

#VOGUE Promotions retail level
vogue <- data.frame(Chain = deo$Chain, Sales = deo$VOGUESales, Promotion = deo$VOGUE)
retail_vogue <- vogue %>% group_by(Chain) %>% summarize(correlation = cor(Sales, Promotion), test = cor.test(Sales, Promotion)$p.value)


retail_corr <- data.frame(Chain = levels(retail_axe$Chain), "8X4" = retail_8x4$correlation, "AXE" = retail_axe$correlation, "DOVE" = retail_dove$correlation,
                     "FA" = retail_fa$correlation, "NIVEA" = retail_nivea$correlation, "REXONA" = retail_rexona$correlation, "SANEX" = retail_sanex$correlation, "VOGUE" = retail_vogue$correlation)
retail_corr <- cbind(retail_corr$Chain, round(retail_corr[,-1], 4))




#############################
### MARKET RESPONSE MODEL ###
#############################

library(MASS)
library(readxl)
library(scales)

# Preparing model-ready data
deo <- read_xlsx("Deodorant.xlsx", sheet="Model")
deo <- deo[,-c(2,3,4)] #Remove times
deo <- deo[,-c(5,6,7,8,9)] #Remove discounts

deo <- deo[-1,] #Remove first observation

deo$other_price_diff <- as.factor(deo$other_price_diff)


#### MODELS

# Linear
full_model <- lm(market_share ~ ., data=deo)
summary(full_model)
linear_model <- stepAIC(full_model, direction = "both", trace = FALSE)
summary(linear_model) #std=0.03517, R2=0.8363, AIC=-465.6181


# Multiplicative
deo$other_price_diff <- as.integer(deo$other_price_diff)
deo_plus1 <- cbind(deo[,c(1:4)], deo[,c(5:10)]+1)
deo_log <- log(deo_plus1)
full_model_multiplicative <- lm(market_share ~ ., data=deo_log)
summary(full_model_multiplicative)
multiplicative_model <- stepAIC(full_model_multiplicative, direction = "both", trace = FALSE)
summary(multiplicative_model)
AIC(multiplicative_model) #-86.92293


# Constrained (market_share < 1)
full_model_constrained <- glm(market_share~., data=deo, family=binomial(link="logit"))
summary(full_model_constrained)
#AIC = 81.366, #Null = 4.94790, #Residual = 0.76712
constrained_model <- stepAIC(full_model_constrained, direction = "both", trace = FALSE)
summary(constrained_model)
AIC(constrained_model) #68.1306


ccc <- glm(formula=linear_model$call$formula, data=deo, family=binomial(link="logit"))
summary(ccc)



# Time series
deo_time <- read_xlsx("Deodorant.xlsx", sheet="Model")[,1:4]
deo_ts <- ts(deo_time[,1], start=c(2003, 46), end=c(2006, 9), frequency = 52)
plot(diff(log(deo_ts)))
max.lag <- round(sqrt(length(deo_ts)))
CADFtest(deo_ts, type = "drift", criterion = "BIC", max.lag.y = max.lag) #Stationary
acf(deo_ts)  #2
pacf(deo_ts) #2
time_series_model <- arima(deo_ts, order=c(2,0,0))
Box.test(time_series_model$residuals, lag = max.lag, type = "Ljung-Box") #Valid model


#### MULTICOLLINEARITY
correlations <- cor(deo)
library(car)
vif(linear_model)
vif(multiplicative_model)
vif(constrained_model)


#### HETEROSCEDASTICITY (GOLDFELD QUANDT TEST)
plot(linear_model$residuals)
plot(multiplicative_model$residuals)
plot(constrained_model$residuals)
plot(time_series_model$residuals)

plot(linear_model)


#Linear
linear_res <- linear_model$residuals
t <- length(linear_res)
r <- t/3
k <- length(linear_model$coefficients)
sub1 <- linear_res[1:((t-r)/2)]
sub2 <- tail(linear_res, ((t-r)/2))
middle <- linear_res[(((t-r)/2)+1):(t-r)]

sum(linear_res == append(sub1, append(middle, sub2))) == length(linear_res)

sse1 <- sum(sub1**2); sse2 <- sum(sub2**2)
test_stat <- sse2/sse1
dof <- ((t-r)/2)-k
1-pf(test_stat, dof, dof)

#Multiplicative
multiplicative_res <- multiplicative_model$residuals
t <- length(multiplicative_res)
r <- t/3
k <- length(multiplicative_model$coefficients)
sub1 <- multiplicative_res[1:((t-r)/2)]
sub2 <- tail(multiplicative_res, ((t-r)/2))
middle <- multiplicative_res[(((t-r)/2)+1):(t-r)]

sum(multiplicative_res == append(sub1, append(middle, sub2))) == length(multiplicative_res)

sse1 <- sum(sub1**2); sse2 <- sum(sub2**2)
test_stat <- sse2/sse1
dof <- ((t-r)/2)-k
1-pf(test_stat, dof, dof)

#Constrained
constrained_res <- constrained_model$residuals
t <- length(constrained_res)
r <- t/3
k <- length(constrained_model$coefficients)
sub1 <- constrained_res[1:((t-r)/2)]
sub2 <- tail(constrained_res, ((t-r)/2))
middle <- constrained_res[(((t-r)/2)+1):(t-r)]

sum(constrained_res == append(sub1, append(middle, sub2))) == length(constrained_res)

sse1 <- sum(sub1**2); sse2 <- sum(sub2**2)
test_stat <- sse2/sse1
dof <- ((t-r)/2)-k
1-pf(test_stat, dof, dof)

#Time series
time_series_res <- time_series_model$residuals
t <- length(time_series_res)
r <- t/3
k <- length(time_series_model$coefficients)
sub1 <- time_series_res[1:((t-r)/2)]
sub2 <- tail(time_series_res, ((t-r)/2))
middle <- time_series_res[(((t-r)/2)+1):(t-r)]

sum(time_series_res == append(sub1, append(middle, sub2))) == length(time_series_res)

sse1 <- sum(sub1**2); sse2 <- sum(sub2**2)
test_stat <- sse1/sse2
dof <- ((t-r)/2)-k
1-pf(test_stat, dof, dof)

#### STABILITY & OUT-OF-SAMPLE PERFORMANCE
library(caret)

#Linear
set.seed(123) 
train.control <- trainControl(method = "cv", number = 5)
linear_cross_validation <- train(linear_model$call$formula, data = deo, method = "lm", trControl = train.control)
print(linear_cross_validation) #RMSE=0.03657209, Rsquared=0.7838862, MAE=0.02860226     
#Multiplicative
set.seed(123) 
train.control <- trainControl(method = "cv", number = 5)
multiplicative_cross_validation <- train(multiplicative_model$call$formula, data = deo_log, method = "lm", trControl = train.control)
print(multiplicative_cross_validation) #RMSE=0.1802897, Rsquared=0.7531931, MAE=0.1353043  
#Constrained
set.seed(123) 
train.control <- trainControl(method = "cv", number = 5)
constrained_cross_validation <- train(constrained_model$call$formula, data = deo, method = "glm", trControl = train.control)
print(constrained_cross_validation) #RMSE=0.05159391, Rsquared=0.6629266, MAE=0.03810132  

#### ULTIMATE MODEL


#### MANAGERIAL INSIGHTS & DECISIONS







set.seed(1)
data = data.frame(x=1:100)
data$y = 1 / (1 + exp(5-0.1*(data$x) + rnorm(100)))

model = glm(y~x, family = 'binomial', data=data)
summary(model)
plot(data$x, data$y)
lines(data$x, predict(model, data, type = 'response'))






###################################
####### ASSIGNMENT 3 ##############
###################################

library(readxl)
deo <- read_xlsx("Deodorant.xlsx", sheet="Model")

#Time series analysis for seasonality
deo_ts <- ts(deo$market_share, frequency = 52)
plot(decompose(deo_ts))

#Relative price variable
deo$axe_relative_price <- deo$axe_price / deo$other_price
deo$axe_relative_r_price <- deo$axe_r_price / deo$other_r_price

#Dummy for quarters
deo$d_q1 <- ifelse(deo$quarter == 1, 1, 0)
deo$d_q2 <- ifelse(deo$quarter == 2, 1, 0)
deo$d_q3 <- ifelse(deo$quarter == 3, 1, 0)

#Dummy for price war
deo$d_p <- ifelse(deo$other_price > 2.25, 0, 1)


### MODELS ###

#Model ready data
deo <- deo[, -c(2,3,4)] #Remove times
deo <- deo[-1,] #Remove first observation
deo <- deo[, -c(2,3,5,6,8,9)] #Remove price variables (use relatives instead)

# Linear
full_model <- lm(market_share ~ ., data=deo)
summary(full_model)
linear_model <- stepAIC(full_model, direction = "both", trace = FALSE)
summary(linear_model) #std=0.03479, R2=0.8399, F-stat: 107.7 (dof = 6,116), AIC=-468.3227

# Multiplicative
deo$other_price_diff <- as.integer(deo$other_price_diff)
deo_mult <- cbind(deo[,c(1,10,15)], deo[,c(2,3)]+2, deo[,c(4:9, 11:14)]+1)
deo_log <- log(deo_mult)
deo_log <- deo_log[-c(16,26,34,35),] #-Inf logs, relative_r_prices are the same.
full_model_multiplicative <- lm(market_share ~ ., data=deo_log)
summary(full_model_multiplicative)
multiplicative_model <- stepAIC(full_model_multiplicative, direction = "both", trace = FALSE)
summary(multiplicative_model) #std=0.1748, R2=0.6988, F-stat: 46.62 (dof = 6,112), AIC=-68.59807

# Constrained (market_share < 1)
full_model_constrained <- glm(market_share~., data=deo, family=binomial(link="logit"))
summary(full_model_constrained)
constrained_model <- stepAIC(full_model_constrained, direction = "both", trace = FALSE)
summary(constrained_model) #AIC = 66.067, #Null = 4.94790 (dof=122), #Residual = 1.6731 (dof=121)


### DYNAMICS ###

# Static model
static_model <- linear_model
summary(static_model)

# Determining number of lags
# Variables: axe_relative_price, axe_relative_r_price, axe_display, other_display, axe_combined, other_combined



#axe_relative_price
summary(lm(market_share~axe_relative_price, data=deo))

lag1 <- head(c(NA, deo$axe_relative_price), -1)
deo_axe_relative_price <- data.frame(market_share = deo$market_share, axe_relative_price = deo$axe_relative_price,
                                     lag1 = lag1)[-1,]
summary(lm(market_share~axe_relative_price+lag1, data=deo_axe_relative_price))

lag2 <- head(c(NA, NA, deo$axe_relative_price), -2)
deo_axe_relative_price <- data.frame(market_share = deo$market_share, axe_relative_price = deo$axe_relative_price,
                                     lag1 = lag1, lag2 = lag2)[-c(1,2),]
summary(lm(market_share~axe_relative_price+lag1+lag2, data=deo_axe_relative_price))

lag3 <- head(c(NA, NA, NA, deo$axe_relative_price), -3)
deo_axe_relative_price <- data.frame(market_share = deo$market_share, axe_relative_price = deo$axe_relative_price,
                                     lag1 = lag1, lag2 = lag2, lag3 = lag3)[-c(1,2,3),]
summary(lm(market_share~axe_relative_price+lag1+lag2+lag3, data=deo_axe_relative_price))

lag4 <- head(c(NA, NA, NA, NA, deo$axe_relative_price), -4)
deo_axe_relative_price <- data.frame(market_share = deo$market_share, axe_relative_price = deo$axe_relative_price,
                                     lag1 = lag1, lag2 = lag2, lag3 = lag3, lag4 = lag4)[-c(1,2,3,4),]
summary(lm(market_share~axe_relative_price+lag1+lag2+lag3+lag4, data=deo_axe_relative_price))

#RESULT: axe_relative_price >>> lag=2

#axe_relative_r_price
summary(lm(market_share~axe_relative_r_price, data=deo))

lag1 <- head(c(NA, deo$axe_relative_r_price), -1)
deo_axe_relative_r_price <- data.frame(market_share = deo$market_share, axe_relative_r_price = deo$axe_relative_r_price,
                                     lag1 = lag1)[-1,]
summary(lm(market_share~axe_relative_r_price+lag1, data=deo_axe_relative_r_price))

lag2 <- head(c(NA, NA, deo$axe_relative_r_price), -2)
deo_axe_relative_r_price <- data.frame(market_share = deo$market_share, axe_relative_r_price = deo$axe_relative_r_price,
                                       lag1 = lag1, lag2 = lag2)[-c(1,2),]
summary(lm(market_share~axe_relative_r_price+lag1+lag2, data=deo_axe_relative_r_price))

#RESULT: axe_relative_price >>> lag=0

#axe_display_promotion
summary(lm(market_share~axe_display, data=deo))

lag1 <- head(c(NA, deo$axe_display), -1)
deo_axe_display <- data.frame(market_share = deo$market_share, axe_display = deo$axe_display,
                                       lag1 = lag1)[-1,]
summary(lm(market_share~axe_display+lag1, data=deo_axe_display))

lag2 <- head(c(NA, NA, deo$axe_display), -2)
deo_axe_display <- data.frame(market_share = deo$market_share, axe_display = deo$axe_display,
                              lag1 = lag1, lag2 = lag2)[-c(1,2),]
summary(lm(market_share~axe_display+lag1+lag2, data=deo_axe_display))

#RESULT: axe_display >>> lag=0

#other_display_promotion
summary(lm(market_share~other_display, data=deo))

lag1 <- head(c(NA, deo$other_display), -1)
deo_other_display <- data.frame(market_share = deo$market_share, other_display = deo$other_display,
                              lag1 = lag1)[-1,]
summary(lm(market_share~other_display+lag1, data=deo_other_display))

lag2 <- head(c(NA, NA, deo$other_display), -2)
deo_other_display <- data.frame(market_share = deo$market_share, other_display = deo$other_display,
                                lag1 = lag1, lag2 = lag2)[-c(1,2),]
summary(lm(market_share~other_display+lag1+lag2, data=deo_other_display))

#RESULT: other_display >>> lag=0

#axe_combined_promotion
summary(lm(market_share~axe_combined, data=deo))

lag1 <- head(c(NA, deo$axe_combined), -1)
deo_axe_combined <- data.frame(market_share = deo$market_share, axe_combined = deo$axe_combined,
                              lag1 = lag1)[-1,]
summary(lm(market_share~axe_combined+lag1, data=deo_axe_combined))

lag2 <- head(c(NA, NA, deo$axe_combined), -2)
deo_axe_combined <- data.frame(market_share = deo$market_share, axe_combined = deo$axe_combined,
                               lag1 = lag1, lag2 = lag2)[-c(1,2),]
summary(lm(market_share~axe_combined+lag1+lag2, data=deo_axe_combined))

lag3 <- head(c(NA, NA, NA, deo$axe_combined), -3)
deo_axe_combined <- data.frame(market_share = deo$market_share, axe_combined = deo$axe_combined,
                               lag1 = lag1, lag2 = lag2, lag3 = lag3)[-c(1,2,3),]
summary(lm(market_share~axe_combined+lag1+lag2+lag3, data=deo_axe_combined))

#RESULT: axe_combined >>> lag=1

#other_combined_promotion
summary(lm(market_share~other_combined, data=deo))

lag1 <- head(c(NA, deo$other_combined), -1)
deo_other_combined <- data.frame(market_share = deo$market_share, other_combined = deo$other_combined,
                               lag1 = lag1)[-1,]
summary(lm(market_share~other_combined+lag1, data=deo_other_combined))

lag2 <- head(c(NA, NA, deo$other_combined), -2)
deo_other_combined <- data.frame(market_share = deo$market_share, other_combined = deo$other_combined,
                                 lag1 = lag1, lag2 = lag2)[-c(1,2),]
summary(lm(market_share~other_combined+lag1+lag2, data=deo_other_combined))

lag3 <- head(c(NA, NA, NA, deo$other_combined), -3)
deo_other_combined <- data.frame(market_share = deo$market_share, other_combined = deo$other_combined,
                                 lag1 = lag1, lag2 = lag2, lag3 = lag3)[-c(1,2,3),]
summary(lm(market_share~other_combined+lag1+lag2+lag3, data=deo_other_combined))

#RESULT: other_combined >>> lag=0

#INCORPORATING LAGS

axe_relative_price_lag1 <- head(c(NA, deo$axe_relative_price), -1)
axe_relative_price_lag2 <- head(c(NA, NA, deo$axe_relative_price), -2)
axe_combined_lag1 <- head(c(NA, deo$axe_combined), -1)

deo_dynamic <- cbind(deo, axe_relative_price_lag1, axe_relative_price_lag2, axe_combined_lag1)

dynamic_model <- lm(market_share ~ axe_display + other_display + axe_combined + other_combined + 
                            axe_relative_price + axe_relative_r_price + axe_relative_price_lag1, data=deo_dynamic[-1,])
summary(dynamic_model)

dynamic_model <- lm(market_share ~ axe_display + other_display + axe_combined + other_combined + 
                            axe_relative_price + axe_relative_r_price + axe_relative_price_lag1 + axe_combined_lag1, data=deo_dynamic[-1,])
summary(dynamic_model)

# MULTICOLLINEARITY & HETEROSCEDASTICITY

library(car)
vif(dynamic_model)
cor.test(deo_dynamic$axe_combined, deo_dynamic$axe_combined_lag1)

dynamic_res <- dynamic_model$residuals
t <- length(dynamic_res)
r <- t/3
k <- length(dynamic_model$coefficients)
sub1 <- dynamic_res[1:((t-r)/2)]
sub2 <- tail(dynamic_res, ((t-r)/2))
middle <- dynamic_res[(((t-r)/2)+1):(t-r)]

sum(dynamic_res == append(sub1, append(middle, sub2))) == length(dynamic_res)

sse1 <- sum(sub1**2); sse2 <- sum(sub2**2)
test_stat <- sse2/sse1
dof <- ((t-r)/2)-k
1-pf(test_stat, dof, dof)

# CROSS VALIDATION, PERFORMANCE COMPARISON

library(caret)
library(mlr)

#Static
set.seed(123) 
train.control <- trainControl(method = "cv", number = 5)
static_cross_validation <- train(static_model$call$formula, data = deo, method = "lm", trControl = train.control)
print(static_cross_validation)    

#Dynamic
set.seed(123) 
train.control <- trainControl(method = "cv", number = 5)
dynamic_cross_validation <- train(market_share ~ axe_display + other_display + axe_combined + other_combined + 
                                          axe_relative_price + axe_relative_r_price + axe_relative_price_lag1 + 
                                          axe_combined_lag1, data = deo_dynamic[-1,], method = "lm", trControl = train.control)
print(dynamic_cross_validation)


dynamic_model$call$formula


### FOLLOWER - LEADER ###
library(readxl)
follower <- read_xlsx("Deodorant.xlsx", sheet="Follower")

#Adding 1-lags
follower$axe_relative_price_lag <- head(c(NA, follower$axe_relative_price), -1)
follower$rexona_relative_price_lag <- head(c(NA, follower$rexona_relative_price), -1)
follower$eightfour_relative_price_lag <- head(c(NA, follower$eightfour_relative_price), -1)

follower$axe_combined_lag <- head(c(NA, follower$axe_combined), -1)
follower$rexona_combined_lag <- head(c(NA, follower$rexona_combined), -1)
follower$eightfour_combined_lag <- head(c(NA, follower$eightfour_combined), -1)

# AXE vs. REXONA
summary(lm(rexona_relative_price~axe_relative_price, data=follower))
summary(lm(log(rexona_relative_price)~log(axe_relative_price), data=follower))

summary(lm(rexona_combined~axe_combined, data=follower))
summary(lm(log(rexona_combined+1)~log(axe_combined+1), data=follower))

#Reaction matrix
axe_rexona_reaction_matrix <- matrix(data=c(1, cor(follower$rexona_relative_price,follower$axe_relative_price), cor(follower$axe_combined,follower$axe_relative_price), cor(follower$rexona_combined,follower$axe_relative_price),
                                           cor(follower$axe_relative_price, follower$rexona_relative_price), 1, cor(follower$axe_combined, follower$rexona_relative_price), cor(follower$rexona_combined, follower$rexona_relative_price),
                                           cor(follower$axe_relative_price, follower$axe_combined), cor(follower$rexona_relative_price, follower$axe_combined), 1, cor(follower$rexona_combined, follower$axe_combined),
                                           cor(follower$axe_relative_price, follower$rexona_combined), cor(follower$rexona_relative_price,follower$rexona_combined), cor(follower$axe_combined, follower$rexona_combined), 1),
                                     nrow=4, ncol=4, byrow=T)


# AXE vs. 8X4
summary(lm(eightfour_relative_price~axe_relative_price, data=follower))
summary(lm(log(eightfour_relative_price)~log(axe_relative_price), data=follower))

summary(lm(eightfour_relative_price~axe_relative_price, data=follower[c(1:60),]))
summary(lm(log(eightfour_relative_price)~log(axe_relative_price), data=follower[c(1:60),]))

summary(lm(eightfour_relative_price~axe_relative_price, data=follower[c(61:124),]))
summary(lm(log(eightfour_relative_price)~log(axe_relative_price), data=follower[c(61:124),]))


summary(lm(eightfour_combined~axe_combined, data=follower))
summary(lm(log(eightfour_combined+1)~log(axe_combined+1), data=follower))




















#Reaction matrix
axe_8x4_reaction_matrix <- matrix(data=c(1, cor(follower$eightfour_relative_price,follower$axe_relative_price), cor(follower$axe_combined,follower$axe_relative_price), cor(follower$eightfour_combined,follower$axe_relative_price),
                                         cor(follower$axe_relative_price, follower$eightfour_relative_price), 1, cor(follower$axe_combined, follower$eightfour_relative_price), cor(follower$eightfour_combined, follower$eightfour_relative_price),
                                         cor(follower$axe_relative_price, follower$axe_combined), cor(follower$eightfour_relative_price, follower$axe_combined), 1, cor(follower$eightfour_combined, follower$axe_combined),
                                         cor(follower$axe_relative_price, follower$eightfour_combined), cor(follower$eightfour_relative_price,follower$eightfour_combined), cor(follower$axe_combined, follower$eightfour_combined), 1),
                                  nrow=4, ncol=4, byrow=T)

eightfour_axe_reaction_matrix <- matrix(data=c(1, cor(follower$axe_relative_price,follower$eightfour_relative_price), cor(follower$eightfour_combined,follower$eightfour_relative_price), cor(follower$axe_combined,follower$eightfour_relative_price),
                                         cor(follower$eightfour_relative_price, follower$axe_relative_price), 1, cor(follower$eightfour_combined, follower$axe_relative_price), cor(follower$axe_combined, follower$axe_relative_price),
                                         cor(follower$eightfour_relative_price, follower$eightfour_combined), cor(follower$axe_relative_price, follower$eightfour_combined), 1, cor(follower$axe_combined, follower$eightfour_combined),
                                         cor(follower$eightfour_relative_price, follower$axe_combined), cor(follower$axe_relative_price,follower$axe_combined), cor(follower$eightfour_combined, follower$axe_combined), 1),
                                  nrow=4, ncol=4, byrow=T)
