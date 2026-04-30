library(forecast)
library(ggplot2)
library(BTSR)
source("Funcoes.R")
source("FitRUR.R")
source("best.RUR.R")
######################
## Data preparation ##
######################
data <- readr::read_delim("https://raw.githubusercontent.com/emanueleg/lora-rssi/master/vineyard-2021_data/combined_hourly_data.csv", 
                          delim = ";", escape_double = FALSE, trim_ws = TRUE) |> 
  dplyr::mutate(timestamp=as.POSIXct(timestamp, tz="GMT",
                                     origin="1970-01-01 00:00:00"),
                hum=hum/100) # available at https://github.com/emanueleg/lora-rssi/blob/master/vineyard-2021_data/combined_hourly_data.csv

data<-data[830:1870,] # only winter
n<-round(dim(data)[1]*.8)

#########################
## Train and test sets ##
#########################
datatrain<-cbind(data[1:n,])  
colnames(datatrain)<-paste0(colnames(datatrain),"_train")

datatest<-cbind(data[(n+1):(dim(data)[1]),])  
colnames(datatest)<-paste0(colnames(datatest),"_test")

suppressMessages(attach(datatrain))
suppressMessages(attach(datatest))
suppressMessages(attach(data))
X<-as.matrix(datatrain[,2])
Xtest<-as.matrix(datatest[,2])
X0<-rbind(X,Xtest)
nX<-dim(X)[2]
######################
## stacionary tests ##
######################
truncation<-round(12*(n/100)^(.25)) # Schwert rule
adf.level<-tseries::adf.test(hum_train,k=truncation) # stacionary
adf.diff<-tseries::adf.test(diff(hum_train),k=truncation) # stacionary
kpss.level<-tseries::kpss.test(hum_train,lshort = F) # non-stacionary
kpss.diff<-tseries::kpss.test(diff(hum_train),lshort = F) # stacionary

table.stationarity<-data.frame(
  Series=c("In Level","1st difference"),
  ADF=c(adf.level$statistic,adf.diff$statistic),
  `p-value ADF`=c(adf.level$p.value,adf.diff$p.value),
  KPSS=c(kpss.level$statistic,kpss.diff$statistic),
  `p-value KPSS`=c(kpss.level$p.value,kpss.diff$p.value)
)

# ACF and PACF of the differenced time series
acf((hum_train), main ="")
pacf((hum_train), main ="")
acf(diff(hum_train), main ="")
pacf(diff(hum_train), main ="")

########################
## Fitting the models ##
########################
# fiting the ARIMA
a01<-auto.arima(hum_train, allowdrift = F)
new1<-Arima(hum_test,model=a01) #one-step-ahead
a02<-auto.arima(hum_train, xreg = X,allowdrift = F)
new2<-Arima(hum_test,xreg = Xtest,model=a02) #one-step-ahead

# fiting the SARIMA
s01<-auto.arima(ts(hum_train,frequency = 24), allowdrift = F)
snew1<-Arima(hum_test,model=s01) #one-step-ahead
s02<-auto.arima(ts(hum_train,frequency = 24), xreg = X,allowdrift = F)
snew2<-Arima(hum_test,xreg = Xtest,model=s02) #one-step-ahead

# fiting the unit ARMA
quant<-.5
order<-matrix(NA,16,8)
cont<-1
for(i in 0:3){
  for(j in 0:3){
    
    #tive que setar yt = hum_train, o pacote nao estava chamando
    
    # ajustes não convergiam adicionado um try para os erros
    
    barma<-try(summary(BARFIMA.fit(yt=hum_train,p=i,d=F,q=j,info=T,
                                   report=F)),silent=T)
    if(inherits(barma,"try-error") || is.null(barma$aic) || is.na(barma$aic)) barma<-list(aic=Inf)
    
    karma1<-try(suppressWarnings(KARFIMA.fit(yt=hum_train,p=i,d=F,q=j,info=T,
                                             rho=quant,
                                             control = list(method="Nelder-Mead",stopcr=1e-2),
                                             report=F)),silent=T)
    if(inherits(karma1,"try-error")){
      karma<-list(aic=Inf)
    }else{
      karma<-try(summary(karma1),silent=T)
      if(inherits(karma,"try-error") || is.null(karma1$convergence) ||
         is.na(karma1$convergence) || karma1$convergence!=0 ||
         is.null(karma$aic) || is.na(karma$aic) || is.nan(karma$aic)) karma$aic=Inf
    }
    
    uwarma1<-try(suppressWarnings(UWARFIMA.fit(yt=hum_train,p=i,d=F,q=j,info=T,rho=quant,
                                               report=F)),silent=T)
    if(inherits(uwarma1,"try-error")){
      uwarma<-list(aic=Inf)
    }else{
      uwarma<-try(summary(uwarma1),silent=T)
      if(inherits(uwarma,"try-error") || is.null(uwarma1$convergence) ||
         is.na(uwarma1$convergence) || uwarma1$convergence!=0 ||
         is.null(uwarma$aic) || is.na(uwarma$aic) || is.nan(uwarma$aic)) uwarma$aic=Inf
    }
    
    barmax<-try(summary(BARFIMA.fit(yt=hum_train,p=i,d=F,q=j,info=T,
                                    xreg = X,
                                    report=F)),silent=T)
    if(inherits(barmax,"try-error") || is.null(barmax$aic) || is.na(barmax$aic)) barmax<-list(aic=Inf)
    
    karmax1<-try(suppressWarnings(KARFIMA.fit(yt=hum_train,p=i,d=F,q=j,info=T,
                                              xreg = X,rho=quant,
                                              control = list(method="Nelder-Mead",stopcr=1e-2),
                                              report=F)),silent=T)
    if(inherits(karmax1,"try-error")){
      karmax<-list(aic=Inf)
    }else{
      karmax<-try(summary(karmax1),silent=T)
      if(inherits(karmax,"try-error") || is.null(karmax1$convergence) ||
         is.na(karmax1$convergence) || karmax1$convergence!=0 ||
         is.null(karmax$aic) || is.na(karmax$aic) || is.nan(karmax$aic)) karmax$aic=Inf
    }
    
    uwarmax1<-try(suppressWarnings(UWARFIMA.fit(yt=hum_train,p=i,d=F,q=j,info=T,
                                                rho=quant, xreg = X,
                                                report=F)),silent=T)
    if(inherits(uwarmax1,"try-error")){
      uwarmax<-list(aic=Inf)
    }else{
      uwarmax<-try(summary(uwarmax1),silent=T)
      if(inherits(uwarmax,"try-error") || is.null(uwarmax1$convergence) ||
         is.na(uwarmax1$convergence) || uwarmax1$convergence!=0 ||
         is.null(uwarmax$aic) || is.na(uwarmax$aic) || is.nan(uwarmax$aic)) uwarmax$aic=Inf
    }
    
    #   print(c(karma1$convergence,karma$aic))
    order[cont,]<-c(i,j,barma$aic,karma$aic,uwarma$aic,
                    barmax$aic,karmax$aic,uwarmax$aic)
    cont<-cont+1
  }
}
order<-order[-1,]

rur_best <- best.RUR(hum_train,
                     sf = c(start = c(2020,356*24), frequency = 24*366),
                     pmax = 3, qmax = 3,
                     h=length(hum_test),
                     nbest = 1)

# rurx_best <- best.RUR(hum_train,
#                      sf = c(start = c(2020,356*24), frequency = 24*366),
#                      pmax = 3, qmax = 3,
#                      h=length(hum_test),
#                      nbest = 1,X=X,X_hat = 0)

# selecting the order of each class
orbarma<-order[which.min(order[,3]),c(1:2)]
orkarma<-order[which.min(order[,4]),c(1:2)]
oruwarma<-order[which.min(order[,5]),c(1:2)]
orbarmax<-order[which.min(order[,6]),c(1:2)]
orkarmax<-order[which.min(order[,7]),c(1:2)]
oruwarmax<-order[which.min(order[,8]),c(1:2)]
orrurarma<-rur_best$best_order

names_rows<-c("BARMAX","KARMAX","UWARMAX",
              "ARIMAX",
              "BARMA","KARMA","UWARMA",
              "RUR-ARMA","SARIMA")

barmax<-BARFIMA.fit(yt=hum_train,p=orbarmax[1],d=F,q=orbarmax[2],
                    xreg=X,info=T,report=F)
karmax<-KARFIMA.fit(yt=hum_train,p=orkarmax[1],d=F,q=orkarmax[2],rho=quant,
                    control = list(method="Nelder-Mead",stopcr=1e-2),
                    xreg=X,info=T,report=F)
uwarmax<-UWARFIMA.fit(yt=hum_train,p=oruwarmax[1],d=F,q=oruwarmax[2],rho=quant,
                      xreg=X,info=T,report=F)

if(orrurarma[1]==0) ar_rur<-NA else ar_rur<-1:orrurarma[1]
if(orrurarma[2]==0) ma_rur<-NA else ma_rur<-1:orrurarma[2]
rurarma<-RURarma.fit(ts(hum_train,start=c(2020,356*24),frequency=24*366),
                     ar=ar_rur,ma=ma_rur,
                     tau=quant,link="logit",h=length(hum_test))

barma<-BARFIMA.fit(yt=hum_train,p=orbarma[1],d=F,q=orbarma[2],
                   info=T,report=F)
karma<-KARFIMA.fit(yt=hum_train,p=orkarma[1],d=F,q=orkarma[2],rho=quant,
                   control = list(method="Nelder-Mead",stopcr=1e-2),
                   info=T,report=F)
uwarma<-UWARFIMA.fit(yt=hum_train,p=oruwarma[1],d=F,q=oruwarma[2],rho=quant,
                     info=T,report=F)

barmax_coeff<-t(summary(barmax)$coefficients[,c(1,4)])
karmax_coeff<-t(summary(karmax)$coefficients[,c(1,4)])
uwarmax_coeff<-t(summary(uwarmax)$coefficients[,c(1,4)])
sarimax_coeff<-t(lmtest::coeftest(s02)[,c(1,4)])
arimax_coeff<-t(lmtest::coeftest(a02)[,c(1,4)])
barma_coeff<-t(summary(barma)$coefficients[,c(1,4)])
karma_coeff<-t(summary(karma)$coefficients[,c(1,4)])
uwarma_coeff<-t(summary(uwarma)$coefficients[,c(1,4)])
rurarma_coeff<-t(rurarma$model[,c(1,4)])
sarima_coeff<-t(lmtest::coeftest(s01)[,c(1,4)])
arima_coeff<-t(lmtest::coeftest(a01)[,c(1,4)])

(round(barmax_coeff,4))
(round(karmax_coeff,4))
(round(uwarmax_coeff,4))
(round(sarimax_coeff[,c(ncol(sarimax_coeff),1:(ncol(sarimax_coeff)-1))],2))
(round(arimax_coeff[,c(ncol(arimax_coeff),1:(ncol(arimax_coeff)-1))],4))
(round(barma_coeff,4))
(round(karma_coeff,4))
(round(uwarma_coeff, 4))
(round(rurarma_coeff,4))
(round(sarima_coeff,2))
(round(arima_coeff,4))

# Box.test BARMA
Box.test(barmax$residuals, lag=20, fitdf = sum(orbarmax))
# Box.test KARMA
Box.test(karmax$residuals, lag=20, fitdf = sum(orkarmax))
# Box.test BARMA
Box.test(uwarmax$residuals, lag=20, fitdf = sum(oruwarmax))
# Box.test ARIMAX
Box.test(a02$residuals, lag=20, fitdf = 4)
# Box.test BARMA
Box.test(barma$residuals, lag=20, fitdf = sum(orbarma))
# Box.test KARMA
Box.test(karma$residuals, lag=20, fitdf = sum(orkarma))
# Box.test BARMA
Box.test(uwarma$residuals, lag=20, fitdf = sum(oruwarma))
# Box.test RUR-ARMA
Box.test(rurarma$residuals, lag=20, fitdf = sum(orrurarma))
# Box.test ARIMA
Box.test(s01$residuals, lag=20, fitdf = 4)

# ACF and PACF of the residuals
acf(barmax$residuals, main ="")
pacf(barmax$residuals, main ="")
acf(karmax$residuals, main ="")
pacf(karmax$residuals, main ="")
acf(uwarmax$residuals, main ="")
pacf(uwarmax$residuals, main ="")
acf(rurarma$residuals, main ="")
pacf(rurarma$residuals, main ="")
acf(a02$residuals, main ="")
pacf(a02$residuals, main ="")

RURARMA.extract<-function(yt, alpha, phi=NULL, theta=NULL,
                          ar=NA, ma=NA, link="logit")
{
  yt<-as.numeric(yt)
  if(min(yt)<=0 || max(yt)>=1) stop("OUT OF RANGE (0,1)!")
  
  linktemp <- substitute(link)
  if(!is.character(linktemp)){
    linktemp <- deparse(linktemp)
    if(linktemp == "link") linktemp <- eval(link)
  }
  if(linktemp %in% c("logit","probit","cloglog")) stats <- make.link(linktemp)
  else stop("link not available")
  
  linkfun<-stats$linkfun
  linkinv<-stats$linkinv
  nfull<-length(yt)
  ynew<-linkfun(yt)
  
  if(any(is.na(ar))) {p1<-0; ar<-1; phi<-0} else {p1<-length(ar); phi<-as.numeric(phi)}
  if(any(is.na(ma))) {q1<-0; ma<-1; theta<-0} else {q1<-length(ma); theta<-as.numeric(theta)}
  
  m<-max(ifelse(p1==0,0,max(ar)),ifelse(q1==0,0,max(ma)))
  etahat<-rep(NA,nfull)
  errorhat<-rep(0,nfull)
  
  for(i in (m+1):nfull)
  {
    etahat[i]<-alpha + as.numeric(phi%*%ynew[i-ar]) + as.numeric(theta%*%errorhat[i-ma])
    errorhat[i]<-ynew[i]-etahat[i]
  }
  
  mut<-linkinv(etahat)
  return(list(mut=mut, eta=etahat, residuals=errorhat))
}

barmax_out<-BARFIMA.extract(yt=hum, xreg = X0,  
                            coefs = list(alpha = barmax$coefficients[1], 
                                         beta = barmax$coefficients[2:(nX+1)],
                                         phi= if(orbarmax[1]==0) NULL else barmax$coefficients[(nX+2):(orbarmax[1]+nX+1)], 
                                         theta = if(orbarmax[2]==0) {NULL} else{
                                           barmax$coefficients[(orbarmax[1]+(nX+2)):(orbarmax[1]+nX+1+orbarmax[2])]},
                                         nu = barmax$coefficients[(orbarmax[1]+(nX+2)+orbarmax[2])])
)
karmax_out<-KARFIMA.extract(yt=hum,xreg = X0,rho=quant,  
                            coefs = list(alpha = karmax$coefficients[1], 
                                         beta = karmax$coefficients[2:(nX+1)],
                                         phi= if(orkarmax[1]==0) NULL else karmax$coefficients[(nX+2):(orkarmax[1]+nX+1)], 
                                         theta = if(orkarmax[2]==0) {NULL} else{
                                           karmax$coefficients[(orkarmax[1]+nX+2):(orkarmax[1]+nX+1+orkarmax[2])]},
                                         nu = karmax$coefficients[(orkarmax[1]+nX+2+orkarmax[2])])
)
uwarmax_out<-UWARFIMA.extract(yt=hum,xreg = X0,rho=quant,  
                              coefs = list(alpha = uwarmax$coefficients[1],
                                           beta = uwarmax$coefficients[2:(nX+1)],
                                           phi= if(oruwarmax[1]==0) NULL else uwarmax$coefficients[(nX+2):(oruwarmax[1]+nX+1)],
                                           theta = if(oruwarmax[2]==0) {NULL} else{
                                             uwarmax$coefficients[(oruwarmax[1]+nX+2):(oruwarmax[1]+nX+1+oruwarmax[2])]},
                                           nu = uwarmax$coefficients[(oruwarmax[1]+nX+2+oruwarmax[2])])
)
barma_out<-BARFIMA.extract(yt=hum,
                           coefs = list(alpha = barma$coefficients[1], 
                                        phi= if(orbarma[1]==0) NULL else barma$coefficients[2:(orbarma[1]+1)], 
                                        theta = if(orbarma[2]==0) {NULL} else{
                                          barma$coefficients[(orbarma[1]+2):(orbarma[1]+1+orbarma[2])]},
                                        nu = barma$coefficients[(orbarma[1]+2+orbarma[2])])
)
karma_out<-KARFIMA.extract(yt=hum,
                           coefs = list(alpha = karma$coefficients[1], 
                                        phi= if(orkarma[1]==0) NULL else karma$coefficients[2:(orkarma[1]+1)],
                                        theta = if(orkarma[2]==0) {NULL} else{
                                          karma$coefficients[(orkarma[1]+2):(orkarma[1]+1+orkarma[2])]},
                                        nu = karma$coefficients[(orkarma[1]+2+orkarma[2])])
)
uwarma_out<-UWARFIMA.extract(yt=hum,
                             coefs = list(alpha = uwarma$coefficients[1], 
                                          phi= if(oruwarma[1]==0) NULL else uwarma$coefficients[2:(oruwarma[1]+1)],
                                          theta = if(oruwarma[2]==0) {NULL} else{
                                            uwarma$coefficients[(oruwarma[1]+2):(oruwarma[1]+1+oruwarma[2])]},
                                          nu = uwarma$coefficients[(oruwarma[1]+2+oruwarma[2])])
)

rurarma_out<-RURARMA.extract(yt=hum,
                             alpha=rurarma$alpha,
                             phi=if(orrurarma[1]==0) NULL else rurarma$phi,
                             theta=if(orrurarma[2]==0) NULL else rurarma$theta,
                             ar=if(orrurarma[1]==0) NA else 1:orrurarma[1],
                             ma=if(orrurarma[2]==0) NA else 1:orrurarma[2],
                             link="logit")

results_outsample<-rbind(
  forecast::accuracy(barmax_out$mut[(n+1):(dim(data)[1])], hum_test),
  forecast::accuracy(karmax_out$mut[(n+1):(dim(data)[1])], hum_test),
  forecast::accuracy(uwarmax_out$mut[(n+1):(dim(data)[1])], hum_test),
  forecast::accuracy(new2$fitted, hum_test),
  forecast::accuracy(barma_out$mut[(n+1):(dim(data)[1])], hum_test),
  forecast::accuracy(karma_out$mut[(n+1):(dim(data)[1])], hum_test),
  forecast::accuracy(uwarma_out$mut[(n+1):(dim(data)[1])], hum_test),
  forecast::accuracy(rurarma_out$mut[(n+1):(dim(data)[1])], hum_test),
  forecast::accuracy(snew1$fitted, hum_test)
)[,c(3,2,5)]

row.names(results_outsample)<-
  names_rows

# Round results
rounded_result <- round(results_outsample, 4)

# Percentage differences
df_percent_diff <- sweep(
  sweep(rounded_result, 2, rounded_result[3, ], FUN = "-"),
  2,
  rounded_result[3, ],
  FUN = "/"
)[-3, ] * 100

# Data frame
df <- data.frame(
  values = as.vector(df_percent_diff),
  model = rep(names_rows[-3], 3),
  measure = c(rep("MAE", 8), rep("MAPE", 8), rep("RMSE", 8))
)

model_order <- c("SARIMA", "ARIMAX","BARMA", "BARMAX",
                 "KARMA", "KARMAX","RUR-ARMA", "UWARMA")

df$model <- factor(df$model, levels = model_order)

# Colors
colors <- c(
  "#00008B", "gray20", "#2B36F8", "grey45", "#566CF2",
  "grey65", "#81A2EC", "grey80"
)

# w1=4.5
# h11=4
# setEPS()
# postscript("comparision1.eps",width = w1, height = h11,family = "Times")
ggplot(df, aes(y = values, x = measure, fill = factor(model))) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9) +
  labs(fill = "", y = "Percentage differences", x = "") +
  scale_fill_manual(
    values = colors,
    labels = c("SARIMA", "ARIMAX",expression(bold(beta * ARMA)),
               expression(bold(beta * ARMAX)),
               "KARMA", "KARMAX","RUR-ARMA", "UWARMA"),
    guide = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  geom_text(
    aes(x = measure, y = values, label = round(values, 2)),
    fontface = "bold",
    hjust = 0,
    vjust = ifelse(df$values >= 0, -0.3, -0.9),
    angle = 55,
    position = position_dodge(width = 0.9),
    color = "black",
    size = 2
  ) +
  coord_cartesian(ylim=c(min(-1,min(df$values,na.rm=T)*1.05),
                         max(df$values,na.rm=T)*1.15)) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 8),
    plot.title = element_text(face = "bold", size = 8),
    legend.text = element_text(face = "bold", size = 8),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.x = unit(0.3, "cm"),
    legend.margin = margin(t = -13, l = -30, unit = "pt"),
    axis.title.y = element_text(face = "bold", color = "black", size = 8),
    axis.title.x = element_text(face = "bold", color = "black", size = 8),
    axis.text.x = element_text(face = "bold", color = "black", size = 8),
    axis.text.y = element_text(face = "bold", color = "black", size = 8),
    panel.background = element_rect(fill = "white", colour = "white")
  )
# dev.off()

hum1=xts::xts(hum_test, order.by=data$timestamp[(n+1):(dim(data)[1])])
hum2<-xts::xts(new2$fitted, order.by=data$timestamp[(n+1):(dim(data)[1])])#
hum3<-xts::xts(uwarmax_out$mut[(n+1):(dim(data)[1])], order.by=data$timestamp[(n+1):(dim(data)[1])])
hum4<-xts::xts(rurarma_out$mut[(n+1):(dim(data)[1])], order.by=data$timestamp[(n+1):(dim(data)[1])])
# w1=4.5
# h11=2.5
# setEPS()
# postscript("comparision2.eps",width = w1, height = h11,family = "Times")
plot(hum1,
     main="", yaxis.right=FALSE, grid.col = "white",
     format.labels="%b-%Y", main.timespan = FALSE,
     lwd=.5,
     ylim=c(min(hum1,hum2,hum3,hum4),
            max(hum1,hum2,hum3,hum4)+.1)
)
lines(hum2,lty=2,lwd=1,col=2)
lines(hum3,lwd=.5,col=4)
lines(hum4,lwd=.5,col=3)
xts::addLegend("topleft",
               legend.names =c("Original data","ARIMAX","UWARMAX","RUR-ARMA"),
               col=c(1,2,4,3), cex=.8,lty=c(1,2,1,1),
               lwd=c(.5,1,.5,.5),
               ncol=4
)
# dev.off()
