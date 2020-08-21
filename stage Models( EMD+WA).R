#####################################################################################
               3-STAGE     (EMD-CEEMDAN-MM,EMD-EEMD-MM,WA-CEEMDAN-MM)
               2-STAGE     (WA-ARIMA,WA-RGMDH,WA-RBFN,EMD-ARIMA,EMD-RGMDH,EMD-RBFN)
               1 STAGE     (ARIMA)
#####################################################################################
library(wmtsa)
library(wavelets)
library(waveslim)
library(wavethresh)
library(Rwave)
library(ggplot2)
library(Rlibeemd)
library(forecast)
library(TSA)
library(EbayesThresh)
library(Rlibeemd)
Time=seq(as.Date("2005/07/1"), as.Date("2019/6/30"), by = "months")
length(Time)

Data_used_in_analysis <- read.csv(file.choose())

gas1=Data_used_in_analysis$gas

w=gas1
length(w)
set.seed(321)
#########################################################
                EMD-CEEMDAN-MM   3-Stage   
#########################################################
set.seed(321)
gas_emd=emd(w, num_imfs = 0, S_number = 40, num_siftings = 50L)



#####    Soft thresholding   #########

IMF1=gas_emd[,1]
E1=mad(IMF1,center=0,constant=0.6745)
T1=0.6*sqrt(2*E1*log(168))
thresh_soft_IMF1=ifelse(abs(IMF1)>T1,(sign(IMF1)*(abs(IMF1)-T1)),0)
#OR
value=(sign(IMF1)*(abs(IMF1)-T1))
thresh_soft_IMF1=ifelse(abs(IMF1)>T1,value,0)

IMF2=gas_emd[,2]
E2=mad(IMF2,center=0,constant=0.6745)
T2=0.6*sqrt(2*E2*log(168))
thresh_soft_IMF2=ifelse(abs(IMF2)>T2,(sign(IMF2)*(abs(IMF2)-T2)),0)

IMF3=gas_emd[,3]
E3=mad(IMF3,center=0,constant=0.6745)
T3=0.6*sqrt(2*E3*log(168))
thresh_soft_IMF3=ifelse(abs(IMF3)>T3,(sign(IMF3)*(abs(IMF3)-T3)),0)

IMF4=gas_emd[,4]
E4=mad(IMF4,center=0,constant=0.6745)
T4=0.6*sqrt(2*E4*log(168))
thresh_soft_IMF4=ifelse(abs(IMF4)>T4,(sign(IMF4)*(abs(IMF4)-T4)),0)

IMF5=gas_emd[,5]
E5=mad(IMF5,center=0,constant=0.6745)
T5=0.6*sqrt(2*E5*log(168))
thresh_soft_IMF5=ifelse(abs(IMF5)>T5,(sign(IMF5)*(abs(IMF5)-T5)),0)

IMF6=gas_emd[,6]
E6=mad(IMF6,center=0,constant=0.6745)
T6=0.6*sqrt(2*E6*log(168))
thresh_soft_IMF6=ifelse(abs(IMF6)>T6,(sign(IMF6)*(abs(IMF6)-T6)),0)

RES=gas_emd[,7]



imfs=cbind(IMF1,IMF2,thresh_soft_IMF3,thresh_soft_IMF4,thresh_soft_IMF5,IMF6,RES)## no need to apply thresghold on last 2 imfs bcz they already smooth
               
imfs_mat_soft=as.matrix(imfs)
gas_rec_imfs_soft=rowSums(imfs_mat_soft)

plot(Time,w,type="o",col=1)
length(w)
length(Time)
lines(Time,gas_rec_imfs_soft,col=4)
legend("topright", c("gas", "Denoised by EMD"), cex=0.75, fill=c(1,4))

###########performance measures###################
##Mean Absolute  Deviation (MAD)##
MAD_gas_emd_2=1/168*(sum(abs(w-gas_rec_imfs_soft)))
MAD_gas_emd_2  	### 64.56036
##mean absolute percentage error (MAPE)##
MAPE_gas_emd_2=1/168*(sum(abs((w-gas_rec_imfs_soft)/w)))
MAPE_gas_emd_2	###  0.0005273928
##MS error (MSE)##  
MSE_gas_emd_2=1/168*(sum(abs(w-gas_rec_imfs_soft)^2))
MSE_gas_emd_2   ### 5849.804

###denoiseeffectiveness###
mean_o=mean(w)
mean_o  ##  122355
sd_o=sd(w)
sd_o    ##5477.39
mean_denoised_soft=mean(gas_rec_imfs_soft)
mean_denoised_soft ##122352.4
sd_denoised_soft=sd(gas_rec_imfs_soft)
sd_denoised_soft  ## 5452.545
aa=w^2
bb=gas_rec_imfs_soft^2
STN_soft=10*log10(sum(aa)/sum(aa-bb)^2)
STN_soft ## -39.49446
write.csv(gas_rec_imfs_soft,file="denoised.semd.gas.gas.csv")


######################CEEMDAN##################
gas_flow=ceemdan(gas_rec_imfs_soft, num_imfs = 7, ensemble_size = 900, noise_strength = 2.5,
                           S_number =40 , num_siftings = 100, rng_seed = 0L, threads = 0L)
plot(gas_flow)
gas_ceemd=gas_flow

IMFs_ceemd=gas_ceemd[,1:4]
plot(IMFs_ceemd)
IMFs_ceemd_irr=gas_ceemd[,5:6]
plot(IMFs_ceemd_irr)
Trend_ceemd=gas_ceemd[,7]
plot(Trend_ceemd)


######################
IMF11=gas_ceemd[,1]
IMF21=gas_ceemd[,2]
IMF31=gas_ceemd[,3]
IMF41=gas_ceemd[,4]
IMF51=gas_ceemd[,5]
IMF61=gas_ceemd[,6]
IMF71=gas_ceemd[,7]
######################IMF1 prediction########################
library(GMDH)
gas_gm_IMF11 = fcast(IMF11, method = "GMDH", input = 5, layer = 1, f.number = 5,
                    level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                      0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
                                                                      
length(gas_gm_IMF11)          
length(IMF11-gas_gm_IMF11$fitted)
##Mean Absolute  Deviation (MAD)##
MAD_gas_ceemd_1=1/163*(sum(abs(IMF11-gas_gm_IMF11$fitted)))
MAD_gas_ceemd_1  	###  1622.949
##mean absolute percentage error (MAPE)##
MAPE_gas_ceemd_2=1/163*(sum(abs((IMF11-gas_gm_IMF11$fitted)/IMF11)))
MAPE_gas_ceemd_2	###  0.8326521
##MS error (MSE)##  
MSE_gas_ceemd_3=1/163*(sum(abs(IMF11-gas_gm_IMF11$fitted)^2))
MSE_gas_ceemd_3       ###  4882711

#####################################RGMDH IMF11###############
library(GMDH)
gas_IMF11 = fcast(IMF11, method = "RGMDH", input = 3, layer = 1, f.number = 5,
                    level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                      0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_ceemd_1=1/165*(sum(abs(IMF11-gas_IMF11$fitted)))
MAD_gas_ceemd_1  	### 1666.069
##mean absolute percentage error (MAPE)
MAPE_gas_ceemd_2=1/165*(sum(abs((IMF11-gas_IMF11$fitted)/IMF11)))
MAPE_gas_ceemd_2	### 0.9020908
##MS error (MSE)##  
MSE_gas_ceemd_3=1/165*(sum(abs(IMF11-gas_IMF11$fitted)^2))
MSE_gas_ceemd_3       ###  5308985  



############################ Neural network ###############
library(neuralnet)
library(RSNNS)
gas_ts=ts(IMF11,frequency=12,start=c(2005,6,1))
axis11=embed(IMF11,2)
axis1=axis11[,-1]
response1=axis11[,1]
set.seed(21)
rbf.model <- rbf(response1, axis1, size = 5, maxit = 1000, linOut = FALSE)
error=response1-fitted(rbf.model)

##Mean Absolute  Deviation (MAD)
MAD_gas_rbf_1=1/168*(sum(abs(error)))
MAD_gas_rbf_1  	### 2950.58
 ##mean absolute percentage error (MAPE)
MAPE_gas_rbf_2=1/168*(sum(abs((error)/response1)))
MAPE_gas_rbf_2	###0.9939282
##MS error (MSE)##  
MSE_gas_rbf_3=1/168*(sum(abs(error)^2))
MSE_gas_rbf_3   ### 10995874

#########################   ARMA  ####################
gas_IMF11=auto.arima(IMF11)
summary(gas_IMF11)
plot(forecast(gas_IMF11,h=30))
points(1:168,fitted(gas_IMF11),type="l",col="green")

##Mean Absolute  Deviation (MAD)##
MAD_gas_arma1=1/168*(sum(abs(IMF11- gas_IMF11$fitted)))
MAD_gas_arma1  	 ###  1596.854
##mean absolute percentage error (MAPE)##
MAPE_gas_arma12=1/168*(sum(abs((IMF11-gas_IMF11$fitted)/IMF11)))
MAPE_gas_arma12	 ### 0.7765151
 ##MS error (MSE)##
MSE_gas_arma13=1/168*(sum(abs(IMF11- gas_IMF11$fitted)^2))
MSE_gas_arma13   ## 4547034

pred_1=gas_IMF11$fitted
##############################IMF21 prediction########################
gas_gm_IMF21 = fcast(IMF21, method = "GMDH", input = 5, layer = 1, f.number = 5,
                    level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                      0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)##
MAD_gas_gm_21=1/163*(sum(abs(IMF21-gas_gm_IMF21$fitted)))
MAD_gas_gm_21  	## 395.967
##mean absolute percentage error (MAPE)##
MAPE_gas_gm_22=1/163*(sum(abs((IMF21-gas_gm_IMF21$fitted)/IMF21)))
MAPE_gas_gm_22	### 1.262854
##MS error (MSE)##  
MSE_gas_gm_32=1/163*(sum(abs(IMF21-gas_gm_IMF21$fitted)^2))
MSE_gas_gm_32       ### 253273.2

#####################################RGMDH IMF21###############
library(GMDH)
gas_rgm_IMF21 = fcast(IMF21, method = "RGMDH", input = 3, layer = 1, f.number = 5,
                        level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                          0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_rgm_21=1/165*(sum(abs(IMF21-gas_rgm_IMF21$fitted)))
MAD_gas_rgm_21  	###  418.5598
##mean absolute percentage error (MAPE)
MAPE_gas_rgm_22=1/165*(sum(abs((IMF21-gas_rgm_IMF21$fitted)/IMF21)))
MAPE_gas_rgm_22	      ### 1.376165
##MS error (MSE)##  
MSE_gas_rgm_23=1/165*(sum(abs(IMF21-gas_rgm_IMF21$fitted)^2))
MSE_gas_rgm_23       ### 287091.7 


###########################Neural network###############
library(neuralnet)
library(RSNNS)
gas_ts=ts(IMF21,frequency=12,start=c(2005,7,1))
axis211=embed(IMF21,2)
axis21=axis211[,-1]
response21=axis211[,1]
set.seed(21)
rbf.model_21 <- rbf(response21, axis21, size = 5, maxit = 1000, linOut = TRUE)
error21=response21-fitted(rbf.model_21)
##Mean Absolute  Deviation (MAD)
MAD_gas_rbf_21=1/168*(sum(abs(error21)))
MAD_gas_rbf_21  	### 794.8186
##mean absolute percentage error (MAPE)
MAPE_gas_rbf_22=1/168*(sum(abs((error21)/response21)))
MAPE_gas_rbf_22	      ### 1.048927
##MS error (MSE)##  
MSE_gas_rbf_22=1/168*(sum(abs(error21)^2))
MSE_gas_rbf_22        ###  988459.5

#########################ARMA####################
gas_arma_IMF21=auto.arima(IMF21)
summary(gas_arma_IMF21)
plot(forecast(gas_arma_IMF21,h=30))
points(1:length(IMF21),fitted(gas_arma_IMF21),type="l",col="green")

##Mean Absolute  Deviation (MAD)##
MAD_gas_arma21=1/168*(sum(abs(IMF21- gas_arma_IMF21$fitted)))
MAD_gas_arma21  	 ### 366.961
##mean absolute percentage error (MAPE)##
MAPE_gas_arma22=1/168*(sum(abs((IMF21- gas_arma_IMF21$fitted)/IMF21)))
MAPE_gas_arma22	       ###  1.178967
##MS error (MSE)##
MSE_gas_arma23=1/168*(sum(abs(IMF21- gas_arma_IMF21$fitted)^2))
MSE_gas_arma23         ### 227199
pred_2=gas_arma_IMF21$fitted
##################################IM31 PREDICTION#################
library(GMDH)
gas_gm_IMF31 = fcast(IMF31, method = "GMDH", input = 5, layer = 1, f.number = 5,
                    level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                      0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_gm_31=1/163*(sum(abs(IMF31-gas_gm_IMF31$fitted)))
MAD_gas_gm_31  	    ### 272.2637
##mean absolute percentage error (MAPE)
MAPE_gas_gm_32=1/163*(sum(abs((IMF31-gas_gm_IMF31$fitted)/IMF31)))
MAPE_gas_gm_32	    ### 1.394959
##MS error (MSE)  
MSE_gas_gm_33=1/163*(sum(abs(IMF31-gas_gm_IMF31$fitted)^2))

MSE_gas_gm_33       ### 130925.3

#####################################RGMDH IMF31###############
library(GMDH)
gas_rgm_IMF31 = fcast(IMF31, method = "RGMDH", input =3, layer = 1, f.number = 5,
                        level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                          0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_rgm_31=1/165*(sum(abs(IMF31-gas_rgm_IMF31$fitted)))
MAD_gas_rgm_31  	### 241.332
##mean absolute percentage error (MAPE)
MAPE_gas_rgm_32=1/165*(sum(abs((IMF31-gas_rgm_IMF31$fitted)/IMF31)))
MAPE_gas_rgm_32	      ### 1.488133
##MS error (MSE)##  
MSE_gas_rgm_33=1/165*(sum(abs(IMF31-gas_rgm_IMF31$fitted)^2))
MSE_gas_rgm_33       ###  98360.26

############################Neural network###############
library(neuralnet)
library(RSNNS)
gas_ts=ts(IMF31,frequency=12,start=c(2005,7,1))
axis311=embed(IMF31,2)
axis31=axis311[,-1]
response31=axis311[,1]
set.seed(21)
rbf.model_31 <- rbf(response31, axis31, size = 5, maxit = 1000, linOut = TRUE)
error31=response31-fitted(rbf.model_31)
##Mean Absolute  Deviation (MAD)
MAD_gas_rbf_31=1/168*(sum(abs(error31)))
MAD_gas_rbf_31  	### 720.5983
##mean absolute percentage error (MAPE)
MAPE_gas_rbf_32=1/168*(sum(abs((error31)/response31)))
MAPE_gas_rbf_32	      ### 1.122901
##MS error (MSE)##  
MSE_gas_rbf_32=1/168*(sum(abs(error31)^2))
MSE_gas_rbf_32        ###   806456.2
#########################ARMA####################
gas_arma_IMF31=auto.arima(IMF31)
summary(gas_arma_IMF31)
plot(forecast(gas_arma_IMF31,h=30))
points(1:length(IMF31),fitted(gas_arma_IMF31),type="l",col="green")

##Mean Absolute  Deviation (MAD)##
MAD_gas_arma31=1/168*(sum(abs(IMF31- gas_arma_IMF31$fitted)))
MAD_gas_arma31  	 ### 187.4474
##mean absolute percentage error (MAPE)##
MAPE_gas_arma32=1/168*(sum(abs((IMF31- gas_arma_IMF31$fitted)/IMF31)))
MAPE_gas_arma32	       ###  1.133286
##MS error (MSE)##
MSE_gas_arma33=1/168*(sum(abs(IMF31- gas_arma_IMF31$fitted)^2))
MSE_gas_arma33         ###  63268.81
pred_3=gas_arma_IMF31$fitted################
##################################IM41 PREDICTION#################
library(GMDH)
gas_gm_IMF41 = fcast(IMF41, method = "GMDH", input = 5, layer = 1, f.number = 5,
                    level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                      0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_gm_41=1/163*(sum(abs(IMF41-gas_gm_IMF41$fitted)))
MAD_gas_gm_41  	    ###  110.8701
##mean absolute percentage error (MAPE)
MAPE_gas_gm_42=1/163*(sum(abs((IMF41-gas_gm_IMF41$fitted)/IMF41)))
MAPE_gas_gm_42	    ### 0.2781001
##MS error (MSE)  
MSE_gas_gm_43=1/163*(sum(abs(IMF41-gas_gm_IMF41$fitted)^2))
MSE_gas_gm_43       ###    17489.48


#####################################RGMDH IMF41###############
library(GMDH)
gas_rgm_IMF41 = fcast(IMF41, method = "RGMDH", input = 3, layer = 1, f.number = 5,
                        level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                          0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_rgm_41=1/165*(sum(abs(IMF41-gas_rgm_IMF41$fitted)))
MAD_gas_rgm_41  	###  84.70092
##mean absolute percentage error (MAPE)
MAPE_gas_rgm_42=1/165*(sum(abs((IMF41-gas_rgm_IMF41$fitted)/IMF41)))
MAPE_gas_rgm_42	      ### 0.209808
##MS error (MSE)##  
MSE_gas_rgm_43=1/165*(sum(abs(IMF41-gas_rgm_IMF41$fitted)^2))
MSE_gas_rgm_43        ###  9949.821

############################Neural network###############
library(neuralnet)
library(RSNNS)
gas_ts=ts(IMF41,frequency=12,start=c(2005,7,1))
axis411=embed(IMF41,2)
axis41=axis411[,-1]
response41=axis411[,1]
set.seed(21)
rbf.model_41 <- rbf(response41, axis41, size = 5, maxit = 1000, linOut = TRUE)


error41=response41-fitted(rbf.model_41)
##Mean Absolute  Deviation (MAD)
MAD_gas_rbf_41=1/168*(sum(abs(error41)))
MAD_gas_rbf_41  	###  975.213
##mean absolute percentage error (MAPE)
MAPE_gas_rbf_42=1/168*(sum(abs((error41)/response41)))
MAPE_gas_rbf_42	      ###1.011031
##MS error (MSE)##  
MSE_gas_rbf_42=1/168*(sum(abs(error41)^2))
MSE_gas_rbf_42        ###   1525335
#########################ARMA####################
gas_arma_IMF41=auto.arima(IMF41)
summary(gas_arma_IMF41)
plot(forecast(gas_arma_IMF41,h=30))
points(1:length(IMF41),fitted(gas_arma_IMF41),type="l",col="green")

##Mean Absolute  Deviation (MAD)##
MAD_gas_arma41=1/168*(sum(abs(IMF41- gas_arma_IMF41$fitted)))

MAD_gas_arma41  	 ### 45.65331
##mean absolute percentage error (MAPE)##
MAPE_gas_arma42=1/168*(sum(abs((IMF41- gas_arma_IMF41$fitted)/IMF41)))
MAPE_gas_arma42	       ###  0.1139251
##MS error (MSE)##
MSE_gas_arma43=1/168*(sum(abs(IMF41-gas_arma_IMF41$fitted)^2))
MSE_gas_arma43         ### 3443.68
pred_4=gas_arma_IMF41$fitted
##################################IM51 PREDICTION#################
library(GMDH)
gas_gm_IMF51 = fcast(IMF51, method = "GMDH", input = 5, layer = 1, f.number = 5,
                    level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                      0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_gm_51=1/163*(sum(abs(IMF51-gas_gm_IMF51$fitted)))
MAD_gas_gm_51  	    ### 33.33644
##mean absolute percentage error (MAPE)
MAPE_gas_gm_52=1/163*(sum(abs((IMF51-gas_gm_IMF51$fitted)/IMF51)))
MAPE_gas_gm_52	    ###0.266158
##MS error (MSE)  
MSE_gas_gm_53=1/163*(sum(abs(IMF51-gas_gm_IMF51$fitted)^2))
MSE_gas_gm_53       ###  1799.853
#####################################RGMDH IMF51###############
library(GMDH)
gas_rgm_IMF51 = fcast(IMF51, method = "RGMDH", input = 3, layer = 1, f.number = 5,
                        level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                          0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_rgm_51=1/165*(sum(abs(IMF51-gas_rgm_IMF51$fitted)))
MAD_gas_rgm_51  	### 14.47725
##mean absolute percentage error (MAPE)
MAPE_Indus_rgm_52=1/165*(sum(abs((IMF51-gas_rgm_IMF51$fitted)/IMF51)))
MAPE_Indus_rgm_52	      ###  0.07457351
##MS error (MSE)##  
MSE_gas_rgm_53=1/165*(sum(abs(IMF51-gas_rgm_IMF51$fitted)^2))
MSE_gas_rgm_53        ### 349.9352

############################Neural network###############
library(neuralnet)
library(RSNNS)
axis511=embed(IMF51,2)
axis51=axis511[,-1]
response51=axis511[,1]
set.seed(21)
rbf.model_51 <- rbf(response51, axis51, size = 5, maxit = 1000, linOut = TRUE)
error51=response51-fitted(rbf.model_51)
##Mean Absolute  Deviation (MAD)
MAD_gas_rbf_51=1/168*(sum(abs(error51)))
MAD_gas_rbf_51  	###  532.5535

##mean absolute percentage error (MAPE)
MAPE_gas_rbf_52=1/168*(sum(abs((error51)/response51)))
MAPE_gas_rbf_52	      ### 0.9756661
##MS error (MSE)##  
MSE_gas_rbf_52=1/168*(sum(abs(error51)^2))
MSE_gas_rbf_52        ### 599190.3
#########################ARMA####################
gas_arma_IMF51=auto.arima(IMF51)
summary(gas_arma_IMF51)
plot(forecast(gas_arma_IMF51,h=30))
points(1:length(IMF51),fitted(gas_arma_IMF51),type="l",col="green")

##Mean Absolute  Deviation (MAD)##
MAD_gas_arma51=1/168*(sum(abs(IMF51- gas_arma_IMF51$fitted)))
MAD_gas_arma51  	 ### 6.806156

##mean absolute percentage error (MAPE)##
MAPE_gas_arma52=1/168*(sum(abs((IMF51- gas_arma_IMF51$fitted)/IMF51)))
MAPE_gas_arma52	       ###   0.04222504

##MS error (MSE)##
MSE_gas_arma53=1/168*(sum(abs(IMF51- gas_arma_IMF51$fitted)^2))
MSE_gas_arma53         ###  75.87722
pred_5=gas_arma_IMF51$fitted##################
##################################IM61 PREDICTION#################
library(GMDH)
gas_IMF61 = fcast(IMF61, method = "GMDH", input = 5, layer = 1, f.number = 5,
                  level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                    0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_gm_61=1/163*(sum(abs(IMF61-gas_IMF61$fitted)))
MAD_gas_gm_61  	    ### 16.1942
##mean absolute percentage error (MAPE)
MAPE_gas_gm_62=1/163*(sum(abs((IMF61-gas_IMF61$fitted)/IMF61)))

MAPE_gas_gm_62	    ### 0.0839542
##MS error (MSE)  
MSE_gas_gm_63=1/163*(sum(abs(IMF61-gas_IMF61$fitted)^2))

MSE_gas_gm_63       ###  442.9472
#####################################RGMDH IMF61###############
library(GMDH)
gas_rgm_IMF61 = fcast(IMF61, method = "RGMDH", input = 3, layer = 1, f.number = 5,
                      level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                        0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_rgm_61=1/163*(sum(abs(IMF61-gas_rgm_IMF61$fitted)))
MAD_gas_rgm_61  	### 4.097449
##mean absolute percentage error (MAPE)
MAPE_gas_rgm_62=1/163*(sum(abs((IMF61-gas_rgm_IMF61$fitted)/IMF61)))
MAPE_gas_rgm_62	      ### 0.03492298
##MS error (MSE)##  
MSE_gas_rgm_63=1/163*(sum(abs(IMF61-gas_rgm_IMF61$fitted)^2))
MSE_gas_rgm_63        ### 24.8032
pred_6=gas_rgm_IMF61$fitted########
############################Neural network###############
library(neuralnet)
library(RSNNS)
gas_ts=ts(IMF61,frequency=12,start=c(2005,7,1))
axis611=embed(IMF61,2)
axis61=axis611[,-1]
response61=axis611[,1]
set.seed(21)
rbf.model_61 <- rbf(response61, axis61, size = 20, maxit = 1000, linOut = TRUE)
error61=response61-fitted(rbf.model_61)
##Mean Absolute  Deviation (MAD)
MAD_gas_rbf_61=1/168*(sum(abs(error61)))
MAD_gas_rbf_61  	###   477.388
##mean absolute percentage error (MAPE)
MAPE_gas_rbf_62=1/168*(sum(abs((error61)/response61)))
MAPE_gas_rbf_62	      ###  1.389289
##MS error (MSE)##  
MSE_gas_rbf_62=1/168*(sum(abs(error61)^2))
MSE_gas_rbf_62        ###  440144.1 

#########################ARMA####################

gas_arma_IMF61=auto.arima(IMF61)
summary(gas_arma_IMF61)
plot(forecast(gas_arma_IMF61,h=30))
points(1:length(IMF61),fitted(gas_arma_IMF61),type="l",col="green")

##Mean Absolute  Deviation (MAD)
MAD_gas_arma61=1/168*(sum(abs(IMF61- gas_arma_IMF61$fitted)))
MAD_gas_arma61  	 ### 1.541011
##mean absolute percentage error (MAPE)
MAPE_gas_arma62=1/168*(sum(abs((IMF61- gas_arma_IMF61$fitted)/IMF61)))
MAPE_gas_arma62	       ### 1.330071e-05
##MS error (MSE)##
MSE_gas_arma63=1/168*(sum(abs(IMF61- gas_arma_IMF61$fitted)^2))
MSE_gas_arma63         ### 154.1213


##################################IM71 PREDICTION#################
library(GMDH)
gas_IMF71 = fcast(IMF71, method = "GMDH", input = 5, layer = 1, f.number = 5,
                  level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                    0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_gm_71=1/163*(sum(abs(IMF71-gas_IMF71$fitted)))
MAD_gas_gm_71  	    ###11.12131
##mean absolute percentage error (MAPE)
MAPE_gas_gm_72=1/163*(sum(abs((IMF71-gas_IMF71$fitted)/IMF71)))

MAPE_gas_gm_72	    ###  9.097413e-05
##MS error (MSE)  
MSE_gas_gm_73=1/163*(sum(abs(IMF71-gas_IMF71$fitted)^2))

MSE_gas_gm_73       ###  195.4106
#####################################RGMDH IMF71###############
library(GMDH)
gas_rgm_IMF71 = fcast(IMF71, method = "RGMDH", input = 3, layer = 1, f.number = 5,
                      level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                        0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_rgm_71=1/165*(sum(abs(IMF71-gas_rgm_IMF71$fitted)))
MAD_gas_rgm_71  	### 9.663011
##mean absolute percentage error (MAPE)
MAPE_gas_rgm_72=1/165*(sum(abs((IMF71-gas_rgm_IMF71$fitted)/IMF71)))
MAPE_gas_rgm_72	      ### 7.836102e-05
##MS error (MSE)##  
MSE_gas_rgm_73=1/165*(sum(abs(IMF71-gas_rgm_IMF71$fitted)^2))
MSE_gas_rgm_73        ### 167.6817
############################Neural network###############
library(neuralnet)
library(RSNNS)
gas_ts=ts(IMF71,frequency=12,start=c(2005,7,1))
axis711=embed(IMF71,2)
axis71=axis711[,-1]
response71=axis711[,1]
set.seed(21)
rbf.model_71 <- rbf(response71, axis71, size = 20, maxit = 1000, linOut = TRUE)
error71=response71-fitted(rbf.model_71)
##Mean Absolute  Deviation (MAD)
MAD_gas_rbf_71=1/168*(sum(abs(error71)))
MAD_gas_rbf_71  	###   22792.35
##mean absolute percentage error (MAPE)
MAPE_gas_rbf_72=1/168*(sum(abs((error71)/response71)))
MAPE_gas_rbf_72	      ### 0.1863759
##MS error (MSE)##  
MSE_gas_rbf_72=1/168*(sum(abs(error71)^2))
MSE_gas_rbf_72        ### 533663591

#########################ARMA####################

gas_arma_IMF71=auto.arima(IMF71)
summary(gas_arma_IMF71)
plot(forecast(gas_arma_IMF71,h=30))
points(1:length(IMF71),fitted(gas_arma_IMF71),type="l",col="green")

##Mean Absolute  Deviation (MAD)
MAD_gas_arma71=1/168*(sum(abs(IMF71- gas_arma_IMF71$fitted)))
MAD_gas_arma71  	 ### 1.541011
##mean absolute percentage error (MAPE)
MAPE_gas_arma72=1/168*(sum(abs((IMF71- gas_arma_IMF71$fitted)/IMF71)))
MAPE_gas_arma72	       ###  1.330071e-05
##MS error (MSE)##
MSE_gas_arma73=1/168*(sum(abs(IMF71- gas_arma_IMF71$fitted)^2))
MSE_gas_arma73         ### 154.1213
pred_7=gas_arma_IMF71$fitted########

##############################################END PREDICTION#########################
pred_gas_1=cbind(pred_1,pred_2,pred_3,pred_4,pred_5,pred_6,pred_7)
pred_gas_sum=rowSums(pred_gas_1)
plot(Time,gas_rec_imfs_soft,type="o",col=1)
lines(Time,pred_gas_sum,col=4)
legend("topright", c("Denoised gas by EMD-CEEMDAN-MM", "Predicted gas"), cex=0.75, fill=c(1,4))
error_final=gas_rec_imfs_soft-pred_gas_sum
length(error_final)
##Mean Absolute  Deviation (MAD)
MAD_gas1=1/165*(sum(abs(error_final[-c(1:3)])))
MAD_gas1 	       ### 1194.999
##mean absolute percentage error (MAPE)
MAPE_gas2=1/165*(sum(abs((error_final[-c(1:3)])/gas_rec_imfs_soft[-c(1:3)])))
MAPE_gas2	       ###0.009735269
##MS error (MSE)
MSE_gas3=1/165*(sum(abs(error_final[-c(1:3)])^2))
MSE_gas3         ###  2336537

#################################################
                      ARIMA 1-Stage
#################################################
g2=auto.arima(w)
summary(g2)   ###(2,1,1)
error=w-g2$fitted
##Mean Absolute  Deviation (MAD)
MAD_g1=1/168*(sum(abs(error)))
MAD_g1 	       ### 3774.267
##mean absolute percentage error (MAPE)
MAPE_g2=1/168*(sum(abs((error)/w)))
MAPE_g2	       ### 0.03094788
##MS error (MSE)
MSE_g3=1/168*(sum(abs(error)^2))
MSE_g3  ##  23538742

###############################################
              EMD-ARIMA
###############################################

smooth_gas=auto.arima(gas_rec_imfs_soft)
summary(smooth_gas)   ###(1,1,1)
error_smooth=gas_rec_imfs_soft-smooth_gas$fitted
head(smooth_gas$fitted)
length(smooth_gas$fitted)
##Mean Absolute  Deviation (MAD)
MAD_gas1=1/168*(sum(abs(error_smooth)))            ##MRE
MAD_gas1 	       ###  3750.389
##mean absolute percentage error (MAPE)              ### MAE
MAPE_gas2=1/168*(sum(abs((error_smooth)/gas_rec_imfs_soft)))
MAPE_gas2	       ### 0.03074798
##MS error (MSE)
MSE_gas3=1/168*(sum(abs(error_smooth)^2))
MSE_gas3         ###23240715

####################################################
                  EMD-RGMDH
####################################################

library(GMDH)
gas_rgm_smooth = fcast(ts(gas_rec_imfs_soft), method = "RGMDH", input = 3, layer = 1, f.number = 5,
                         level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                           0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
head(gas_rgm_smooth$fitted)
##Mean Absolute  Deviation (MAD)
error_rgm_sm=ts(gas_rec_imfs_soft)-gas_rgm_smooth$fitted
MAD_gas_rgm_sm=1/165*(sum(abs(error_rgm_sm)))
MAD_gas_rgm_sm  	### 3753.873
##mean absolute percentage error (MAPE)
MAPE_gas_rgm_sm=1/165*(sum(abs((error_rgm_sm)/ts(gas_rec_imfs_soft))))
MAPE_gas_rgm_sm	      ### 0.03096583
##MS error (MSE)##  
MSE_gas_rgm_sm=1/165*(sum(abs(error_rgm_sm)^2))
MSE_gas_rgm_sm        ### 23340680

##########################################
              EMD-RBFNN
##########################################

library(neuralnet)
library(RSNNS)

axis_smooth=embed(gas_rec_imfs_soft,2)
axis_sm=axis_smooth[,-1]
response_sm=axis_smooth[,1]
set.seed(21)
rbf.model_sm <- rbf(response_sm, axis_sm, size = 20, maxit = 1000, linOut = TRUE)
error_sm=response_sm-fitted(rbf.model_sm)
head(fitted(rbf.model_sm))
length(fitted(rbf.model_sm))
##Mean Absolute  Deviation (MAD)
MAD_gas_rbf_sm=1/168*(sum(abs(error_sm)))
MAD_gas_rbf_sm  	### 22916.16
##mean absolute percentage error (MAPE)
MAPE_gas_rbf_sm=1/168*(sum(abs((error_sm)/response_sm)))
MAPE_gas_rbf_sm	      ### 0.1856141
##MS error (MSE)##  
MSE_gas_rbf_sm=1/168*(sum(abs(error_sm)^2))

MSE_gas_rbf_sm        ###556913324

##########################################################
                WA-CEEMDAN-MM
##########################################################

n.level=log10(length(w))
gas.modwt <- modwt(w, boundary="reflection", n.levels=2)
gassmooth.modwt <- ebayesthresh.wavelet(gas.modwt, vscale="level")
gassmooth.emodwt <- imodwt(gassmooth.modwt)[1:168]
error_wa=w-gassmooth.emodwt 
###denoiseeffectiveness###
mean_o=mean(w)
mean_o # 122355
sd_o=sd(w)
sd_o 	#5477.39
mean_denoised_emodwt=mean(gassmooth.emodwt)
mean_denoised_emodwt # 122355
sd_denoised_emodwt=sd(gassmooth.emodwt)
sd_denoised_emodwt   # 3757.383
a=w^2
b=gassmooth.emodwt^2
STN=10*log10(sum(w^2)/sum(a-b)^2)
STN  ## -64.45929
write.csv(gassmooth.emodwt,file="denoised.wav.gas.csv")

###########performance measures###################

##Mean Absolute  Deviation (MAD)##
MAD_gas_wa_1=1/168*(sum(abs(error_wa)))
MAD_gas_wa_1  	### 3258.046
##mean absolute percentage error (MAPE)##
MAPE_gas_wa_2=1/168*(sum(abs((error_wa)/w)))
MAPE_gas_wa_2	### 0.0267866
##MS error (MSE)  
MSE_gas_wa_3=1/168*(sum((error_wa)^2))
MSE_gas_wa_3   ### 14792170

######################CEEMDAN##################
gas_rec_imfs1=gassmooth.emodwt
gas_flow=ceemdan(gas_rec_imfs1, num_imfs = 7, ensemble_size = 900, noise_strength = 0.3,
                           S_number = 40, num_siftings = 100, rng_seed = 0L, threads = 0L)
plot(gas_flow)
gas_ceemd=gas_flow
IMFs_ceemd=gas_ceemd[,1:4]
plot(IMFs_ceemd)
IMFs_ceemd_irr=gas_ceemd[,5:6]
plot(IMFs_ceemd_irr)
Trend_ceemd=gas_ceemd[,7]
plot(Trend_ceemd)

######################
IMF11=gas_ceemd[,1]
IMF21=gas_ceemd[,2]
IMF31=gas_ceemd[,3]
IMF41=gas_ceemd[,4]
IMF51=gas_ceemd[,5]
IMF61=gas_ceemd[,6]
IMF71=gas_ceemd[,7]
######################IMF11 prediction########################
library(GMDH)
gas_gm_IMF11 = fcast(IMF11, method = "GMDH", input = 5, layer = 1, f.number = 5,
                    level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                      0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)##
MAD_gas_gm_ceemd_1=1/163*(sum(abs(IMF11-gas_gm_IMF11$fitted)))
MAD_gas_gm_ceemd_1  	###  37.20776

##mean absolute percentage error (MAPE)##
MAPE_gas_gm_ceemd_2=1/163*(sum(abs((IMF11-gas_gm_IMF11$fitted)/IMF11)))
MAPE_gas_gm_ceemd_2	###  2.166595
##MS error (MSE)##  
MSE_gas_gm_ceemd_3=1/163*(sum(abs(IMF11-gas_gm_IMF11$fitted)^2))
MSE_gas_gm_ceemd_3       ###   2176.231

#####################################RGMDH IMF11###############
library(GMDH)
gas_rgm_IMF11 = fcast(IMF11, method = "RGMDH", input = 3, layer = 1, f.number = 5,
                    level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                      0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_ceemd_1=1/65*(sum(abs(IMF11-gas_rgm_IMF11$fitted)))
MAD_gas_ceemd_1  	###  99.5001
##mean absolute percentage error (MAPE)
MAPE_gas_ceemd_2=1/165*(sum(abs((IMF11-gas_rgm_IMF11$fitted)/IMF11)))
MAPE_gas_ceemd_2	### 1.292929
##MS error (MSE)##  
MSE_gas_ceemd_3=1/165*(sum(abs(IMF11-gas_rgm_IMF11$fitted)^2))
MSE_gas_ceemd_3       ### 2863.394

############################Neural network###############
library(neuralnet)
library(RSNNS)
gas_ts=ts(IMF11,frequency=12,start=c(2005,7,1))
axis=embed(IMF11,2)
axis1=axis[,-1]
response1=axis[,1]
set.seed(21)
rbf.model <- rbf(response1, axis1, size = 5, maxit = 1000, linOut = FALSE)
error=response1-fitted(rbf.model)
##Mean Absolute  Deviation (MAD)
MAD_gas_rbf_1=1/168*(sum(abs(error)))
MAD_gas_rbf_1  	###   46.78049
##mean absolute percentage error (MAPE)
MAPE_gas_rbf_2=1/168*(sum(abs((error)/response1)))
MAPE_gas_rbf_2	###  0.9962112
##MS error (MSE)##  
MSE_gas_rbf_3=1/168*(sum(abs(error)^2))
MSE_gas_rbf_3 ###   5750.531
########################ARMA####################
gas_IMF11=auto.arima(IMF11)
summary(gas_IMF11)
plot(forecast(gas_IMF11,h=30))
points(1:length(IMF11),fitted(gas_IMF11),type="l",col="green")

##Mean Absolute  Deviation (MAD)##
MAD_gas_arma1=1/168*(sum(abs(IMF11- gas_IMF11$fitted)))
MAD_gas_arma1  	 ###   2087.989
##mean absolute percentage error (MAPE)##
MAPE_gas_arma12=1/168*(sum(abs((IMF11- gas_IMF11$fitted)/IMF11)))
MAPE_gas_arma12	 ###  186.1527
##MS error (MSE)##
MSE_gas_arma13=1/168*(sum(abs(IMF11- gas_IMF11$fitted)^2))
MSE_gas_arma13   ### 6528912
pred1=gas_IMF11$fitted
##############################IMF21 prediction########################
gas_IMF21 = fcast(IMF21, method = "GMDH", input = 5, layer = 1, f.number = 5,
                    level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                      0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)##
MAD_gas_gm_21=1/163*(sum(abs(IMF21-gas_IMF21$fitted)))
MAD_gas_gm_21  	### 57.66629
##mean absolute percentage error (MAPE)##
MAPE_gas_gm_22=1/163*(sum(abs((IMF21-gas_IMF21$fitted)/IMF21)))
MAPE_gas_gm_22	###  1.764813
##MS error (MSE)##  
MSE_gas_gm_32=1/163*(sum(abs(IMF21-gas_IMF21$fitted)^2))
MSE_gas_gm_32       ###4980.681
#####################################RGMDH IMF21 ###############
library(GMDH)
gas_rgm_IMF21 = fcast(IMF21, method = "RGMDH", input = 3, layer = 1, f.number = 5,
                        level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                          0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_rgm_21=1/165*(sum(abs(IMF21-gas_rgm_IMF21$fitted)))
MAD_gas_rgm_21  	### 57.88134
##mean absolute percentage error (MAPE)
MAPE_gas_rgm_22=1/165*(sum(abs((IMF21-gas_rgm_IMF21$fitted)/IMF21)))
MAPE_gas_rgm_22	      ###  1.734524
##MS error (MSE)##  
MSE_gas_rgm_23=1/165*(sum(abs(IMF21-gas_rgm_IMF21$fitted)^2))
MSE_gas_rgm_23       ### 5016.591

############################Neural network###############
library(neuralnet)
library(RSNNS)
gas_ts=ts(IMF21,frequency=12,start=c(2005,7,1))
axis211=embed(IMF21,2)
axis21=axis211[,-1]
response21=axis211[,1]
set.seed(21)
rbf.model_21 <- rbf(response21, axis21, size = 5, maxit = 1000, linOut = FALSE)
error21=response21-fitted(rbf.model_21)
##Mean Absolute  Deviation (MAD)
MAD_gas_rbf_21=1/168*(sum(abs(error21)))
MAD_gas_rbf_21  	### 111.3154
##mean absolute percentage error (MAPE)
MAPE_gas_rbf_22=1/168*(sum(abs((error21)/response21)))
MAPE_gas_rbf_22	      ### 0.9906184
##MS error (MSE)##  
MSE_gas_rbf_23=1/168*(sum(abs(error21)^2))
MSE_gas_rbf_23        ###  19073.35
#########################ARMA####################
gas_arma_IMF21=auto.arima(IMF21)
summary(gas_arma_IMF21)
plot(forecast(gas_arma_IMF21,h=30))
points(1:length(IMF21),fitted(gas_arma_IMF21),type="l",col="green")

##Mean Absolute  Deviation (MAD)##
MAD_gas_arma21=1/168*(sum(abs(IMF21- gas_arma_IMF21$fitted)))
MAD_gas_arma21  	 ### 40.48732
##mean absolute percentage error (MAPE)##
MAPE_gas_arma22=1/168*(sum(abs((IMF21- gas_arma_IMF21$fitted)/IMF21)))
MAPE_gas_arma22	       ### 1.041543
##MS error (MSE)##
MSE_gas_arma23=1/168*(sum(abs(IMF21- gas_arma_IMF21$fitted)^2))
MSE_gas_arma23         ###  2620.582
pred2=gas_arma_IMF21$fitted
##################################IM31 PREDICTION#################
library(GMDH)
gas_IMF31 = fcast(IMF31, method = "GMDH", input = 5, layer = 1, f.number = 5,
                    level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                      0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_gm_31=1/163*(sum(abs(IMF31-gas_IMF31$fitted)))
MAD_gas_gm_31  	    ### 259.7587
##mean absolute percentage error (MAPE)
MAPE_gas_gm_32=1/163*(sum(abs((IMF31-gas_IMF31$fitted)/IMF31)))
MAPE_gas_gm_32	    ###   3.506361
##MS error (MSE)  
MSE_gas_gm_33=1/163*(sum(abs(IMF31-gas_IMF31$fitted)^2))
MSE_gas_gm_33       ###   108923.4
#####################################RGMDH IMF31 #############
library(GMDH)
gas_rgm_IMF31 = fcast(IMF31, method = "RGMDH", input = 3, layer = 1, f.number = 5,
                        level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                          0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_rgm_31=1/165*(sum(abs(IMF31-gas_rgm_IMF31$fitted)))
MAD_gas_rgm_31  	###233.076
##mean absolute percentage error (MAPE)
MAPE_gas_rgm_32=1/165*(sum(abs((IMF31-gas_rgm_IMF31$fitted)/IMF31)))
MAPE_gas_rgm_32	      ###   3.424744
##MS error (MSE)##  
MSE_gas_rgm_33=1/165*(sum(abs(IMF31-gas_rgm_IMF31$fitted)^2))
MSE_gas_rgm_33       ###  83375.3

############################Neural network###############
library(neuralnet)
library(RSNNS)
gas_ts=ts(IMF31,frequency=12,start=c(2005,7,1))
axis311=embed(IMF31,2)
axis31=axis311[,-1]
response31=axis311[,1]
set.seed(21)
rbf.model_31 <- rbf(response31, axis31, size = 5, maxit = 1000, linOut = TRUE)
error31=response31-fitted(rbf.model_31)
##Mean Absolute  Deviation (MAD)
MAD_gas_rbf_31=1/168*(sum(abs(error31)))
MAD_gas_rbf_31  	### 795.9788
##mean absolute percentage error (MAPE)
MAPE_gas_rbf_32=1/168*(sum(abs((error31)/response31)))
MAPE_gas_rbf_32	      ### 2.160593
##MS error (MSE)##  
MSE_gas_rbf_32=1/168*(sum(abs(error31)^2))
MSE_gas_rbf_32        ### 985595.2

#########################ARMA####################
gas_arma_IMF31=auto.arima(IMF31)
summary(gas_arma_IMF31)
plot(forecast(gas_arma_IMF31,h=30))
points(1:length(IMF31),fitted(gas_arma_IMF31),type="l",col="green")

##Mean Absolute  Deviation (MAD)##
MAD_gas_arma31=1/168*(sum(abs(IMF31- gas_arma_IMF31$fitted)))
MAD_gas_arma31  	 ### 226.8366
##mean absolute percentage error (MAPE)##
MAPE_gas_arma32=1/168*(sum(abs((IMF31- gas_arma_IMF31$fitted)/IMF31)))
MAPE_gas_arma32	       ###  3.007421
##MS error (MSE)##
MSE_gas_arma33=1/168*(sum(abs(IMF31- gas_arma_IMF31$fitted)^2))
MSE_gas_arma33         ### 80334.51
pred3=gas_arma_IMF31$fitted


##################################IM41 PREDICTION#################
library(GMDH)
gas_IMF41 = fcast(IMF41, method = "GMDH", input = 5, layer = 1, f.number = 5,
                    level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                      0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_gm_41=1/163*(sum(abs(IMF41-gas_IMF41$fitted)))
MAD_gas_gm_41  	    ###  145.4758
 ##mean absolute percentage error (MAPE)
MAPE_gas_gm_42=1/163*(sum(abs((IMF41-gas_IMF41$fitted)/IMF41)))
MAPE_gas_gm_42	    ###  0.4170381
##MS error (MSE)  
MSE_gas_gm_43=1/163*(sum(abs(IMF41-gas_IMF41$fitted)^2))
MSE_gas_gm_43       ###  32928.94

#####################################RGMDH IMF41###############
library(GMDH)
gas_rgm_IMF41 = fcast(IMF41, method = "RGMDH", input = 3, layer = 1, f.number = 5,
                        level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                          0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_rgm_41=1/165*(sum(abs(IMF41-gas_rgm_IMF41$fitted)))
MAD_gas_rgm_41  	### 98.42367
##mean absolute percentage error (MAPE)
MAPE_gas_rgm_42=1/165*(sum(abs((IMF41-gas_rgm_IMF41$fitted)/IMF41)))
MAPE_gas_rgm_42	      ###  0.3162278
##MS error (MSE)##  
MSE_gas_rgm_43=1/165*(sum(abs(IMF41-gas_rgm_IMF41$fitted)^2))
MSE_gas_rgm_43        ### 14103.05
############################Neural network###############
library(neuralnet)
library(RSNNS)
gas_ts=ts(IMF41,frequency=12,start=c(2005,1,1))
axis411=embed(IMF41,2)
axis41=axis411[,-1]
response41=axis411[,1]
set.seed(21)
rbf.model_41 <- rbf(response41, axis41, size = 5, maxit = 1000, linOut = TRUE)
error41=response41-fitted(rbf.model_41)
##Mean Absolute  Deviation (MAD)
MAD_gas_rbf_41=1/168*(sum(abs(error41)))
MAD_gas_rbf_41  	### 966.1017
##mean absolute percentage error (MAPE)
MAPE_gas_rbf_42=1/168*(sum(abs((error41)/response41)))
MAPE_gas_rbf_42	      ###   1.349518
##MS error (MSE)##  
MSE_gas_rbf_42=1/168*(sum(abs(error41)^2))
MSE_gas_rbf_42        ### 1431643

#########################ARMA####################
gas_arma_IMF41=auto.arima(IMF41)
summary(gas_arma_IMF41)
plot(forecast(gas_arma_IMF41,h=30))
points(1:length(IMF41),fitted(gas_arma_IMF41),type="l",col="green")

##Mean Absolute  Deviation (MAD)##
MAD_gas_arma41=1/168*(sum(abs(IMF41- gas_arma_IMF41$fitted)))
MAD_gas_arma41  	 ### 55.92283
##mean absolute percentage error (MAPE)##
MAPE_gas_arma42=1/168*(sum(abs((IMF41- gas_arma_IMF41$fitted)/IMF41)))
MAPE_gas_arma42	       ### 0.1971141
##MS error (MSE)##
MSE_gas_arma43=1/168*(sum(abs(IMF41- gas_arma_IMF41$fitted)^2))
MSE_gas_arma43         ### 4996.979
pred4=gas_arma_IMF41$fitted

##################################IM51 PREDICTION#################
library(GMDH)
gas_IMF51 = fcast(IMF51, method = "GMDH", input = 5, layer = 1, f.number = 5,
                    level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                      0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_gm_51=1/163*(sum(abs(IMF51-gas_IMF51$fitted)))
MAD_gas_gm_51  	    ###    48.02102
##mean absolute percentage error (MAPE)
MAPE_gas_gm_52=1/163*(sum(abs((IMF51-gas_IMF51$fitted)/IMF51)))
MAPE_gas_gm_52	    ###0.2162084
##MS error (MSE)  
MSE_gas_gm_53=1/163*(sum(abs(IMF51-gas_IMF51$fitted)^2))
MSE_gas_gm_53       ###  3463.1375

#####################################RGMDH IMF51###############
library(GMDH)
gas_rgm_IMF51 = fcast(IMF51, method = "RGMDH", input =3 , layer = 1, f.number = 5,
                        level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                          0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_rgm_51=1/165*(sum(abs(IMF51-gas_rgm_IMF51$fitted)))
MAD_gas_rgm_51  	### 19.11052
##mean absolute percentage error (MAPE)
MAPE_gas_rgm_52=1/165*(sum(abs((IMF51-gas_rgm_IMF51$fitted)/IMF51)))
MAPE_gas_rgm_52	      ###   0.06981058
##MS error (MSE)##  
MSE_gas_rgm_53=1/165*(sum(abs(IMF51-gas_rgm_IMF51$fitted)^2))
MSE_gas_rgm_53        ### 616.9022
############################Neural network###############
library(neuralnet)
library(RSNNS)
axis511=embed(IMF51,2)
axis51=axis511[,-1]
response51=axis511[,1]
set.seed(22)
rbf.model_51 <- rbf(response51, axis51, size = 20, maxit = 1000, linOut = TRUE)
error51=response51-fitted(rbf.model_51)
##Mean Absolute  Deviation (MAD)
MAD_gas_rbf_51=1/168*(sum(abs(error51)))
MAD_gas_rbf_51  	### 687.6102
##mean absolute percentage error (MAPE)
MAPE_gas_rbf_52=1/168*(sum(abs((error51)/response51)))
MAPE_gas_rbf_52	      ### 1.048248
##MS error (MSE)##  
MSE_gas_rbf_52=1/168*(sum(abs(error51)^2))
MSE_gas_rbf_52        ### 1004072
set.seed(21)
#########################ARMA####################
gas_arma_IMF51=auto.arima(IMF51)
summary(gas_arma_IMF51)
plot(forecast(gas_arma_IMF51,h=30))
points(1:length(IMF51),fitted(gas_arma_IMF51),type="l",col="green")

##Mean Absolute  Deviation (MAD)##
MAD_gas_arma51=1/168*(sum(abs(IMF51- gas_arma_IMF51$fitted)))
MAD_gas_arma51  	 ### 8.50624
##mean absolute percentage error (MAPE)##
MAPE_gas_arma52=1/168*(sum(abs((IMF51- gas_arma_IMF51$fitted)/IMF51)))
MAPE_gas_arma52	       ### 0.03633302
##MS error (MSE)##
MSE_gas_arma53=1/168*(sum(abs(IMF51- gas_arma_IMF51$fitted)^2))
MSE_gas_arma53         ###  120.0723 
pred5=gas_arma_IMF51$fitted

##################################IM61 PREDICTION#################
library(GMDH)
gas_IMF61 = fcast(IMF61, method = "GMDH", input = 5, layer = 1, f.number = 5,
                    level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                      0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_gm_61=1/163*(sum(abs(IMF61-gas_IMF61$fitted)))
MAD_gas_gm_61  	    ###  23.50361
##mean absolute percentage error (MAPE)
MAPE_gas_gm_62=1/163*(sum(abs((IMF61-gas_IMF61$fitted)/IMF61)))
MAPE_gas_gm_62	    ### 0.0001895834
##MS error (MSE)  
MSE_gas_gm_63=1/163*(sum(abs(IMF61-gas_IMF61$fitted)^2))
MSE_gas_gm_63       ###  850.0206

#####################################RGMDH IMF61###############
library(GMDH)
gas_rgm_IMF61 = fcast(IMF61, method = "RGMDH", input = 3, layer = 1, f.number = 5,
                        level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                          0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_rgm_61=1/165*(sum(abs(IMF61-gas_rgm_IMF61$fitted)))
MAD_gas_rgm_61  	###   23.10029
##mean absolute percentage error (MAPE)
MAPE_gas_rgm_62=1/165*(sum(abs((IMF61-gas_rgm_IMF61$fitted)/IMF61)))
MAPE_gas_rgm_62	      ###  0.0001863473
##MS error (MSE)##  
MSE_gas_rgm_63=1/165*(sum(abs(IMF61-gas_rgm_IMF61$fitted)^2))
MSE_gas_rgm_63        ### 842.5707

############################Neural network###############
library(neuralnet)
library(RSNNS)
gas_ts=ts(IMF61,frequency=12,start=c(2005,1,1))
axis611=embed(IMF61,2)
axis61=axis611[,-1]
response61=axis611[,1]
set.seed(21)
rbf.model_61 <- rbf(response61, axis61, size = 20, maxit = 1000, linOut = TRUE)
error61=response61-fitted(rbf.model_61)
##Mean Absolute  Deviation (MAD)
MAD_gas_rbf_61=1/168*(sum(abs(error61)))
MAD_gas_rbf_61  	### 39306.91
##mean absolute percentage error (MAPE)
MAPE_gas_rbf_62=1/168*(sum(abs((error61)/response61)))
MAPE_gas_rbf_62	      ###  0.3197175
##MS error (MSE)##  
MSE_gas_rbf_62=1/168*(sum(abs(error61)^2))
MSE_gas_rbf_62        ### 1573028310

#########################ARMA####################
gas_arma_IMF61=auto.arima(IMF61)
summary(gas_arma_IMF61)
plot(forecast(gas_arma_IMF61,h=30))
points(1:length(IMF61),fitted(gas_arma_IMF61),type="l",col="green")

##Mean Absolute  Deviation (MAD)
MAD_gas_arma61=1/168*(sum(abs(IMF61- gas_arma_IMF61$fitted)))
MAD_gas_arma61  	 ###2.931013
##mean absolute percentage error (MAPE)
MAPE_gas_arma62=1/168*(sum(abs((IMF61-gas_arma_IMF61$fitted)/IMF61)))
MAPE_gas_arma62	       ### 2.498287e-05
##MS error (MSE)##
MSE_gas_arma63=1/168*(sum(abs(IMF61- gas_arma_IMF61$fitted)^2))
MSE_gas_arma63         ### 227.469
pred6=gas_arma_IMF61$fitted
##################################IM71 PREDICTION#################
library(GMDH)
gas_IMF71 = fcast(IMF71, method = "GMDH", input = 5, layer = 1, f.number = 5,
                  level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                    0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_gm_71=1/163*(sum(abs(IMF71-gas_IMF71$fitted)))
MAD_gas_gm_71  	    ### 6.75654
##mean absolute percentage error (MAPE)
MAPE_gas_gm_72=1/163*(sum(abs((IMF71-gas_IMF71$fitted)/IMF71)))

MAPE_gas_gm_72	    ###  0.01265406
##MS error (MSE)  
MSE_gas_gm_73=1/163*(sum(abs(IMF71-gas_IMF71$fitted)^2))

MSE_gas_gm_73       ###  112.0789
#####################################RGMDH IMF71###############
library(GMDH)
gas_rgm_IMF71 = fcast(IMF71, method = "RGMDH", input = 3, layer = 1, f.number = 5,
                      level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                        0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_rgm_71=1/165*(sum(abs(IMF71-gas_rgm_IMF71$fitted)))
MAD_gas_rgm_71  	### 7.657317
##mean absolute percentage error (MAPE)
MAPE_gas_rgm_72=1/165*(sum(abs((IMF71-gas_rgm_IMF71$fitted)/IMF71)))
MAPE_gas_rgm_72	      ###0.01407951
##MS error (MSE)##  
MSE_gas_rgm_73=1/165*(sum(abs(IMF71-gas_rgm_IMF71$fitted)^2))
MSE_gas_rgm_73        ### 144.7766
############################Neural network###############
library(neuralnet)
library(RSNNS)
gas_ts=ts(IMF71,frequency=12,start=c(2005,7,1))
axis711=embed(IMF71,2)
axis71=axis711[,-1]
response71=axis711[,1]
set.seed(21)
rbf.model_71 <- rbf(response71, axis71, size = 20, maxit = 1000, linOut = TRUE)
error71=response71-fitted(rbf.model_71)
##Mean Absolute  Deviation (MAD)
MAD_gas_rbf_71=1/168*(sum(abs(error71)))
MAD_gas_rbf_71  	###   477.1664
##mean absolute percentage error (MAPE)
MAPE_gas_rbf_72=1/168*(sum(abs((error71)/response71)))
MAPE_gas_rbf_72	      ### 0.4960022
##MS error (MSE)##  
MSE_gas_rbf_72=1/168*(sum(abs(error71)^2))
MSE_gas_rbf_72        ### 584087.1

#########################ARMA####################

gas_arma_IMF71=auto.arima(IMF71)
summary(gas_arma_IMF71)
plot(forecast(gas_arma_IMF71,h=30))
points(1:length(IMF71),fitted(gas_arma_IMF71),type="l",col="green")

##Mean Absolute  Deviation (MAD)
MAD_gas_arma71=1/168*(sum(abs(IMF71- gas_arma_IMF71$fitted)))
MAD_gas_arma71  	 ### 6.349155
##mean absolute percentage error (MAPE)
MAPE_gas_arma72=1/168*(sum(abs((IMF71- gas_arma_IMF71$fitted)/IMF71)))
MAPE_gas_arma72	       ###  0.01204036
##MS error (MSE)##
MSE_gas_arma73=1/168*(sum(abs(IMF71- gas_arma_IMF71$fitted)^2))
MSE_gas_arma73         ### 94.92911
pred7=gas_arma_IMF71$fitted########

##############################################END PREDICTION#########################
 
pred_gas=cbind(pred1,pred2,pred3,pred4,pred5,pred6,pred7)
pred_gas_sum=rowSums(pred_gas)
plot(Time,gas_rec_imfs1,type="o",col=1)
lines(Time,pred_gas_sum,col=4)
legend("topright", c("Denoised gas", "Predicted gas"), cex=0.75, fill=c(1,4))
error_final=gas_rec_imfs1-pred_gas_sum
length(error_final)
##Mean Absolute  Deviation (MAD)
MAD_gas1=1/168*(sum(abs(error_final)))
MAD_gas1 	       ###  158.4641
##mean absolute percentage error (MAPE)
MAPE_gas2=1/168*(sum(abs((error_final)/gas_rec_imfs1)))
MAPE_gas2	       ### 0.001297095
##MS error (MSE)
MSE_gas3=1/168*(sum(abs(error_final)^2))
MSE_gas3         ### 38472.38
 

#################################################################
                                 3- stage Model 
                                   EMD-EEMD-MM
#################################################################

Time=seq(as.Date("2005/07/1"), as.Date("2019/6/30"), by = "months")
length(Time)

Data_used_in_analysis <- read.csv(file.choose())

gas1=Data_used_in_analysis$gas

w=gas1

length(w)
set.seed(321)

##################Hard thresholding######################

gas_emd=emd(w, num_imfs = 0, S_number = 40, num_siftings = 50L)

plot(gas_emd)
IMFs=gas_emd[,1:5]
plot(IMFs)
Trend=gas_emd[,6]
length(Trend)
plot(Trend)
length(Trend)
length(IMFs)
gas=IMFs+Trend
R_emd_1=w-gas
##Mean Absolute  Deviation (MAD)##
MAD_gas_emd_1=1/168*(sum(abs(w-gas)))
MAD_gas_emd_1  	### 610832.9
##mean absolute percentage error (MAPE)##
MAPE_gas_emd_1=1/168*(sum(abs((w-gas)/w)))
MAPE_gas_emd_1	### 4.993811
##MS error (MSE)##
MSE_gas_emd_1=1/168*(sum(abs(w-gas)^2))
MSE_gas_emd_1   ###74755081289

####Hard thresholding####
IMF1=gas_emd[,1]
E1=mad(IMF1,center=0,constant=0.6745)
T1=0.6*sqrt(2*E1*log(168))
T12=0.7*sqrt(2*E1*log(168))
T13=0.8*sqrt(2*E1*log(168))
thresh_IMF1=ifelse(abs(IMF1)>T1,IMF1,0)

IMF2=gas_emd[,2]
E2=mad(IMF2,center=0,constant=0.6745)
T2=0.6*sqrt(2*E2*log(168))
#T22=0.7*sqrt(2*E2*log(168))
#T23=0.8*sqrt(2*E2*log(168))
thresh_IMF2=ifelse(abs(IMF2)>T2,IMF2,0)

IMF3=gas_emd[,3]
E3=mad(IMF3,center=0,constant=0.6745)
T3=0.6*sqrt(2*E3*log(168))
#T32=0.7*sqrt(2*E3*log(168))
#T33=0.8*sqrt(2*E3*log(168))
thresh_IMF3=ifelse(abs(IMF3)>T3,IMF3,0)

IMF4=gas_emd[,4]
E4=mad(IMF4,center=0,constant=0.6745)
T4=0.6*sqrt(2*E4*log(168))
#T42=0.7*sqrt(2*E4*log(168))
#T43=0.8*sqrt(2*E4*log(168))
thresh_IMF4=ifelse(abs(IMF4)>T4,IMF4,0)

IMF5=gas_emd[,5]
E5=mad(IMF5,center=0,constant=0.6745)
T5=0.6*sqrt(2*E5*log(168))
#T52=0.7*sqrt(2*E5*log(168))
#T53=0.8*sqrt(2*E5*log(168))
thresh_IMF5=ifelse(abs(IMF5)>T5,IMF5,0)

IMF6=gas_emd[,6]
E6=mad(IMF6,center=0,constant=0.6745)
T6=0.6*sqrt(2*E6*log(168))
#T62=0.7*sqrt(2*E6*log(168))
#T63=0.8*sqrt(2*E6*log(168))
thresh_IMF6=ifelse(abs(IMF5)>T6,IMF6,0)

RES=gas_emd[,7]

imfs=cbind(IMF1,IMF2,thresh_IMF3,thresh_IMF4,thresh_IMF5,IMF6,RES)
           
imfs_mat=as.matrix(imfs)
gas_rec_imfs=rowSums(imfs_mat)
par("mar"=c(1,1,1,1))
plot(Time,w,type="o",col=1)
length(gas_rec_imfs)
length(Time)
lines(Time,gas_rec_imfs,col=4)
legend("topright", c("gas ", "Denoised by EMD"), cex=0.75, fill=c(1,4))
###########performance measures###################
##Mean Absolute  Deviation (MAD)##
MAD_gas_emd_1=1/168*(sum(abs(w-gas_rec_imfs)))
MAD_gas_emd_1  	###  2.293892
##mean absolute percentage error (MAPE)##
MAPE_gas_emd_1=1/168*(sum(abs((w-gas_rec_imfs)/w)))
MAPE_gas_emd_1	### 
##MS error (MSE)##  1.915927e-05
MSE_gas_emd_1=1/168*(sum(abs(w-gas_rec_imfs)^2))
MSE_gas_emd_1   ### 71.5795
###denoiseeffectiveness###
mean_o=mean(w)
mean_o #122355
sd_o=sd(w)
sd_o 	#5477.39
mean_denoised=mean(gas_rec_imfs)
mean_denoised #122355.2
sd_denoised=sd(gas_rec_imfs)
sd_denoised   # 5475.99
a=w^2
b=gas_rec_imfs^2
STN=10*log10(sum(w^2)/sum(a-b)^2)
STN  ##-9.643599
################now work on denoised series################

gas_ceemd=eemd(gas_rec_imfs, num_imfs = 7, ensemble_size = 900, noise_strength = 2.5,
                 S_number = 8, num_siftings = 50L, rng_seed = 0L, threads = 0L)
IMFs_ceemd=gas_ceemd[,1:4]
plot(IMFs_ceemd)
IMFs_ceemd_irr=gas_ceemd[,5:6]
plot(IMFs_ceemd_irr)
Trend_ceemd=gas_ceemd[,7]
plot(Trend_ceemd)

######################
IMF11=gas_ceemd[,1]
IMF21=gas_ceemd[,2]
IMF31=gas_ceemd[,3]
IMF41=gas_ceemd[,4]
IMF51=gas_ceemd[,5]
IMF61=gas_ceemd[,6]
IMF71=gas_ceemd[,7]
#####################
######################IMF1 prediction########################
library(GMDH)
gas_IMF11 = fcast(IMF11, method = "GMDH", input = 5, layer = 1, f.number = 5,
                    level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                      0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)##
MAD_gas_ceemd_1=1/163*(sum(abs(IMF11-gas_IMF11$fitted)))
MAD_gas_ceemd_1  	### 1563.32
##mean absolute percentage error (MAPE)##
MAPE_gas_ceemd_2=1/163*(sum(abs((IMF11-gas_IMF11$fitted)/IMF11)))
MAPE_gas_ceemd_2	### 1.025371
##MS error (MSE)##  
MSE_gas_ceemd_3=1/163*(sum(abs(IMF11-gas_IMF11$fitted)^2))
MSE_gas_ceemd_3       ###  4255894
#####################################RGMDH IMF11###############
library(GMDH)
gas_IMF11 = fcast(IMF11, method = "RGMDH", input = 3, layer = 1, f.number = 5,
                    level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                      0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_ceemd_1=1/165*(sum(abs(IMF11-gas_IMF11$fitted)))
MAD_gas_ceemd_1  	### 1527.932
##mean absolute percentage error (MAPE)
MAPE_gas_ceemd_2=1/165*(sum(abs((IMF11-gas_IMF11$fitted)/IMF11)))
MAPE_gas_ceemd_2	### 0.9231777
##MS error (MSE)##  
MSE_gas_ceemd_3=1/165*(sum(abs(IMF11-gas_IMF11$fitted)^2))
MSE_gas_ceemd_3       ### 4030778
pred11=gas_IMF11$fitted#################
############################Neural network###############
library(neuralnet)
library(RSNNS)
gas_ts=ts(IMF11,frequency=12,start=c(2005,7,1))
axis=embed(IMF11,2)
axis1=axis[,-1]
response=axis[,1]
set.seed(21)
rbf.model <- rbf(response, axis1, size = 5, maxit = 1000, linOut = FALSE)
error=response-fitted(rbf.model)
##Mean Absolute  Deviation (MAD)
MAD_gas_rbf_1=1/168*(sum(abs(error)))
MAD_gas_rbf_1  	### 2907.041
##mean absolute percentage error (MAPE)
MAPE_gas_rbf_2=1/168*(sum(abs((error)/response)))
MAPE_gas_rbf_2	### 0.9940712
##MS error (MSE)##  
MSE_gas_rbf_3=1/168*(sum(abs(error)^2))
MSE_gas_rbf_3   ### 11756843

#########################ARMA####################
gas_IMF11=auto.arima(IMF11)
summary(gas_IMF11)

##Mean Absolute  Deviation (MAD)##
MAD_gas_arma1=1/168*(sum(abs(IMF11- gas_IMF11$fitted)))
MAD_gas_arma1  	 ### 1335.33
##mean absolute percentage error (MAPE)##
MAPE_gas_arma12=1/168*(sum(abs((IMF11- gas_IMF11$fitted)/IMF11)))
MAPE_gas_arma12	 ### 0.8309913
##MS error (MSE)##
MSE_gas_arma13=1/168*(sum(abs(IMF11- gas_IMF11$fitted)^2))
MSE_gas_arma13   ### 2676793




##############################IMF21 prediction########################
gas_IMF21 = fcast(IMF21, method = "GMDH", input = 5, layer = 1, f.number = 5,
                    level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                      0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)##
MAD_gas_gm_21=1/163*(sum(abs(IMF21-gas_IMF21$fitted)))
MAD_gas_gm_21  	    ### 429.6846
##mean absolute percentage error (MAPE)##
MAPE_gas_gm_22=1/163*(sum(abs((IMF21-gas_IMF21$fitted)/IMF21)))
MAPE_gas_gm_22	    ### 1.482818
##MS error (MSE)##  
MSE_gas_gm_32=1/163*(sum(abs(IMF21-gas_IMF21$fitted)^2))
MSE_gas_gm_32       ###  286325.8
#####################################RGMDH IMF21###############
gas_rgm_IMF21 = fcast(IMF21, method = "RGMDH", input = 3, layer = 1, f.number = 5,
                        level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                          0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_rgm_21=1/165*(sum(abs(IMF21-gas_rgm_IMF21$fitted)))
MAD_gas_rgm_21  	###  400.0308
##mean absolute percentage error (MAPE)
MAPE_gas_rgm_22=1/165*(sum(abs((IMF21-gas_rgm_IMF21$fitted)/IMF21)))
MAPE_gas_rgm_22	      ###  1.49854
##MS error (MSE)##  
MSE_gas_rgm_23=1/165*(sum(abs(IMF21-gas_rgm_IMF21$fitted)^2))
MSE_gas_rgm_23       ### 241241.9
#####
############################Neural network###############
library(neuralnet)
library(RSNNS)
axis211=embed(IMF21,2)
axis21=axis211[,-1]
response21=axis211[,1]
set.seed(21)
rbf.model_21 <- rbf(response21, axis21, size = 5, maxit = 1000, linOut = FALSE)
error21=response21-fitted(rbf.model_21)
##Mean Absolute  Deviation (MAD)
MAD_gas_rbf_21=1/168*(sum(abs(error21)))
MAD_gas_rbf_21  	### 1118.269
##mean absolute percentage error (MAPE)
MAPE_gas_rbf_22=1/168*(sum(abs((error21)/response21)))
MAPE_gas_rbf_22	      ### 1118.269
##MS error (MSE)##  
MSE_gas_rbf_23=1/168*(sum(abs(error21)^2))
MSE_gas_rbf_23        ### 1893641

#########################ARMA 21####################

gas_arma_IMF21=auto.arima(IMF21)
summary(gas_arma_IMF21)
##Mean Absolute  Deviation (MAD)
MAD_gas_arma21=1/168*(sum(abs(IMF21- gas_arma_IMF21$fitted)))
MAD_gas_arma21  	 ### 1.469901
##mean absolute percentage error (MAPE)
MAPE_gas_arma22=1/168*(sum(abs((IMF21- gas_arma_IMF21$fitted)/IMF21)))
MAPE_gas_arma22	       ### 244.6553
##MS error (MSE)##244.6553
MSE_gas_arma23=1/168*(sum(abs(IMF21- gas_arma_IMF21$fitted)^2))
MSE_gas_arma23         ### 99122.57
pred22=gas_arma_IMF21$fitted##########################

##################################IM31 PREDICTION#################
gas_IMF31 = fcast(IMF31, method = "GMDH", input = 5, layer = 1, f.number = 5,
                    level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                      0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_gm_31=1/163*(sum(abs(IMF31-gas_IMF31$fitted)))
MAD_gas_gm_31  	    ### 81.17751
##mean absolute percentage error (MAPE)
MAPE_gas_gm_32=1/163*(sum(abs((IMF31-gas_IMF31$fitted)/IMF31)))
MAPE_gas_gm_32	    ### 1.362759
##MS error (MSE)  
MSE_gas_gm_33=1/163*(sum(abs(IMF31-gas_IMF31$fitted)^2))
MSE_gas_gm_33       ### 10876.17

#####################################RGMDH IMF31###############
gas_rgm_IMF31 = fcast(IMF31, method = "RGMDH", input = 3, layer = 1, f.number = 5,
                        level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                          0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_rgm_31=1/165*(sum(abs(IMF31-gas_rgm_IMF31$fitted)))
MAD_gas_rgm_31  	### 49.58508
##mean absolute percentage error (MAPE)
MAPE_gas_rgm_32=1/165*(sum(abs((IMF31-gas_rgm_IMF31$fitted)/IMF31)))
MAPE_gas_rgm_32	      ### 1.048442
##MS error (MSE)##  
MSE_gas_rgm_33=1/165*(sum(abs(IMF31-gas_rgm_IMF31$fitted)^2))
MSE_gas_rgm_33        ###  4120.564

############################Neural network###############
library(neuralnet)
library(RSNNS)
gas_ts=ts(IMF31,frequency=12,start=c(2005,7,1))
axis311=embed(IMF31,2)
axis31=axis311[,-1]
response31=axis311[,1]
set.seed(21)
rbf.model_31 <- rbf(response31, axis31, size = 5, maxit = 1000, linOut = TRUE)
error31=response31-fitted(rbf.model_31)
##Mean Absolute  Deviation (MAD)
MAD_gas_rbf_31=1/168*(sum(abs(error31)))
MAD_gas_rbf_31  	###  1039.67
##mean absolute percentage error (MAPE)
MAPE_gas_rbf_32=1/168*(sum(abs((error31)/response31)))
MAPE_gas_rbf_32	      ### 1.9682
##MS error (MSE)##  
MSE_gas_rbf_32=1/168*(sum(abs(error31)^2))
MSE_gas_rbf_32        ### 1781044

#########################ARMA####################
gas_arma_IMF31=auto.arima(IMF31)
summary(gas_arma_IMF31)

##Mean Absolute  Deviation (MAD)##
MAD_gas_arma31=1/168*(sum(abs(IMF31- gas_arma_IMF31$fitted)))
MAD_gas_arma31  	 ### 24.94387
##mean absolute percentage error (MAPE)##
MAPE_gas_arma32=1/168*(sum(abs((IMF31- gas_arma_IMF31$fitted)/IMF31)))
MAPE_gas_arma32	       ### 1.076262
##MS error (MSE)##
MSE_gas_arma33=1/168*(sum(abs(IMF31- gas_arma_IMF31$fitted)^2))
MSE_gas_arma33         ### 1085.691

pred33=gas_arma_IMF31$fitted  ################

##################################IM41 PREDICTION#################
gas_IMF41 = fcast(IMF41, method = "GMDH", input = 5, layer = 1, f.number = 5,
                    level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                      0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_gm_41=1/163*(sum(abs(IMF41-gas_IMF41$fitted)))
MAD_gas_gm_41  	    ### 17.81153
##mean absolute percentage error (MAPE)
MAPE_gas_gm_42=1/163*(sum(abs((IMF41-gas_IMF41$fitted)/IMF41)))
MAPE_gas_gm_42	    ###  0.09671955
##MS error (MSE)  
MSE_gas_gm_43=1/163*(sum(abs(IMF41-gas_IMF41$fitted)^2))
MSE_gas_gm_43       ###  519.266

#####################################RGMDH IMF41###############
library(GMDH)
gas_rgm_IMF41 = fcast(IMF41, method = "RGMDH", input = 3, layer = 1, f.number = 5,
                        level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                          0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_rgm_41=1/165*(sum(abs(IMF41-gas_rgm_IMF41$fitted)))
MAD_gas_rgm_41  	### 5.77092
##mean absolute percentage error (MAPE)
MAPE_gas_rgm_42=1/165*(sum(abs((IMF41-gas_rgm_IMF41$fitted)/IMF41)))
MAPE_gas_rgm_42	      ### 0.03783595
##MS error (MSE)##  
MSE_gas_rgm_43=1/165*(sum(abs(IMF41-gas_rgm_IMF41$fitted)^2))
MSE_gas_rgm_43        ### 49.91772

############################Neural network###############
library(neuralnet)
library(RSNNS)
gas_ts=ts(IMF31,frequency=12,start=c(2005,7,1))
axis411=embed(IMF41,2)
axis41=axis411[,-1]
response41=axis411[,1]
set.seed(21)
rbf.model_41 <- rbf(response41, axis41, size = 5, maxit = 1000, linOut = TRUE)
error41=response41-fitted(rbf.model_41)
##Mean Absolute  Deviation (MAD)
MAD_gas_rbf_41=1/168*(sum(abs(error41)))
MAD_gas_rbf_41  	### 550.682
##mean absolute percentage error (MAPE)
MAPE_gas_rbf_42=1/168*(sum(abs((error41)/response41)))
MAPE_gas_rbf_42	      ###0.9388052
##MS error (MSE)##  
MSE_gas_rbf_42=1/168*(sum(abs(error41)^2))
MSE_gas_rbf_42        ###  705580.3

#########################ARMA####################
gas_arma_IMF41=auto.arima(IMF41)
summary(gas_arma_IMF41)
##Mean Absolute  Deviation (MAD)##
MAD_gas_arma41=1/168*(sum(abs(IMF41- gas_arma_IMF41$fitted)))
MAD_gas_arma41  	 ### 149.5262
##mean absolute percentage error (MAPE)##
MAPE_gas_arma42=1/168*(sum(abs((IMF41- gas_arma_IMF41$fitted)/IMF41)))
MAPE_gas_arma42	       ### 0.7617945
##MS error (MSE)##
MSE_gas_arma43=1/168*(sum(abs(IMF41- gas_arma_IMF41$fitted)^2))
MSE_gas_arma43         ### 45547.69
pred44=gas_arma_IMF41$fitted################

##################################IM51 PREDICTION#################
gas_IMF51 = fcast(IMF51, method = "GMDH", input = 5, layer = 1, f.number = 5,
                    level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                      0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_gm_51=1/163*(sum(abs(IMF51-gas_IMF51$fitted)))
MAD_gas_gm_51  	    ### 3.605029
##mean absolute percentage error (MAPE)
MAPE_gas_gm_52=1/163*(sum(abs((IMF51-gas_IMF51$fitted)/IMF51)))
MAPE_gas_gm_52	    ### 0.02102941
##MS error (MSE)  
MSE_gas_gm_53=1/163*(sum(abs(IMF51-gas_IMF51$fitted)^2))
MSE_gas_gm_53       ###   22.90622

#####################################RGMDH IMF51###############
gas_rgm_IMF51 = fcast(IMF51, method = "RGMDH", input = 3, layer = 1, f.number = 5,
                        level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                          0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_rgm_51=1/165*(sum(abs(IMF51-gas_rgm_IMF51$fitted)))
MAD_gas_rgm_51  	### 0.7001262
##mean absolute percentage error (MAPE)
MAPE_gas_rgm_52=1/165*(sum(abs((IMF51-gas_rgm_IMF51$fitted)/IMF51)))
MAPE_gas_rgm_52	      ### 0.003315759
##MS error (MSE)##  
MSE_gas_rgm_53=1/165*(sum(abs(IMF51-gas_rgm_IMF51$fitted)^2))
MSE_gas_rgm_53        ### 0.7692934
pred55=gas_rgm_IMF51$fitted###################
############################Neural network###############
library(neuralnet)
library(RSNNS)
axis511=embed(IMF51,2)
axis51=axis511[,-1]
response51=axis511[,1]
set.seed(21)
rbf.model_51 <- rbf(response51, axis51, size = 20, maxit = 1000, linOut = TRUE)
error51=response51-fitted(rbf.model_51)
##Mean Absolute  Deviation (MAD)
MAD_gas_rbf_51=1/168*(sum(abs(error51)))
MAD_gas_rbf_51  	### 549.6609
##mean absolute percentage error (MAPE)
MAPE_gas_rbf_52=1/168*(sum(abs((error51)/response51)))
MAPE_gas_rbf_52	      ###0.8568592
##MS error (MSE)##  
MSE_gas_rbf_52=1/168*(sum(abs(error51)^2))
MSE_gas_rbf_52        ### 542965.8

#########################ARMA####################
gas_arma_IMF51=auto.arima(IMF51)
summary(gas_arma_IMF51)

##Mean Absolute  Deviation (MAD)##
MAD_gas_arma51=1/168*(sum(abs(IMF51- gas_arma_IMF51$fitted)))
MAD_gas_arma51  	 ###659.7642
##mean absolute percentage error (MAPE)##
MAPE_gas_arma52=1/168*(sum(abs((IMF51- gas_arma_IMF51$fitted)/IMF51)))
MAPE_gas_arma52	       ### 1
##MS error (MSE)##
MSE_gas_arma53=1/168*(sum(abs(IMF51- gas_arma_IMF51$fitted)^2))
MSE_gas_arma53         ### 657500.8



##################################IM61 PREDICTION#################
gas_IMF61 = fcast(IMF61, method = "GMDH", input = 5, layer = 1, f.number = 5,
                    level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                      0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_gm_61=1/163*(sum(abs(IMF61-gas_IMF61$fitted)))
MAD_gas_gm_61  	    ###  0.154451
##mean absolute percentage error (MAPE)
MAPE_gas_gm_62=1/163*(sum(abs((IMF61-gas_IMF61$fitted)/IMF61)))
MAPE_gas_gm_62	    ### 0.008648799
##MS error (MSE)  
MSE_gas_gm_63=1/163*(sum(abs(IMF61-gas_IMF61$fitted)^2))
MSE_gas_gm_63       ###0.05905085

#####################################RGMDH IMF61###############
gas_rgm_IMF61 = fcast(IMF61, method = "RGMDH", input = 3, layer = 1, f.number = 5,
                        level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                          0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_rgm_61=1/165*(sum(abs(IMF61-gas_rgm_IMF61$fitted)))
MAD_gas_rgm_61  	###  0.2102163
##mean absolute percentage error (MAPE)
MAPE_gas_rgm_62=1/165*(sum(abs((IMF61-gas_rgm_IMF61$fitted)/IMF61)))
MAPE_gas_rgm_62	      ### 0.006403599
##MS error (MSE)##  
MSE_gas_rgm_63=1/165*(sum(abs(IMF61-gas_rgm_IMF61$fitted)^2))
MSE_gas_rgm_63        ### 0.07239722

############################Neural network###############
library(neuralnet)
library(RSNNS)
gas_ts=ts(IMF61,frequency=12,start=c(2005,7,1))
axis611=embed(IMF61,2)
axis61=axis611[,-1]
response61=axis611[,1]
set.seed(21)
rbf.model_61 <- rbf(response61, axis61, size = 20, maxit = 1000, linOut = TRUE)
error61=response61-fitted(rbf.model_61)
##Mean Absolute  Deviation (MAD)
MAD_gas_rbf_61=1/168*(sum(abs(error61)))
MAD_gas_rbf_61  	### 430.0549
##mean absolute percentage error (MAPE)
MAPE_gas_rbf_62=1/168*(sum(abs((error61)/response61)))
MAPE_gas_rbf_62	      ### 5.183348
##MS error (MSE)##  
MSE_gas_rbf_63=1/168*(sum(abs(error61)^2))
MSE_gas_rbf_63        ###  308170.4

#########################ARMA####################
gas_arma_IMF61=auto.arima(IMF61)
summary(gas_arma_IMF61)

##Mean Absolute  Deviation (MAD)
MAD_gas_arma61=1/168*(sum(abs(IMF61- gas_arma_IMF61$fitted)))
MAD_gas_arma61  	 ### 1.565191
##mean absolute percentage error (MAPE)
MAPE_gas_arma62=1/168*(sum(abs((IMF61- gas_arma_IMF61$fitted)/IMF61)))
MAPE_gas_arma62	       ### 0.03226629
##MS error (MSE)##
MSE_gas_arma63=1/168*(sum(abs(IMF61- gas_arma_IMF61$fitted)^2))
MSE_gas_arma63         ### 3.234601
pred66=gas_arma_IMF61$fitted########
##################################IM71 PREDICTION#################
library(GMDH)
gas_IMF71 = fcast(IMF71, method = "GMDH", input = 5, layer = 1, f.number = 5,
                  level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                    0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_gm_71=1/163*(sum(abs(IMF71-gas_IMF71$fitted)))
MAD_gas_gm_71  	    ###  9.422851
##mean absolute percentage error (MAPE)
MAPE_gas_gm_72=1/163*(sum(abs((IMF71-gas_IMF71$fitted)/IMF71)))

MAPE_gas_gm_72	    ###  7.749626e-05
##MS error (MSE)  
MSE_gas_gm_73=1/163*(sum(abs(IMF71-gas_IMF71$fitted)^2))

MSE_gas_gm_73       ###   135.9448
#####################################RGMDH IMF71###############
library(GMDH)
gas_rgm_IMF71 = fcast(IMF71, method = "RGMDH", input = 3, layer = 1, f.number = 5,
                      level = 95, tf = "all", weight = 0.70, lambda = c(0, 0.01, 0.02, 0.04, 0.08, 0.16,
                                                                        0.32, 0.64, 1.28, 2.56, 5.12, 10.24))
##Mean Absolute  Deviation (MAD)
MAD_gas_rgm_71=1/165*(sum(abs(IMF71-gas_rgm_IMF71$fitted)))
MAD_gas_rgm_71  	### 0.002968392
##mean absolute percentage error (MAPE)
MAPE_gas_rgm_72=1/165*(sum(abs((IMF71-gas_rgm_IMF71$fitted)/IMF71)))
MAPE_gas_rgm_72	      ###2.423086e-08
##MS error (MSE)##  
MSE_gas_rgm_73=1/165*(sum(abs(IMF71-gas_rgm_IMF71$fitted)^2))
MSE_gas_rgm_73        ### 1.251892e-05
pred77=gas_rgm_IMF71$fitted########
############################Neural network###############
library(neuralnet)
library(RSNNS)
gas_ts=ts(IMF71,frequency=12,start=c(2005,7,1))
axis711=embed(IMF71,2)
axis71=axis711[,-1]
response71=axis711[,1]
set.seed(21)
rbf.model_71 <- rbf(response71, axis71, size = 20, maxit = 1000, linOut = TRUE)
error71=response71-fitted(rbf.model_71)
##Mean Absolute  Deviation (MAD)
MAD_gas_rbf_71=1/168*(sum(abs(error71)))
MAD_gas_rbf_71  	###   22850.69
##mean absolute percentage error (MAPE)
MAPE_gas_rbf_72=1/168*(sum(abs((error71)/response71)))
MAPE_gas_rbf_72	      ###  0.1866985
##MS error (MSE)##  
MSE_gas_rbf_72=1/168*(sum(abs(error71)^2))
MSE_gas_rbf_72        ### 531156229

#########################ARMA####################

gas_arma_IMF71=auto.arima(IMF71)
summary(gas_arma_IMF71)
plot(forecast(gas_arma_IMF71,h=30))
points(1:length(IMF71),fitted(gas_arma_IMF71),type="l",col="green")

##Mean Absolute  Deviation (MAD)
MAD_gas_arma71=1/168*(sum(abs(IMF71- gas_arma_IMF71$fitted)))
MAD_gas_arma71  	 ### 1.420566
##mean absolute percentage error (MAPE)
MAPE_gas_arma72=1/168*(sum(abs((IMF71- gas_arma_IMF71$fitted)/IMF71)))
MAPE_gas_arma72	       ### 1.22089e-05
##MS error (MSE)##
MSE_gas_arma73=1/168*(sum(abs(IMF71- gas_arma_IMF71$fitted)^2))
MSE_gas_arma73         ### 157.7217


############################################################################
pred_gas_ac=cbind(pred11,pred22,pred33,pred44,pred55,pred66,pred77)
pred_gas_sum_ac=rowSums(pred_gas_ac)
plot(Time,gas_rec_imfs,type="o",col=1)
lines(Time,pred_gas_sum_ac,col=4)
legend("topright", c("Denoised gas through EMD-EEMD-MM", "Predicted gas"), cex=0.7, fill=c(1,4))
error_final_ac=gas_rec_imfs-pred_gas_sum_ac
##Mean Absolute  Deviation (MAD)
MAD_gas12=1/165*(sum(abs(error_final_ac[-c(1:3)])))
MAD_gas12 	       ### 2020.751
##mean absolute percentage error (MAPE)
MAPE_gas22=1/165*(sum(abs((error_final_ac[-c(1:3)])/gas_rec_imfs[-c(1:3)])))
MAPE_gas22	       ### 0.01644532
##MS error (MSE)
MSE_gas32=1/165*(sum(abs(error_final_ac[-c(1:3)])^2))
MSE_gas32              ###  7766789