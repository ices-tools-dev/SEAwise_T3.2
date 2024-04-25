
DF <- read.table("Rec-SSB_GSA20_dep1__env_yearly.csv", sep=",",header = T)


#add a column with sum chla
DF$ChlTot<-DF$chl1+DF$chl2+DF$chl3+DF$chl4

names(DF)
#[1] "year"        "HKE_rec"     "HKE_ssb"     "HKE_ssb_1yl" "MUT_rec"     "MUT_ssb"    
#[7] "MUT_ssb_1yl" "Temp"        "S"           "ssh"         "grossPP"     "netPP"      
#[13] "P1c"         "P2c"         "P3c"         "P4c"         "Z4c"         "Z5c"        
#[19] "Z6c"         "B1c"         "chl1"        "chl2"        "chl3"        "chl4"       
#[25] "pH"          "O"           "ComProd"     "MLD"         "ben_DIC"     "ben_N"      
#[31] "ben_P"       "ben_Si"      "ben_O"       "ben_dep"     "ben_sus"     "ben_DOC"    
#[37] "optChl"      "optDepth"    "PC1"         "PC2"         "HS_HKE_large" "HS_HKE_small"
#[43] "HS_HKE_tot"   "HS_MUT_large" "HS_MUT_small" "HS_MUT_tot" "ChlTot"     



# Apply GAMs ----

library(mgcv)
library(fmsb)
library(Hmisc)
library(dplyr)
library(dsm)
library(tidyverse)
library(rio) 


##  data filtering----
DF1<- DF %>% 
  filter(year>2002)


## correlations ---- 

#correlation matrix
#exclude columns
DF2 = subset(DF1, select = -c(year,HKE_rec, HKE_ssb,HKE_ssb_1yl,MUT_rec,MUT_ssb,MUT_ssb_1yl) )
#estimate correlation matrix
M<-cor(DF2)
head(round(M,2))

library(corrplot)
corrplot(M, method="number", number.cex = 0.5)

jpeg (filename="cor2_d3.3.jpeg", pointsize = 16, width=1300, height=1100)
corrplot(M, method="number", number.cex = 0.8)
graphics.off()


#correlation table
rcorr(as.matrix(DF1)) 
corr.mat<-rcorr(as.matrix(DF1)) 
df.corr.mat.r=data.frame(corr.mat$r)
write.csv(df.corr.mat.r, file ="corr.csv")



     #------------------------#
     #                        #
     ##          HKE           ----     
     #                        #
     #------------------------#



#_______
### CONSIDERING ALL AVAILABLE VARIABLES ----
#_______

## Stepwise forward inclusion for GAM modeling


#step 1
mod01 <- gam(HKE_rec ~ s(HKE_ssb, k=5), family= gaussian (link="log"),data=DF1, select=T) 
mod02 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(Temp, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod03 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(S, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod04 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(ssh, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod05 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(grossPP, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod06 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(netPP, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod07 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(P1c, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod08 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(P2c, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod09 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(P3c, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod10 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(P4c, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod11 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(Z4c, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod12 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(Z5c, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod13 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(Z6c, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod14 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(B1c, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod15 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(chl1, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod16 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(chl2, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod17 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(chl3, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod18 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(chl4, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod19 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(pH, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod20 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(O, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod21 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(ComProd, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod22 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(MLD, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod23 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(ben_DIC, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod24 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(ben_N, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod25 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(ben_P, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod26 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(ben_Si, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod27 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(ben_O, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod28 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(ben_dep, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod29 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(ben_sus, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod30 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(ben_DOC, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod31 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(optChl, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod32 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(optDepth, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod33 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(ChlTot, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod34 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(PC1, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod35 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(PC2, k=5), family= gaussian (link="log"),data=DF1, select=T)

#w lags
mod34b <- gam(HKE_rec ~ s(HKE_ssb_1yl, k=5)+s(PC1, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod35b <- gam(HKE_rec ~ s(HKE_ssb_1yl, k=5)+s(PC2, k=5), family= gaussian (link="log"),data=DF1, select=T)

#additional for D3.3
mod36 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(HS_HKE_large, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod37 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(HS_HKE_small, k=5), family= gaussian (link="log"),data=DF1, select=T)

s01 <- summary(mod01)
s02 <- summary(mod02)
s03 <- summary(mod03)
s04 <- summary(mod04)
s05 <- summary(mod05)
s06 <- summary(mod06)
s07 <- summary(mod07)
s08 <- summary(mod08)
s09 <- summary(mod09)
s10 <- summary(mod10)
s11 <- summary(mod11)
s12 <- summary(mod12)
s13 <- summary(mod13)
s14 <- summary(mod14)
s15 <- summary(mod15)
s16 <- summary(mod16)
s17 <- summary(mod17)
s18 <- summary(mod18)
s19 <- summary(mod19)
s20 <- summary(mod20)
s21 <- summary(mod21)
s22 <- summary(mod22)
s23 <- summary(mod23)
s24 <- summary(mod24)
s25 <- summary(mod25)
s26 <- summary(mod26)
s27 <- summary(mod27)
s28 <- summary(mod28)
s29 <- summary(mod29)
s30 <- summary(mod30)
s31 <- summary(mod31)
s32 <- summary(mod32)
s33 <- summary(mod33)
s34 <- summary(mod34)
s35 <- summary(mod35)
s36 <- summary(mod36)
s37 <- summary(mod37)


res37 <- data.frame(mod = c("mod01","mod02","mod03","mod04","mod05","mod06","mod07","mod08","mod09","mod10","mod11","mod12", "mod13",
                            "mod14", "mod15", "mod16", "mod17", "mod18", "mod19", "mod20", "mod21", "mod22", "mod23", "mod24", "mod25", 
                            "mod26", "mod27", "mod28", "mod29", "mod30", "mod31", "mod32", "mod33", "mod34", "mod35", "mod36", "mod37"),
                    dev.expl = c(s01$dev.expl,  s02$dev.expl,  s03$dev.expl, s04$dev.expl, s05$dev.expl,s06$dev.expl, s07$dev.expl, 
                                 s08$dev.expl, s09$dev.expl, s10$dev.expl, s11$dev.expl, s12$dev.expl,  s13$dev.expl,  
                                 s14$dev.expl,  s15$dev.expl,  s16$dev.expl,  s17$dev.expl,  s18$dev.expl,  s19$dev.expl,  s20$dev.expl,
                                 s21$dev.expl,  s22$dev.expl,  s23$dev.expl,  s24$dev.expl,  s25$dev.expl,  s26$dev.expl,  s27$dev.expl,  
                                 s28$dev.expl,  s29$dev.expl,  s30$dev.expl,  s31$dev.expl,  s32$dev.expl, s33$dev.expl, s34$dev.expl, s35$dev.expl, s36$dev.expl,
                                 s37$dev.expll), 
                    AIC = AIC(mod01,mod02,mod03,mod04,mod05,mod06, mod07, mod08, mod09, mod10, mod11, mod12, mod13, mod14, mod15, 
                              mod16, mod17, mod18, mod19, mod20, mod21, mod22, mod23, mod24, mod25, mod26, mod27, mod28, mod29,
                              mod30, mod31, mod32, mod33, mod34, mod35, mod36, mod37)[,2],
                    GCV=c(s01$ sp.criterion, s02$sp.criterion, s03$sp.criterion, s04$sp.criterion, s05$sp.criterion, s06$sp.criterion,
                          s07$sp.criterion, s08$sp.criterion, s09$sp.criterion, s10$sp.criterion, s11$sp.criterion, s12$sp.criterion, 
                          s13$sp.criterion, s14$sp.criterion, s15$sp.criterion, s16$sp.criterion, s17$sp.criterion, s18$sp.criterion, 
                          s19$sp.criterion, s20$sp.criterion, s21$sp.criterion, s22$sp.criterion, s23$sp.criterion, s24$sp.criterion, 
                          s25$sp.criterion, s26$sp.criterion, s27$sp.criterion, s28$sp.criterion, s29$sp.criterion, s30$sp.criterion, 
                          s31$sp.criterion, s32$sp.criterion, s33$sp.criterion, s34$sp.criterion, s35$sp.criterion, s36$sp.criterion, s37$sp.criterion))
res37


summary(mod25)
AIC(mod25)
gam.check(mod25)
plot(mod01, page=1, all.terms=T)




#________________________________
### ONLY WITH SELECTED VARIABLES ----
#________________________________

#step 1
mod01 <- gam(HKE_rec ~ s(HKE_ssb, k=5), family= gaussian (link="log"),data=DF1, select=T) 
mod02 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(Temp, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod03 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(S, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod04 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(ssh, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod05 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(grossPP, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod06 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(netPP, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod07 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(Z4c, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod08 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(chl4, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod09 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(O, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod10 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(ComProd, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod11 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(MLD, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod12 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(optChl, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod13 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(ChlTot, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod14 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(HS_HKE_small, k=5), family= gaussian (link="log"),data=DF1, select=T)


s01 <- summary(mod01)
s02 <- summary(mod02)
s03 <- summary(mod03)
s04 <- summary(mod04)
s05 <- summary(mod05)
s06 <- summary(mod06)
s07 <- summary(mod07)
s08 <- summary(mod08)
s09 <- summary(mod09)
s10 <- summary(mod10)
s11 <- summary(mod11)
s12 <- summary(mod12)
s13 <- summary(mod13)
s14 <- summary(mod14)

res14 <- data.frame(mod = c("mod01","mod02","mod03","mod04","mod05","mod06","mod07","mod08","mod09","mod10","mod11","mod12", "mod13", "mod14"),
                    dev.expl = c(s01$dev.expl,  s02$dev.expl,  s03$dev.expl, s04$dev.expl, s05$dev.expl,s06$dev.expl, s07$dev.expl, 
                                 s08$dev.expl, s09$dev.expl, s10$dev.expl, s11$dev.expl, s12$dev.expl,  s13$dev.expl,  s14$dev.expl), 
                    AIC = AIC(mod01,mod02,mod03,mod04,mod05,mod06, mod07, mod08, mod09, mod10, mod11, mod12, mod13, mod14)[,2],
                    GCV=c(s01$ sp.criterion, s02$sp.criterion, s03$sp.criterion, s04$sp.criterion, s05$sp.criterion, s06$sp.criterion,
                          s07$sp.criterion, s08$sp.criterion, s09$sp.criterion, s10$sp.criterion, s11$sp.criterion, s12$sp.criterion, 
                          s13$sp.criterion, s14$sp.criterion))
res14

#mod  dev.expl      AIC        GCV
#1  mod01 0.2053747 436.5086 1815393773
#2  mod02 0.4144748 431.9283 1419369010
#3  mod03 0.3858528 433.6200 1573611011
#4  mod04 0.4171955 432.5716 1482729358
#5  mod05 0.4268925 433.7368 1613361124
#6  mod06 0.2053782 436.5086 1815393988
#7  mod07 0.2053749 436.5086 1815393138
#8  mod08 0.5628492 429.4466 1283116425#
#9  mod09 0.2053748 436.5086 1815393191
#10 mod10 0.2053756 436.5086 1815393641
#11 mod11 0.2053748 436.5086 1815393141
#12 mod12 0.3706238 433.8600 1590979100
#13 mod13 0.2953191 435.3879 1722361513
#14 mod14 0.4568324 433.0973 1565003414

summary(mod08)
AIC(mod25)
gam.check(mod08)
plot(mod08, page=1, all.terms=T)


#step2
#exclude correlated variable
##grossPP, optChl, ChlTot

mod01 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(chl4), family= gaussian (link="log"),data=DF1, select=T) 
mod02 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(chl4)+s(Temp, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod03 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(chl4)+s(S, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod04 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(chl4)+s(ssh, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod05 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(chl4)+s(netPP, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod06 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(chl4)+s(Z4c, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod07 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(chl4)+s(O, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod08 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(chl4)+s(ComProd, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod09 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(chl4)+s(MLD, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod10 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(chl4)+s(HS_HKE_small, k=5), family= gaussian (link="log"),data=DF1, select=T)


s01 <- summary(mod01)
s02 <- summary(mod02)
s03 <- summary(mod03)
s04 <- summary(mod04)
s05 <- summary(mod05)
s06 <- summary(mod06)
s07 <- summary(mod07)
s08 <- summary(mod08)
s09 <- summary(mod09)
s10 <- summary(mod10)

res10 <- data.frame(mod = c("mod01","mod02","mod03","mod04","mod05","mod06","mod07","mod08","mod09","mod10"),
                    dev.expl = c(s01$dev.expl,  s02$dev.expl,  s03$dev.expl, s04$dev.expl, s05$dev.expl,s06$dev.expl, s07$dev.expl, 
                                 s08$dev.expl, s09$dev.expl, s10$dev.expl), 
                    AIC = AIC(mod01,mod02,mod03,mod04,mod05,mod06, mod07, mod08, mod09, mod10)[,2],
                    GCV=c(s01$ sp.criterion, s02$sp.criterion, s03$sp.criterion, s04$sp.criterion, s05$sp.criterion, s06$sp.criterion,
                          s07$sp.criterion, s08$sp.criterion, s09$sp.criterion, s10$sp.criterion))
res10

#    mod  dev.expl      AIC        GCV
#1 mod01 0.5515163 429.7151 1298291201
#2 mod02 0.6816634 423.4114  912740514@
#3 mod03 0.5590353 429.6179 1295710111
#4 mod04 0.5541485 429.6448 1293979212
#5 mod05 0.5515164 429.7151 1298290945
#6 mod06 0.6388484 427.9813 1228944237
#7 mod07 0.5729755 429.4577 1293541651
#8 mod08 0.5515164 429.7151 1298291090
#9 mod09 0.8886240 412.7274  639768882# but MLD excluded because values don't seem to correspond to Upper mixed layer depth
#10 mod10 0.7975124 420.7321  896181437# lower AIC but ??AIC very small: mod01 was selected as the most parsimonious



#_______________________
### Best model ----
#_______________________

mod01 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(chl4), family= gaussian (link="log"),data=DF1, select=T) 


summary(mod01)
gam.check(mod01)
plot(mod01, page=1, all.terms=T)

png (filename="plot_HKE_D3.3.png", width=800, height=800, pointsize = 20)
par(mfrow=c(2,2))
plot(mod01, all.terms=T, scale=0, shade=T, scheme=2,too.far=0.1,residuals=FALSE)#, page=1
graphics.off()

par(mfrow=c(1,1))

#save gam.check
png (filename="check_HKE_D3.3.png", width=800, height=800, pointsize = 20)
par(mfrow=c(2,2))
gam.check(mod01)
graphics.off()


### evaluate prediction skill of a model ---------------------------------
#Script info 

# SEAwise Task 3.3.
# example code to evaluate prediction skill of a model
# following the approach applied in NS (Brooks, Taylor et al)
# sequential hindcasting with rolling origin (up to 15 years)

# last changes: 23/01/2023

library(ggplot2)
library(tidyr)
library(dplyr)

theme_set(theme_bw(base_size=16))


# short name of the case study and species. it will be used to store objects
cs <- "Eastern Ionian Sea" 
sp <- "HKE"


#### Prediction skill --------------------------------------------------------

# selected model 
dat=DF1
m0 <- gam(HKE_rec ~ s(HKE_ssb, k=5)+s(chl4), family= gaussian (link="log"),data=DF1, select=T)
#m0 <- lm(formula=lw0 ~ lssb, data=dat) # CHANGE THE MODEL FOR EACH CASE STUDY
summary(m0)

# check there are no missing values
nrow(dat) == nrow(na.omit(dat[, all.vars(formula(m0))]))

# number of data points used in the model fit
n <- nrow(dat) 

# number of years to be removed (15 at most), ensure we have at last half of the data for the model fit
 #reduced 15 to 4 because "Model has more coefficients than data"
x <- min(floor(n/2), 4) 
 
# loop

out <- NULL 
for (i in 1:x){
  m <- gam(formula(m0), data=dat[1:(n-i), ]) # adjust the function accordingly. it may need for arguments, e.g. family in case of glm/gam
  pred.dat <- expand.grid(ntrain=(n-i), h=1:i)
  pred.dat$obs <- dat[(n-i+1):n, all.vars(formula(m0))[1], drop=T]
  pred.dat$fit.mod <- predict(m, newdata=dat[(n-i+1):n, ])
  pred.dat$fit.naive <- mean(dat[(n-i-2):(n-i), all.vars(formula(m0))[1], drop=T], na.rm=T)
  pred.dat$e.mod <- pred.dat$obs-pred.dat$fit.mod
  pred.dat$e.naive <- pred.dat$obs-pred.dat$fit.naive
  out <- rbind(out, pred.dat)
}
out <- out %>% 
  mutate(hclass= cut(h, breaks=c(0, 3, 6, 15)),
         hterm=cut(h, breaks=c(0, 3, 6, 15), labels=c("short (1-3)", "medium (4-6)","long (7-15)")))

# summarise by ntrain (nb years for training)


avg <- out %>% 
    summarise(relmae=mean(abs(e.mod))/mean(abs(e.naive)) )

res.ntrain <- out %>% 
  group_by(ntrain) %>% 
  summarise(relmae=mean(abs(e.mod))/mean(abs(e.naive)) )

p <- ggplot(res.ntrain, aes(ntrain, relmae))+
  geom_point(cex=2)+
  geom_line()+
  expand_limits(y=0)+
  geom_hline(yintercept=1, col=2)+
  xlab("Nb years training dataset")+
  ylab(expression(rMAE=MAE[mod]/MAE[naive]))
print(p)

png (filename="MAE_HKE_ntrain.png", width=600, height=600, pointsize = 15)
plot(p)#, page=1
graphics.off()


# summarise by h (prediction time horizon) 

res.h <- out %>% 
  group_by(h) %>% 
  summarise(relmae=mean(abs(e.mod))/mean(abs(e.naive)) )

p <- ggplot(res.h, aes(h, relmae))+
  geom_point(cex=2)+
  geom_line()+
  expand_limits(y=0)+
  geom_hline(yintercept=1, col=2)+
  scale_x_continuous(breaks=seq(1, 15))+
  xlab("Time horizon (years)")+
  ylab(expression(rMAE=MAE[mod]/MAE[naive]))
print(p)

png (filename="MAE_HKE_h.png", width=600, height=600, pointsize = 15)
plot(p)#, page=1
graphics.off()

# summarise by hclass or hterm (grouped time horizons: short, medium, long)
# PLEASE report the values from res.hterm into the excel file

res.hterm <- out %>% 
  group_by(hterm) %>% 
  summarise(relmae=mean(abs(e.mod))/mean(abs(e.naive)) )

p <- ggplot(res.hterm, aes(hterm, relmae))+
  geom_bar(stat="identity", fill="grey")+
  expand_limits(y=0)+
  geom_hline(yintercept=1, col=2)+
  xlab("Time horizon classes")+
  ylab(expression(rMAE=MAE[mod]/MAE[naive]))
print(p)

png (filename="MAE_HKE_hclass.png", width=600, height=600, pointsize = 15)
plot(p)#, page=1
graphics.off()


# save out results
# PLEASE store these objects into the sharepoint

saveRDS(out, file=paste0("output_predskill_",cs,"_",sp,".rds"))
saveRDS(m0, file=paste0("model_predskill_",cs,"_",sp,".rds"))




#### Other packages

# note that the function tsCV() in the forecast package
# and the function ro() in the greybox package have some options to
# run rolling origin approaches to evaluate prediction skills
# These options have not been attempted here





#------------------------#
#                        #
##          MUT           ----
#                        #
#------------------------#



#_______________________
### CONSIDERING ALL AVAILABLE VARIABLES ----
#_______________________

## Stepwise forward inclusion for GAM modeling


#step 1
mod01 <- gam(MUT_rec ~ s(MUT_ssb, k=5), family= gaussian (link="log"),data=DF1, select=T) 
mod02 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(Temp, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod03 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(S, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod04 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(ssh, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod05 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(grossPP, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod06 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(netPP, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod07 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(P1c, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod08 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(P2c, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod09 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(P3c, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod10 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(P4c, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod11 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(Z4c, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod12 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(Z5c, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod13 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(Z6c, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod14 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(B1c, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod15 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(chl1, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod16 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(chl2, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod17 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(chl3, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod18 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(chl4, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod19 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(pH, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod20 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(O, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod21 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(ComProd, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod22 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(MLD, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod23 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(ben_DIC, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod24 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(ben_N, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod25 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(ben_P, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod26 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(ben_Si, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod27 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(ben_O, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod28 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(ben_dep, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod29 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(ben_sus, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod30 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(ben_DOC, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod31 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(optChl, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod32 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(optDepth, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod33 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(ChlTot, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod34 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(PC1, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod35 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(PC2, k=5), family= gaussian (link="log"),data=DF1, select=T)

mod34b <- gam(MUT_rec ~ s(MUT_ssb_1yl, k=5)+s(PC1, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod35b <- gam(MUT_rec ~ s(MUT_ssb_1yl, k=5)+s(PC2, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod01b <- gam(MUT_rec ~ s(MUT_ssb_1yl, k=5), family= gaussian (link="log"),data=DF1, select=T) 

s01 <- summary(mod01)
s02 <- summary(mod02)
s03 <- summary(mod03)
s04 <- summary(mod04)
s05 <- summary(mod05)
s06 <- summary(mod06)
s07 <- summary(mod07)
s08 <- summary(mod08)
s09 <- summary(mod09)
s10 <- summary(mod10)
s11 <- summary(mod11)
s12 <- summary(mod12)
s13 <- summary(mod13)
s14 <- summary(mod14)
s15 <- summary(mod15)
s16 <- summary(mod16)
s17 <- summary(mod17)
s18 <- summary(mod18)
s19 <- summary(mod19)
s20 <- summary(mod20)
s21 <- summary(mod21)
s22 <- summary(mod22)
s23 <- summary(mod23)
s24 <- summary(mod24)
s25 <- summary(mod25)
s26 <- summary(mod26)
s27 <- summary(mod27)
s28 <- summary(mod28)
s29 <- summary(mod29)
s30 <- summary(mod30)
s31 <- summary(mod31)
s32 <- summary(mod32)
s33 <- summary(mod33)
s34 <- summary(mod34)
s35 <- summary(mod35)

res35 <- data.frame(mod = c("mod01","mod02","mod03","mod04","mod05","mod06","mod07","mod08","mod09","mod10","mod11","mod12", "mod13",
                            "mod14", "mod15", "mod16", "mod17", "mod18", "mod19", "mod20", "mod21", "mod22", "mod23", "mod24", "mod25", 
                            "mod26", "mod27", "mod28", "mod29", "mod30", "mod31", "mod32", "mod33", "mod34", "mod35"),
                    dev.expl = c(s01$dev.expl,  s02$dev.expl,  s03$dev.expl, s04$dev.expl, s05$dev.expl,s06$dev.expl, s07$dev.expl, 
                                 s08$dev.expl, s09$dev.expl, s10$dev.expl, s11$dev.expl, s12$dev.expl,  s13$dev.expl,  
                                 s14$dev.expl,  s15$dev.expl,  s16$dev.expl,  s17$dev.expl,  s18$dev.expl,  s19$dev.expl,  s20$dev.expl,
                                 s21$dev.expl,  s22$dev.expl,  s23$dev.expl,  s24$dev.expl,  s25$dev.expl,  s26$dev.expl,  s27$dev.expl,  
                                 s28$dev.expl,  s29$dev.expl,  s30$dev.expl,  s31$dev.expl,  s32$dev.expl, s33$dev.expl, s34$dev.expl, s35$dev.expl), 
                    AIC = AIC(mod01,mod02,mod03,mod04,mod05,mod06, mod07, mod08, mod09, mod10, mod11, mod12, mod13, mod14, mod15, 
                              mod16, mod17, mod18, mod19, mod20, mod21, mod22, mod23, mod24, mod25, mod26, mod27, mod28, mod29,
                              mod30, mod31, mod32, mod33, mod34, mod35)[,2],
                    GCV=c(s01$ sp.criterion, s02$sp.criterion, s03$sp.criterion, s04$sp.criterion, s05$sp.criterion, s06$sp.criterion,
                          s07$sp.criterion, s08$sp.criterion, s09$sp.criterion, s10$sp.criterion, s11$sp.criterion, s12$sp.criterion, 
                          s13$sp.criterion, s14$sp.criterion, s15$sp.criterion, s16$sp.criterion, s17$sp.criterion, s18$sp.criterion, 
                          s19$sp.criterion, s20$sp.criterion, s21$sp.criterion, s22$sp.criterion, s23$sp.criterion, s24$sp.criterion, 
                          s25$sp.criterion, s26$sp.criterion, s27$sp.criterion, s28$sp.criterion, s29$sp.criterion, s30$sp.criterion, 
                          s31$sp.criterion, s32$sp.criterion, s33$sp.criterion, s34$sp.criterion, s35$sp.criterion))
res35

summary(mod25)
AIC(mod25)
gam.check(mod25)
plot(mod25, page=1, all.terms=T)


#mod  dev.expl      AIC       GCV
#1  mod01 0.7620110 307.5853 1440076.9
#2  mod02 0.7907919 306.8647 1416988.2
#3  mod03 0.7699375 307.4484 1438347.1
#4  mod04 0.8908117 300.4492 1132042.7
#5  mod05 0.8844279 298.8879  962255.9
#6  mod06 0.8787158 300.1402 1041717.9
#7  mod07 0.8662151 301.5249 1114155.3
#8  mod08 0.7620112 307.5853 1440076.7
#9  mod09 0.7620110 307.5853 1440077.2
#10 mod10 0.7978239 305.7566 1321731.3
#11 mod11 0.7620112 307.5853 1440077.3
#12 mod12 0.8866816 298.3581  930327.5
#13 mod13 0.7828325 307.1330 1428769.2
#14 mod14 0.7620112 307.5853 1440076.8
#15 mod15 0.8938984 297.3163  881088.0
#16 mod16 0.7620115 307.5853 1440076.8
#17 mod17 0.7620109 307.5853 1440077.5
#18 mod18 0.8227870 304.2293 1231532.4
#19 mod19 0.8655800 302.8313 1238000.2
#20 mod20 0.7620109 307.5853 1440077.2
#21 mod21 0.8719457 300.5309 1048983.3
#22 mod22 0.8231596 304.5301 1260091.3
#23 mod23 0.9795211 270.5579  216867.5
#24 mod24 0.9029203 295.0772  766319.4
#25 mod25 0.9719327 275.1028  268993.3
#26 mod26 0.8970281 295.7392  788047.7
#27 mod27 0.9440166 287.9277  555381.4
#28 mod28 0.7620114 307.5853 1440076.8
#29 mod29 0.7620110 307.5853 1440076.8
#30 mod30 0.9471585 287.0865  533459.7
#31 mod31 0.9242564 292.2519  682833.9#
#32 mod32 0.8626832 300.0053  980536.5
#33 mod33 0.8727751 298.7184  914392.9
#34 mod34 0.8957779 297.0349  868284.0
#35 mod35 0.7620110 307.5853 1440077.0

summary(mod01b)
summary(mod34b)
summary(mod35b)

#   _______________________
#---> no lag |
#   _______________________




### ONLY WITH SELECTED VARIABLES ----


#step 1
mod01 <- gam(MUT_rec ~ s(MUT_ssb, k=5), family= gaussian (link="log"),data=DF1, select=T) 
mod02 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(Temp, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod03 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(S, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod04 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(ssh, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod05 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(grossPP, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod06 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(netPP, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod07 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(Z4c, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod08 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(chl4, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod09 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(O, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod10 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(ComProd, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod11 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(optChl, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod12 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(ChlTot, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod13 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(PC1, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod14 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(PC2, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod15 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(HS_MUT_small, k=5), family= gaussian (link="log"),data=DF1, select=T)


s01 <- summary(mod01)
s02 <- summary(mod02)
s03 <- summary(mod03)
s04 <- summary(mod04)
s05 <- summary(mod05)
s06 <- summary(mod06)
s07 <- summary(mod07)
s08 <- summary(mod08)
s09 <- summary(mod09)
s10 <- summary(mod10)
s11 <- summary(mod11)
s12 <- summary(mod12)
s13 <- summary(mod13)
s14 <- summary(mod14)
s15 <- summary(mod15)


res15 <- data.frame(mod = c("mod01","mod02","mod03","mod04","mod05","mod06","mod07","mod08","mod09","mod10","mod11","mod12", "mod13", "mod14", "mod15"),
                    dev.expl = c(s01$dev.expl,  s02$dev.expl,  s03$dev.expl, s04$dev.expl, s05$dev.expl,s06$dev.expl, s07$dev.expl, 
                                 s08$dev.expl, s09$dev.expl, s10$dev.expl, s11$dev.expl, s12$dev.expl,  s13$dev.expl, s14$dev.expl, s15$dev.expl), 
                    AIC = AIC(mod01,mod02,mod03,mod04,mod05,mod06, mod07, mod08, mod09, mod10, mod11, mod12, mod13, mod14, mod15)[,2],
                    GCV=c(s01$ sp.criterion, s02$sp.criterion, s03$sp.criterion, s04$sp.criterion, s05$sp.criterion, s06$sp.criterion,
                          s07$sp.criterion, s08$sp.criterion, s09$sp.criterion, s10$sp.criterion, s11$sp.criterion, s12$sp.criterion, 
                          s13$sp.criterion, s14$sp.criterion, s15$sp.criterion))
res15


#     mod  dev.expl      AIC        GCV
#1  mod01 0.7620110 307.5853 1440076.9
#2  mod02 0.7907919 306.8647 1416988.2
#3  mod03 0.7699375 307.4484 1438347.1
#4  mod04 0.8908117 300.4492 1132042.7
#5  mod05 0.8844279 298.8879  962255.9
#6  mod06 0.8787158 300.1402 1041717.9
#7  mod07 0.7620112 307.5853 1440077.3
#8  mod08 0.8227870 304.2293 1231532.4
#9  mod09 0.7620109 307.5853 1440077.2
#10 mod10 0.8719457 300.5309 1048983.3
#11 mod11 0.9242564 292.2519  682833.9##
#12 mod12 0.8727751 298.7184  914392.9
#13 mod13 0.8957779 297.0349  868284.0
#14 mod14 0.7620110 307.5853 1440077.0
#15 mod15 0.8218366 305.3603 1337876.3


# explore correlations 
DF3 = subset(DF1, select = c(MUT_ssb,Temp,S,ssh,grossPP,netPP,Z4c,chl4,O,ComProd,optChl,ChlTot,PC1,PC2,HS_MUT_large,HS_MUT_small) )
#estimate correlation matrix
M<-cor(DF3)
head(round(M,2))

library(corrplot)
corrplot(M, method="number", number.cex = 0.5)

jpeg (filename="cor2_d3.3_selected_MUT.jpeg", pointsize = 16, width=800, height=800)
corrplot(M, method="number", number.cex = 0.6)
graphics.off()

#step2
#exclude correlated variables
#grossPP, chl4, comprod,PC1, PC2, ChlTot

mod01 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(optChl, k=5), family= gaussian (link="log"),data=DF1, select=T) 
mod02 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(optChl, k=5)+s(Temp, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod03 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(optChl, k=5)+s(S, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod04 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(optChl, k=5)+s(Z4c, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod05 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(optChl, k=5)+s(O, k=5), family= gaussian (link="log"),data=DF1, select=T)
mod06 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(optChl, k=5)+s(HS_MUT_small, k=5), family= gaussian (link="log"),data=DF1, select=T)




s01 <- summary(mod01)
s02 <- summary(mod02)
s03 <- summary(mod03)
s04 <- summary(mod04)
s05 <- summary(mod05)
s06 <- summary(mod06)



res6 <- data.frame(mod = c("mod01","mod02","mod03","mod04","mod05", "mod06"),
                    dev.expl = c(s01$dev.expl,  s02$dev.expl,  s03$dev.expl, s04$dev.expl, s05$dev.expl, s06$dev.expl), 
                    AIC = AIC(mod01,mod02,mod03,mod04,mod05, mod06)[,2],
                    GCV=c(s01$ sp.criterion, s02$sp.criterion, s03$sp.criterion, s04$sp.criterion, s05$sp.criterion, s06$sp.criterion))
res6



#    mod  dev.expl      AIC      GCV
#1 mod01 0.9242564 292.2519 682833.9
#2 mod02 0.9242565 292.2519 682834.1
#3 mod03 0.9567412 284.8115 492905.8#
#4 mod04 0.9242564 292.2519 682833.9
#5 mod05 0.9242565 292.2519 682834.1
#6 mod06 0.9473958 288.7151 621567.0

par(mfrow=c(2,2))
summary(mod03)
gam.check(mod03)
plot(mod17, page=1, all.terms=T)
AIC(mod03)
summary(mod10)
gam.check(mod10)
plot(mod06, page=1, all.terms=T)


#_______________________
### Best model ----
#_______________________

mod01 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(optChl, k=5)+s(S, k=5), family= gaussian (link="log"),data=DF1, select=T)


summary(mod01)
gam.check(mod01)
plot(mod01, page=1, all.terms=T)

png (filename="plot_MUT_D3.3.png", width=800, height=800, pointsize = 20)
par(mfrow=c(2,2))
plot(mod01, all.terms=T, scale=0, shade=T, scheme=2,too.far=0.1,residuals=FALSE)#, page=1
graphics.off()

par(mfrow=c(1,1))

#save gam.check
png (filename="check_MUT_D3.3.png", width=800, height=800, pointsize = 20)
par(mfrow=c(2,2))
gam.check(mod01)
graphics.off()


### evaluate prediction skill of a model ---------------------------------
#Script info 

# SEAwise Task 3.3.
# example code to evaluate prediction skill of a model
# following the approach applied in NS (Brooks, Taylor et al)
# sequential hindcasting with rolling origin (up to 15 years)

# last changes: 23/01/2023

library(ggplot2)
library(tidyr)
library(dplyr)

theme_set(theme_bw(base_size=16))


# short name of the case study and species. it will be used to store objects
cs <- "Eastern Ionian Sea" 
sp <- "MUT"



#### Prediction skill --------------------------------------------------------

# selected model 
dat=DF1
m0 <- gam(MUT_rec ~ s(MUT_ssb, k=5)+s(optChl, k=5)+s(S, k=5), family= gaussian (link="log"),data=DF1, select=T)
#m0 <- lm(formula=lw0 ~ lssb, data=dat) # CHANGE THE MODEL FOR EACH CASE STUDY
summary(m0)

# check there are no missing values
nrow(dat) == nrow(na.omit(dat[, all.vars(formula(m0))]))

# number of data points used in the model fit
n <- nrow(dat) 

# number of years to be removed (15 at most), ensure we have at last half of the data for the model fit
#reduced 15 to 4 because "Model has more coefficients than data"
x <- min(floor(n/2), 5) 

# loop

out <- NULL 
for (i in 1:x){
  m <- gam(formula(m0), data=dat[1:(n-i), ]) # adjust the function accordingly. it may need for arguments, e.g. family in case of glm/gam
  pred.dat <- expand.grid(ntrain=(n-i), h=1:i)
  pred.dat$obs <- dat[(n-i+1):n, all.vars(formula(m0))[1], drop=T]
  pred.dat$fit.mod <- predict(m, newdata=dat[(n-i+1):n, ])
  pred.dat$fit.naive <- mean(dat[(n-i-2):(n-i), all.vars(formula(m0))[1], drop=T], na.rm=T)
  pred.dat$e.mod <- pred.dat$obs-pred.dat$fit.mod
  pred.dat$e.naive <- pred.dat$obs-pred.dat$fit.naive
  out <- rbind(out, pred.dat)
}
out <- out %>% 
  mutate(hclass= cut(h, breaks=c(0, 3, 6, 15)),
         hterm=cut(h, breaks=c(0, 3, 6, 15), labels=c("short (1-3)", "medium (4-6)","long (7-15)")))

avg <- out %>% 
  summarise(relmae=mean(abs(e.mod))/mean(abs(e.naive)) )

# summarise by ntrain (nb years for training)

res.ntrain <- out %>% 
  group_by(ntrain) %>% 
  summarise(relmae=mean(abs(e.mod))/mean(abs(e.naive)) )

p <- ggplot(res.ntrain, aes(ntrain, relmae))+
  geom_point(cex=2)+
  geom_line()+
  expand_limits(y=0)+
  geom_hline(yintercept=1, col=2)+
  xlab("Nb years training dataset")+
  ylab(expression(rMAE=MAE[mod]/MAE[naive]))
print(p)

png (filename="MAE_MUT_ntrain.png", width=600, height=600, pointsize = 15)
plot(p)#, page=1
graphics.off()


# summarise by h (prediction time horizon) 

res.h <- out %>% 
  group_by(h) %>% 
  summarise(relmae=mean(abs(e.mod))/mean(abs(e.naive)) )

p <- ggplot(res.h, aes(h, relmae))+
  geom_point(cex=2)+
  geom_line()+
  expand_limits(y=0)+
  geom_hline(yintercept=1, col=2)+
  scale_x_continuous(breaks=seq(1, 15))+
  xlab("Time horizon (years)")+
  ylab(expression(rMAE=MAE[mod]/MAE[naive]))
print(p)

png (filename="MAE_MUT_h.png", width=600, height=600, pointsize = 15)
plot(p)#, page=1
graphics.off()

# summarise by hclass or hterm (grouped time horizons: short, medium, long)
# PLEASE report the values from res.hterm into the excel file

res.hterm <- out %>% 
  group_by(hterm) %>% 
  summarise(relmae=mean(abs(e.mod))/mean(abs(e.naive)) )

p <- ggplot(res.hterm, aes(hterm, relmae))+
  geom_bar(stat="identity", fill="grey")+
  expand_limits(y=0)+
  geom_hline(yintercept=1, col=2)+
  xlab("Time horizon classes")+
  ylab(expression(rMAE=MAE[mod]/MAE[naive]))
print(p)

png (filename="MAE_MUT_hclass.png", width=600, height=600, pointsize = 15)
plot(p)#, page=1
graphics.off()


# save out results
# PLEASE store these objects into the sharepoint

saveRDS(out, file=paste0("output_predskill_",cs,"_",sp,".rds"))
saveRDS(m0, file=paste0("model_predskill_",cs,"_",sp,".rds"))
