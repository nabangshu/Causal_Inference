library(RMySQL)
library(reshape2)
library(plyr)
library(igraph)
library(ggplot2)
library(lattice)
library(lubridate)
library(CausalImpact)
library(stringr)
library(zoo)
library(ergm)
library(igraph)
library(network)
library(intergraph)
library(readxl)
library(dplyr)
library(beepr)
library(lmtest)
library(sandwich)
library(clubSandwich)
library(lubridate)

gcenter <- function(df1,group) {
  x <-ncol(df1)
  variables <- paste(rep("C", x), colnames(df1), sep=".")
  copydf <- df1
  for (i in 1:ncol(df1)) {
    copydf[,i] <- df1[,i] - ave(df1[,i], group,FUN=mean)}
  colnames(copydf) <- variables
  return(cbind(df1,copydf))
}

clx <- function(fm, dfcw, cluster){
  M <- length(unique(cluster))
  N <- length(cluster)
  dfc <- (M/(M-1))*((N-1)/(N-fm$rank))
  u <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum))
  vcovCL <- dfc*sandwich(fm, meat=crossprod(u)/N)*dfcw
  return(coeftest(fm, vcovCL))
}


data <- read.csv("./data.csv")
mem <- read.csv("./mem.csv")

# Data preprocessing
mem$age<- round(as.numeric((as.Date("2013-01-01")-as.Date(mem$birthday))/365.25), digits = 0)
mem$reg_yr <- year(mem$reg_date)
data$age <- round(as.numeric((as.Date("2013-01-01")-as.Date(data$birthday))/365.25), digits = 0)
data$reg_yr <- year(data$reg_date)
data$in_nw <- ifelse(is.na(data$follower_reg_date) & is.na(data$followed_reg_date), 0,1)
data$leave <- ifelse(is.na(data$leave_date),0,1)
data$count <- ifelse(data$order_cnt>4,0,1)
mem$in_nw <- ifelse(is.na(mem$follower_reg_date) & is.na(mem$followed_reg_date), 0,1)
mem$leave <- ifelse(is.na(mem$leave_date),0,1)
mem$count <- ifelse(mem$order_cnt>4,0,1)

members=list()
for(var in 2013:2019)
{
  members[[var-2012]]<-unique(data[which(data$yr==var),] %>% select(mem_no))
}

filtered=list()
filtered[[1]]<- members[[7]] %>% inner_join(members[[6]], 
                                            by="mem_no")%>% inner_join(members[[5]],by="mem_no")%>% inner_join(members[[4]], by="mem_no")%>% inner_join(members[[3]], 
                                                                                                                                                        by="mem_no")%>% inner_join(members[[2]], by="mem_no")%>% inner_join(members[[1]], by="mem_no")
for(var in 2:7)
{
  x<- mem[which(mem$reg_yr==2012+var),]
  y<-x
  y<-inner_join(y,x)
  for(i in var:7)
  {
    y<-inner_join(y,members[[i]], by="mem_no")
  }
  filtered[[var]]<-select(y,mem_no)
}

active_members_list <- rbind(filtered[[7]],filtered[[6]])%>%rbind(filtered[[5]])%>%rbind(filtered[[4]])%>%rbind(filtered[[3]])%>%rbind(filtered[[2]])%>%rbind(filtered[[1]])%>%merge(mem[which(as.Date(mem$reg_date)<as.Date("2016-01-25")),])%>%select(mem_no)

active_members <- merge(active_members_list, mem, by= "mem_no")

purchse_data <- merge(active_members_list, data, by="mem_no")
purchse_data <- purchse_data[which(as.Date(purchse_data$order_date)>=as.Date("2013-01-01")&as.Date(purchse_data$reg_date)<as.Date("2016-11-25")),]

###Start of Causal Inference
pre <- as.Date(c("2013-01-01","2016-01-25"))
post <- as.Date(c("2016-01-26","2019-01-01"))

# Using Google's Causal Impact package
# Synthetic control group
p <- group_by(purchse_data[which(purchse_data$in_nw==1&purchse_data$age>=40&
                                   purchse_data$age<50&purchse_data$gender=="W"&purchse_data$married=="Y"),],order_date)%>%summarise(innwsales = sum(rev_krw_sum)/100000,num=length(unique(mem_no)),count=sum(rev_krw_cnt))
innw <- cbind( p[,1], p[,4])


pp <- group_by(purchse_data[which(purchse_data$in_nw==0&purchse_data$age>=40&
                                    purchse_data$age<50&purchse_data$gender=="W"&purchse_data$married=="Y"),],order_date)%>%summarise(notinnwsales = sum(rev_krw_sum)/100000,num=length(unique(mem_no)),count=sum(rev_krw_cnt))

pp2 <- group_by(purchse_data[which(purchse_data$in_nw==0&purchse_data$age>=40&
                                     purchse_data$age<50&purchse_data$gender=="W"&purchse_data$married=="N"),],order_date)%>%summarise(notinnwsales = sum(rev_krw_sum)/100000,num=length(unique(mem_no)),count=sum(rev_krw_cnt))

pp3 <- group_by(purchse_data[which(purchse_data$in_nw==0&purchse_data$age>=30&
                                     purchse_data$age<40&purchse_data$gender=="W"&purchse_data$married=="Y"),],order_date)%>%summarise(notinnwsales = sum(rev_krw_sum)/100000,num=length(unique(mem_no)),count=sum(rev_krw_cnt))

pp4 <- group_by(purchse_data[which(purchse_data$in_nw==0&purchse_data$age>=30&
                                     purchse_data$age<40&purchse_data$gender=="W"&purchse_data$married=="N"),],order_date)%>%summarise(notinnwsales = sum(rev_krw_sum)/100000,num=length(unique(mem_no)),count=sum(rev_krw_cnt))
pp5 <- group_by(purchse_data[which(purchse_data$in_nw==0&purchse_data$age>=50&
                                     purchse_data$age<60&purchse_data$gender=="W"&purchse_data$married=="Y"),],order_date)%>%summarise(notinnwsales = sum(rev_krw_sum)/100000,num=length(unique(mem_no)),count=sum(rev_krw_cnt))
pp6 <- group_by(purchse_data[which(purchse_data$in_nw==0&purchse_data$age>=50&
                                     purchse_data$age<60&purchse_data$gender=="W"&purchse_data$married=="N"),],order_date)%>%summarise(notinnwsales = sum(rev_krw_sum)/100000,num=length(unique(mem_no)),count=sum(rev_krw_cnt))

pp11 <- group_by(purchse_data[which(purchse_data$in_nw==0&purchse_data$age>=50&
                                      purchse_data$age<60&purchse_data$gender=="M"&purchse_data$married=="Y"),],order_date)%>%summarise(notinnwsales = sum(rev_krw_sum)/100000,num=length(unique(mem_no)),count=sum(rev_krw_cnt))

notinnw <- cbind( pp[,1], pp[,4])
notinnw2 <- cbind( pp2[,1], pp2[,4])
notinnw3 <- cbind( pp3[,1], pp3[,4])
notinnw4 <- cbind( pp4[,1], pp4[,4])
notinnw5 <- cbind( pp5[,1], pp5[,4])
notinnw6 <- cbind( pp6[,1], pp6[,4])
notinnw11 <- cbind( pp11[,1], pp11[,4])
df1<- inner_join(notinnw3, notinnw,by="order_date")
df2<- inner_join(notinnw3, notinnw2,by="order_date")
df3<- inner_join(notinnw3, notinnw4,by="order_date")
df4<- inner_join(notinnw3, notinnw5,by="order_date")
df5<- inner_join(notinnw3, notinnw6,by="order_date")
df10<- inner_join(notinnw3, notinnw11,by="order_date")
df6 <- inner_join(notinnw3, innw,by="order_date")

df.m <- cbind(df1[,1],df6[,3],df1[,2],df2[,3],df4[,3],df5[,3],df10[,3])
df.m <- as.data.frame(df.m)
colnames(df.m) <- c("order_date","y","x1","x2","x3","x4","x5")
df.m$order_date <- as.Date(df.m$order_date)
df.m$y <- as.numeric(df.m$y)
df.m$x1 <- as.numeric(df.m$x1)
df.m$x2 <- as.numeric(df.m$x2)
df.m$x3 <- as.numeric(df.m$x3)
df.m$x4 <- as.numeric(df.m$x4)
df.m$x5 <- as.numeric(df.m$x5)
pre <- as.Date(c("2013-01-01","2016-01-25"))
post <- as.Date(c("2016-01-26","2019-01-01"))
dates <- as.Date(df.m$order_date)
data.impact <- zoo(cbind(df.m$y,df.m$x1,df.m$x2,df.m$x3,df.m$x4,df.m$x5),dates)
impact <- CausalImpact(data.impact,pre,post,alpha=0.1,model.args = list(niter = 1000, standardize.data=TRUE))
summary(impact)
plot(impact)

## Comparing with vanilla DD with checks for heteroskedasticity

did <- purchse_data[which(purchse_data$in_nw==1&purchse_data$age>=40&purchse_data$age<50&purchse_data$gender=="W"&purchse_data$married=="Y"),]
data_did <- cbind(did[,c(1,24,25,30)])
data_did$order_date <- as.Date(data_did$order_date)
data_did <- gcenter(data_did, data_did$mem_no)
data_did$did <- data_did$time*data_did$in_nw
fit_did <- lm(rev_krw_sum/100000~order_date, data = data_did)
M <- length(unique(data_did$mem_no))
dfcw <- fit_did$df / (fit_did$df - (M - 1))
print(summary(fit_did))
clx(fit_did,dfcw,data_did$mem_no)


#checking with count of sales

fit_did <- lm(rev_krw_cnt~order_date, data = data_did)
M <- length(unique(data_did$mem_no))
dfcw <- fit_did$df / (fit_did$df - (M - 1))
print(summary(fit_did))
clx(fit_did,dfcw,data_did$mem_no)