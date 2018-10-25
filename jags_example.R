## Checks to see if package rjags is installed. if it is not, package is installed and loaded
rm(list = ls())
if (!"rjags" %in% installed.packages()) install.packages("rjags")
library(rjags)

## Changes working directory and reads in data
working_directory = '/john/home/'
setwd(working_directory)

#--Observed Data---------------------------------
observations = read.csv(paste(working_directory,'/plot_data/plot_observations.csv',sep=''))

#CONVERT TO METRIC UNITS
observations$BA = observations$BA*0.2295
observations$stem_count = observations$stem_count/0.404686
observations$height = observations$height * 0.3048
observations$SI = observations$SI * 0.3048

#FIGURE OUT HOW MANY UNIQUE PLOTS AND ASSIGN ID
observations$PlotID[1]=1
plot_index = 1
for(i in 2:length(observations$PlotID)){
  if(observations$PLOT[i] == observations$PLOT[i-1])
    observations$PlotID[i]=plot_index
  else{
    observations$PlotID[i]=plot_index+1
    plot_index = plot_index + 1
  }
}
## Creates vector to store plotlist
plotlist = unique(observations$PlotID)
plotlist_multi_measure = rep(0,length(plotlist))

## Checks to see which plotlists have multiple measurements, stores binary response in plotlist_multi_measure
for(i in 1:length(plotlist)){
  tmp = length(observations$PlotID[which(observations$PlotID == plotlist[i])])
  if(tmp > 1){
    plotlist_multi_measure[i]=1
  }
}

#PROCESS DATA TO FORMAT FOR BUGS

## Sets plots to be the PlotIDs that have multiple measurements
plots = plotlist[which(plotlist_multi_measure == 1)]
## Sets observations to be the observations from the csv that are have PlotIDs that have multiple measurements
observations =  observations[observations$PlotID %in% plots, ]
## Sets n.plots to be the number of unique plots w multiple measurements
n.plots = length(unique(plots))
## Determines earliest and latest years from the (now subsetted) observations, creates year sequence and determines length
earliest = min(observations$YearMeas)
latest = max(observations$YearMeas)
years = seq(earliest,latest,1)
n.years=length(years)
## Initializes basal area vectors, age, nobs, start and end year
BA.int = rep(NA,n.plots)
BA.int.meas.error = rep(NA,n.plots)
BA.int.sample.error = rep(NA,n.plots)
age = array(0,dim=c(n.plots,n.years))
nobs = rep(0,n.plots)
start.year = rep(0,n.plots)
end.year = rep(0,n.plots)
for(plotnum in 1:n.plots){
  tmp = observations[which(observations$PlotID == plots[plotnum]),]
  plot_youngest = min(tmp$AgeMeas)
  plot_oldest = max(tmp$AgeMeas)
  plot_earliest = min(tmp$YearMeas)
  plot_latest = max(tmp$YearMeas)
  plot_earliest_index = (plot_earliest - earliest)+1
  plot_latest_index = plot_earliest_index + (plot_latest-plot_earliest)
  start.year[plotnum] = plot_earliest_index
  end.year[plotnum] = plot_latest_index
  plot_age = seq(plot_youngest,plot_oldest,1)
  if(plot_latest_index == 24){
  age[plotnum,plot_earliest_index:24] = c(plot_age)
  }else{
    age[plotnum,plot_earliest_index:24] = c(plot_age,seq(from=(plot_oldest+1),to=(plot_oldest)+(24-plot_latest_index),1))    
  }
  nobs[plotnum] = length(tmp$YearMeas)-1
  BA.int[plotnum] = tmp$BA[1]
  BA.int.meas.error[plotnum] = BA.int[plotnum]*0.01
  BA.int.sample.error[plotnum] = BA.int[plotnum]*0.01
  if(plotnum == 4){
    BA.int.sample.error[plotnum] = BA.int[plotnum]*0.1
  }

}

largest_obs = max(nobs)
BA.obs = array(0,dim=c(n.plots,largest_obs))
year = array(0,dim=c(n.plots,largest_obs))
meas.error = array(0,dim=c(n.plots,largest_obs))
sample.error = array(0,dim=c(n.plots,largest_obs))
for(plotnum in 1:n.plots){
  tmp = observations[which(observations$PlotID == plots[plotnum]),]
  BA.obs[plotnum,1:nobs[plotnum]] = tmp$BA[2:(1+nobs[plotnum])]
  meas.error[plotnum,1:nobs[plotnum]] = BA.obs[plotnum,1:nobs[plotnum]]*0.01
  sample.error[plotnum,1:nobs[plotnum]] = BA.obs[plotnum,1:nobs[plotnum]]*0.05
  if(plotnum == 4){
    sample.error[plotnum,1:nobs[plotnum]] = BA.obs[plotnum,1:nobs[plotnum]]*0.05
  }
  for(i in 2:(length(tmp$PlotID))){
    year[plotnum,i-1] = which(years == tmp$YearMeas[i])
  }
}

## Sets sample error for plot4 to be 10 times greater
sample.error[4,nobs[4]] = sample.error[4,nobs[4]]*10
#nobs[4]=0

## Reads in bugs model using sink function
sink("gy.bug")
cat('model {
      #loop through the plots
      for(plotnum in 1:n.plots){
        #yr is the total number of years covered by all plots
        #only model the years with 
        for(yr in (start.year[plotnum]+1):end.year[plotnum]){
          #PROCESS MODEL
          BA.pred[plotnum,yr] <-  BA[plotnum,yr-1]*pow(((1 - exp(-k*(age[plotnum,yr])))/(1-exp(-k*(age[plotnum,yr-1])))),(1/b))
          BA[plotnum,yr] ~ dnorm(BA.pred[plotnum,yr],res.prec)
        }
        #MEASUREMENT AND SAMPLING MODEL
        for(obs in 1:nobs[plotnum]){
          BA.sample[plotnum,obs]~dnorm(BA[plotnum,year[plotnum,obs]],sample.prec[plotnum,obs])
          BA.obs[plotnum,obs]~dnorm(BA.sample[plotnum,obs],obs.prec[plotnum,obs])

          #Each mesurement has its own sample and measurement error
          sample.prec[plotnum,obs] <- pow(sample.sd[plotnum,obs],-2)
          sample.sd[plotnum,obs] <- sample.error[plotnum,obs] 
    
          obs.prec[plotnum,obs] <- pow(obs.sd[plotnum,obs],-2)
          obs.sd[plotnum,obs] <- meas.error[plotnum,obs] 
        }

        #INITIAL CONDITIONS MODEL (required to be separate in a state-space model)
        BA[plotnum,start.year[plotnum]]~dnorm(BA.int.sample[plotnum],int.sample.prec[plotnum])
        BA.int.sample[plotnum] ~ dnorm(BA.int[plotnum],int.obs.prec[plotnum])

        int.sample.prec[plotnum] <- pow(int.sample.sd[plotnum],-2)
        int.sample.sd[plotnum] <- BA.int.sample.error[plotnum] 

        int.obs.prec[plotnum] <- pow(int.obs.sd[plotnum],-2)
        int.obs.sd[plotnum] <- BA.int.meas.error[plotnum] 
      }
      #FORECASTING STEP
      for(yr in (end.year[plotum_for]+1):n.years){
        BA.pred[plotum_for,yr] <- BA[plotum_for,yr-1]*pow(((1 - exp(-k*(age[plotum_for,yr])))/(1-exp(-k*(age[plotum_for,yr-1])))),(1/b))
        BA[plotum_for,yr] ~ dnorm(BA.pred[plotum_for,yr],res.prec)
      }

  #PRIORS ON PARAMETERS
  k ~ dnorm(0.1,1/1000)
  b ~ dnorm(0.1,1/1000)

  #PRIORS ON PROCESS ERRORS
  res.prec<- pow(res.sd,-2)
  res.sd ~ dunif(0, 100)
}'
)
sink()

jags <- jags.model('gy.bug',
                   data = list('n.plots' = n.plots,
                               'start.year' = start.year,
                               'end.year' = end.year,
                               'BA.obs' = BA.obs,
                               'age' = age,
                               'nobs' = nobs,
                               'BA.int' = BA.int,
                               'year' = year,
                               'n.years' = n.years,
                               'plotum_for' = 4,
                               'sample.error' = sample.error,
                               'meas.error' = meas.error,
                               'BA.int.meas.error' = BA.int.meas.error,
                               'BA.int.sample.error' = BA.int.meas.error
                                ),
                   n.chains = 4,
                   n.adapt = 1000)

#burn in, this updates the jags$state()
update(jags,n.iter = 1000)

#sample from posterier starting from the end of the burn in.  The coda.samples track samples for a trace plot
samples = coda.samples(model = jags,
                       variable.names = c('k', 'b','res.sd',
                                          'BA[4,2]','BA[4,3]','BA[4,4]','BA[4,5]','BA[4,6]','BA[4,7]'
                                          ,'BA[4,8]','BA[4,9]','BA[4,10]','BA[4,11]','BA[4,12]','BA[4,13]','BA[4,14]'
                                          ,'BA[4,15]','BA[4,16]','BA[4,17]','BA[4,18]','BA[4,19]','BA[4,20]','BA[4,21]'
                                          ,'BA[4,22]','BA[4,23]','BA[4,24]'
                                          ),
                       n.iter = 1000)

#plot(samples)
#summary(samples)


#PLOT RESULTS FOR PLOT ID = plotum_for

obs_x = observations[8:11,3]
obs_y = observations[8:11,4]

pred_x = seq(min(obs_x),min(obs_x)+22,1)
pred_y = rep(NA,length(pred_x))
pred_y_upper = rep(NA,length(pred_x))
pred_y_lower = rep(NA,length(pred_x))
index=0
for(i in 16:23){
  index = index + 1
  tmp = unlist(samples[1][,i])
  pred_y[index] = median(tmp)
  pred_y_upper[index] = quantile(tmp,0.975)
  pred_y_lower[index] = quantile(tmp,0.0275)
}
for(i in 1:15){
  index = index + 1
  tmp = unlist(samples[1][,i])
  pred_y[index] = median(tmp)
  pred_y_upper[index] = quantile(tmp,0.975)
  pred_y_lower[index] = quantile(tmp,0.0275)
}

k = median(unlist(samples[1][,25]))
b = median(unlist(samples[1][,24]))

model_x = rep(NA,length(pred_x))
model_y= rep(NA,length(pred_x))
model_y_up= rep(NA,length(pred_x))
model_x[1] = obs_x[1]
model_y[1] = obs_y[1]
model_y_up[1] = obs_y[1]
for(i in 2:length(pred_x)){
  model_x[i] = pred_x[i]
  model_y[i] = model_y[i-1]*(((1 - exp(-k*(pred_x[i])))/(1-exp(-k*(pred_x[i-1]))))^(1/b))
}

x_p1 = c(min(pred_x),max(pred_x))
y_p1 = c(obs_y[1],obs_y[1]*(((1 - exp(-k*(x_p1[2])))/(1-exp(-k*(x_p1[1]))))^(1/b)))

(((1 - exp(-k*(41)))/(1-exp(-k*(40))))^(1/b))
(((1 - exp(-k*(21)))/(1-exp(-k*(20))))^(1/b))

plot(obs_x,obs_y,ylim=c(0,60),xlim=range(pred_x),xlab='age',ylab='Basal Area (m2/ha)',pch=20)
points(pred_x,pred_y,type='l',col='red')
points(pred_x,pred_y_upper,type='l',col='blue')
points(pred_x,pred_y_lower,type='l',col='blue')
#points(model_x,model_y,type='l',col='green')
#points(x_p1[2],y_p1[2],type='p',col='green')





