

library(rstan)
library(dplyr)
library(tidyverse)
library(lubridate)
library(data.table)
library(here)
library(geosphere)
library(MASS)
library(spacetime)
library(gstat)

output_path = ""

load(paste0(output_path,"beijing_data.RData"))

beijing_data <- beijing_data %>% mutate(int_id=1:nrow(beijing_data)) |> 
  dplyr::filter(pollutant %in% c("PM25","PM10","NO2","O3")) %>% 
  dplyr::filter(!is.na(value)) |> 
  group_by(pollutant) %>%
  filter(value < quantile(value, 0.98)) %>%
  ungroup()

### set up missing ####
### single pollutant missing, take 50 days out #######
set.seed(518)
valid_data=beijing_data %>% 
  group_by(pollutant,station) %>% 
  sample_n(50, replace=FALSE)
test_data=beijing_data %>% 
  anti_join(valid_data, by=c("pollutant","station","date","hour"))


#################### glm to get residuals #############################
beijing_data$date = strptime(
  as.character(beijing_data$date), format = '%Y%m%d',
  tz='UTC'
)

# get the covariates ready for stan model ---------------------------------

beijing_data$day = trunc(beijing_data$date, units='days')


beijing_data$dow = factor(wday(beijing_data$day), levels = 
                            wday(seq(ISOdate(2000,1,3), len=7, by='days')))

beijing_data$pollutant <- as.numeric(factor(beijing_data$pollutant, levels = c("PM25","PM10","NO2","O3")))
beijing_data$station <- as.numeric(factor(beijing_data$station))

# get the distance matrix for the stations, be carefull about the order of stations

station_matrix <- beijing_data %>% dplyr::select(station, lon, lat) %>% unique() %>% arrange(station)



test_data <- beijing_data |> filter(int_id %in% test_data$int_id) 
valid_data <- beijing_data |> filter(int_id %in% valid_data$int_id) 

for(pollutant_index in 1:4)
{
  
  
  ## two step to estimate ####
  model <- glm(value~as.factor(hour)+as.factor(dow), family  = Gamma(link = "log"),data = test_data[test_data$pollutant==pollutant_index,])
  print(gamma.shape(model))
  
  test_data$residuals[test_data$pollutant==pollutant_index] <- residuals(model,type = "response")
  test_data$fitted[test_data$pollutant==pollutant_index] <- model$fitted.values
  
  
  valid_pred <- predict(model, newdata = valid_data[valid_data$pollutant==pollutant_index,], type = "response")
  valid_data$fitted[valid_data$pollutant==pollutant_index] <- valid_pred
  valid_data$residuals[valid_data$pollutant==pollutant_index] <- valid_data$value[valid_data$pollutant==pollutant_index]- valid_pred
  
}
############################

for(pollutant_index in 1:1)
{
  for(hour_index in 0:23)
  {
    
    print(pollutant_index)
    print(hour_index)
    temp_data <- test_data |> filter(pollutant==pollutant_index, hour==hour_index)
    temp_data <- temp_data |> dplyr::select(value, residuals, date, lon, lat) |> as.data.frame()
    temp_data$date <- as.Date(temp_data$date)
    temp_data <- as.data.frame(temp_data)
    
    temp_data_2 <- valid_data |> filter(pollutant==pollutant_index, hour==hour_index) 
    temp_data_2 <- temp_data_2 |> dplyr::select(value, residuals, date, lon, lat) |> as.data.frame()
    temp_data_2$date <- as.Date(temp_data_2$date)
    temp_data_2 <- as.data.frame(temp_data_2)
    
    
    STObj <- stConstruct(x = temp_data , # data set
                         space = c("lon","lat"), # spatial fields
                         time = "date") # time field
    
    STObj_pred <- stConstruct(x = temp_data_2 , # data set
                              space = c("lon","lat"), # spatial fields
                              time = "date") # time field
    
    vv <- variogram(object = residuals ~ 1 , # fixed effect component
                    data = STObj, # 
                    width = 10, # spatial bin (80 km)
                    cutoff = 100, # consider pts < 100 km apart
                    tlags = 0.01:6.01) # 0 days to 6 day
    
    
    metricVgm <- vgmST(stModel = "metric",
                       joint = vgm(10, "Exp", 400, nugget = 0.1),
                       sill = 10,
                       stAni = 100)
    metricVgm <- fit.StVariogram(vv, metricVgm)
    
    
    pred_kriged <- krigeST(residuals ~ 1 , # latitude trend
                           data = STObj, # data set w/o 14 July
                           newdata = STObj_pred, # prediction grid
                           modelList = metricVgm, # semivariogram
                           computeVar = FALSE) # compute variances
    
    valid_data$STkriging[valid_data$pollutant==pollutant_index & valid_data$hour==hour_index]=pred_kriged[[1]]
  }
}


valid_data |> 
  mutate(bias=STkriging-residuals) |> 
  mutate(se=(STkriging-residuals)^2) |> 
  group_by(pollutant) |> 
  summarise(prediction_bias=mean(bias, na.rm = T), prediction_rmse=sqrt(mean(se, na.rm = T)))

