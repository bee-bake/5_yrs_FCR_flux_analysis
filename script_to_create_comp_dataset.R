#Script to create a comprehensive dataset of daily methane and CO2 fluxes with other environmental and lake drivers
#created by: Bibek Kandel on 2025 Feb 02

#Load libraries
# Load in libraries
pacman::p_load(tidyverse,ncdf4,ggplot2,ggpubr,LakeMetabolizer,zoo,scales,lubridate,
               lognorm,MuMIn,rsq,Metrics,astsa,DescTools,kSamples,rLakeAnalyzer)

###########
# Read in Met file from EDI
met_all <- read_csv("https://pasta.lternet.edu/package/data/eml/edi/389/9/62647ecf8525cdfc069b8aaee14c0478",
                    col_select=c("DateTime","PAR_umolm2s_Average","Rain_Total_mm",
                                  "WindSpeed_Average_m_s")) |> 
  
  mutate(DateTime = force_tz(DateTime, tzone="EST")) |> 
  # Start time series on the 00:15:00 to facilitate 30-min averages
  filter(DateTime >= ymd_hms("2020-04-04 00:15:00"))

# Bind files together if need to use current file
# met_curr <- read_csv("https://raw.githubusercontent.com/FLARE-forecast/FCRE-data/fcre-metstation-data-qaqc/FCRmet_L1.csv",
#                      col_select=c("DateTime","PAR_umolm2s_Average","PAR_Total_mmol_m2"  ,"BP_Average_kPa",
#                                   "AirTemp_C_Average","RH_percent","Rain_Total_mm",
#                                   "ShortwaveRadiationUp_Average_W_m2", "ShortwaveRadiationDown_Average_W_m2",
#                                   "InfraredRadiationUp_Average_W_m2","InfraredRadiationDown_Average_W_m2",
#                                   "Albedo_Average_W_m2","WindSpeed_Average_m_s","WindDir_degrees"))%>%
#   mutate(DateTime = force_tz(DateTime, tzone="EST"))

#met_all <- dplyr::bind_rows(met_curr, met_all) # bind everything together


# Start time series on the 00:15:00 to facilitate 30-min averages

# Select data every 30 minutes from Jan 2020 to end of met data
met_all$Breaks <- cut(met_all$DateTime,breaks = "30 mins",right=FALSE)
met_all$Breaks <- ymd_hms(as.character(met_all$Breaks))

# Average met data to the 30 min mark (excluding Total Rain and Total PAR)
met_30 <- met_all |> 
  group_by(Breaks) |> 
  summarise_all(mean,na.rm=TRUE)

met_daily <- met_30 |> 
  mutate(date=as.Date(DateTime)) |> 
  group_by(date) |> 
  summarise_all(mean,na.rm=TRUE) |> 
  select(-Breaks, -DateTime)
############


#EC fluxes
#Load in EC data
# Can't load this file because it is on your local drive. Put it in the repo. 
ec <- read_csv("C:/Users/13188/Desktop/Reservoirs/Data/DataNotYetUploadedToEDI/EddyFlux_Processing/Eddy_fcr_footprint_full.csv") 

ec2 <- ec |> 
  mutate(DateTime = as.POSIXct(datetime, "%Y-%m-%dT%H:%M:%SZ", tz = "EST"))

#In the following steps, we only included the days where there are more than or equal to 20 half hourly fluxes to calculate daily methane fluxes.
ec_selected <- ec2 |> 
  select(DateTime, ch4_flux_umolm2s) |> 
  filter(DateTime >= "2020-05-01" & DateTime <= "2024-12-31") |> 
  mutate(date = as.Date(DateTime))

# Count the number of data points per day for each variable
  ch4_counts <- ec_selected |> 
  select(date, ch4_flux_umolm2s) |> 
  na.omit() |> 
  group_by(date) |> 
  summarise(count = n(), .groups = "drop")

#Inner join counts with ec_selected
ec_more_than_20 <- inner_join(ch4_counts, ec_selected, by = "date")

#calculate the daily flux values
ec_daily <- ec_more_than_20 |> 
  filter(count >= 20) |>  #this step filters days with half hourly flux values more than or equal to 20
  group_by(date) |> 
  summarise(daily_CH4 = sum(ch4_flux_umolm2s, na.rm = TRUE))


######
#GHG dataset
##Filter the date and depth needed and convert to fluxes
ghg_EDI <- read_csv("https://pasta.lternet.edu/package/data/eml/edi/551/9/98f19e7acae8bea7d127c463b1bb5fbc") |> 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M", tz="EST"))) |>  
  filter(Reservoir == "FCR" & Site == 50) |> 
  filter(DateTime >= "2020-05-01") |>  
  mutate(DateTime = round_date(DateTime, "30 mins")) 

## Will also want to include GHG data from 2024 - forthcoming!
# ghg_current <- read.csv("C:/Users/13188/Desktop/Reservoirs/Data/DataNotYetUploadedToEDI/Raw_GHG/L1_manual_GHG.csv") %>%
#   mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M", tz="EST"))) %>% 
#   filter(Reservoir == "FCR" & Site == 50) %>% 
#   filter(DateTime >= "2020-01-01") %>% 
#   mutate(DateTime = round_date(DateTime, "30 mins")) 

#Combine the datasets
#ghg_FCR <- rbind(ghg_EDI, ghg_current)

## Plot to check
## Methane
ggplot(ghg_EDI,mapping=aes(x=DateTime,y=CH4_umolL))+
  geom_point()+
  geom_line()+
  theme_classic(base_size = 15)

## Carbon dioxide
ggplot(ghg_EDI,mapping=aes(x=DateTime,y=CO2_umolL))+
  geom_point()+
  geom_line()+
  theme_classic(base_size = 15)

## Pivot wider
ghg_fluxes_FCR <- ghg_EDI |> 
  mutate(date = as.Date(DateTime)) |> 
  select(date,Rep,CH4_umolL,CO2_umolL, Depth_m) |> 
  summarize(mean_CH4 = mean(CH4_umolL, na.rm = TRUE),
            mean_CO2 = mean(CO2_umolL, na.rm = TRUE),
            .by = c("date", "Depth_m"))


############
#Load in Catwalk data
catwalk_edi <- read_csv("https://pasta.lternet.edu/package/data/eml/edi/271/9/f23d27b67f71c25cb8e6232af739f986") |>
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="EST"))) |> 
  filter(DateTime >= "2020-05-01" & DateTime <= "2024-12-31")

# catwalk_2024 <- read_csv("https://raw.githubusercontent.com/FLARE-forecast/FCRE-data/fcre-catwalk-data-qaqc/fcre-waterquality_L1.csv") %>%
#   filter(DateTime >= "2024-01-01 00:00:00") %>%
#   mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="EST")))
# 
# catwalk_all <- bind_rows(catwalk_edi, catwalk_2024)%>%
#   filter(DateTime >= as.POSIXct("2020-01-01 01:00:00") & DateTime < as.POSIXct("2024-05-01 01:00:00"))

# Calculate average half-hourly values
# Select data every 30 minutes from Jan 2020 to end of met data
catwalk_edi$Breaks <- cut(catwalk_edi$DateTime,breaks = "30 mins",right=FALSE)
catwalk_edi$Breaks <- ymd_hms(as.character(catwalk_edi$Breaks))

# Average to half-hourly catwalk measurements
catwalk_30 <- catwalk_edi |>
  select(DateTime,RDO_mgL_5_adjusted, RDO_mgL_9_adjusted, EXODO_mgL_1, Breaks, ThermistorTemp_C_1,
         ThermistorTemp_C_5, ThermistorTemp_C_9, EXOfDOM_QSU_1, EXOChla_ugL_1, EXOTemp_C_1) |> 
  group_by(Breaks) |>
  summarise_all(mean,na.rm=TRUE)

catwalk_30 <- catwalk_30 |> 
  select(Breaks,RDO_mgL_5_adjusted, RDO_mgL_9_adjusted, EXODO_mgL_1, ThermistorTemp_C_1,
         ThermistorTemp_C_5, ThermistorTemp_C_9, EXOfDOM_QSU_1, EXOChla_ugL_1, EXOTemp_C_1) |> 
  mutate(Breaks = force_tz(Breaks,"EST"))

names(catwalk_30)[names(catwalk_30) == 'Breaks'] <- 'DateTime'

water_temp_9m <- catwalk_edi |>
  mutate(date = as.Date(DateTime),
         Month = month(date)) |>
  select(date, ThermistorTemp_C_surface, ThermistorTemp_C_1, ThermistorTemp_C_2, ThermistorTemp_C_3,
         ThermistorTemp_C_4, ThermistorTemp_C_5, ThermistorTemp_C_6, ThermistorTemp_C_7, 
         ThermistorTemp_C_8, ThermistorTemp_C_9, EXOfDOM_QSU_1, EXOChla_ugL_1, EXODOsat_percent_1, RDOsat_percent_5_adjusted, RDOsat_percent_9_adjusted) |>
  na.omit() |>
  group_by(date) |>
  summarise(temp_surface = mean(ThermistorTemp_C_surface),
            d1 = mean(ThermistorTemp_C_1),
            d2 = mean(ThermistorTemp_C_2),
            d3 = mean(ThermistorTemp_C_3),
            temp_4m = mean(ThermistorTemp_C_4),
            temp_5m = mean(ThermistorTemp_C_5),
            d6 = mean(ThermistorTemp_C_6),
            d7 = mean(ThermistorTemp_C_7),
            d8 = mean(ThermistorTemp_C_8),
            temp_1.6m = (d1 + d2)/2, #mean of temp at 1m and 2m
            temp_9m = mean(ThermistorTemp_C_9),
            DOsat_surface = mean(EXODOsat_percent_1),
            DOsat_5m = mean(RDOsat_percent_5_adjusted),
            DOsat_9m = mean(RDOsat_percent_9_adjusted),
            mean_Chla = mean(EXOChla_ugL_1),
            mean_fDOM = mean(EXOfDOM_QSU_1)) |>
  select(date, temp_surface, temp_9m, temp_5m, mean_fDOM, mean_Chla, DOsat_surface, DOsat_5m, DOsat_9m, temp_4m, temp_1.6m)


#############
#PREPARE thermocline depth and buoyancy dataset
# Format catwalk temp data for use in LakeAnalyzer in Matlab
catwalk_temp <- catwalk_edi |> 
  select(DateTime,ThermistorTemp_C_surface, ThermistorTemp_C_1, ThermistorTemp_C_2, ThermistorTemp_C_3,
         ThermistorTemp_C_4, ThermistorTemp_C_5, ThermistorTemp_C_6, ThermistorTemp_C_7, ThermistorTemp_C_8,
         ThermistorTemp_C_9) |> 
  group_by(by30 = cut(DateTime, "30 min")) |> 
  summarise_all(mean,na.rm=TRUE) |>
  select(-DateTime) |>
  rename(dateTime = by30, wtr_0.1 = ThermistorTemp_C_surface, wtr_1.0 = ThermistorTemp_C_1, 
         wtr_2.0 = ThermistorTemp_C_2, wtr_3.0 = ThermistorTemp_C_3, wtr_4.0 = ThermistorTemp_C_4,
         wtr_5.0 = ThermistorTemp_C_5, wtr_6.0 = ThermistorTemp_C_6, wtr_7.0 = ThermistorTemp_C_7,
         wtr_8.0 = ThermistorTemp_C_8, wtr_9.0 = ThermistorTemp_C_9)

catwalk_temp$dateTime <- ymd_hms(catwalk_temp$dateTime)

catwalk_temp = data.frame(catwalk_temp)

# Export out for LakeAnalyzer
write.table(catwalk_temp, "C:/Users/13188/Desktop/Data_repository/fcr.wtr", sep='\t', row.names=FALSE)

## Load in data for Lake Analyzer - move path to rLakeAnalyzer folder
#wtr.path = system.file('extdata', 'fcr.wtr', package="rLakeAnalyzer")
fcr_temp = load.ts("C:/Users/13188/Desktop/Data_repository/fcr.wtr")

thermo_depth = ts.thermo.depth(fcr_temp,seasonal=TRUE)

# Replace NaNs with NA
thermo_depth <- thermo_depth |>
  mutate(thermo.depth = ifelse(thermo.depth == "NaN", NA, thermo.depth))

# Plot to look at results
ggplot(thermo_depth,mapping=aes(x=datetime,y=-thermo.depth)) +
  geom_line()+
  xlim(as.POSIXct("2020-04-01"), as.POSIXct("2024-12-31"))

## Calculate N2
n2_freq <- ts.buoyancy.freq(fcr_temp,seasonal=TRUE)

# Remove NA rows in the key column before joining
n2_freq <- n2_freq |> filter(!is.na(datetime))

fcr_results_la <- full_join(thermo_depth,n2_freq,by="datetime")
ggplot()+
  geom_line(data=fcr_results_la, aes(x=datetime,y=thermo.depth))

#Calculate schmidt stability
# Create sample bathymetric data
date = seq(as.POSIXct("2020-05-01 00:10:00"), as.POSIXct("2024-12-31 23:59:00"), by = "30 mins")

bathy <- data.frame(
  depths = c(0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0),
  areas = c(46101,46101,18700,17201,17200,6900,7292,7291,7290)  # Areas for each depth adapted from Hounshell et al. 2020
)
# Export out for LakeAnalyzer
write.table(bathy, "C:/Users/13188/Desktop/Data_repository/bathy.wtr", sep='\t', row.names=FALSE)

## Load in data for Lake Analyzer - move path to rLakeAnalyzer folder
#wtr.path = system.file('extdata', 'fcr.wtr', package="rLakeAnalyzer")
fcr_bathy = load.bathy("C:/Users/13188/Desktop/Data_repository/bathy.wtr")


#bathy <- bathy[bathy_wide$datetime %in% fcr_temp$datetime, ]
sch_stab <- ts.schmidt.stability(fcr_temp,fcr_bathy, na.rm = TRUE)

# Remove NA rows in the key column before joining
sch_stab <- sch_stab |> filter(!is.na(datetime))

fcr_results_la <- full_join(thermo_depth,n2_freq,by="datetime")
fcr_results_la_sch_stab <- full_join(fcr_results_la,sch_stab,by="datetime")

#daily mean
fcr_results_la_sch_stab <- fcr_results_la_sch_stab |>
  mutate(doy = yday(datetime),
         Year= year(datetime)) |>
  na.omit() |>
  group_by(doy,Year) |>
  summarise(mean_sch = mean(schmidt.stability),
            mean_buo = mean(n2),
            mean_thermodepth = mean(thermo.depth))

#plot
plot <- ggplot()+
  geom_line(data = fcr_results_la_sch_stab,aes(x=doy,y=mean_thermodepth,color=as.factor(Year)),size=1)
plot+scale_y_reverse()

write.csv(fcr_results_la_sch_stab,"C:/Users/13188/Desktop/Data_repository/FCR_results_LA.csv")

#############
# Load in buoyancy frequency: currently daily; will need to update if we want to include for hourly!
la <- read_csv("C:/Users/13188/Desktop/Data_repository/FCR_results_LA.csv")

la_daily <- la |>
  group_by(Year, doy) |>
  summarise(thermo.depth = mean(mean_thermodepth,na.rm=TRUE),
            n2 = mean(mean_buo,na.rm=TRUE),
            sch_stab = mean(mean_sch,na.rm=TRUE))

la_daily$date <- as.Date(la_daily$doy - 1, origin = paste0(la_daily$Year, "-01-01"))


####################
#PREPARE HYPOLIMNETIC OXYGEN ADDED DATASET
#O2 data (needs update for 2024)
FCR_O2 <- read_csv("C:/Users/13188/Desktop/Data_repository/O2_data_FCR.csv")

#Change date format
FCR_O2$Date <- strptime(as.character(FCR_O2$time), "%m/%d/%Y")
FCR_O2$Date = as.POSIXct(FCR_O2$Date, "%Y-%m-%dT%H:%M:%SZ", tz = "EST")

#O2 added per day
FCR_O2_day <- FCR_O2 |>
  select(Date, o2_kg_d) |>
  rename(date = Date) |>
  filter(date >= "2020-05-01" & date < "2025-01-01")


##########

#FINALLY MERGE ALL DATASETS TO CREATE A COMPREHENSIVE DAILY DATASET
#Select the columns of interest to avoid unnecessary columns
#we want ghg with depths 
ghg_selected <- ghg_fluxes_FCR |>
  select(date,Depth_m, mean_CO2,mean_CH4) |>
  filter(Depth_m == 0.1 | Depth_m == 9) |> #we only want ghgs at 0.1m and 9m!
  pivot_longer(cols = c(mean_CO2, mean_CH4), names_to = "Category", values_to = "Value") |>
  pivot_wider(names_from = c(Category,Depth_m), values_from = Value)


#we want single daily value of thermocline depth, buoyancy freq and schmidt stability
la_daily_selected <- la_daily |>
  select(date, thermo.depth, n2, sch_stab)


# List of datasets to join
datasets <- list(la_daily_selected, water_temp_9m, FCR_O2_day, ec_daily, ghg_selected, met_daily)


# Join multiple datasets using reduce
comp_data <- reduce(datasets, left_join, by = "date")


# Convert the Date column to Date class
comp_data$date <- as.Date(comp_data$date, format = "%Y-%m-%d")


#write the dataset to csv
write.csv(comp_data,"C:/Users/13188/Desktop/Data_repository/comp_data.csv")



  









