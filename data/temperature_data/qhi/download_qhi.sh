for year in `seq 2001 2017`
do wget --content-disposition "http://climate.weather.gc.ca/climate_data/bulk_data_e.html?format=csv&stationID=1560&Year=${year}&Day=14&timeframe=2&submit= Download+Data"
done
