# download Barrow climate data from ERSL long term observatory
# data website https://www.esrl.noaa.gov/gmd/obop/brw/
for year in $(seq 1994 2016)
do
  echo ${year}
  wget ftp://aftp.cmdl.noaa.gov/data/meteorology/in-situ/brw/met_brw_insitu_1_obop_hour_${year}.txt ;
done
