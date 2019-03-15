### Script to extract the Environment Canada daily weather station data
### Jakob Assmann 2 November 2016

require(dplyr)

# First lets create a vector of all the weather station files contained in the folder of interest
# for this we use a variable containg the folder path of interest and the Sys.glob command for
# globbing the file paths (wildcard expansion on file paths)
folder_path <- "scripts/users/jassmann/phenology/temperature_data/komakuk/"
files <- Sys.glob(paste0(folder_path, "eng*"))

# check that it has worked:
files
# nice!

# read in all files and bind them into one data frame bind_rows() is super helpful!
komakuk_weather <- files %>% lapply(read.csv, skip = 25, stringsAsFactors = FALSE) %>% bind_rows()

# excellent now we have one big data frame of all the years. time to get familiar with it.

# collumn names
names(komakuk_weather)

# actually we are only interested in the mean temeprature so lets get rid of the other stuff
komakuk_weather <- komakuk_weather %>% select(Date.Time, Year, Month, Day, Data.Quality, Mean.Temp...C., Mean.Temp.Flag)

# a quick check on rows with estimated values:
# function to extract rows with context (thanks to stackexchange user flodel)
extract.with.context <- function(x, rows, after = 0, before = 0) {
  
  match.idx  <- which(rownames(x) %in% rows)
  span       <- seq(from = -before, to = after)
  extend.idx <- c(outer(match.idx, span, `+`))
  extend.idx <- Filter(function(i) i > 0 & i <= nrow(x), extend.idx)
  extend.idx <- sort(unique(extend.idx))
  
  return(x[extend.idx, , drop = FALSE])
}

# Let's have a look at the flagged rows
flagged_rows <- komakuk_weather[which(komakuk_weather$Mean.Temp.Flag != ""),]
# Get rows with estimated values
estimated <- which(komakuk_weather$Mean.Temp.Flag == "E")
# extract rows with surrounding two columns
estimated_context <- extract.with.context(komakuk_weather, estimated, 2, 2)
# the following dates seem a bit off:
dates_to_be_removed <- c("2005-01-08", "2006-12-01","2009-12-20","2010-06-10","2011-10-05","2012-05-28","2014-02-11","2014-10-10","2016-01-11","2016-07-02")
# set them to be NA
komakuk_weather[komakuk_weather$Date.Time %in% dates_to_be_removed,6] <- NA

# Get rid of all but the date and mean temp columns
komakuk_weather <- komakuk_weather %>% select(-Data.Quality, -Mean.Temp.Flag)

# Now, let's find out how many values are missing
nrow(komakuk_weather[is.na(komakuk_weather$Mean.Temp...C.),])
# of total
nrow(komakuk_weather)
# 698 of 5844 = 12% quite a lot. we will have to do some gap filling.
# maybe another day...

# let's clean up first:
rm(estimated_context, flagged_rows, weather_data, dates_to_be_removed, estimated, files, extract.with.context)

