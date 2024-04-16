library("xml2")
library("httr")
library("rvest")
library("ncdf4")
library("dplyr")


#extract all the names of the links from the parent directory
url.main <- "https://www.ncei.noaa.gov/data/land-normalized-difference-vegetation-index/access/"
resource <- GET(url.main)
parse <- htmlParse(resource)
years <- xpathSApply(parse, path = "//a", xmlGetAttr, "href")

yearnames <- as.character(c(1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 
                            1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010,
                            2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023))

#convert each name to a working link (separate by 1900s and 2000s due to structure of link names)
yearlinks1 <- sprintf("%s%s", url.main, sub("^.1*", "1", years[6:24]))
yearlinks2 <- sprintf("%s%s", url.main, sub("^.2*", "2", years[25:48]))
yearlinks <- c(yearlinks1, yearlinks2) 

#nested for loop to download all files for each year
for(i in 1:length(yearlinks)) {
  
  #output provides all links to files
  url <- yearlinks[i]
  res <- GET(url)
  pars <- htmlParse(res)
  daysfiles <- xpathSApply(pars, path = "//a", xmlGetAttr, "href")
  dayslinks <- sprintf("%s%s", url, sub("^.*AVHRR-Land", "AVHRR-Land", daysfiles))
 
   pg <- read_html(url)
  links <- html_attr(html_nodes(pg, "a"), "href")
  #Then identify and subset all of the links with ".nc" in the name
  filenames <- daysfiles[6:length(daysfiles)]
  LINK.nc <- dayslinks[12:length(dayslinks)]
  
  #download the files for each year
  for(j in 1:length(LINK.nc)) {
      
      download.file(LINK.nc[j], destfile = paste0("Canada/NDVI/NOAA_Files/",paste(filenames[i])))

    Sys.sleep(2)  #pause between downloads to not overwhelm server
    
    }
  
  Sys.sleep(2)  #pause between downloads
  
}

#test one file to make sure files have downloaded correctly
file <- "Canada/NDVI/NOAA_Files/AVHRR-Land_v005_AVH13C1_NOAA-07_19810624_c20170610041337.nc"

test <- terra::rast(file)




