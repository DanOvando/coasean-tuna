

MakeBoundingBoxShapefile <- function(path, data, file_name){
  
  # require the necessary packages
  require(sp)
  require(maptools)
  require(rgdal)
  
  # read in file
  A <- data
  
  Q <- unique(A[,c('latitude_north', 'latitude_south', 'longitude_west', 'longitude_east')])
  
  # set name for shapefile
  shp_name <- paste(file_name, '.shp', sep = '')
  name <- file_name
  
  # check to see if shapefile already exists
  if(file.exists(paste(path, '/', name, '/', shp_name, sep = '') )){
    cat(shp_name, " already exists...manually remove before proceeding\n\n")
  }
  
  
  # index vector to loop over coordinates with
  ind <- 1:nrow(Q)
  
  cat("generating shapefile..\n\n")
  
  # loop over index to Build shapefile
  for(i in ind) {
    # i <- ind[1]
    x <- Q[i,]
    
    # generate polygons for any occurrences  with coordinates
    if(!any(is.na(c( x["latitude_north"], 
                     x["latitude_south"], 
                     x["longitude_west"], 
                     x["longitude_east"])))){
      
      # create matrix with the bounding box coordinates, fucked up way, but only way?
      
      bb.coords <- matrix(unlist(list(c(x["longitude_west"], 
                                        x["longitude_east"], 
                                        x["longitude_east"], 
                                        x["longitude_west"], 
                                        x["longitude_west"],
                                        x["latitude_north"], 
                                        x["latitude_north"], 
                                        x["latitude_south"], 
                                        x["latitude_south"], 
                                        x["latitude_north"]))) , ncol=2)
      
      # make the spatial bounding box
      sp <-  SpatialPolygons(list(Polygons(list(Polygon(bb.coords)), i)), proj4string=CRS("+proj=longlat +datum=WGS84"))
      spProj <- spTransform(sp, CRS("+proj=longlat +datum=WGS84") )
      
      # # add attributes
      sp.data <- SpatialPolygonsDataFrame(spProj, data.frame(x, 
                                                             row.names=i)) 
      
      # bind the shapefiles into one unit
      if(!exists("master.sp") ) {
        master.sp <- sp.data
       } else{
         master.sp <- spRbind(master.sp, sp.data)
      }
    }
  }
  cat("saving new shapefile...\n\n")
  # write the results to a shapefile
  writeOGR(master.sp, paste(path, '/', name, '/', sep = ''), name, "ESRI Shapefile")
  
}