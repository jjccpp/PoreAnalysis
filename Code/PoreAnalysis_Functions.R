# Function to get asc MLT files to tiff 
# Author: Julio Pachon
#


#Motivation ####
# Obtain porosity information from a from the Multistripe Laser Triangulation.  This is an Rmarkdown HEAVILY based on Daniel Hirmas's 20220225_ProximalSensingProject_Functions.R and made for the purpose of cutting down the time and parallelizing the function
# 

#modification to ensure that the x, y grid includes all rows and there is no N/As is adding +resol to xdim and y dim
#xdim <- max(nmlt$x)+resol # width
#ydim <- max(nmlt$y)+resol # length





library(compiler)
library(tiff)
library(matlab)



#
#
#
Aggregating_asc <- function(ascname, resol=0.18, dname = NULL, oname = NULL, date = FALSE) {
  
  # ascname: character; file name of the point cloud data in the .asc format
  # outputted from RapidWorks; this file contains columns of data, the 
  # first three of which, refer to the x, y, z coordinates of the points;
  # additional columns refer to the RGB values associated with these 
  # points and the normals associated with these points
  #
  # resol: numeric; value that represents the resolution of the resulting
  # matrix grid in mm
  #
  # dname: character; folder where the asc file is stored if not in the
  # the default folder
  #
  # oname: character; folder where the tiff file is to be outputted if 
  # not in the the default folder
  #
  # date: logical; TRUE means output the system date in the tiff file
  # name (default); FALSE means not to output the date in the file name
  # 
  # NOTE: There are two large loops in this function whose progress are
  # displayed as two progress bars
  #
  #
  # Require the tiff package:
  require(tiff)
  #
  # Save the file name without the extension for later use:
  mltfile <- sub('.asc', '', ascname)
  #
  ### Read in the asc file: ####
  if (length(dname)==0) {
    mlt <- read.table(ascname,header=FALSE,sep="\t")
    if( length(names(mlt)) == 9){
      names(mlt) <- c("x","y","z","R","G","B","a","b","c")}else{
        names(mlt) <- c("x","y","z","R","G","B")
      }
  } else { # file is stored in another folder
    currentwd <- getwd()
    setwd(dname)
    mlt <- read.table(ascname,header=FALSE,sep="\t")
    if( length(names(mlt)) == 9){
      names(mlt) <- c("x","y","z","R","G","B","a","b","c")}else{
        names(mlt) <- c("x","y","z","R","G","B")
      }
    setwd(currentwd) # resets to the original working directory
  }
  #
  ### Define the origin based on the minimum data in each column: ####
  nmlt <- mlt # save the "normalized" mlt values into a new dataframe
  nmlt$x <- nmlt$x-min(mlt$x)
  nmlt$y <- nmlt$y-min(mlt$y)
  nmlt$z <- nmlt$z-min(mlt$z)
  #
  ### Calculate the spatial extent of the data (in mm): ####
  xdim <- max(nmlt$x)+resol # width
  ydim <- max(nmlt$y)+resol # length
  #
  # Setup a grid based on the extent data:
  xbound <- seq(0,xdim,by=resol) # provides the x coordinate for the left 
  # grid cell boundary
  ybound <- seq(0,ydim,by=resol) # provides the y coordinate for the upper 
  # grid cell boundary
  
  # Assign the column values of mltmat to a new dataframe, xmlt:
  xmlt <- data.frame(nmlt,row=NA,col=NA) # switch of variables
  xmlt <- xmlt[order(xmlt$x),] # sort the dataframe on the x values
  row.names(xmlt) <- NULL
  x <- round(xmlt$x,round(-log10(resol/10))) # rounds the x data to a tenth
  # of the specified resolution; if the resolution is greater than 10 this
  # will cap the rounding to a single mm
  
  dummy_mat <- matrix(0,nrow=length(ybound),ncol=length(xbound)) # matrix to
  # store the answer to the question: Is there a point in this grid cell? A:
  # Yes (1) or No (0); Note: columns = x, rows = y
  
  pb <- txtProgressBar(min = 0, max = (ncol(dummy_mat)-1), style = 3) # create
  # progress bar
  for (i in 1:(ncol(dummy_mat)-1)) {
    xmlt$col[which((x > xbound[i]) & (x <= xbound[i+1]))] <- i 
    Sys.sleep(0.1)
    setTxtProgressBar(pb, i) # update progress bar
  }
  close(pb)
  #
  # Assign the row values of mltmat to xmlt:
  xmlt <- xmlt[order(xmlt$y),] # sort the dataframe on the y values
  row.names(xmlt) <- NULL
  y <- round(xmlt$y,round(-log10(resol/10)))
  pb <- txtProgressBar(min = 0, max = (nrow(dummy_mat)-1), style = 3) # create
  # progress bar
  for (i in 1:(nrow(dummy_mat)-1)) {
    xmlt$row[which((y > ybound[i]) & (y <= ybound[i+1]))] <- i
    Sys.sleep(0.1)
    setTxtProgressBar(pb, i) # update progress bar
  }
  close(pb)
  #
  # Create a new rowcol dataframe and order by rows then cols:
  rowcol <- xmlt[,c("row","col")]
  rowcol <- rowcol[order(rowcol$row, rowcol$col),]
  row.names(rowcol) <- NULL
  
  
  rowcol <- unique(rowcol) #rowcol[!duplicated(rowcol),] # removes any duplicated cells
  row.names(rowcol) <- NULL
  
  #Aggregating data ####
  
  # Calculates the mean height (z) and R, G, B for each cell:
  xmlt <- xmlt[order(xmlt$row, xmlt$col),]
  row.names(xmlt) <- NULL
  xmlt <- data.frame(xmlt, # used to aggregate the data below
                     matid = as.factor(paste(xmlt$row,xmlt$col,sep="_")))
  xmlta <- aggregate(cbind(z, R, G, B) ~ matid, xmlt, FUN = "mean")
  xmlta$R <- round(xmlta$R,0) # rounds to the nearest integer
  xmlta$G <- round(xmlta$G,0) # rounds to the nearest integer
  xmlta$B <- round(xmlta$B,0) # rounds to the nearest integer
  xmltm <- merge(xmlt[,c("row","col","matid")],xmlta,by="matid") # holds
  
  
  
  ### the row and col indices and the mean values for z, R, G, and B ####
  xmltm <- xmltm[,c("row","col","z","R","G","B")] # removes the matid col
  
  
  xmltm <- xmltm[-which(is.na(xmltm$col) | is.na(xmltm$row)),] # removes the rows or columns with NA values
  
  
  ### removes any duplicated cells####
  xmltm <- xmltm[!duplicated(xmltm),] 
  
  
  
  xmltm <- xmltm[order(xmltm$row, xmltm$col),] # order the data.frame
  row.names(xmltm) <- NULL
  
  
  if (length(oname)==0) {
    if (date == FALSE) {
      save(xmltm, ybound, xbound , rowcol,mltfile,
           file= paste(mltfile,'_',
                       resol,
                       'mmResolution',
                       '_aggregatedASCData.Rdata',sep=''))
    } else { # include the date in the file name
      save(xmltm, ybound, xbound , rowcol,mltfile,
           file= paste(gsub("-", "", as.character(Sys.Date())),
                       '_',mltfile,'_',resol,'mmResolution',
                       '_aggregatedASCData.Rdata',sep=''))
    }
  } else { # store the tiff in a different folder
    currentwd <- getwd()
    setwd(oname)
    if (date == FALSE) {
      save(xmltm, ybound, xbound , rowcol,mltfile,
           file=paste(mltfile,'_',resol,'mmResolution',
                      '_aggregatedASCData.Rdata',sep=''))
    } else { # include the date in the file name
      save(xmltm, ybound, xbound , rowcol,mltfile,
           file=paste(gsub("-", "", as.character(Sys.Date())),
                      '_',mltfile,'_',resol,'mmResolution',
                      '_aggregatedASCData.Rdata',sep=''))
    }
    setwd(currentwd) # resets to the original working directory
  }
  
  #return(xmltm, xbound, ybound, rowcol) 
  
}


asctotiffs <- function(rdata_file,resol=0.18, tray_depth_cm=4, depth_from_tray_bottom_mm=15, dname = NULL, oname = NULL, date = FALSE) {
  
  load(rdata_file)
  
  mltfile <- substr(rdata_file, 25, nchar(rdata_file)-13 )
  
  library(matlab)
  library(compiler)
  require(tiff)
  
  
  # Function for rotating a matrix #####
  rotate <- function(x,direction="clockwise") {
    # Rotates a matrix either clockwise or counter clockwise as specified
    # by direction
    # x = matrix
    # direction = either 'clockwise' or 'counter' (any value other than
    # 'clockwise' will rotate counter clockwise)
    if (direction == "clockwise") {
      t(apply(x, 2, rev))
    } else {
      apply(t(x), 2, rev)
    }
  }
  
  cmp_rotate <- cmpfun(rotate)
  
  
  # Place the mean values into the correct matrices: ####
  ##RGB matrix####
  rgbmat <- array(NA,c(length(ybound),length(xbound),3)) # array to
  # store the answer to the question: What is the color of this grid cell?
  # Note: columns = x, rows = y
  
  # the single index notation for matrices
  rgbmat <- rgbmat[1:(nrow(rgbmat[,,1])-1),1:(ncol(rgbmat[,,1])-1),]
  rgbmat[,,1][xmltm$row + nrow(rgbmat) * (xmltm$col - 1)] <- xmltm$R # uses 
  # the single index notation for matrices
  rgbmat[,,2][xmltm$row + nrow(rgbmat) * (xmltm$col - 1)] <- xmltm$G # uses 
  # the single index notation for matrices
  rgbmat[,,3][xmltm$row + nrow(rgbmat) * (xmltm$col - 1)] <- xmltm$B # uses 
  # the single index notation for matrices
  
  
  dimhold <- dim(cmp_rotate(t(rgbmat[,,1]), direction = 'counter'))
  colormat <- array(NA,c(dimhold[1],dimhold[2],3)) # creates a new array to
  # hold the correctly oriented rgb values
  colormat[,,1] <- cmp_rotate(t(rgbmat[,,1]), direction = 'counter')/255
  colormat[,,2] <- cmp_rotate(t(rgbmat[,,2]), direction = 'counter')/255
  colormat[,,3] <- cmp_rotate(t(rgbmat[,,3]), direction = 'counter')/255
  
  
  
  ### Save the resulting matrices as tiff files: ####
  if (length(oname)==0) {
    if (date == FALSE) {
      writeTIFF(colormat,paste(mltfile,'_',resol,'mmResolution',
                               '_RGB.tiff',sep=''))
    } else { # include the date in the file name
      writeTIFF(colormat,paste(gsub("-", "", as.character(Sys.Date())),
                               '_',mltfile,'_',resol,'mmResolution',
                               '_RGB.tiff',sep=''))    }
  } else { # store the tiff in a different folder
    currentwd <- getwd()
    setwd(currentwd)
    if (date == FALSE) {
      writeTIFF(colormat,paste(mltfile,'_',resol,'mmResolution',
                               '_RGB.tiff',sep=''))
    } else { # include the date in the file name
      writeTIFF(colormat,paste(gsub("-", "", as.character(Sys.Date())),
                               '_',mltfile,'_',resol,'mmResolution',
                               '_RGB.tiff',sep=''))
    }
    setwd(currentwd) # resets to the original working directory
  }
  
  ##Macropore matrix ####
  # matrix to store the answer to the question: Is there a point in this grid cell? A:
  # Yes (1) or No (0); Note: columns = x, rows = y
  mltmat <- matrix(0,nrow=length(ybound),ncol=length(xbound)) 
  
  
  # Assign zero values in mltmat for the appropriate missing indices in 
  # rowcol:
  mltmat <- mltmat[1:(nrow(mltmat)-1),
                   1:(ncol(mltmat)-1)] # remove the last row and col 
  # in mltmat; Note: these are removed because they represent the right 
  # or bottom extent of the last full cells
  
  # uses the single index notation for matrices
  mltmat[rowcol$row + nrow(mltmat) * (rowcol$col - 1)] <- 1 
  #
  # Fix the orientation of the matrix for the tiff below:
  mltmat <- cmp_rotate(t(mltmat), direction = 'counter')
  #
  
  
  
  #DEM matrix####
  
  #view the 
  
  demmat <- matrix(NA,nrow=length(ybound),ncol=length(xbound)) # matrix to
  # store the answer to the question: What is the z of this grid cell?
  # Note: columns = x, rows = y
  
  demmat <- demmat[1:(nrow(demmat)-1),
                   1:(ncol(demmat)-1)] # remove the last row and col 
  # in demmat; Note: these are removed because they represent the right 
  # or bottom extent of the last full cells
  demmat[xmltm$row + nrow(demmat) * (xmltm$col - 1)] <- xmltm$z # uses 
  
  
  # Fix the orientation of the following matrices to match mltmat:
  demmat <- cmp_rotate(t(demmat), direction = 'counter')
  dimhold <- dim(cmp_rotate(t(rgbmat[,,1]), direction = 'counter'))
  colormat <- array(NA,c(dimhold[1],dimhold[2],3)) # creates a new array to
  # hold the correctly oriented rgb values
  colormat[,,1] <- cmp_rotate(t(rgbmat[,,1]), direction = 'counter')/255
  colormat[,,2] <- cmp_rotate(t(rgbmat[,,2]), direction = 'counter')/255
  colormat[,,3] <- cmp_rotate(t(rgbmat[,,3]), direction = 'counter')/255
  #
  ### Assigns color values to the DEM matrix to output it as a tiff: ####
  zRGB <- col2rgb(jet.colors(100)[as.numeric(cut(xmltm$z,breaks=100))])
  zcol <- array(NA,c(dim(colormat)[1],dim(colormat)[2],dim(colormat)[3]))
  zcol[,,1][xmltm$row + nrow(zcol) * (xmltm$col - 1)] <- zRGB[1,] # uses 
  # the single index notation for matrices (R)
  zcol[,,2][xmltm$row + nrow(zcol) * (xmltm$col - 1)] <- zRGB[2,] # uses 
  # the single index notation for matrices (G)
  zcol[,,3][xmltm$row + nrow(zcol) * (xmltm$col - 1)] <- zRGB[3,] # uses 
  # the single index notation for matrices (B)
  dimhold <- dim(cmp_rotate(t(zcol[,,1]), direction = 'counter'))
  
  
  
  
  demcol <- array(NA,c(dimhold[1],dimhold[2],3)) # creates a new array to
  
  
  
  # hold the correctly oriented rgb values
  demcol[,,1] <- cmp_rotate(t(zcol[,,1]), direction = 'counter')/255
  demcol[,,2] <- cmp_rotate(t(zcol[,,2]), direction = 'counter')/255
  demcol[,,3] <- cmp_rotate(t(zcol[,,3]), direction = 'counter')/255
  #
  
  
  ###view original DEM#### 
  par(mfrow = c(1,2))
  image(as.matrix(demmat), useRaster=TRUE, axes=FALSE, col = hcl.colors(500, "spectral"))
  
  
  #Correct DEM ####
  #With the image above, you can see there is a gradient, use code below to correct it.
  
  
  rightMost <- dim(demmat)[2]
  bottomMost <- dim(demmat)[1]
  Top_left <- max(demmat[1:300,1:300], na.rm = TRUE)
  Top_right <- max(demmat[1:300,(rightMost-300):rightMost], na.rm = TRUE)
  Bottom_left <- max(demmat[(bottomMost-300):bottomMost,1:300], na.rm = TRUE)
  Bottom_right <- max(demmat[(bottomMost-300):bottomMost,(rightMost-300):rightMost], na.rm = TRUE)
  
  Reg_table <- data.frame()
  
  Reg_table[1,1:3] <-   xmltm[xmltm$z==Top_left,1:3]
  Reg_table[2,1:3] <- xmltm[xmltm$z==Top_right,1:3]
  Reg_table[3,1:3] <- xmltm[xmltm$z==Bottom_left,1:3]
  Reg_table[4,1:3] <- xmltm[xmltm$z==Bottom_right,1:3]
  
  reg_eq <- lm(z~row+col,Reg_table)
  
  Modified_Z <- xmltm$z-predict(reg_eq,xmltm)+reg_eq[[1]][1]
  
  
  mod_demmat <- matrix(NA,nrow=length(ybound),ncol=length(xbound)) # matrix to
  # store the answer to the question: What is the z of this grid cell?
  # Note: columns = x, rows = y
  
  mod_demmat <- mod_demmat[1:(nrow(mod_demmat)-1),
                           1:(ncol(mod_demmat)-1)] # remove the last row and col 
  # in mod_demmat; Note: these are removed because they represent the right 
  # or bottom extent of the last full cells
  
  
  mod_demmat[xmltm$row + nrow(mod_demmat) * (xmltm$col - 1)] <- Modified_Z  
  mod_demmat <- cmp_rotate(t(mod_demmat), direction = 'counter')
  
  ###view corrected DEM####
  #A sligh gradient can still be seen but keep in mind that the rulers are not perfectly straight and there are other sources of uncertainty.
  image(as.matrix(mod_demmat), useRaster=TRUE, axes=FALSE, col = hcl.colors(500, "spectral"))
  
  
  ###Create Corrected DEM tiff file ####
  
  colormat <- array(NA,c(dimhold[1],dimhold[2],3)) # creates a new array to
  # hold the correctly oriented rgb values
  colormat[,,1] <- cmp_rotate(t(rgbmat[,,1]), direction = 'counter')/255
  colormat[,,2] <- cmp_rotate(t(rgbmat[,,2]), direction = 'counter')/255
  colormat[,,3] <- cmp_rotate(t(rgbmat[,,3]), direction = 'counter')/255
  #
  ### Assigns color values to the mod_dem matrix to output it as a tiff: ####
  zRGB <- col2rgb(jet.colors(100)[as.numeric(cut(Modified_Z,breaks=100))])
  zcol <- array(NA,c(dim(colormat)[1],dim(colormat)[2],dim(colormat)[3]))
  zcol[,,1][xmltm$row + nrow(zcol) * (xmltm$col - 1)] <- zRGB[1,] # uses 
  # the single index notation for matrices (R)
  zcol[,,2][xmltm$row + nrow(zcol) * (xmltm$col - 1)] <- zRGB[2,] # uses 
  # the single index notation for matrices (G)
  zcol[,,3][xmltm$row + nrow(zcol) * (xmltm$col - 1)] <- zRGB[3,] # uses 
  # the single index notation for matrices (B)
  dimhold <- dim(cmp_rotate(t(zcol[,,1]), direction = 'counter'))
  
  
  
  
  mod_demcol <- array(NA,c(dimhold[1],dimhold[2],3)) # creates a new array to
  
  
  
  # hold the correctly oriented rgb values
  mod_demcol[,,1] <- cmp_rotate(t(zcol[,,1]), direction = 'counter')/255
  mod_demcol[,,2] <- cmp_rotate(t(zcol[,,2]), direction = 'counter')/255
  mod_demcol[,,3] <- cmp_rotate(t(zcol[,,3]), direction = 'counter')/255
  #
  
  
  
  
  ### Save the resulting Corrected DEM as tiff files: ####
  if (length(oname)==0) {
    if (date == FALSE) {
      writeTIFF(mod_demcol,paste(mltfile,'_',
                                 resol,
                                 'mmResolution',
                                 '_DEM.tiff',sep=''))
    } else { # include the date in the file name
      writeTIFF(mod_demcol,paste(gsub("-", "", as.character(Sys.Date())),
                                 '_',mltfile,'_',resol,'mmResolution',
                                 '_DEM.tiff',sep=''))
    }
  } else { # store the tiff in a different folder
    currentwd <- getwd()
    setwd(oname)
    if (date == FALSE) {
      writeTIFF(mod_demcol,paste(mltfile,'_',resol,'mmResolution',
                                 '_DEM.tiff',sep=''))
    } else { # include the date in the file name
      writeTIFF(mod_demcol,paste(gsub("-", "", as.character(Sys.Date())),
                                 '_',mltfile,'_',resol,'mmResolution',
                                 '_DEM.tiff',sep=''))
    }
    setwd(currentwd) # resets to the original working directory
  }
  
  
  
  
  #Macroporosity tiff corrected using DEM corrected DEM ####
  # matrix to store the answer to the question: Is there a point in this grid cell? A:
  # Yes (1) or No (0); Note: columns = x, rows = y
  mltmat <- matrix(0,nrow=length(ybound),ncol=length(xbound)) 
  
  
  # Assign zero values in mltmat for the appropriate missing indices in 
  # rowcol:
  mltmat <- mltmat[1:(nrow(mltmat)-1),
                   1:(ncol(mltmat)-1)] # remove the last row and col 
  # in mltmat_corrected; Note: these are removed because they represent the right 
  # or bottom extent of the last full cells
  
  # uses the single index notation for matrices
  
  mltmat[rowcol$row + nrow(mltmat) * (rowcol$col - 1)] <- 1 
  
  
  
  
  ##Change the 1 to 0s where the DEM corrected shows low values.  
  
  
  
  rightMost <- dim(mod_demmat)[2]
  bottomMost <- dim(mod_demmat)[1]
  Top_left <- max(mod_demmat[1:300,1:300], na.rm = TRUE)
  Top_right <- max(mod_demmat[1:300,(rightMost-300):rightMost], na.rm = TRUE)
  Bottom_left <- max(mod_demmat[(bottomMost-300):bottomMost,1:300], na.rm = TRUE)
  Bottom_right <- max(mod_demmat[(bottomMost-300):bottomMost,(rightMost-300):rightMost], na.rm = TRUE)
  
  Average_depth_rulers <- mean(c(Top_left,Top_right, Bottom_left, Bottom_right))
  
  ##Create a threshold in which we want to exclude any pixels coming from depths 20 mm from te bottom of the tray. the "*10" converts cm to mm, the *20 makes it the z units from th bottom to 20 mm 
  DEM_threshold <-depth_from_tray_bottom_mm*Average_depth_rulers/(tray_depth_cm*10)
  
  Correction_sites <- which(mod_demmat<DEM_threshold)
  
  
  mltmat <- cmp_rotate(t(mltmat), direction = 'counter')
  mltmat_corrected <- mltmat
  mltmat_corrected[Correction_sites] <- 0 
  
  
  #mltmat_corrected <- cmp_rotate(t(mltmat_corrected), direction = 'counter')
  
  
  par(mfrow=c(1,3))
  image(as.matrix(mltmat), useRaster=TRUE, axes=FALSE, col = hcl.colors(500, "grays"))
  image(as.matrix(mltmat_corrected), useRaster=TRUE, axes=FALSE, col = hcl.colors(500, "grays"))
  image(as.matrix(mod_demmat), useRaster=TRUE, axes=FALSE, col = hcl.colors(500, "spectral"))
  
  #
  # Fix the orientation of the matrix for the tiff below:
  
  #
  ### Save the resulting macropore as tiff files: ####
  if (length(oname)==0) {
    if (date == FALSE) {
      writeTIFF(mltmat_corrected,paste(mltfile,'_',
                                       resol,
                                       'mmResolution',
                                       '_Macropores.tiff',sep=''))
    } else { # include the date in the file name
      writeTIFF(mltmat_corrected,paste(gsub("-", "", as.character(Sys.Date())),
                                       '_',mltfile,'_',resol,'mmResolution',
                                       '_Macropores.tiff',sep=''))
    }
  } else { # store the tiff in a different folder
    currentwd <- getwd()
    setwd(oname)
    if (date == FALSE) {
      writeTIFF(mltmat_corrected,paste(mltfile,'_',resol,'mmResolution',
                                       '_Macropores.tiff',sep=''))
    } else { # include the date in the file name
      writeTIFF(mltmat_corrected,paste(gsub("-", "", as.character(Sys.Date())),
                                       '_',mltfile,'_',resol,'mmResolution',
                                       '_Macropores.tiff',sep=''))
    }
    setwd(currentwd) # resets to the original working directory
  }
  
  #
}




# Function for rotating MLT point data (xyz) into a matrix and outputting a tiff image ####
Rotating_ascFile <- function(ascname, resol=0.18, dname = NULL, oname = NULL, date = FALSE, angle) {
  
  # ascname: character; file name of the point cloud data in the .asc format
  # outputted from RapidWorks; this file contains columns of data, the 
  # first three of which, refer to the x, y, z coordinates of the points;
  # additional columns refer to the RGB values associated with these 
  # points and the normals associated with these points
  #
  # resol: numeric; value that represents the resolution of the resulting
  # matrix grid in mm
  #
  # dname: character; folder where the asc file is stored if not in the
  # the default folder
  #
  # oname: character; folder where the tiff file is to be outputted if 
  # not in the the default folder
  #
  # date: logical; TRUE means output the system date in the tiff file
  # name (default); FALSE means not to output the date in the file name
  # 
  # NOTE: There are two large loops in this function whose progress are
  # displayed as two progress bars
  #
  #
  # Require the tiff package:
  require(tiff)
  #
  # Save the file name without the extension for later use:
  mltfile <- sub('.asc', '', ascname)
  #
  ### Read in the asc file: ####
  if (length(dname)==0) {
    mlt <- read.table(ascname,header=FALSE,sep="\t")
    if( length(names(mlt)) == 9){
      names(mlt) <- c("x","y","z","R","G","B","a","b","c")}else{
        names(mlt) <- c("x","y","z","R","G","B")
      }
  } else { # file is stored in another folder
    currentwd <- getwd()
    setwd(dname)
    mlt <- read.table(ascname,header=FALSE,sep="\t")
    if( length(names(mlt)) == 9){
      names(mlt) <- c("x","y","z","R","G","B","a","b","c")}else{
        names(mlt) <- c("x","y","z","R","G","B")
      }
    setwd(currentwd) # resets to the original working directory
  }
  
  
  
  ### Define the origin based on the minimum data in each column: ####
  nmlt <- mlt # save the "normalized" mlt values into a new dataframe
  nmlt$x <- nmlt$x-min(mlt$x)
  nmlt$y <- nmlt$y-min(mlt$y)
  nmlt$z <- nmlt$z-min(mlt$z)
  
  
  ###Change to Polar Coordinates ####
  
  nmlt$r <- sqrt((nmlt$y)^2+(nmlt$x)^2)
  
  nmlt$theta <- atan((nmlt$y)/(nmlt$x))
  
  
  ### Perform polar rotation ####
  angle_rad <- angle*(pi/180)
  nmlt$theta_rotated <- (nmlt$theta-angle_rad)
  
  ### Change back to cartesian ####
  nmlt$rotated_y <- nmlt$r*sin(nmlt$theta_rotated)
  nmlt$rotated_x <- nmlt$r*cos(nmlt$theta_rotated)
  
  
  nmlt$rotated_x <- nmlt$rotated_x-min(nmlt$rotated_x)
  nmlt$rotated_y <- nmlt$rotated_y-min(nmlt$rotated_y)
  
  
  
  
  
  #
  ### Calculate the spatial extent of the data (in mm): ####
  xdim <- max(nmlt$rotated_x)+resol # width
  ydim <- max(nmlt$rotated_y)+resol # length
  #
  # Setup a grid based on the extent data:
  xbound <- seq(0,xdim,by=resol) # provides the x coordinate for the left 
  # grid cell boundary
  ybound <- seq(0,ydim,by=resol) # provides the y coordinate for the upper 
  # grid cell boundary
  
  # Assign the column values of mltmat to a new dataframe, xmlt:
  xmlt <- data.frame(nmlt,row=NA,col=NA) # switch of variables
  xmlt <- xmlt[order(xmlt$rotated_x),] # sort the dataframe on the x values
  row.names(xmlt) <- NULL
  x <- round(xmlt$rotated_x,round(-log10(resol/10))) # rounds the x data to a tenth
  # of the specified resolution; if the resolution is greater than 10 this
  # will cap the rounding to a single mm
  
  dummy_mat <- matrix(0,nrow=length(ybound),ncol=length(xbound)) # matrix to
  # store the answer to the question: Is there a point in this grid cell? A:
  # Yes (1) or No (0); Note: columns = x, rows = y
  
  pb <- txtProgressBar(min = 0, max = (ncol(dummy_mat)-1), style = 3) # create
  # progress bar
  for (i in 1:(ncol(dummy_mat)-1)) {
    xmlt$col[which((x > xbound[i]) & (x <= xbound[i+1]))] <- i 
    Sys.sleep(0.1)
    setTxtProgressBar(pb, i) # update progress bar
  }
  close(pb)
  #
  # Assign the row values of mltmat to xmlt:
  xmlt <- xmlt[order(xmlt$rotated_y),] # sort the dataframe on the y values
  row.names(xmlt) <- NULL
  y <- round(xmlt$rotated_y,round(-log10(resol/10)))
  pb <- txtProgressBar(min = 0, max = (nrow(dummy_mat)-1), style = 3) # create
  # progress bar
  for (i in 1:(nrow(dummy_mat)-1)) {
    xmlt$row[which((y > ybound[i]) & (y <= ybound[i+1]))] <- i
    Sys.sleep(0.1)
    setTxtProgressBar(pb, i) # update progress bar
  }
  close(pb)
  #
  # Create a new rowcol dataframe and order by rows then cols:
  rowcol <- xmlt[,c("row","col")]
  rowcol <- rowcol[order(rowcol$row, rowcol$col),]
  row.names(rowcol) <- NULL
  
  
  rowcol <- unique(rowcol) #rowcol[!duplicated(rowcol),] # removes any duplicated cells
  row.names(rowcol) <- NULL
  
  #Aggregating data ####
  
  # Calculates the mean height (z) and R, G, B for each cell:
  xmlt <- xmlt[order(xmlt$row, xmlt$col),]
  row.names(xmlt) <- NULL
  xmlt <- data.frame(xmlt, # used to aggregate the data below
                     matid = as.factor(paste(xmlt$row,xmlt$col,sep="_")))
  xmlta <- aggregate(cbind(z, R, G, B) ~ matid, xmlt, FUN = "mean")
  xmlta$R <- round(xmlta$R,0) # rounds to the nearest integer
  xmlta$G <- round(xmlta$G,0) # rounds to the nearest integer
  xmlta$B <- round(xmlta$B,0) # rounds to the nearest integer
  xmltm <- merge(xmlt[,c("row","col","matid")],xmlta,by="matid") # holds
  
  
  
  ### the row and col indices and the mean values for z, R, G, and B ####
  xmltm <- xmltm[,c("row","col","z","R","G","B")] # removes the matid col
  
  
  xmltm <- xmltm[-which(is.na(xmltm$col) | is.na(xmltm$row)),] # removes the rows or columns with NA values
  
  
  ### removes any duplicated cells####
  xmltm <- xmltm[!duplicated(xmltm),] 
  
  
  
  xmltm <- xmltm[order(xmltm$row, xmltm$col),] # order the data.frame
  row.names(xmltm) <- NULL
  
  
  if (length(oname)==0) {
    if (date == FALSE) {
      save(xmltm, ybound, xbound , rowcol,mltfile,
           file= paste(mltfile,'_',
                       resol,
                       'mmResolution',
                       '_rotated',angle,'_ASCData.Rdata',sep=''))
    } else { # include the date in the file name
      save(xmltm, ybound, xbound , rowcol,mltfile,
           file= paste(gsub("-", "", as.character(Sys.Date())),
                       '_',mltfile,'_',resol,'mmResolution',
                       '_rotated',angle,'_ASCData.Rdata',sep=''))
    }
  } else { # store the tiff in a different folder
    currentwd <- getwd()
    setwd(oname)
    if (date == FALSE) {
      save(xmltm, ybound, xbound , rowcol,mltfile,
           file=paste(mltfile,'_',resol,'mmResolution',
                      '_rotated',angle,'_ASCData.Rdata',sep=''))
    } else { # include the date in the file name
      save(xmltm, ybound, xbound , rowcol,mltfile,
           file=paste(gsub("-", "", as.character(Sys.Date())),
                      '_',mltfile,'_',resol,'mmResolution',
                      '_rotated',angle,'_ASCData.Rdata',sep=''))
    }
    setwd(currentwd) # resets to the original working directory
  }
  
  return(xmltm) 
  
}


