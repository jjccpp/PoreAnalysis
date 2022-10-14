
#-----------------------------------------#

# Processing soil asc files from MLT to obtain Consistent Pore data ####   
# Author: Julio Pachon

#-----------------------------------------#

##Motivation ####
# Obtain porosity information from a from the Multistripe Laser Triangulation. This code is meant to increase recording of image processing once it has been converted from .asc files to .tiff files.

# An advantage of this code for rotation vs using the rotate function in ImageJ is that we have three different tiff files and here we can ensure the same rotation if they are ever needed to be worked as a raster brick, and ensure that it lines up.  

# The advantage of cropping here is having control and understanding of how much was cropped from each side and keeping consistence across all three tiff files.

##Sample used to create code: KONZA_APT_S128148_T5

##Using EBImage package: 
#https://bioconductor.org/packages/devel/bioc/vignettes/EBImage/inst/doc/EBImage-introduction.html#9_Image_segmentation
#install.packages("BiocManager")
#BiocManager::install("EBImage", force = TRUE)



library(compiler)

library(tiff)
library(raster)
library("EBImage")





source("./Code/PoreAnalysis_Functions.R")

#Compile the functions to make them run faster
Aggregating_asc <- cmpfun(Aggregating_asc)
asctotiffs <- cmpfun(asctotiffs)
Rotating_ascFile <- cmpfun(Rotating_ascFile)




#Aggregating_asc  ####


param_files <- list.files("./ASC Files/", pattern = "*.asc")


##Running files in line
for  (i in 1:length(param_files)){  
  
  ascname <- paste0("./ASC Files/",param_files[i])
  

  Aggregating_asc(ascname)
  }


##parallel Run ####
library(parallel)
ncores <- detectCores(logical = FALSE)/2
n <- ncores:1

# Create a cluster
cl <- makeCluster(ncores)

# Use clusterApply to call rnorm for each n in parallel,
# again setting mean to 10 and sd to 2 

start <- Sys.time()
start
clusterApply(cl, ascname , Aggregating_asc)

# Stop the cluster
stopCluster(cl)
(end <- Sys.time())
end-start




#asctotiffs  ####
rdata_files <- list.files("./Rdata Files/Unrotated/", pattern = "*.Rdata")


for  (i in 1:length(rdata_files)){  
  
  rdata_file <- paste0("./Rdata Files/Unrotated/",rdata_files[i])
  
  asctotiffs(rdata_file, depth_from_tray_bottom_mm = 7, tray_depth_cm = 4)
}














 #Rotating_asc  ####
# EBIimage allows for the manipulation of the tiff in R to zoom in (scrolling with the mouse or the top right buttons) and view the location of the top and left boundaries of the trays by hovering with the mouse over their location and reading the coordiantes in the lower portion.


##Angle determination####

RGB_unrotated_RGBfile <- list.files("./tiff files/Unrotated/RGB/", pattern = "*.tiff")

RGB_unrotated_RGB <- readImage(paste0("./tiff files/Unrotated/RGB/",RGB_unrotated_RGBfile[1]))


RGB_unrotated_DEMfile <- list.files("./tiff files/Unrotated/DEM/", pattern = "*.tiff")

RGB_unrotated_DEM <- readImage(paste0("./tiff files/Unrotated/DEM/",RGB_unrotated_DEMfile[1]))


display(RGB_unrotated_RGB)
display(RGB_unrotated_DEM)



#find the angle by looking at two points that should be straight and doing basic math:
(angle <- atan((144-121)/(1233-168))*180/pi)



param_files <- list.files("./ASC Files/", pattern = "*.asc")


ascname <- paste0("./ASC Files/",param_files[1])

Rotating_ascFile(ascname,  angle=angle)




#asctotiffs  ####
rdata_files <- list.files("./Rdata Files/Rotated/", pattern = "*.Rdata")


  rdata_file <- paste0("./Rdata Files/Rotated/",rdata_files[j])
  
  asctotiffs(rdata_file, depth_from_tray_bottom_mm = 15, tray_depth_cm = 4)

  

  
  
  #Processing the horizons ####
  
  RBGBfile <- list.files("./tiff files/Rotated/RBG/", pattern = "*.tiff")
  RBG <- readImage(paste0("./tiff files/Rotated/RBG/",RRGBfile[1]))
  
  
  DEMfile <- list.files("./tiff files/Rotated/DEM/", pattern = "*.tiff")
  DEM <- readImage(paste0("./tiff files/Rotated/DEM/",DEMfile[1]))
  
  Macroporesfile <- list.files("./tiff files/Rotated/Macropore/", pattern = "*.tiff")
  Macropores <- readImage(paste0("./tiff files/Rotated/Macropore/",Macroporesfile[1]))
  
  
  display(RBG)
  display(DEM)
  display(Macropores)
  
  
  ### Correcting rotation of the image
  
  # If the image is not done correctly (check that the ruler shows T at the top, L on the left, B on the bottom, R on the right), rotate to correct that.
  
  
  Orientation_change= -90
  
  RBG = EBImage::rotate(RBG, Orientation_change)
  
  DEM = EBImage::rotate(DEM, Orientation_change)
  
  Macropores = EBImage::rotate(Macropores, Orientation_change)
  
  display(RBG)
  display(DEM)
  display(Macropores)
  
  
  #Basic Info- check by using the RBG file
  if(dim(RBG)[2]>2200){
    Length_tray_mm = 400
  } else {
    Length_tray_mm = 200
  }
  
  Length_tray_mm 
  Width_tray_mm <-300
  Resolution_mm <- 0.18 #resolution of the camera, which is also the size used in the aggregation when transforming from .asc to .tiff files.
  Top_tray <-      135  #row
  Left_tray <-     129  #col
  Bot_crop_mm <-   10 
  Top_crop_mm <-   10
  
  
  
  length_tray_pixel <- Length_tray_mm/Resolution_mm #300 mm 
  Width_tray_pixel <- Width_tray_mm/Resolution_mm
  
  
  
  
  ##Cropping the outer edges ####
  #pick column where it is easy to tell where the lower limit of the metal tray is
  
  #Use the browser and plot above to choose the uppercut
  
  
  
  
  #so the bottom of the tray should be 1666.667+140
  
  RBG_crop <- RBG[,(Top_tray):(Top_tray+length_tray_pixel),]
  display(RBG_crop)
  
  
  
  #However both the top and bottom have gaps so take away Buffer cm from each which is 56 pixels
  
  
  Top_crop_pixels <- Top_crop_mm/Resolution_mm
  Bot_crop_pixels <- Bot_crop_mm/Resolution_mm
  
  
  RBG_crop <- RBG[,(Top_tray+Top_crop_pixels):(Top_tray+length_tray_pixel-Bot_crop_pixels),]
  display(RBG_crop)
  
  
  #If it looks good, crop DEM and Macropore too
  
  Macropore_crop <- Macropores[,
                               (Top_tray+Top_crop_pixels):(Top_tray+length_tray_pixel-Bot_crop_pixels)]
  
  DEM_crop<- DEM[,
                 (Top_tray+Top_crop_pixels):(Top_tray+length_tray_pixel-Bot_crop_pixels),]
  
  display(Macropore_crop)
  display(DEM_crop)
  
 
  
##Look at file name for information on horizon
  RBGBfile
  
  tray_top_depth <-    74 
  (Tray_bottom_depth <- tray_top_depth+(Length_tray_mm/10))
  
  Horizon_depth_top <-    68 
  Horizon_depth_bottom <- 103
  
  Bot_crop_mm
  Top_crop_mm
  
  
  if(Horizon_depth_top>tray_top_depth+ (Top_crop_pixels*Resolution_mm/10)){
    Top_cut_pixels <- Top_tray+((Horizon_depth_top-tray_top_depth)*10/Resolution_mm)
  } else {
    Top_cut_pixels <- Top_tray+Top_crop_pixels 
  }
  
  
  
  if(Tray_bottom_depth>Horizon_depth_bottom- (Bot_crop_pixels*Resolution_mm/10)){
    Bottom_cut_pixels <- Top_tray+length_tray_pixel-((Tray_bottom_depth-Horizon_depth_bottom)*10/Resolution_mm)
    
  } else {
    Bottom_cut_pixels <- Top_tray+length_tray_pixel-Bot_crop_pixels
  }
  
  
  Horizon_RBG_test <- RBG[,Top_cut_pixels:Bottom_cut_pixels,]
  
  Horizon_Macropore_test <- Macropores[,Top_cut_pixels:Bottom_cut_pixels]
  
  Horizon_DEM_test <- DEM[,Top_cut_pixels:Bottom_cut_pixels,]
  
  
  display(Horizon_RBG_test)
  display(Horizon_Macropore_test)
  display(Horizon_DEM_test)
  
  
  
  #Select the left and right boundaries for the section selected:
  Left_crop_mm <-  20 
  Right_crop_mm <- 10 
  
  #find left tray inner edge
  
  Left_crop_pixels <- Left_crop_mm/Resolution_mm
  Righ_crop_pixels <- Right_crop_mm/Resolution_mm
  
  
  #Check and create the images that will become tiffs
  
  Horizon_RBG <- RBG[(Left_tray+Left_crop_pixels):(Left_tray+Width_tray_pixel-Righ_crop_pixels),
                     Top_cut_pixels:Bottom_cut_pixels,]
  display(Horizon_RBG)
  
  Horizon_Macropores <- Macropores[(Left_tray+Left_crop_pixels):(Left_tray+Width_tray_pixel-Righ_crop_pixels),
                                   Top_cut_pixels:Bottom_cut_pixels]
  display(Horizon_Macropores)
  
  Horizon_DEM <- DEM[(Left_tray+Left_crop_pixels):(Left_tray+Width_tray_pixel-Righ_crop_pixels),
                     Top_cut_pixels:Bottom_cut_pixels,]
  display(Horizon_DEM)
  
  
  ##Image for metadata
  
  top_ruler <- 1:(Top_tray-25)
  bottom_ruler <- (Top_tray+length_tray_pixel+25):dim(RBG)[2]
  
  (Left_tray+Left_crop_pixels):(Left_tray+Width_tray_pixel-Righ_crop_pixels)
  
  Horizon_RBG_test <- RBG[c(1:(Left_tray-25),
                            rep(1,20), #creates some space between ruler and the horizon
                            (Left_tray+Left_crop_pixels):(Left_tray+Width_tray_pixel-Righ_crop_pixels),
                            rep(1,20),
                            (Left_tray+Width_tray_pixel+25):dim(RBG)[1]),
                          c(top_ruler,
                            rep(1,20),  
                            Top_cut_pixels:Bottom_cut_pixels,
                            rep(1,20),
                            bottom_ruler),]
  
  display(Horizon_RBG_test)
  
  
  
  
  ##Create Metadata for the horizons ####
  
  RBGBfile
  
  (tray_top_depth+(Top_cut_pixels-Top_tray)*Resolution_mm/10)
  (Tray_bottom_depth-(Top_tray+length_tray_pixel-Bottom_cut_pixels)*Resolution_mm/10)
  
  # horizon_name <- paste0("mlt_", "KONZA_PPT_Bt_0.18mmResolution", "_T1_22-38cm")
  
  
  (horizon_name <- paste0("mlt_KONZA_PPT_Btss1_0.18mmResolution_T3_",
                         (tray_top_depth+(Top_cut_pixels-Top_tray)*Resolution_mm/10),
                         "-",
                         (Tray_bottom_depth-(Top_tray+length_tray_pixel-
                                               Bottom_cut_pixels)*Resolution_mm/10),
                         "cm"))
  
  angle
  
  Horizon_metadata <- data.frame(
    file=horizon_name, 
    hor="KONZA_PPT_T4_91-131cm_0.18mmResolution",
    tray=RBGBfile,
    Width_tray_cm=Width_tray_mm/10,
    Length_tray_cm=Length_tray_mm/10,
    rotation_angle_degrees=angle, 
    bottom_crop_fromBottomTray_cm=(Top_tray+length_tray_pixel-Bottom_cut_pixels)*Resolution_mm/10, 
    top_crop_toptray_cm=(Top_cut_pixels-Top_tray)*Resolution_mm/10,
    left_crop_fromLeftTray_cm=Left_crop_mm/10, 
    Right_crop_fromRightTray_cm=Right_crop_mm/10, 
    left_tray_col=Left_tray, 
    Top_tray_row=Top_tray)
  
  write.table(Horizon_metadata , "clipboard-16384", sep="\t", row.names=FALSE)
  
  
  # #Save Horizon cropped tiffs ####
  
  
  writeImage(Horizon_RBG_test,paste0("./Cropped Horizons/",horizon_name,'_Rulers_CroppedHorizon.tiff'))
  
  writeImage(Horizon_RBG,paste0("./Cropped Horizons/RBG/",horizon_name,'_RBG_CroppedHorizon.tiff'))
  
  writeImage(Horizon_DEM,paste0("./Cropped Horizons/DEM/",horizon_name,'_DEM_CroppedHorizon.tiff'))
  
  writeImage(Horizon_Macropores,paste0("./Cropped Horizons/Macropore/",horizon_name,'_Macropores_CroppedHorizon.tiff'))
  
  
  
  
  
  
  
  