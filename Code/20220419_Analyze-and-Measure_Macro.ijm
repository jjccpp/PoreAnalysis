input=getDirectory("Choose Source Directory");
list = getFileList(input);
dir=getDirectory("Choose Results Directory");

for (i = 0; i < list.length; i++)
        Particle_analyze(input, list[i]);


function Particle_analyze(input,filename){
	open (input + filename);
	
// Assign raw opened image the name "currentImage"
currentImage=getImageID();


run("Set Scale...", "distance=5.5556 known=1 pixel=1 unit=mm global");
run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction stack limit redirect=None decimal=3");
//since we output a binary file from the R script, white is 255 and black is 0; the threshold below should recognize all the pores
setThreshold(0, 250);



run("Analyze Particles...", "display clear overlay");

originalName = getTitle();
originalNameWithoutExt = replace( originalName , ".tiff" , "" ); 
resultName = originalNameWithoutExt + "_Analyze.csv";
   saveAs("Results", dir+resultName); 
run("Clear Results");


run("Measure");

   //Change the file name here
resultName = originalNameWithoutExt + "_Measure.csv";
   saveAs("Results", dir+resultName); 
run("Clear Results");


run("*");

}

