// @File(label = "Input folder:", style = "directory") inputDir
// @File(label = "Output folder:", style = "directory") outputDir
// @String (label = "File suffix", value = ".nd2") fileSuffix
// @int(label= "Channel to analyze", style = "spinner", val = 1) channelNum


// analyze_protrusions_confocal.ijm
// by Theresa Swayne for Xu Zhang, 2024
// identifies processes in 2d fluorescence images and determines length and other parameters
// Thanks to Ignacio Argando Carreras for the beanshell script prunebysize_.bsh 
// Note that prunebysize_.bsh must be in the fiji plugins/scripts folder!

// INPUT: a folder of 2D fluorescence images (can be multichannel)
// OUTPUT: skeleton information (length in PIXELS), cell counts, and overlay of detected processes and cells.

// TO USE: Run the macro and specify folders for input and output.
// If images are multichannel, select the channel to use for analysis.

// TIPS: This code is designed for confocal images with pixel size ~0.3 µm. 
//		 It works best if cell bodies are NOT saturated.


// ---- Setup ----

while (nImages>0) { // clean up open images
	selectImage(nImages);
	close();
}

print("\\Clear"); // clear Log window
roiManager("reset"); // clear ROI Mgr
run("Clear Results"); // clear results window
setOption("ScaleConversions", true); // ensures good pixel value?
run("Set Measurements...", "area centroid shape feret's display redirect=None decimal=3"); // get shape measurements

setBatchMode(true); // faster performance
run("Bio-Formats Macro Extensions"); // support native microscope files

// set up a table for batch results
summaryTable = "SkeletonLengths";
Table.create(summaryTable);
Table.setLocationAndSize(200, 200, 250, 600);

// ---- Run ----
processFolder(inputDir, outputDir, fileSuffix, channelNum);
print("Finished");

// ---- Functions ----


function processFolder(input, output, suffix, channel) {
	fileNum = -1;
	// function to scan folder tree to find files with correct suffix
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i])) {
			processFolder(input + File.separator + list[i], output, suffix, channel);
		}
		if(endsWith(list[i], suffix)) {
			fileNum = fileNum + 1;
			processFile(input, output, list[i], channel, fileNum);
		}
	}
}

function processFile(input, output, file, channel, filenumber) {

	// function to process a single image
	
	print("Processing channel",  channel, "of" , input + File.separator + file);
	
	path = input + File.separator + file;

	// ---- Open image and get information ----
	run("Bio-Formats", "open=&path");
	id = getImageID();
	title = getTitle();
	dotIndex = indexOf(title, ".");
	basename = substring(title, 0, dotIndex);
	extension = substring(title, dotIndex);
	getDimensions(width, height, channels, slices, frames);
	print("Processing",title, "with basename",basename);
	
	// ---- Prepare images ----
	if (channels > 1) {
	
		run("Split Channels");
		img = "C"+channel+"-"+title;
		
	}
	
	else {
		img = title;
	}	
	
	selectImage(img);
	resetMinAndMax();
	run("8-bit"); // enables local threshold to work
	run("Duplicate...", "title=orig"); // make a copy for overlaying later
	
	// ---- Segmentation of cell bodies
	
	// this count is used to normalize the skeleton length later
	selectImage("orig");
	run("Duplicate...", "title=Cells");
	selectImage("Cells");
	setAutoThreshold("Default dark");
	run("Convert to Mask");
	run("Options...", "iterations=5 count=1 do=Open");
	run("Watershed");
	run("Analyze Particles...", "size=100-Infinity show=Nothing add clear summarize");
	selectWindow("Summary");
	IJ.renameResults("Summary","Results");
	cellNum = getResult("Count", 0);
	//roiManager("Save", output + basename + "_cellROIs.zip");

	// ---- Segmentation of processes ----
	
	// apply a local threshold to identify cell area including processes	
	selectImage(img);
	run("Auto Local Threshold", "method=Phansalkar radius=15 parameter_1=0 parameter_2=0 white");
	
	// apply a filter to enhance tube-like structures
	run("Frangi Vesselness", "input=[&img] dogauss=true spacingstring=[1, 1] scalestring=1");
	
	selectImage("result");
	//rename("vesselness");
	
	// apply a global threshold to preserve the most tube-like structures
	setAutoThreshold("Percentile dark");
	setOption("BlackBackground", false);
	run("Convert to Mask");
	
	// identify objects, removing small and more circular objects
	run("Analyze Particles...", "size=150-Infinity circularity=0.00-0.50 show=[Masks] clear");
	
	// generate a skeleton from the remaining objects
	selectImage("Mask of result");
	run("Fill Holes"); // make cell bodies into solid objects that don't contribute much to the skeleton
	run("Skeletonize (2D/3D)");
	
	// remove short (usually artifactual) segments from the skeleton -- more accurate than pruning end segments in Analyze
	selectImage("Mask of result");
	run("Analyze Particles...", "size=10-Infinity circularity=0.00-1 show=Masks");
	selectImage("Mask of Mask of result");
	rename("Skeleton");
	
	// script by igancio arganda, image.sc forum https://forum.image.sc/t/analyzeskeleton-gui-prune-by-length/3657/18?u=iarganda
	// the beanshell script prunebysize_.bsh must be in the fiji plugins/scripts folder!
	// threshold is length in pixels -- smaller segments will be eliminated
	run("prunebysize ", "image=Skeleton threshold=30.0");
	
	
	// measure the segments of the skeleton
	selectImage("Skeleton-pruned");
	//rename("Skeleton");
	run("Analyze Skeleton (2D/3D)", "prune=none show");

		
	// Overlay skeleton and orig image
	//run("Merge Channels...", "c1=Skeleton c2=orig create keep");
	run("Merge Channels...", "c1=Skeleton-pruned c2=orig create keep");
	selectImage("Composite");
	rename(basename + "_overlay");
	// Property.set("CompositeProjection", "null");
	Stack.setDisplayMode("color");
	Stack.setChannel(1);
	run("Red");
	Property.set("CompositeProjection", "Sum");
	Stack.setDisplayMode("composite");
	Property.set("CompositeProjection", "null");
	// add the overlay of cells
	run("From ROI Manager");

	

	// ---- Calculate total branch length per image ----

	selectWindow("Results");
	totalLength = 0;
	for (row=0; row<nResults; row++) {
		
		// total length for one skeleton
		skelLength = getResult("# Branches", row) * getResult("Average Branch Length", row);
		setResult("Total Length", row, skelLength);
		updateResults();
		
		// total for all skeletons in the image
		totalLength = totalLength + skelLength;
	    
	}

	// ---- Calculate median branch length per image ---

	selectWindow("Branch information");
	lengths = Table.getColumn("Branch length", "Branch information");
	medianLength = median(lengths);
	
	// ---- Update skeleton data ---- 
	
	Table.set("Image", filenumber, basename, summaryTable);
	Table.set("Total Length", filenumber, totalLength, summaryTable);
	Table.set("Number of cells", filenumber, cellNum, summaryTable);
	Table.set("Normalized length", filenumber, totalLength/cellNum, summaryTable);
	Table.set("Median branch length", filenumber, medianLength, summaryTable);
	Table.update(summaryTable);

	
	// ---- Save skeleton results ---- 
	
	print("Saving to " + output);

	//vesselName = basename+"_vesselness.tif";
	//selectImage(basename + "_vesselness");
	//saveAs("tiff", output + File.separator + vesselName);
	
	overlayName = basename+"_overlay.tif";
	selectImage(basename + "_overlay");
	saveAs("tiff", output + File.separator + overlayName);
	
	skelName = basename+"_skeleton.tif";
	selectImage("Skeleton-pruned");
	saveAs("tiff", output + File.separator + skelName);
	
	skelDataName = basename + "_skel_info.csv";
	selectWindow("Results");
	saveAs("Results", output + File.separator + skelDataName);
	
	branchDataName = basename + "_branch_info.csv";
	selectWindow("Branch information");
	saveAs("Results", output + File.separator + branchDataName);
	
	summaryName = "batch_summary.csv";
	Table.save(output + File.separator + summaryName, summaryTable);

	// clean up open images and tables
	while (nImages>0) {
	selectImage(nImages);
	close();
	}
	
	// selectWindow("Results");
	// close();
	selectWindow(branchDataName); // close branch results
	run("Close");
	run("Clear Results"); // clear skeleton results
} // processFile function


function median(x){
	// determine the median of an array
	x=Array.sort(x);
	if (x.length%2>0.5) {
		m=x[floor(x.length/2)];
	}else{
		m=(x[x.length/2]+x[x.length/2-1])/2;
	};
	return m
}


