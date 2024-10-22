// @ File(label="File directory", style="directory") dir
// @ File(label="Pixel classifier", style="file") pixel_class
// @ String (label="What working station are you using?", choices={"PC7_Paris", "PC_Nicolas", "NAS_Tokyo", "PC2_Mexico"}, style="listBox") pc
// @ Float (label="Image scalling", min=0.0, max=1.0, value=0.20) scaling
// @ String (label="Perform analysis", choices={"yes", "no"}, style="listBox") analysis_param
// @ Integer (label="Bin size for local thickness[µm]", min=0.1, max=10000, value=20) binesize_thick
// @ String (label="File suffix", choices={".nd2", ".tif"}, style="listBox") suffix

/* 7/06/2024
Nicolas Peredo
VIB BioImaging Core Leuven - Center for Brain and Disease Research
Nikon Center of Excellence
Campus Gasthuisberg - ON5 - room 04.367
Herestraat 49 - box 62
3000 Leuven
Belgium
phone +32 (0)16/37.70.03

When you publish data analyzed with this script please add the references of the used plug-ins:

Legland, D., Arganda-Carreras, I., & Andrey, P. (2016). MorphoLibJ: integrated library and plugins for mathematical morphology with ImageJ. Bioinformatics, 32(22), 3532–3534. doi:10.1093/bioinformatics/btw413
*/

//Get directories and lists of files
setOption("ExpandableArrays", true);

//Blood vessel directory
fileList = getFilesList(dir, suffix);
Array.sort(fileList);

//Create the different folders with results
File.makeDirectory(dir + "/Analysis");
File.makeDirectory(dir + "/Segmented");

//GPU parameters
computerarray = newArray("PC7_Paris", "PC_Nicolas", "NAS_Tokyo", "PC2_Mexico");
gpuparameters = newArray("[NVIDIA GeForce RTX 3090]", "[NVIDIA GeForce RTX 3070]", "[Quadro RTX 8000]", "[Quadro K2200]");
gpu = "";
for (i = 0; i < computerarray.length; i++) {
	if (computerarray[i] == pc) {
		gpu = gpuparameters[i];
}

//Start of the arrays for the different measurements
filename_column = newArray();
volume_organoid_column = newArray();
label_cysts_column = newArray();
volume_cysts_column = newArray();
cystnb_column = newArray(); 

for (files = 0; files < fileList.length; files++) {
//for (files = 5; files < 6; files++) {
	
	//File and ROI names
	file = fileList[files];
	name = getBasename(file, suffix);
	
	//Open image
	run("Bio-Formats Importer", "open=[" + dir + File.separator + file + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");

	rename("raw_image");
	run("Duplicate...", "title=actin duplicate channels=1-1");
	
	selectWindow("raw_image");
	run("Duplicate...", "title=nuclei duplicate channels=2-2");
	
	selectWindow("actin");
	run("Enhance Contrast...", "saturated=0 equalize process_all use");
	
	selectWindow("nuclei");
	run("Enhance Contrast...", "saturated=0 equalize process_all use");
	
	imageCalculator("Add create 32-bit stack", "actin","nuclei");
	setSlice(nSlices/2);
	resetMinAndMax();
	run("8-bit");
	dimensions = scaling2D(scaling, true);
	
	//This step is only necessary for dataset2
	run("Subtract Background...", "rolling=50 stack");
	rename("scaled");
	
	//Get parameters from image
	getDimensions(width, height, channels, slices, frames);
	getVoxelSize(width, height, depth, unit);
	
	run("Run Pixel Classification Prediction", "projectfilename=" + pixel_class + " inputimage=scaled pixelclassificationtype=Segmentation");
	
	//Reassign pixel size since the generated filtered image is not calibrated anymore
	Stack.setXUnit(unit);
	run("Properties...", "channels=1 slices=" + slices + " frames=1 pixel_width=" + width + " pixel_height=" + width + " voxel_depth=" + depth);
	
	//Classified image
	rename("segmented");
	
	//Isolate cysts
	selectWindow("segmented");
	setThreshold(255, 255);
	setOption("BlackBackground", true);
	run("Convert to Mask", "background=Dark black create");
	rename("cysts");
	
	//Get the organoid
	selectWindow("segmented");
	setThreshold(253, 253);
	setOption("BlackBackground", true);
	run("Convert to Mask", "background=Dark black create");
	rename("organoid");
	
	imageCalculator("Add create stack", "organoid","cysts");
	saveAs("Tiff", dir + "/Segmented/" + name + "_organoid");
	rename("organoid_filled");
	
	//Filter cysts out of the organoid
	run("Morphological Filters (3D)", "operation=Erosion element=Ball x-radius=10 y-radius=10 z-radius=10");
	rename("temp-erode");
	
	run("Keep Largest Region");
	rename("temp-erode-filter");
	run("Morphological Filters (3D)", "operation=Dilation element=Ball x-radius=6 y-radius=6 z-radius=6");
	
	rename("organoid_eroded");
	
	imageCalculator("AND create stack", "cysts","organoid_eroded");
	rename("cysts_inside");
	
	//Median filter
	run("CLIJ2 Macro Extensions", "cl_device=" + gpu);
	
	// median
	image1 = "cysts_inside";
	Ext.CLIJ2_push(image1);
	image2 = "cysts_inside_median";
	radius_x = 4.0;
	radius_y = 4.0;
	radius_z = 4.0;
	Ext.CLIJ2_median3DSphere(image1, image2, radius_x, radius_y, radius_z);
	Ext.CLIJ2_pull(image2);
	
	run("Size Opening 2D/3D", "min=100");
	
	run("3D Fill Holes");
	rename("cysts_binary");
	
	run("Connected Components Labeling", "connectivity=6 type=[16 bits]");
	//run("Kill Borders");
	run("Remap Labels");
	resetMinAndMax();
	rename("cysts_inside_median_labelled");
	
	//Reassign pixel size since the generated filtered image is not calibrated anymore
	Stack.setXUnit(unit);
	run("Properties...", "channels=1 slices=" + slices + " frames=1 pixel_width=" + width + " pixel_height=" + width + " voxel_depth=" + depth);
	
	//Saving cyst images
	saveAs("Tiff", dir + "/Segmented/" + name + "_cysts");
	rename("cysts_final");
	
	//Analysis
	if (analysis_param == "yes") {
		//Starting arrays
		filename_array = newArray();
		volume_organoid_array = newArray();
		cystnb_array = newArray();
		label_cysts_array = newArray();
		volume_cysts_array = newArray();
		
		//Calculate organoid volume
		selectWindow("organoid_filled");
		run("Analyze Regions 3D", "volume surface_area_method=[Crofton (13 dirs.)] euler_connectivity=6");
		volume_organoid = getResult("Volume", 0);

		//Calculate the cyst volume and store it in an array		
		selectWindow("cysts_final");
		run("Analyze Regions 3D", "volume surface_area_method=[Crofton (13 dirs.)] euler_connectivity=6");
		Table.rename("cysts_final-morpho", "Results");
		
		if (nResults > 0) {
			selectWindow("Results");
			label_cysts_array = Table.getColumn("Label");
			volume_cysts_array = Table.getColumn("Volume");
			for (label = 0; label < label_cysts_array.length; label++) {
				filename_array[label] = name;
				volume_organoid_array[label] = volume_organoid;
				cystnb_array[label] = label_cysts_array.length;
			}
		}
		else {
			
			filename_array[0] = name;
			volume_organoid_array[0] = volume_organoid;
			cystnb_array[0] = 0;
			label_cysts_array[0] = 0;
			volume_cysts_array[0] = 0;
			
		}
		
		//Concatenate the arrays into the full column content
		filename_column = Array.concat(filename_column,filename_array);
		volume_organoid_column = Array.concat(volume_organoid_column,volume_organoid_array);
		label_cysts_column = Array.concat(label_cysts_column,label_cysts_array);
		volume_cysts_column = Array.concat(volume_cysts_column,volume_cysts_array);
		cystnb_column = Array.concat(cystnb_column,cystnb_array);
		close("Results");
		
		//Local thickness analysis
		selectWindow("cysts_binary");
		run("Local Thickness (masked, calibrated, silent)");
		rename("cysts_thickness");
		getstackhisto("cysts_thickness", binesize_thick);
		Table.rename("Results", "Localthickness");
		
		saveAs("Results", dir + File.separator + "Analysis" + File.separator + name + "_LocalThickness.csv");
		close(name + "_LocalThickness.csv");
	}
	
	//Close non important windows
	close("*");
	run("Collect Garbage");
	//close("\\Others");
}

Table.create("Pooled_Results");
Table.setColumn("Filename", filename_column);
Table.setColumn("Cyst_Nb", cystnb_column);
Table.setColumn("Organoid volume (um^3)", volume_organoid_column);
Table.setColumn("Cyst_ID", label_cysts_column);
Table.setColumn("Cyst volume (um^3)", volume_cysts_column);

saveAs("Results", dir + File.separator + "Analysis" + File.separator + "Pooled_Results.csv");


//Extract a string from another string at the given input smaller string (eg ".")
function getBasename(filename, SubString){
  dotIndex = indexOf(filename, SubString);
  basename = substring(filename, 0, dotIndex);
  return basename;
}

//Return a file list contain in the directory dir filtered by extension.
function getFilesList(dir, fileExtension) {  
  tmplist=getFileList(dir);
  list = newArray(0);
  imageNr=0;
  for (i=0; i<tmplist.length; i++)
  {
    if (endsWith(tmplist[i], fileExtension)==true)
    {
      list[imageNr]=tmplist[i];
      imageNr=imageNr+1;
      //print(tmplist[i]);
    }
  }
  Array.sort(list);
  return list;
}


// this scaling only occurs for XY and leaves the Z calibration the same
function scaling2D(ScaleNumber, interp){
	//Preprocessing
	getDimensions(width, height, channels, slices, frames);
	originaldimensions = newArray(width, height, channels, slices, frames);
	NumberSlices = nSlices/channels;
	Scaled_X = width*ScaleNumber;
	Scaled_Y = height*ScaleNumber;
	if (interp) {
		run("Scale...", "x=" + ScaleNumber + " y=" + ScaleNumber + " z=1.0 width=" + Scaled_X + " height=" + Scaled_Y + " depth=" + NumberSlices + " interpolation=Bilinear average process create");
	}
	else {
		run("Scale...", "x=" + ScaleNumber + " y=" + ScaleNumber + " z=1.0 width=" + Scaled_X + " height=" + Scaled_Y + " depth=" + NumberSlices + " interpolation=Bilinear average process create");	
	}
	return originaldimensions;
}

//Get the histogram of a stack with the counts in fractions of 1 additionally
function getstackhisto(image, binsize) {
	setOption("ExpandableArrays", true);
	//Select image of interest
	selectWindow(image);
	
	//Calculate the different values to feed into the histogram function
	getMinAndMax(min, max);
	max = (binsize * Math.ceil(max / binsize)) + binsize; // Round up max to the nearest multiple of binsize + binsize to round up the range
	
	bins = max / binsize;
	
	// Initialize the resulting array and the counts_bis array
	result_fraction = newArray(bins);
	result_count = newArray(bins);
	
	// Process each slice in the stack
	for (slice = 1; slice <= nSlices; slice++) {
		setSlice(slice);

		// Calculate the histogram for the current slice
		getHistogram(values, counts, bins, 0, max);

		// Update the resulting arrays with counts and fractions
		for (i = 0; i < values.length; i++) {
			result_count[i] += counts[i];
		}
	}
	                
	//Get the fraction of the volume numbers
	//Get the total amount of pixel counts
	totalpx = 0;
	for (i = 0; i < values.length; i++) {
		totalpx += result_count[i];
	}
	
	//Get the fraction array
	for (i = 0; i < values.length; i++) {
		result_fraction[i] = result_count[i]/totalpx;
	}
	
	//Generate a results table containing the stack histogram
	Table.create("Results");
	Table.setColumn("Values", values);
	Table.setColumn("Counts", result_count);
	Table.setColumn("Counts_fraction", result_fraction);
	close("Montage");
}