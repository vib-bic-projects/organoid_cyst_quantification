// @ Float (label="Image scalling", min=0.0, max=1.0, value=0.20) scaling

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
rename("processed");
close("\\Others");


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
