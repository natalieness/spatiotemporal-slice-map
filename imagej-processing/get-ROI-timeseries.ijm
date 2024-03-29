/*This script measures intensity over time in a grid overlaid 
 on the sample based on an ROI
 Inputs:
 	1. A stack image, can be multichannel. 
 	2. An associated array of ROIs of the area to be mapped. 

 Parameters to set: 
 	1. loc - file path location to save outputs 
 	2. array_of_nam is a list of all time periods to be considered, e.g. BSL, EXP, etc. Needs to be the name of each respective file.
*/


//get directory to save results to 
loc = getDirectory("Choose directory to save your data to");
print(loc);

//Set parameters using dialog

openwins=getList("image.titles");
L = lengthOf(openwins);
print(L);
//dialog 
Dialog.create("Set images and ROIs");


for (i = 0; i < L; i++) {
	Dialog.addChoice("Next image: ", openwins);
	Dialog.addNumber("ROI number, starting at 0: ", 0);
	
}
Dialog.show();
//get dialog input
order = newArray();
rois = newArray();
for (i = 0; i < L; i++) {
	o = Dialog.getChoice();
	r = Dialog.getNumber();
	order = Array.concat(order, newArray(o));
	rois = Array.concat(rois, r);	
}

array_of_nam = order;
array_of_roi = rois;
//print for reference
Array.print(array_of_nam);
Array.print(array_of_roi);
// setting of parameters ends here 

//get number of channels
Stack.getDimensions(width, height, c, slices, frames);


//if (endsWith(loc, "/") == 0){
//	exit("File path is not formatted correctly");


c = c+1;
for (i=0; i<array_of_nam.length;i+=1) {
	nam = array_of_nam[i];
	roi_no = array_of_roi[i];
	print(nam);
	selectWindow(nam);
	roiManager("deselect");
	run("Duplicate...", "title="+nam+"_duplicate duplicate" );
	selectWindow(nam);
	run("Split Channels");
	for (j=1; j<c;j+=1) {
		selectWindow("C"+j+"-"+nam);
		if (j>0){	
			run("stack detrend linear jru v1", "segments=1 maintain");
			selectWindow("Detrended");
        	rename("C"+j+"-"+nam+"_dtr");
        	roiManager("Select", roi_no);
        	run("Duplicate...", "title=temp duplicate" );
        	setBackgroundColor(0, 0, 0);
        	run("Clear Outside", "stack");
       		selectWindow("C"+j+"-"+nam+"_dtr");
        	close();
        	selectWindow("temp");
        	rename("C"+j+"-"+nam+"_dtr");
		}
		selectWindow("C"+j+"-"+nam);
		roiManager("Select", roi_no);
		run("Duplicate...", "title=temp duplicate" );
		setBackgroundColor(0, 0, 0);
		run("Clear Outside", "stack");	
		selectWindow("C"+j+"-"+nam);
		close();
		selectWindow("temp");
		rename("C"+j+"-"+nam);
	}
}

k=1;
selectWindow("C"+k+"-"+nam);
h = getHeight();
w = getWidth();
rows = Math.round(h/3);
cols = Math.round(w/3);
print(rows, cols);

roiManager("reset"); //clear all ROIs from manager to make space for grid ROI

x = 0;
y = 0;
width = 3;
height = 3;
spacing = 0;
numCol = cols; //note cols and rows are reversed in this macro
numRow = rows;


for(i = 0; i < numRow; i++)
{
                for(j = 0; j < numCol; j++)
                {
                                xOffset = j * (width + spacing);
                                print(xOffset);
                                yOffset = i * (height + spacing);
                                print(yOffset);
                                makeRectangle(x + xOffset, y + yOffset, width, height);
                                roiManager("Add");
                                if (roiManager("count") > 100000)
                                                {
                                                print("Maximum reached: 100000 entries have been created.");
                                                exit;
                                                }
                } 
}                             


for (i=1; i<c;i+=1) {
    nam = "C"+i;
    for (j=0; j<array_of_nam.length;j+=1) {
    	period = array_of_nam[j];
    
	    new_nam = nam+"-"+period;
	    selectWindow(new_nam);
	    roiManager("Select", 5);
	    roiManager("Deselect");
	    roiManager("multi-measure measure_all");
	    saveAs("Results", loc+new_nam+"_measure.csv");
	    
	    new_nam = nam+"-"+period+"_dtr";
	    selectWindow(new_nam);
	    roiManager("Select", 5);
	    roiManager("Deselect");
	    roiManager("multi-measure measure_all");
	    saveAs("Results", loc+new_nam+"_measure.csv");
    }
}

print(rows, cols);
f= File.open(loc+"dimensions.csv");
print(f, rows+","+cols);
File.close(f);

exit("Finished :\)");
