This folder provides the files necessary to take the frames from a reflection high energy electron diffraction 
(RHEED) video and extract the length, width, and intensity of each diffraction spot for each frame.  The 
features can then be used as parameters to build and test models.

The R script, “FunctionsForWidthLengthIntensity_OSF.R”, extracts the features: height, width, and intensity, 
from RHEED videos.  The file then lets you save the final data frame as a CSV file.  

The folder “Example_Frames” contains a few frames of a RHEED video as an example to run through 
“FunctionsForWidthLengthIntensity_OSF.R”.  These also show you the optimal type of RHEED patterns to use.  
For instance, only the first order Laue zone is showing and the diffraction spots are in a line horizontally 
with the direct spot above them.  

Three example result .csv files are also contained in this folder.  The file “Example_X_OSF.csv” is just the 
results from the data obtained about the width of the diffraction spots.  The file “Example_Y_OSF.csv” contains 
information about the length of the diffraction spots.  Finally, “Example_Result.csv” contains the final 
dataframe.  

If “FunctionsForWidthLengthIntensity_OSF.R” is used on several videos of different samples, the data frames 
(like “Example_Result.csv) could be combined together and then used in the final file, the R markdown file 
“Classification_OSF.Rmd”.  This file contains methods for separating each video into either testing or training.  
You can then build models using each training set.  The models in the script are linear discriminant analysis, 
quadratic discriminant analysis, support vector machine, and logarithmic regression.  The models can then each 
be evaluated using the testing data and the resulting confusion matrices.
