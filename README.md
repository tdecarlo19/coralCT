# coralCT
Software tool to analyze computerized tomography scans of coral skeletal cores for annual calcification and bioerosion rates

This project was developed in the Cohen Lab at Woods Hole Oceanographic Institution: http://www.whoi.edu/cohenlab/
with funding from the National Science Foundation (NSF) grants OCE 1041106 and OCE 1220529, The Nature Conservancy award PNA/WHOI061810, and an NSF Graduate Research Fellowship award to Thomas DeCarlo.

coralCT is a MATLAB code for coral calcification and bioerosion analyses. Please see the license, tutorial, and example analysis within the GitHub repository.


INSTRUCTIONS TO RUN IN MATLAB:

coralCT_v1_1(action, coralDirectory, dirLim, consisSTDdir, consisSTDdens, stdDirectory, saveFileName, genus, ...)

REQUIRED INPUTS

ACTION is an integer 1-10 that specifies which function is performed:

1 ==> 'Identify Bands' ==> A series of 2D slabs are displayed and the user can identify the locations of annual low- OR high- density bands by clicking on the images. Before each band is identified, the user can zoom and pan around the image. Once an appropriate field of view is set, press any key. Begin clicking on the first band (5-10 clicks per band per image is usually sufficient), and press 'enter' when done. Repeat for all bands visible in the image. If a mistake is made while clicking on a band, press '1' and then 'enter', and then re-define that band. When all bands on an image have been identified, press 'space bar' while defining the final band. After pressing 'enter' the next image will be displayed. If a mistake is made so that all of the bands on the current image need to be re-defined, press '2' and then 'space', and then 'enter' while defining a band, and a new image will appear.

2 ==> 'Edit Bands' ==> Allows the user to edit the bands on a coral for which the bands have already been defined with the 'Identify Bands' function. Bands can be deleted, added in the middle, or added to the bottom.

3 ==> 'Process Calcification' ==> Processing density, extension, and calcification for corals for which the annual bands have been defined by the 'Identify Bands? function.

4 ==> 'Process Annual Density' ==> Processes annual density for corals for which the annual bands have been defined by the 'Identify Bands function.

5 ==> 'Process Whole-Core Density' ==> Processes the mean density of an entire core. This does not require annual bands to be defined.

6 ==> 'Identify Borings' ==> A series of 2D slabs are displayed and the user can identify the locations of bioeroded regions can be identified. Each bioeroded region needs to be clicked just once.

7 ==> 'Process Identified Borings' ==> Processes the volume of the core that is bioeroded using the regions of defined bioerosion from the 'Identify Borings' function.

8 ==> 'Process All Borings' ==> Processes all regions that are missing from a cylinder in the core as bioerosion. It is recommended to use 'Identify Borings' and 'Process Identified Borings' rather than this function because small anomalies around the outside of the core are likely to be identified as bioerosion with this function.

9 ==> 'Plot Density Standard Curve' ==> Creates a standard curve to convert from Hounsfield units to density (g cm^-3), and plots the curve. If a consistency standard was used, the intercept of the calibration is fit to the consistency standard, while the slope is determined from the 'main' standards. This function does not contribute to any analyses, and only serves to check the density calibration.

10 ==> 'Split Band Files' ==> divides the annual density bands  identified using coralCT into 2 parts. Use this function to break long  cores into smaller pieces for processing, if processing all the bands  at once does not work.


CORALDIRECTORY is a structure that specifies the path names of CT scans. If the paths are saved as a text file 'coralDirectory.txt', then CORALDIRECTORY can be input as: importdata('coralDirectory.txt'). Alternatively, enter one scan as follows: {'/path...'}

DIRLIMIT specifies which entries in the directory to use. DIRLIMIT must be 2 positive integers with the second integer greater than or equal to the first. For example, [3 5] will process the third to fifth entries in the directory. If DIRLIMIT is empty by inputing as [], then all entries in the directory will be processed.

STDDIRECTORY specifies the path to where the CT scans of coral density standards can be found. This path should identify the location of the folder than contains the scans of multiple standards, each in their own subfolder. These subfolders must be named exactly as the 3-digit name of each standard. For example, '267', '169', or '235'.

CONSISSTDDIR is the directory of the consistency standard. The consistency standard is a coral standard that was scanned alongside the sample cores, and is used to adjust the intercept term of the standard curve. If a consistency standard was not used, enter as [].

CONSISSTDDENS is the known density (g/cm^3) of the consistency standard.  POR 225 is 1.0887 g/cm^3.

SAVEFILENAME is a string that specifies the name of the file to be saved or loaded. For the density band and bioerosion identification functions, files are saved with this name in the folder containing the images of a particular core. For other functions, files containing defined annual density bands or bioerosion will be loaded using this name.

GENUS is a string that specifies the coral genus. The first letter must be capitalized. The following genera are currently available: 'Porites', 'Favia', 'Diploastrea', 'Siderastrea', 'Montastrea', or 'Diploria' (only for whole-core density).

OPTIONAL INPUTS

*** Enter optional inputs in typical MATLAB double entry style. The input is first specified as one of the options below, followed by a comma, and then the input value ***

'densityCalibration' allows the user to define the equation to convert HU to density. This can be used if coral skeletal density standards are not available and the calibration equation for a particular CT scanner is known. The input is [a b] where a is the slope and b is the intercept in the equation HU = a*density+b. Based on the calibration presented in DeCarlo et al., (2015), this is [1485.5 -768.9].

'inputStandards' allows use of any set of coral skeletal standards. By default, coralCT uses the WHOI set of standards described in DeCarlo et al., 2015. A new set of standards can be used with the option input 'inputStandards' followed by a 2-column vector in which the first column in a numerical name for each standard and the second column in the acnpucepted density (g/cm^3) of each standard. For example, coralCT_v1_1(...,'inputStandards',[269,0.809; 167,1.165; 221,1.537]). If this optional input is used, the required 'STDDIRECTORY' input must point to a folder containing subfolders, each with a NUMERICAL name corresponding to the input standards information.

'bandAdjustTol' specifies how many mm above or below the user-defined low-density bands to search for true local density minima. Use 1 mm as a default. Currently, this only acts on Porites.

'bandType' specifies whether the user identifies high-density or low-density bands as annual bands. The default value is 1, which indicates low-density bands. Enter a value of 2 for high-density bands.

'dispMov' specifies whether movies of core and corallite identification are displayed. An input of 1 will display movies, 2 will display and save movies, anything else does nothing. These movies will significantly slow down processing, but can help check that the code is working properly for a particular core.

'distTol' specifies the minimum length (mm) of a corallite to be considered. This can be used to filter out short corallites and improve processing speed. Use 1 mm as a default.

'cracksTol' specifies the minimum  cross-sectional area of a core that must be missing to begin considering if there is a crack. Use 1 as a default.

'labInt' specifies the axis label interval (cm) on printed figures

'nSlabs' specifies how many images are displayed when identifying annual density bands and bioerosion. The actual number of images will be the input integer plus 1. It is recommended that NSLABS is at least 3 for band identification and at least 6 for bioerosion identification.

'threshOnOff' specifies whether to use Otsu's thresholding method. This affects how the core is identified. For certain genus, thresholding is toggled automatically. This input allows the user to force a certain threshold for genus that normally use Otsu's method (e.g. Porites). This should only be toggled if regions of the core that contain skeleton are not being defined as core. Use 0 as a default. Enter 1 to turn forced thresholding on.

'coreTolHU' specifies the threshold Hounsfeld Units to use if 'THRESHONOFF' is set to 1. If 'THRESHONOFF' is 0, CORETOLHU does nothing. Use -300 as a default. Use a lower number (e.g. -400) to include more regions are skeleton.

'maxHdist' specifies to maximum horizontal distance (mm) a corallite is allowed to move between successive images. Use 1 mm as a default.

==============

EXAMPLE USEAGE (these examples are not exhaustive of all functions of coralCT, but rather serve as a guide to become accustomed with how to interact with the program)

* Set the directories first
corDir = importdata('coralDirectory.txt');
stdDir = '/Users/tomdecarlo/Documents/WHOI/CT Scans/Taiwan/Standards 18-Oct-13';
consisStd = '/Users/tomdecarlo/Documents/WHOI/CT scans/Taiwan/std 225 scans/756 scan';

* Without using consistency standard
coralCT_v1_1(1,corDir,[],stdDir,[],[],'bands_TMD_May_2013','Porites')

* Analyze bands with a pre-determined density calibration
coralCT_v1_1(3,corDir,[],[],[],[],'bands_TMD_May_2013','Porites','densityCalibration',[1485.5 -768.9])

* Using consistency standard
coralCT_v1_1(1,corDir,[],stdDir,consisStd,1.0887,'bands_TMD_May_2013','Porites')

* With some optional inputs
coralCT_v1_1(1,corDir,[],stdDir,[],[],'bands_TMD_May_2013','Porites','nSlabs',4)

* Processing bands that were already defined
coralCT_v1_1(3,corDir,[],stdDir,[],[],'bands_TMD_May_2013','Porites')

* Processing bands that were already defined, selecting only the 4th-6th corals in the directory
coralCT_v1_1(3,corDir,[4 6],stdDir,[],[],'bands_TMD_May_2013','Porites')

* Set the coral directory in the function input
coralCT_v1_1(3,importdata('coralDirectory.txt'),[],stdDir,[],[],'bands_TMD_May_2013','Porites')

* Check a standard curve 
coralCT_v1_1(9,[],[],stdDir,consisSTD,1.0887,[],[])

* Check a standard curve without consistency standard
coralCT_v1_1(9,[],[],stdDir,[],[],[],[])

* use a different set of standards besides the WHOI ones
coralCT_v1_1(3,corDir,[],stdDir,[],[],'bands_TMD_May_2013','Porites','inputStandards',[269 0.8095; 167 1.1655; 221 1.5374])



version 1.0 (initial public release)
Last Modified March 23, 2016, Thomas M DeCarlo (WHOI)

Please contact Thomas DeCarlo (tdecarlo@whoi.edu) with any problems, questions, or concerns.

PLEASE CITE AS:
	
GitHub/Zenodo code DOI:

AND

Peer-reviewed publication:
DeCarlo T.M., Cohen A.L., Barkley H., Cobban Q., Young C., Shamberger
K.E., Brainard R.E., Golbuu Y. (2015) Coral macrobioerosion is
accelerated by ocean acidification and nutrients. Geology 43 (1) 7-10.


