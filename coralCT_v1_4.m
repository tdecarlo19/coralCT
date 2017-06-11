function coralCT_v1_4(action,coralDirectory,dirLim,stdDirectory,consisSTDdir,consisSTDdens,saveFileName,genus,varargin)
%
% CORALCT VERSION 1.4
%
% CORALCT analysis tool for computerized tomography (CT) scans of coral
% skeleton cores
%
%   This project was developed in the Cohen Lab at Woods Hole Oceanographic
%   Institution: http://www.whoi.edu/cohenlab/ with funding from the 
%   National Science Foundation (NSF) grants OCE 1041106 and OCE 1220529, 
%   The Nature Conservancy award PNA/WHOI061810, and an NSF Graduate 
%   Research Fellowship award to Thomas DeCarlo.
%
%   coralCT_v1_4(action,coralDirectory,dirLim,consisSTDdir,consisSTDdens,stdDirectory,saveFileName,genus,...)
%
%
%
%   ===============
%   REQUIRED INPUTS
%   ===============
%
%
%   ACTION is an integer 1-10 that specifies which function is performed.
%
%   ============
%
%   1 ==> 'Identify Bands' ==> A series of 2D slabs are displayed and the
%   user can identify the locations of annual low- OR high- density bands
%   by clicking on the images. Before each band is identified, the user can
%   zoom and pan around the image. Once an appropriate field of view is
%   set, press any key. Begin clicking on the first band (5-10 clicks per
%   band per image is usually sufficient), and press 'enter' when done.
%   Repeat for all bands visible in the image. If a mistake is made while
%   clicking on a band, press '1' and then 'enter', and then re-define that
%   band. When all bands on an image have been identified, press 'space
%   bar' while defining the final band. After pressing 'enter' the next
%   image will be displayed. If a mistake is made so that all of the bands
%   on the current image need to be re-defined, press '2' and then 'space',
%   and then 'enter' while defining a band, and a new image will appear.
%
%   2 ==> 'Edit Bands' ==> Allows the user to edit the bands on a coral for
%   which the bands have already been defined with the 'Identify Bands'
%   function. Bands can be deleted, added in the middle, or added to the
%   bottom.
%   
%   3 ==> 'Process Calcification' ==> Processing density, extension, and
%   calcification for corals for which the annual bands have been defined
%   by the 'Identify Bands function.
%   
%   4 ==> 'Process Annual Density' ==> Processes annual density for corals 
%   for which the annual bands have been defined by the 'Identify Bands' 
%   function.
%
%   5 ==> 'Process Whole-Core Density' ==> Processes the mean density of an
%   entire core. Does not require annual bands to be defined.
%
%   6 ==> 'Identify Borings' ==> A series of 2D slabs are displayed and the
%   user can identify the locations of bioeroded regions can be identified.
%   Each bioeroded region needs to be clicked just once.
%
%   7 ==> 'Process Identified Borings' ==> Processes the % volume of the
%   core that is bioeroded using the regions of defined bioerosion from the
%   'Identify Borings' function.
%
%   8 ==> 'Process All Borings' ==> Processes all regions that are missing
%   from a cylinder in the core as bioerosion. It is recommended to use
%   'Identify Borings' and 'Process Identified Borings' rather than this
%   function because small anomalies around the outside of the core are
%   likely to be identified as bioerosion with this function.
%
%   9 ==> 'Plot Density Standard Curve' ==> Creates a standard curve to
%   convert from Hounsfield units to density (g cm^-3), and plots the
%   curve. If a consistency standard was used, the intercept of the
%   calibration is fit to the consistency standard, while the slope is
%   determined from the 'main' standards. This function does not contribute
%   to any analyses, and only serves to check the density calibration.
%
%   10 ==> 'Split Band Files' ==> divides the annual density bands 
%   identified using coralCT into 2 parts. Use this function to break long 
%   cores into smaller pieces for processing, if processing all the bands 
%   at once does not work.
%
%   ============
%
%   CORALDIRECTORY is a structure that specifies the path names of CT
%   scans. If the paths are saved as a text file 'coralDirectory.txt', then
%   CORALDIRECTORY can be input as: importdata('coralDirectory.txt').
%   Alternatively, enter one scan as follows: {'/path...'}
%
%   DIRLIMIT specifies which entries in the directory to use. DIRLIMIT must
%   be 2 positive integers with the second integer greater than or equal to
%   the first. For example, [3 5] will process the third to fifth entries
%   in the directory. If DIRLIMIT is empty by inputing as [], then all
%   entries in the directory will be processed.
%
%   STDDIRECTORY specifies the path to where the CT scans of coral density
%   standards can be found. This path should identify the location of the
%   folder than contains the scans of multiple standards, each in their own
%   subfolder. These subfolders must be named exactly as the 3-digit name
%   of each standard. For example, '267', '169', or '235'.
%
%   CONSISSTDDIR is the directory of the consistency standard. The
%   consistency standard is a coral standard that was scanned alongside the
%   sample cores, and is used to adjust the intercept term of the standard 
%   curve. If a consistency standard was not used, enter as [].
%
%   CONSISSTDDENS is the known density (g/cm^3) of the consistency
%   standard.  POR 225 is 1.0887 g/cm^3.
%
%   SAVEFILENAME is a string that specifies the name of the file to be
%   saved or loaded. For the density band and bioerosion identification
%   functions, files are saved with this name in the folder containing the
%   images of a particular core. For other functions, files containing
%   defined annual density bands or bioerosion will be loaded using this
%   name.
%
%   GENUS is a string that specifies the coral genus. The first letter must
%   be capitalized. The following genera are currently available: 'Porites',
%   'Favia', 'Diploastrea', 'Siderastrea', 'Montastrea', 
%   or 'Diploria' (only for whole-core density).
%
%   
%   ===============
%   OPTIONAL INPUTS
%   ===============
%
%   *** Enter optional inputs in typical MATLAB double entry style. The
%   input is first specified as one of the options below, followed by a
%   comma, and then the input value ***
%
%   'densityCalibration' allows the user to define the equation to convert
%   HU to density. This can be used if coral skeletal density standards are
%   not available and the calibration equation for a particular CT scanner
%   is known. The input is [a b] where a is the slope and b is the
%   intercept in the equation HU = a*density+b. Based on the calibration
%   presented in DeCarlo et al., (2015), this is [1485.5 -768.9].
%
%   'inputStandards' allows use of any set of coral skeletal standards. By
%   default, coralCT uses the WHOI set of standards described in DeCarlo et
%   al., 2015. A new set of standards can be used with the option input
%   'inputStandards' followed by a 2-column vector in which the first
%   column in a numerical name for each standard and the second column in
%   the acnpucepted density (g/cm^3) of each standard. For example,
%   coralCT(...,'inputStandards',[269,0.809; 167,1.165; 221,1.537]).
%   If this optional input is used, the required 'STDDIRECTORY' input must
%   point to a folder containing subfolders, each with a NUMERICAL name
%   corresponding to the input standards information.
%
%   'bandAdjustTol' specifies how many mm above or below the user-defined
%   low-density bands to search for true local density minima. Use 1 mm as
%   a default. Currently only acts on Porites.
%
%   'bandType' specifies whether the user identifies high-density or
%   low-density bands as annual bands. The default value is 1, which 
%   indicates low-density bands. Enter a value of 2 for high-density bands.
%
%   'dispMov' specifies whether movies of core and corallite identification
%   are displayed. An input of 1 will display movies, 2 will display and
%   save movies, anything else does nothing. These movies will
%   significantly slow down processing, but can help check that the code is
%   working properly for a particular core.
%
%   'distTol' specifies the minimum length (mm) of a corallite to be
%   considered. This can be used to filter out short corallites and improve
%   processing speed. Use 1 mm as a default.
%
%   'cracksTol' specifies the minimum % cross-sectional area of a core that
%   must be missing to begin considering if there is a crack. Use 1 % as a
%   default.
%
%   'labInt' specifies the axis label interval (cm) on printed figures
%
%   'nSlabs' specifies how many images are displayed when identifying annual
%   density bands and bioerosion. The actual number of images will be the
%   input integer plus 1. It is recommended that NSLABS is at least 3 for
%   band identification and at least 6 for bioerosion identification.
%
%   'threshOnOff' specifies whether to use Otsu's thresholding method. This
%   affects how the core is identified. For certain genus, thresholding
%   is toggled automatically. This input allows the user to force a certain
%   threshold for genus that normally use Otsu's method (e.g. Porites).
%   This should only be toggled if regions of the core that contain
%   skeleton are not being defined as core. Use 0 as a default. Enter 1 to
%   turn forced thresholding on.
%
%   'coreTolHU' specifies the threshold Hounsfeld Units to use if
%   'THRESHONOFF' is set to 1. If 'THRESHONOFF' is 0, CORETOLHU does
%   nothing. Use -300 as a default. Use a lower number (e.g. -400) to
%   include more regions are skeleton.
%
%   'maxHdist' specifies to maximum horizontal distance (mm) a corallite is
%   allowed to move between successive images. Use 1 mm as a default.
%
%   'voxelResize' changes the images size by a specified factor. For
%   example, setting this parameter to 0.5 would reduce the pixel spacing
%   in each direction by half, which makes the 3D image 1/8 its original
%   size. Use this for large datasets that fill your computer's memory.
%   Note that the actual resizing may be slightly different than the
%   specified value, depending on the voxel dimensions. See the output
%   files for the resulting voxel spacings.
%
%   ==============
%
%   EXAMPLE USEAGE (note that examples below are not exhaustive of all
%   functions; see also inputs list above)
%
%   * Set the directories first
%   corDir = importdata('coralDirectory.txt');
%   stdDir = '/Users/tomdecarlo/Documents/WHOI/CT Scans/Taiwan/Standards 18-Oct-13';
%   consisStd = '/Users/tomdecarlo/Documents/WHOI/CT scans/Taiwan/std 225 scans/756 scan';
%
%   * Without using consistency standard
%   coralCT_v1_4(1,corDir,[],stdDir,[],[],'bands_TMD_May_2013','Porites')
%
%   * Analyze bands with a pre-determined density calibration
%   coralCT_v1_4(3,corDir,[],[],[],[],'bands_TMD_May_2013','Porites','densityCalibration',[1485.5 -768.9])
%
%   * Using consistency standard
%   coralCT_v1_4(1,corDir,[],stdDir,consisStd,1.0887,'bands_TMD_May_2013','Porites')
%
%   * With some optional inputs
%   coralCT_v1_4(1,corDir,[],stdDir,[],[],'bands_TMD_May_2013','Porites','nSlabs',4)
%
%   * Processing bands that were already defined
%   coralCT_v1_4(3,corDir,[],stdDir,[],[],'bands_TMD_May_2013','Porites')
%
%   * Processing bands that were already defined, selecting only the 4th-6th corals in the directory
%   coralCT_v1_4(3,corDir,[4 6],stdDir,[],[],'bands_TMD_May_2013','Porites')
%
%   * Set the coral directory in the function input
%   coralCT_v1_4(3,importdata('coralDirectory.txt'),[],stdDir,[],[],'bands_TMD_May_2013','Porites')
%
%   * Check a standard curve 
%   coralCT_v1_4(9,[],[],stdDir,consisSTD,1.0887,[],[])
%
%   * Check a standard curve without consistency standard
%   coralCT_v1_4(9,[],[],stdDir,[],[],[],[])
%
%   * use a different set of standards besides the WHOI ones
%   coralCT_v1_4(3,corDir,[],stdDir,[],[],'bands_TMD_May_2013','Porites','inputStandards',[269 0.8095; 167 1.1655; 221 1.5374])
%
%   ============
%   version 1.1 (initial public release)
%   Last Modified March 23, 2016, Thomas M DeCarlo (WHOI)
%
%   version 1.2 (updated to work with SkyScan1176 microCT)
%   Last Modified June 03, 2017, Thomas M DeCarlo (UWA)
%
%   version 1.3 (updated for surface area output and voxel resizing)
%   *surface area code 'imSurface' by David Legland included here
%   Last Modified June 11, 2017, Thomas M DeCarlo (UWA)
%
%   version 1.4 (updated core filtering)
%   Last Modified June 11, 2017, Thomas M DeCarlo (UWA)
%
%   Please contact Thomas DeCarlo (thomas.decarlo@uwa.edu.au) 
%   with any problems, questions, or concerns.
%
%   PLEASE CITE AS:
%	
%	GitHub/Zenodo code DOI:
%   DeCarlo, Thomas M, & Cohen, Anne L. (2016, July 14). coralCT: software tool 
%   to analyze computerized tomography (CT) scans of coral skeletal cores for 
%   calcification and bioerosion rates. Zenodo. http://doi.org/10.5281/zenodo.595503
%
%	and
%
%	Peer-reviewed publication:
%   DeCarlo T.M., Cohen A.L., Barkley H., Cobban Q., Young C., Shamberger
%   K.E., Brainard R.E., Golbuu Y. (2015) Coral macrobioerosion is
%   accelerated by ocean acidification and nutrients. Geology 43 (1) 7-10.
%
%   DOI for DeCarlo (2017) application of this code in the journal Matters:
%   https://dx.doi.org/10.24433/CO.884ebe95-3c91-41af-a60e-9ba1fd765ffa
%
%   ============================
%
% The MIT License (MIT)
% 
% Copyright (c) 2016 Thomas DeCarlo
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

% ============





% Defaults
bandAdjustTol = 1;
dispMov = 3;
distTol = 1;
cracksTol = 1;
labInt = 5;
nSlabs = 3;
threshOnOff = 0;
coreTolHU = -400;
maxHdist = 1;
bandType = 1; % low density (high-density = 2)
inputStd = 1; % default WHOI standards (any other standards = 2)
nominalStd = [];
voxResize = 1;

% read input request
iInput = 1;
delete = 0;
densEqCheck = 0;
while iInput <= length(varargin)  
    
    delete = 0;
    
    if strcmp('densityCalibration',varargin(iInput))
        densEq = varargin{iInput+1};
        densEqCheck = 1;
        delete = 2;
    end
    if strcmp('bandAdjustTol',varargin(iInput))
        bandAdjustTol = varargin{iInput+1};
        delete = 2;
    end
    if strcmp('dispMov',varargin(iInput))
        dispMov = varargin{iInput+1};
        delete = 2;
    end
    if strcmp('distTol',varargin(iInput))
        distTol = varargin{iInput+1};
        delete = 2;
    end
    if strcmp('cracksTol',varargin(iInput))
        cracksTol = varargin{iInput+1};
        delete = 2;
    end
    if strcmp('labInt',varargin(iInput))
        labInt = varargin{iInput+1};
        delete = 2;
    end
    if strcmp('nSlabs',varargin(iInput))
        nSlabs = varargin{iInput+1};
        delete = 2;
    end
    if strcmp('threshOnOff',varargin(iInput))
        threshOnOff = varargin{iInput+1};
        delete = 2;
    end
    if strcmp('coreTolHU',varargin(iInput))
        coreTolHU = varargin{iInput+1};
        delete = 2;
    end
    if strcmp('maxHdist',varargin(iInput))
        maxHdist = varargin{iInput+1};
        delete = 2;
    end
    if strcmp('bandType',varargin(iInput))
        bandType = varargin{iInput+1};
        delete = 2;
    end
    if strcmp('inputStandards',varargin(iInput))
        inputStd = 2;
        nominalStd = varargin{iInput+1};
        delete = 2;
    end
    if strcmp('voxelResize',varargin(iInput))
        voxResize = varargin{iInput+1};
        delete = 2;
    end
    
    if delete > 0
        varargin(iInput:iInput-1+delete) = [];
    else
        iInput = iInput+1;
    end
end

% initialize and defaults
directory = coralDirectory;
if densEqCheck == 1
    HU2dens = densEq;
end
dirSet = dirLim;
standardOpen = stdDirectory;
cStd = consisSTDdir;
cStdDens = consisSTDdens;
saveName = saveFileName;
gen = genus;
bandTolSet = bandAdjustTol;
inpMov = dispMov;
distanceTolSet = distTol;
crackTol = cracksTol;
labIntSet = labInt;
numSlabs = nSlabs;
thresh = threshOnOff;
coreTol = coreTolHU;
maxJump = maxHdist;
inputStandards = inputStd;
voxRz = voxResize;


if isempty(dirSet)
    dirSet = [1 length(directory)];
end

% do what?
if action == 1
    IDbands
elseif action == 2
    reviseBands
elseif action == 3
    processBands
elseif action == 4
    processAnnualDensity
elseif action == 5
    processWholeCoreDensity
elseif action == 6
    selectBorings
elseif action == 7
    processBoringsID
elseif action == 8
    processBoringsAll
elseif action == 9
    plotStdCurve
elseif action == 10
    splitBands
end


function plotStdCurve
    
    global processTotal processNumber stdCurve consisStdDensityMeasured
    
    processTotal = 1;
    processNumber = 0;
    progressUpdate
    
    coralStandardCurve
    
    fStd = figure;
    
    % check correlation between known density and HU
    r2std = corrcoef(stdCurve(:,1),stdCurve(:,2));
    r2std = r2std(1,2);

    % plot fit line
    xstd = [min(stdCurve(:,1)),max(stdCurve(:,1))];
    hold on
    p1 = plot(xstd,xstd*HU2dens(1)+HU2dens(2),'r','LineWidth',1);
    p2 = plot(stdCurve(:,1),stdCurve(:,2),'ko');
    
    if length(consisStdDensityMeasured) == length(cStdDens)
        p3 = plot(cStdDens,consisStdDensityMeasured,'b*');
    end
    
    set(gca,'FontWeight','bold','FontSize',16)
    xlabel('Actual Density (g cm^{-3})','FontWeight','bold','FontSize',16)
    ylabel('Mean Hounsfeld Units','FontWeight','bold','FontSize',16)
    text(0.9,1100,strcat('r^2 = ^{}',num2str(r2std)),'FontWeight','bold','FontSize',16)
    if length(consisStdDensityMeasured) == length(cStdDens)
        legend([p1 p2 p3],'Standard curve','Main standards','Consistency standard','location','southeast')
    else
        legend([p1 p2],'Standard curve','Main standards','location','southeast')
    end
    
    legend boxoff
        
end



function IDbands
    
    global processTotal processNumber iLoop name fileOpen h3

    processTotal = 1;
    processNumber = 0;
    progressUpdate    
    
    for iLoop = dirSet(1):dirSet(2)
        
        processNumber = 0;
        
        fileOpen = strcat(directory(iLoop)); % read the path name
        
        loadData % load the DICOMs
        
        chooseCoreFilter
        
        manualBandID4loop % allow user to ID bands
        
        h3 = fspecial('gaussian',10, 10);
        
        bands3D4loop % interpolate bands in 3D from user coordinates
        
        displaySlabBands4loop % print a slab with user-defined bands
        
        extractName
        
        print(gcf,'-dpdf',fullfile(fileOpen{:},strcat('band_IDs_1_',name)),sprintf('-r%d',150)) % save a pdf of the slabs
        close
        print(gcf,'-dpdf',fullfile(fileOpen{:},strcat('band_IDs_2_',name)),sprintf('-r%d',150)) % save a pdf of the slabs
        close
    
    end
end


 function splitBands
     
     global iLoop fileOpen

     
     for iLoop = dirSet(1):dirSet(2)
     
        fileOpen = strcat(directory(iLoop)); % read the path name
        
        % make sure the filename ends with a '/'
        lastL = fileOpen{:};
        if lastL(end) ~= '/'
            fileOpen{:} = [fileOpen{:},'/'];
        end
        
        % load the previously defined annual band coordinates
        load0 = load(fullfile(fileOpen{:},saveName), 'userBands');
        userBands= load0.userBands;
        
        directory = coralDirectory;
        fileOpen = strcat(directory(iLoop)); % read the path name
     
        n = length(userBands(1,1,:));
     
        n2 = round(n/2);
     
        userBands1 = userBands(:,:,1:n2+1);
     
        userBands2 = userBands(:,:,n2:end);
     
        userBands = userBands1;
        
        save(fullfile(fileOpen{:},strcat(saveName,'_part_1')), 'userBands');
        
        userBands = userBands2;
     
        save(fullfile(fileOpen{:},strcat(saveName,'_part_2')), 'userBands');
        
     end

end




function processBands
    
    global fileOpen coreMovieName polypMovieName X Xfull hpxS h1 h3 userBands 
    global processTotal processNumber row layers slabDraw1 SA
    global topCore bottomCore pxS volume densityWholeCore iLoop name
    global densityTotal travelTotal calcifTotal densityPerYear polypsPerBand
    
    processTotal = 7;
    processNumber = 0;
    progressUpdate
    
    % open file to save output summary to
    fid = fopen('coralSummary.csv','w');
    
    % print column headers
    fprintf(fid,'%s\n',strcat('Using bands defined by ',saveName,'.m'));
    fprintf(fid,'%s\n',['name, ','path, ','whole-core mean density, ','whole-core mean HU, ','volume (cm^3), ','Surface area (cm^2),','layers, ','horizontal pixel spacing (mm), ','vertical pixel spacing (mm), ']);
    
    % main loop
    for iLoop = dirSet(1):dirSet(2) % loop through each coral in the directory
        
        processNumber = 0;
        
        tic
        
        fileOpen = strcat(directory(iLoop)); % read the path name
        
        % make sure the filename ends with a '/'
        lastL = fileOpen{:};
        if lastL(end) ~= '/'
            fileOpen{:} = [fileOpen{:},'/'];
        end
        
        % load the previously defined annual band coordinates
        load0 = load(fullfile(fileOpen{:},saveName), 'userBands');
        userBands= load0.userBands;
        
        directory = coralDirectory;
        fileOpen = strcat(directory(iLoop)); % read the path name
        
        extractName
        
        h3 = fspecial('gaussian',10, 10);
        
        % filename for movie of core and borings identification
        coreMovieName = {'core.avi'};
        
        % filename for movie of polyp identification
        polypMovieName = {'polyp.avi'};
        
        if ~densEqCheck % if calibration is required
            coralStandardCurve
        end
        
        loadData;
        
        chooseCoreFilter
        
        Xfull = X;
        
        % filter for polyp identification
        % first number is size (in voxels), second number is standard deviation
        % (in voxels). Enter 'hpxS' in workspace after loading a core to
        % see the length of one voxel (in mm). Enter 'surf(h1)' in workspace to
        % visualize what the filter looks like.
        if strcmp(gen,'Porites')
            h1 = fspecial('gaussian',round(10/hpxS*0.1), round(3/hpxS*0.1)); % use this for small porites cores
        elseif strcmp(gen,'Siderastrea')
            h1 = fspecial('gaussian',round(15/hpxS*0.1), round(12/hpxS*0.1)); % use for Siderastrea
        elseif strcmp(gen,'Montastrea')
            h1 = fspecial('gaussian',round(15/hpxS*0.1), round(5/hpxS*0.1)); % use for Montastrea
        elseif strcmp(gen,'Favia')
            h1 = fspecial('disk',round(14/hpxS*0.1)); % use for Favia
        elseif strcmp(gen,'Diploastrea')
            h1 = fspecial('disk',round(30/hpxS*0.1)); % use for Diploastrea
        end
        
        % Build core, check for cracks
        buildCore
        
        bandTol = round(bandTolSet/pxS);
        a = userBands(:);
        a(a==0) = [];
        bottomCore = floor(min(a))-bandTol;
        topCore = ceil(max(a))+bandTol;
        
        if bottomCore <= 0
            bottomCore = 1;
            bandTol = floor(min(a))-1;
        end
        
        if topCore > layers
            topCore = layers;
            bandTol = layers-ceil(max(a));
        end
        
        buildPolyps
        
        Xfull = []; % no longer needed
        
        % Put low-density bands into 3D
        if strcmp(gen,'Porites')
            bands3D4loopNoFilt
        end
        
        % Corallite density tracks
        if strcmp(gen,'Porites')
            firstBands
        end
        
        % Adjust band to find local corallite density minima
        if strcmp(gen,'Porites')
            adjustBands
        end
        
        % Redo bands in 3D
        if strcmp(gen,'Porites')
            bands3D4loop
        elseif strcmp(gen,'Favia') || ...
                strcmp(gen,'Diploastrea') || strcmp(gen,'Siderastrea') || ...
                strcmp(gen,'Montastrea')
            bands3D4loopNoFilt
        end
        
        % create slab
        slabDraw1 = zeros(row,1,layers);
        slabDraw1(:,:,1:layers)  = min(X(row/2:row/2+25,:,:));
        slabDraw1 = permute(slabDraw1,[3,1,2]);
        
        % average density between bands (all voxels between the bands)
        densityBetweenBands
        
        % Compute volume and density
        volumeDensity
        
        X = []; % no longer needed
        
        % Compute calcification rates
        calcification
        
        if strcmp(gen,'Favia') || strcmp(gen,'Diploastrea') % use 'whole-core' density for calcification
            calcifTotal = flipud(densityPerYear).*travelTotal;
        end
        
        % display a slab with polyp tracks and bands printed on top
        displayPolypSlab
        
        print(gcf,'-dpdf',fullfile(fileOpen{:},strcat('polypSlab_band_adjust_',name)),sprintf('-r%d',150))
        
        close
        
        % display all of the bands and 20% of the polyps in 3D
        displayBandsPolyps
        
        print(gcf,'-dpdf',fullfile(fileOpen{:},strcat('polyps3d_band_adjust_',name)),sprintf('-r%d',150))
        
        close
        
        % results plots
        resultsPlot
        print(gcf,'-dpdf',fullfile(fileOpen{:},strcat('results_band_adjust_',name)),sprintf('-r%d',150))
        
        close
        
        % save outputs as .csv file
        fprintf(fid, '%s\n',  strcat(name,', ',directory{iLoop},...
            ', ',num2str(densityWholeCore),', ',num2str(densityWholeCore*HU2dens(1)+HU2dens(2)),...
            ', ',num2str(volume),', ',num2str(SA),', ',num2str(layers),', ',num2str(hpxS),', ',num2str(pxS)));
        fprintf(fid,'%s\n','Top of Core');
        fprintf(fid,'%s\n',strcat('band',', ','Density (g cm^-3)',', ',...
            'Extension (cm)',', ','Calcification (g cm^-2 yr^-1)',', ',...
            'Polyps per band'));
        for iY = length(densityTotal):-1:1
            fprintf(fid,'%s\n', strcat(num2str(length(densityTotal)-iY+1),', ',...
                num2str(densityTotal(iY)),', ',num2str(travelTotal(iY)/10),', ',...
                num2str(calcifTotal(iY)/10),', ',num2str(polypsPerBand(iY))));
        end
        fprintf(fid,'%s\n','Bottom of Core');
        fprintf(fid,'%s\n',' ');
        fprintf(fid,'%s\n',' ');
        
        % save output in coral folder
        fid2 = fopen(fullfile(fileOpen{:},strcat('calcification_output_',saveName,'.csv')),'w');
        % print column headers
        fprintf(fid2,'%s\n',strcat('Using bands defined by ^{}',saveName,'.m'));
        fprintf(fid2,'%s\n',['name, ','path, ','whole-core mean HU, ','whole-core mean density, ','volume (cm^3), ','Surface area (cm^2),','layers, ','horizontal pixel spacing (mm), ','vertical pixel spacing (mm), ']);
        fprintf(fid2, '%s\n',  strcat(name,', ',directory{iLoop},...
            ', ',num2str(densityWholeCore),', ',num2str(densityWholeCore*HU2dens(1)+HU2dens(2)),...
            ', ',num2str(volume),', ',num2str(SA),', ',num2str(layers),', ',num2str(hpxS),', ',num2str(pxS)));
        fprintf(fid2,'%s\n','Top of Core');
        fprintf(fid2,'%s\n',strcat('band',', ','Density (g cm^-3)',', ',...
            'Extension (cm)',', ','Calcification (g cm^-2 yr^-1)',', ',...
            'Polyps per band'));
        for iY = length(densityTotal):-1:1
            fprintf(fid2,'%s\n', strcat(num2str(length(densityTotal)-iY+1),', ',...
                num2str(densityTotal(iY)),', ',num2str(travelTotal(iY)/10),', ',...
                num2str(calcifTotal(iY)/10),', ',num2str(polypsPerBand(iY))));
        end
        fprintf(fid2,'%s\n','Bottom of Core');
        fprintf(fid2,'%s\n',' ');
        fprintf(fid2,'%s\n',' ');
        coralSave
        t = toc;
        fprintf('\nAnalysis of %s\n was completed in %2.0f minutes',name,t/60)
        
    end
end



function processAnnualDensity
    
    global fileOpen coreMovieName polypMovieName X Xfull hpxS h3 
    global processTotal processNumber row layers slabDraw1 SA pxS
    global volume densityWholeCore densityTotal userBands iLoop name
    
    % loop through coral CT scans and process for annual skeletal density
    
    processTotal = 4;
    processNumber = 0;
    progressUpdate
    
    % open file to save output summary to
    fid = fopen('coralSummary.csv','w');
    
    % print column headers
    fprintf(fid,'%s\n',strcat('Using bands defined by ',saveName,'.m'));
    fprintf(fid,'%s\n',['name, ','path, ','whole-core mean density, ','whole-core mean HU, ','volume (cm^3), ','Surface area (cm^2),','layers, ','horizontal pixel spacing (mm), ','vertical pixel spacing (mm), ']);
    
    % main loop
    for iLoop = dirSet(1):dirSet(2) % loop through each coral in the directory
        
        processNumber = 0;
        
        tic
        
        fileOpen = strcat(directory(iLoop)); % read the path name
        
        % make sure the filename ends with a '/'
        lastL = fileOpen{:};
        if lastL(end) ~= '/'
            fileOpen{:} = [fileOpen{:},'/'];
        end
        
        % load the previously defined annual band coordinates
        load0 = load(fullfile(fileOpen{:},saveName), 'userBands');
        userBands = load0.userBands;
        
        directory = coralDirectory;
        fileOpen = strcat(directory(iLoop)); % read the path name
        
        extractName
        
        h3 = fspecial('gaussian',10, 10);
        
        % filename for movie of core and borings identification
        coreMovieName = {'core.avi'};
        
        % filename for movie of polyp identification
        polypMovieName = {'polyp.avi'};
        
        % Make standard curve
        if ~densEqCheck
            coralStandardCurve
        end
        
        % load data
        loadData;
                
        chooseCoreFilter
        
        Xfull = X;
        
        % Build core, check for cracks
        buildCore
         
        % Put low-density bands into 3D
        bands3D4loop
        
        % create slab
        slabDraw1 = zeros(row,1,layers);
        slabDraw1(:,:,1:layers)  = min(X(row/2:row/2+25,:,:));
        slabDraw1 = permute(slabDraw1,[3,1,2]);
        
        % average density between bands (all voxels between the bands)
        densityBetweenBands
        
        % Compute volume and density
        volumeDensity
        
        X = [];
       

        
        % results plots
        resultsPlotDensity
        print(gcf,'-dpdf',fullfile(fileOpen{:},strcat('results_band_adjust_',name)),sprintf('-r%d',150))
        
        close
        
        % save outputs
        
        % save outputs as .csv file
        fprintf(fid, '%s\n',  strcat(name,', ',directory{iLoop},...
            ', ',num2str(densityWholeCore),', ',num2str(densityWholeCore*HU2dens(1)+HU2dens(2)),...
            ', ',num2str(volume),', ',num2str(SA),', ',num2str(layers),', ',num2str(hpxS),', ',num2str(pxS)));
        fprintf(fid,'%s\n','Top of Core');
        fprintf(fid,'%s\n',strcat('band',', ','Density (g cm^-3)'));
        for iY = length(densityTotal):-1:1
            fprintf(fid,'%s\n', strcat(num2str(length(densityTotal)-iY+1),', ',...
                num2str(densityTotal(iY))));
        end
        fprintf(fid,'%s\n','Bottom of Core');
        fprintf(fid,'%s\n',' ');
        fprintf(fid,'%s\n',' ');
        
        % save output in coral folder
        fid2 = fopen(fullfile(fileOpen{:},'density_output.csv'),'w');
        % print column headers
        fprintf(fid2,'%s\n',strcat('Using bands defined by ^{}',saveName,'.m'));
        fprintf(fid2,'%s\n',['name, ','path, ','whole-core mean HU, ','whole-core mean density, ','volume (cm^3), ','Surface area (cm^2),','layers, ','horizontal pixel spacing (mm), ','vertical pixel spacing (mm), ']);
        fprintf(fid2, '%s\n',  strcat(name,', ',directory{iLoop},...
            ', ',num2str(densityWholeCore),', ',num2str(densityWholeCore*HU2dens(1)+HU2dens(2)),...
            ', ',num2str(volume),', ',num2str(SA),', ',num2str(layers),', ',num2str(hpxS),', ',num2str(pxS)));
        fprintf(fid2,'%s\n','Top of Core');
        fprintf(fid2,'%s\n',strcat('band',', ','Density (g cm^-3)'));
        for iY = length(densityTotal):-1:1
            fprintf(fid2,'%s\n', strcat(num2str(length(densityTotal)-iY+1),', ',...
                num2str(densityTotal(iY))));
        end
        fprintf(fid2,'%s\n','Bottom of Core');
        fprintf(fid2,'%s\n',' ');
        fprintf(fid2,'%s\n',' ');
        
        coralSave
        
        t = toc;
        
        fprintf('\nAnalysis of %s\n was completed in %2.0f minutes',name,t/60)
        
    end
end



function processWholeCoreDensity
    
    global fileOpen coreMovieName hpxS pxS
    global processTotal processNumber layers 
    global volume densityWholeCore iLoop name SA

    % open file to save output summary
    fid = fopen('coralSummary.csv','w');
    
    % print column headers
    fprintf(fid,'%s\n',strcat('Using bands defined by ',saveName,'.m'));
    fprintf(fid,'%s\n',['name, ','path, ','whole-core mean density, ','whole-core mean HU, ','volume (cm^3), ','Surface area (cm^2),','layers, ','horizontal pixel spacing (mm), ','vertical pixel spacing (mm), ']);
    
    processTotal = 3;
    processNumber = 0;
    progressUpdate
    
    % main loop
    
    for iLoop = dirSet(1):dirSet(2) % loop through each coral in the directory
        
        processNumber = 0;
        
        tic
        
        fileOpen = strcat(directory(iLoop)); % read the path name
        
        extractName
        
        % filename for movie of core and borings identification
        coreMovieName = {'core.avi'};
        if ~densEqCheck
            coralStandardCurve
        end
        
        % load data
        loadData;
        
        % choose filter
        chooseCoreFilter
        
        % Build core
        buildCore
        
        % Compute volume and density
        volumeDensity
        
        % save outputs
        
        % save outputs as .csv file
        fprintf(fid, '%s\n',  strcat(name,', ',directory{iLoop},...
            ', ',num2str(densityWholeCore),', ',num2str(densityWholeCore*HU2dens(1)+HU2dens(2)),...
            ', ',num2str(volume),', ',num2str(SA),', ',num2str(layers),', ',num2str(hpxS),', ',num2str(pxS)));

        % save output in coral folder
        fid2 = fopen(fullfile(fileOpen{:},'density_output.csv'),'w');
        % print column headers
        fprintf(fid2,'%s\n',strcat('Using bands defined by ^{}',saveName,'.m'));
        fprintf(fid2,'%s\n',['name, ','path, ','whole-core mean density, ','whole-core mean HU, ','volume (cm^3), ','Surface area (cm^2),','layers, ','horizontal pixel spacing (mm), ','vertical pixel spacing (mm), ']);
        %fid2 = fopen(fullfile(fileOpen{:},'density_whole_core_output.csv'),'w');
        fprintf(fid2, '%s\n',  strcat(name,', ',directory{iLoop},...
            ', ',num2str(densityWholeCore),', ',num2str(densityWholeCore*HU2dens(1)+HU2dens(2)),...
            ', ',num2str(volume),', ',num2str(SA),', ',num2str(layers),', ',num2str(hpxS),', ',num2str(pxS)));
        
        coralSave
        
        t = toc;
        
        fprintf('\nAnalysis of %s\n was completed in %2.0f minutes\n',name,t/60)
        
    end
end



function loadData
    
    global hpxS row col layers X distanceTol pxS fileOpen processNumber titleName

    % load data
    
    processNumber = processNumber + 1;
    titleName = 'loading CT data';
    
    [X,metadata,sliceLoc] = read_dcm3(fileOpen{:},1,voxRz);
    
    [row,col,layers] = size(X); % size of the image
    
    % MATLAB reads the core in 'upside down', so let's just reverse it
    X = flipdim(X,3);
    
    % image pixel spacing
    hpxS = metadata.PixelSpacing(1)/voxRz;
    
    % vertical pixel spacing (mm)
    sliceDif = median(sliceLoc(2:end)-sliceLoc(1:end-1));
    pxS = abs(sliceDif);
    if max(abs(min(sliceLoc(2:end)-sliceLoc(1:end-1))-sliceDif)) > 0.0001 || ...
            max(abs(max(sliceLoc(2:end)-sliceLoc(1:end-1))-sliceDif)) > 0.0001
        fprintf('WARNING: UNEVEN DICOM SPACING!')
        
        % sort X by slice location
        sliceLoc = flipdim(sliceLoc,2);
        [b,idx] = sort(sliceLoc);
        x2 = X(:,:,idx);
        X = x2;
        
    end
    
    labInt = labIntSet/pxS*10;
    
    % convert maxJump mm to voxels
    maxJump = maxJump/hpxS; % in voxels
    
    % convert distanceTol mm to voxels
    distanceTol = distanceTolSet/pxS; % convert to voxels

end



function [X,metadata,sliceLoc] = read_dcm3(dIn,p,voxRz) 
    
    % read DCM files from input directory into matrix
    
    inpath = dIn;
    
    % make sure the filename ends with a '/'
    if inpath(end) ~= '/'
        inpath = [inpath '/'];
    end
    
    % directory of subfolders within set path
    folders = dir(inpath);
    
    layerCount = 0; % keep track of where to write files in matrix
    
    allSlice = [];
    
    check1 = 1; % check for whether we have found image size
    
    % initialize
    X = [];
    
    breakCheck = 0;
    
    for j = 1:length(folders)
        
        % directory of DICOM files within subfolders
        D = dir([[inpath folders(j).name '/'] '*.dcm']);
        if length(D) < 1
            D = dir([[inpath folders(j).name '/'] '*.IMA']);
        end
        
        % remove the invisible files added by some USB drives:
        remove = [];
        for jj = 1:length(D)
            if strcmp('._',D(jj).name(1:2))
                remove = [remove jj];
            end
        end
        D(remove) = [];
        
        % check image size
        if length(D) && check1
            metadata = dicominfo([[inpath folders(j).name '/'] filesep D(1).name]);
            ro = ceil(metadata.Height * voxRz);
            co = ceil(metadata.Width * voxRz);
            check1 = 0;
        end
        
        
        % we know each image is roXco, initialize here
        checkX = 0;
        try isempty(X);
            nAvg = round(length(D)/ceil(length(D)*voxRz));
            nLayers = ceil(length(D)*(1/nAvg));
            X(:,:,end+1:end+nLayers) = 0;
            sliceLoc(end+1:end+nLayers) = 0;
            checkX = 1;
        catch
        end
        if checkX == 0 && check1 == 0 & ~isnan(ceil(length(D)*(1/nAvg)))
            nAvg = round(length(D)/ceil(length(D)*voxRz));
            nLayers = ceil(length(D)*(1/nAvg));
            X = zeros(ro,co,nLayers,'double');
            sliceLoc = zeros(1,nLayers);
        end
        
        skip = 0;
        
        if isnan(nLayers)
            nLayers = 0;
        end
        
        % iterating over each file, read the image and populate the appropriate
        % layer in matrix X
        for i1 = 1:nLayers
            
            skipCheck = 0;
            
            % read metadata
            if i1 == 1
                metadata = dicominfo([[inpath folders(j).name '/'] filesep D(1).name]);
            end
            
            metadata1 = dicominfo([[inpath folders(j).name '/'] filesep D((i1-1)*(nAvg)+1).name]);
            
            if min(abs(allSlice-metadata1.SliceLocation)) == 0
                skipCheck = 1;
            end
            
            % read DICOM
            x = NaN(ro,co,nAvg);
            for jj = 1:nAvg
                x(:,:,jj) = imresize(dicomread([[inpath folders(j).name '/'] filesep D((i1-1)*(nAvg)+1).name]),voxRz);
            end
            X(:,:,i1+layerCount-skip) = mean(x,3);
            sliceLoc(i1+layerCount-skip) = metadata1.SliceLocation;
            
            % delete DICOM if this is a repeat
            if skipCheck == 1
                X(:,:,i1+layerCount-skip) = [];
                sliceLoc(i1+layerCount-skip) = [];
            end
            
            if min(abs(allSlice-metadata1.SliceLocation)) == 0
                skip = skip + 1;
            end
            
            allSlice = [allSlice metadata1.SliceLocation];
            
        end
        if length(X)
            layerCount = length(X(1,1,:)); % keep track of size of X
        end
        
        if p == 1
            progressUpdate(j/length(folders))
        end
        
    end
    
    % now rescale all the intensity values in the matrix so that the matrix
    % contains the original intensity values rather than the scaled values that
    % dicomread produces
    X = X.*metadata.RescaleSlope + metadata.RescaleIntercept;
    
    %nAvg = length(X(1,1,:))/ceil(length(X(1,1,:))*voxRz);
    %X2 = NaN(ro,co,ceil(length(X(1,1,:))*voxRz));
    %for jj = 1:ro
    %    for kk = 1:co
    %        X(jj,kk,1:ceil(length(X(1,1,:))*voxRz)) = arrayfun(@(i) mean(X(jj,kk,i:i+nAvg-1),3),1:nAvg:length(X(1,1,:))-nAvg+1)';
    %    end
    %end
    %X(:,:,ceil(length(X(1,1,:))*voxRz)+1:end) = [];
    %X = X2;
    %X2=[];
 
end


function manualBandID4loop
    
    global layers X h2 hashMarks row col pxS hpxS
    global userBands totBands doUpdate slab1 mid slab2 i j ldbDraw1
    global thick contra fileOpen proj iBand f

    % User defined band locations
    
    chooseCoreFilter
    
    % set up slab positions: sample core for general location within image
    samp = round(layers/2);
    filteredXcore(:,:) = imfilter(X(:,:,samp), h2, 'replicate');
    level = 0;
    if thresh == 0
        level = graythresh(((filteredXcore-min(min(filteredXcore)))./max(max(filteredXcore-min(min(filteredXcore))))));%.*max(max(filteredXcore-min(min(filteredXcore))))+min(min(filteredXcore));%graythresh(filteredXcore);
        coralSamp = im2bw(((filteredXcore-min(min(filteredXcore)))./max(max(filteredXcore-min(min(filteredXcore))))),level);
    else
        coralSamp = ind2sub(size(filteredXcore(:,:)),filteredXcore(:,:) > coreTol);
    end
    
    [r,c] = find(coralSamp);
    
    [val,loc] = max(r);
    topMost = [c(loc),r(loc)];
    [val,loc] = min(r);
    bottomMost = [c(loc),r(loc)];
    [val,loc] = max(c);
    rightMost = [c(loc),r(loc)];
    [val,loc] = min(c);
    leftMost = [c(loc),r(loc)];
    
    center = [round((rightMost(1)+leftMost(1))/2),round((topMost(2)+bottomMost(2))/2)];
    
    spacing = linspace((numSlabs+1)/2-numSlabs,(numSlabs-1)/2,numSlabs)...
        .*round((mean(topMost(2))-mean(bottomMost(2)))/numSlabs-1);
    if numSlabs == 3
        spacing = round(spacing*0.9); % if only 3 slabs, bring a little closer
    end
    slab1 = round(mean(center)+spacing);
    slab2 = round(mean(center));
    
    maxBands = 100; % can adjust if needed
    userBands = zeros(row,col,maxBands);
    
    totBands = 0;
    
    ldbDraw1 = zeros(row,col,layers);
        
    contra = [-1000 1800];
    thick = 15;
    proj = 'min';
    mid = round((numSlabs+1)/2);
    ords = [mid, numSlabs+1, 1:mid-1 mid+1:numSlabs];
    numLoop  = numSlabs+1;
    iBand = 0;
    while iBand < numLoop
        setCo = 0;
        iBand = iBand+1;
        i = ords(iBand);
        
        if i < numSlabs+1
            slabDraw3 = zeros(row,1,layers);
            if strcmp(proj,'min')
                slabDraw3(:,:,1:layers)  = min(X(:,slab1(i)-thick:slab1(i)+thick,:),[],2);
            elseif strcmp(proj,'mean')
                slabDraw3(:,:,1:layers)  = mean(X(:,slab1(i)-thick:slab1(i)+thick,:),2);
            elseif strcmp(proj,'max')
                slabDraw3(:,:,1:layers)  = max(X(:,slab1(i)-thick:slab1(i)+thick,:),[],2);
            end
            slabDraw3 = permute(slabDraw3,[3,1,2]);
        else % the last slab is perpendicular to the others
            slabDraw3 = zeros(row,1,layers);
            if strcmp(proj,'min')
                slabDraw3(:,:,1:layers)  = min(X(slab2-thick:slab2+thick,:,:),[],1);
            elseif strcmp(proj,'mean')
                slabDraw3(:,:,1:layers)  = mean(X(slab2-thick:slab2+thick,:,:),1);
            elseif strcmp(proj,'max')
                slabDraw3(:,:,1:layers)  = max(X(slab2-thick:slab2+thick,:,:),[],1);
            end
            slabDraw3 = permute(slabDraw3,[3,1,2]);
        end
        
        
        f(iBand) = figure;
        %imagesc(flipud(slabDraw3))
        pcolor(slabDraw3)
        colormap(bone)
        shading interp
        set(f(iBand),'Units','inches');
        set(f(iBand), 'Position', [1 1 6 10]);
        %axis equal
        set(gca,'CLim',contra)
        
        topC = layers;
        
        int = get(gca,'YTick');
        sp = int(1)-int(2);
        set(gca,'YTick',layers-topC:-sp:layers)
        set(gca,'YTickLabel',(0:-sp:layers)*pxS/10)
        ylabel('Distance Downcore (cm)','FontWeight','bold','FontSize',16)
        set(gca,'FontWeight','bold','FontSize',16)
        
        button = 1; % reset button
        j = 0; % reset band number
        
        
        while setCo == 0;
            prompt = {'Enter min HU:','Enter max HU:','"max", "min", or "mean"',...
                'slab thickness (mm)','slab position (mm)','Enter "1" when done'};
            dlg_title = 'Input';
            num_lines = 1;
            if i ~= numSlabs+1
                slabPos = (slab1(i))*hpxS;
            else
                slabPos = (slab2)*hpxS;
            end
            def = {num2str(contra(1)),num2str(contra(2)),proj,...
                num2str((thick+1)*2*hpxS),num2str(slabPos),'0'};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            contra(1) = str2num(answer{1});
            contra(2) = str2num(answer{2});
            setCo = str2num(answer{end});
            thick = round((str2num(answer{4})/2/hpxS)-1);
            if i ~= numSlabs+1
                slab1(i) = round(str2num(answer{5})/hpxS);
            else
                slab2 = round(str2num(answer{5})/hpxS);
            end
            proj = answer{3};
            
            if i < numSlabs+1
                slabDraw3 = zeros(row,1,layers);
                if strcmp(proj,'min')
                    slabDraw3(:,:,1:layers)  = min(X(:,slab1(i)-thick:slab1(i)+thick,:),[],2);
                elseif strcmp(proj,'mean')
                    slabDraw3(:,:,1:layers)  = mean(X(:,slab1(i)-thick:slab1(i)+thick,:),2);
                elseif strcmp(proj,'max')
                    slabDraw3(:,:,1:layers)  = max(X(:,slab1(i)-thick:slab1(i)+thick,:),[],2);
                end
                slabDraw3 = permute(slabDraw3,[3,1,2]);
            else % the last slab is perpendicular to the others
                slabDraw3 = zeros(row,1,layers);
                if strcmp(proj,'min')
                    slabDraw3(:,:,1:layers)  = min(X(slab2-thick:slab2+thick,:,:),[],1);
                elseif strcmp(proj,'mean')
                    slabDraw3(:,:,1:layers)  = mean(X(slab2-thick:slab2+thick,:,:),1);
                elseif strcmp(proj,'max')
                    slabDraw3(:,:,1:layers)  = max(X(slab2-thick:slab2+thick,:,:),[],1);
                end
                slabDraw3 = permute(slabDraw3,[3,1,2]);
            end
            %imagesc(flipud(slabDraw3))
            pcolor(slabDraw3)
            shading interp
            
            topC = layers;
            
            int = get(gca,'YTick');
            sp = int(1)-int(2);
            set(gca,'YTick',layers-topC:-sp:layers)
            set(gca,'YTickLabel',(0:-sp:layers)*pxS/10)
            ylabel('Distance Downcore (cm)','FontWeight','bold','FontSize',16)
            set(gca,'FontWeight','bold','FontSize',16)
            %axis equal
            set(gca,'CLim',contra)
            drawnow
        end
        
        doUpdate = 0;
        hashMarks = 0;
        x1 = 0;
        while button ~= 32 % spacebar terminates loop
            j = j+1;
            
            pause % allow user to adjust zoom
            
            uistack(f(iBand),'top')
            
            ax1 = get(gca,'XLim');
            ax2 = get(gca,'YLim');
            
            if j == 1 && i ~= mid
                hashMarks = 1;
                updateBands
            else
                hashMarks = 0;
            end
            
            if j > 1 && bands == 1
                updateBands
            end
            
            set(gca,'YLim',ax2);
            set(gca,'XLim',ax1);
            
            drawnow
            
            x = [];
            y = [];
            b = [];
            b1 = 1;
            bands = 0;
            while b1 == 1 || b1 == 32 || b1 == 49 || b1 == 3 || b1 == 50
                figure(f(iBand))
                [x1,y1,b1] = ginput(1);
                x = [x x1];
                y = [y y1];
                b = [b b1];
                
                numClicks = length(b);
                for iClick = 1:numClicks
                    if b(iClick) ~= 1 && b(iClick) ~= 32 && b(iClick) ~= 49 &&...
                            b(iClick) ~= 3 && b(iClick) ~= 50
                        x(iClick) = [];
                        y(iClick) = [];
                        b(iClick) = [];
                        numClicks = numClicks-1;
                    end
                end
                
                if ~isempty(x) > 0
                    bands = 1;
                end
                
                if isempty(b1)
                    break
                end
                
                hold on
                
                if b1 == 1 || b1 == 3
                    scatter(x1,y1,100,3000)
                elseif b1 == 49
                    plot(x,y,'yx','markersize',15)
                elseif b1 == 32
                    doUpdate = 1;
                end
                
            end
            
            [xSort,idx] = sort(x);
            ySort = y(idx);
            plot(xSort,ySort,'-','Color','white')
            drawnow
            save(fullfile(fileOpen{:},saveName), 'userBands');
            
            if min(abs(b-32)) == 0 % finished with slab  - terminate while loop
                x(b==32) = []; % delete null locations
                y(b==32) = []; % delete null locations
                button = 32;
            end
            
            if min(abs(b-49)) == 0 % oops! redo this band
                x = []; % delete null locations
                y = []; % delete null locations
                j = j-1;
            end
            
            % keep track of how many bands identified for whole core
            clicks = length(b(b==1)) + length(b(b==3));
            if j > totBands && clicks > 1
                totBands = j;
            end
            
            % delete clicks off image
            inds = find(round(x)<1 | round(x)>row);
            y(inds) = [];
            x(inds) = [];
            inds2 = find((layers-round(y))<1 | round(y)<1);
            y(inds2) = [];
            x(inds2) = [];
            
            % assign inputs
            for k = 1:length(x)
                if i < numSlabs+1 % for first set of slabs (all but last one)
                    if round(x(k)) > 0 && (layers-round(y(k))) > 0
                        userBands(round(x(k)),slab1(i),j) = layers-round(y(k));
                        %userBands(round(x(k)),slab1(i),j) = round(y(k));
                    end
                else % reverse x and y assignment order
                    if round(x(k)) > 0 && (layers-round(y(k))) > 0
                        userBands(slab2,round(x(k)),j) = layers-round(y(k));
                        %userBands(slab2,round(x(k)),j) = round(y(k));
                    end
                end
            end
            
            skipUpdate = 0;
            if min(abs(b-32)) == 0 % finished with slab  - terminate while loop
                if min(abs(b-50)) == 0
                    skipUpdate = 1;
                    if i < numSlabs+1
                        userBands(:,slab1(i),:) = 0;
                    else
                        userBands(slab2,:,:) = 0;
                    end
                    if ~isempty(ldbDraw1)
                        ldbDraw1 = zeros(size(ldbDraw1));
                    end
                    ord = [mid, numSlabs+1, 1:mid-1 mid+1:numSlabs];
                    ind = find(ord == i);
                    ords = [ords(1:ind),ords(ind),ords(ind+1:end)];
                    numLoop = numLoop+1;
                end
            end
            
            if doUpdate == 1 && skipUpdate == 0;
                updateBands
            end
            
            text(15,mean(round(y)),num2str(j),'Color','white')
            text(470,mean(round(y)),num2str(j),'Color','white')
            text(256,mean(round(y)),num2str(j),'Color','m','FontWeight','bold','FontSize',10)
            text(100,mean(round(y)),num2str(j),'Color','m','FontWeight','bold','FontSize',10)
            text(390,mean(round(y)),num2str(j),'Color','m','FontWeight','bold','FontSize',10)
            
        end
    end
    
    % delete unused layers in userBands
    userBands(:,:,totBands+1:end) = [];
    userBands2 = zeros(size(userBands));
    for i = 1:length(userBands(1,1,:))
        userBands2(:,:,i) = userBands(:,:,length(userBands(1,1,:))-i+1);
    end
    userBands = userBands2;
    
    userBands = layers-userBands;
    userBands(userBands==layers) = 0;
    
    for jjj = 1:nSlabs+1
        close
    end
    
    save(fullfile(fileOpen{:},saveName), 'userBands');

end



function updateBands
    
    global userBands row col hashMarks layers totBands doUpdate slab1 mid 
    global slab2 i j ldbDraw1 iBand f
    
    userBands2 = userBands;
    
    if action ~= 2
        userBands2(:,:,totBands+1:end) = [];
    end
    
    blockLDBstore = userBands2;
    LDBdata = zeros(row,col,length(blockLDBstore(1,1,:)));
    LDBfilter = zeros(row,col,length(blockLDBstore(1,1,:)));
    
    % create row and col attachments here
    [rowMesh ,colMesh] = meshgrid(1:row,1:col);
    
    add2 = [];
    if doUpdate == 1
        add2 = min([j totBands]);
    elseif j > 1
        add2 = j-1;
    else
        add2 = 1;
    end
    
    if hashMarks == 1
        for i4 = 1:length(blockLDBstore(1,1,:))
            if i == numSlabs+1
                [r,c] = find(blockLDBstore(:,:,i4));
                
                v = zeros(1,length(r));
                for j2 = 1:length(r)
                    v(j2) = blockLDBstore(r(j2),c(j2),i4);
                end
                if length(v)==length(unique(v))  && ~isempty(v)
                    LDBfilter(:,slab2,i4) = round(interp1(r,v,slab2));
                    LDBdata = permute(LDBfilter,[1,2,3]);
                end
            else
                [r,c] = find(blockLDBstore(:,slab1(mid),i4));
                if isempty(r)
                    break
                end
                v = zeros(1,length(r));
                for j2 = 1:length(r)
                    v(j2) = blockLDBstore(r(j2),c(j2),i4);
                end
                if length(v)==length(unique(v))  && ~isempty(v)
                    LDBfilter(slab1(mid),:,i4) = round(interp1(r,v,slab2));
                    LDBdata = permute(LDBfilter,[2,1,3]);
                end
            end
        end
    end
    
    if hashMarks == 1
        add2 = 1:length(LDBdata(1,1,:));
    end
    
    % Loop to do all the surface trending
    for i4 = add2
        
        [r,c] = find(blockLDBstore(:,:,i4));
        
        v = zeros(1,length(r));
        for j2 = 1:length(r)
            v(j2) = blockLDBstore(r(j2),c(j2),i4);
        end
        
        warning('off','all');
        
        % 2-d interpolation
        if length(round(griddata(r,c,v,rowMesh,colMesh)))>1
            LDBdata(:,:,i4) = round(griddata(r,c,v,rowMesh,colMesh));
            LDBfilter(:,:,i4) = LDBdata(:,:,i4);
            LDBdata = permute(LDBfilter,[2,1,3]);
            if i == numSlabs+1
                LDBdata = permute(LDBfilter,[2,1,3]);
            end
        elseif i < numSlabs+1
            try round(interp1(r,v,min(r):max(r)));
                LDBfilter(min(r):max(r),c(1),i4) = round(interp1(r,v,min(r):max(r)));
                LDBdata = permute(LDBfilter,[1,2,3]);
            catch
            end
        elseif i == numSlabs+1 && j > 1
            try round(interp1(c,v,min(c):max(c)));
                LDBfilter(r(1),min(c):max(c),i4) = round(interp1(c,v,min(c):max(c)));
                LDBdata = permute(LDBfilter,[1,2,3]);
            catch
            end
        elseif i == numSlabs+1 && j == 1 && action ~= 2
            try round(interp1(r,v,slab2));
                LDBfilter(:,min(c)-2:max(c)+2,i4) = round(interp1(r,v,slab2));
                LDBdata = permute(LDBfilter,[1,2,3]);
            catch
            end
        end
        
        warning('on','all');
        
    end
    
    LDBdata(isnan(LDBdata)) = 0;
    
    % print
    if i < numSlabs+1
        ldbDraw2 = zeros(1,row,layers);
        if hashMarks == 1
            add2 = 1:length(LDBdata(1,1,:));
        end
        for i3 = add2
            [r,c] = find(LDBdata(:,:,i3));
            for j3 = 1:length(r)
                ldbDraw1(r(j3),c(j3),LDBdata(r(j3),c(j3),i3)) = 3000;
            end
            textMark = [];
            textLoc = [];
            if i ~= mid && hashMarks == 1 && LDBdata(slab2,slab1(i),i3)
                [textMark(i3),textLoc(i3)] = max(max(max(ldbDraw1(slab2,slab1(i):slab1(i),LDBdata(slab2,slab1(i),i3)),[],2)));
                if textMark(i3) > 2999
                    figure(f(iBand))
                    hold on
                    text(slab2+15,layers-LDBdata(slab2,slab1(i),i3),num2str(i3),'Color','yellow');
                end
            end
            if action == 2 && hashMarks == 1 && LDBdata(slab2,slab1(i),i3)
                [textMark(i3),textLoc(i3)] = max(max(max(ldbDraw1(slab2,slab1(i):slab1(i),LDBdata(slab2,slab1(i),i3)),[],2)));
                if textMark(i3) > 2999
                    figure(f(iBand))
                    hold on
                    text(slab2+15,layers-LDBdata(slab2,slab1(i),i3),num2str(i3),'Color','yellow');
                end
            end
        end
        ldbDraw2(:,:,1:layers)  = max(ldbDraw1(:,slab1(i):slab1(i),:),[],2);
        ldbDraw2 = permute(ldbDraw2,[3,2,1]);
        if i ~= mid && hashMarks == 1
            ldbDraw2(:,[1:slab2-3, slab2+3:row]) = 0;
            ldbDraw1 = zeros(size(ldbDraw1));
        end
        
    else
        
        LDBdata(isnan(LDBdata)) = 0;
        ldbDraw2 = zeros(1,row,layers);
        
        for i3 = add2
            [r,c] = find(LDBdata(:,:,i3));
            for j3 = 1:length(r)
                ldbDraw1(r(j3),c(j3),LDBdata(r(j3),c(j3),i3)) = 3000;
            end
            textMark = [];
            textLoc = [];
            if i ~= mid && hashMarks == 1 && LDBdata(slab2,slab1(mid),i3)
                [textMark(i3),textLoc(i3)] = max(max(max(ldbDraw1(slab2,slab1(mid),LDBdata(slab2,slab1(mid),i3)),[],2)));
                if textMark(i3) > 2999
                    figure(f(iBand))
                    hold on
                    text(slab2+15,layers-LDBdata(slab2,slab1(mid),i3),num2str(i3),'Color','yellow');
                end
            end
            
            if action == 2 && hashMarks == 1
                [r,c] = find(blockLDBstore(:,:,i3));
        
                v = zeros(1,length(r));
                for j2 = 1:length(r)
                    v(j2) = blockLDBstore(r(j2),c(j2),i3);
                end
                if length(unique(v))>1  && ~isempty(v)
                    try round(interp1(r,v,slab2));
                        LDBfilter(slab1(mid),:,i3) = round(interp1(r,v,slab2));
                        LDBdata = permute(LDBfilter,[2,1,3]);
                        [textMark(i3),textLoc(i3)] = max(max(max(ldbDraw1(slab2,slab1(mid),LDBdata(slab2,slab1(mid),i3)),[],2)));
                        if textMark(i3) > 2999
                            figure(f(iBand))
                            hold on
                            text(slab2+15,layers-LDBdata(slab2,slab1(mid),i3),num2str(i3),'Color','yellow');
                        end
                    catch
                    end
                end
            end
        end
        ldbDraw2(:,:,1:layers)  = max(ldbDraw1(slab2:slab2,:,:),[],1);
        ldbDraw2 = permute(ldbDraw2,[3,2,1]);
    end
    
    
    hold on
    
    [r,c] = find(ldbDraw2);
    a = ldbDraw2(ldbDraw2>0);
    %scatter(c,layers-r,1,a*3000)
    %drawnow

end



function reviseBands
    
    global layers X h2 hashMarks row col pxS hpxS iLoop h3
    global userBands totBandsIn slab1 mid slab2 j i ldbDraw1
    global processTotal processNumber slabDraw1 name LDBdata fileOpen iBand f
    
    processTotal = 1;
    processNumber = 0;
    progressUpdate

    chooseCoreFilter
    
    h3 = fspecial('gaussian',10, 10);

    for iLoop = dirSet(1):dirSet(2)
    
        fileOpen = strcat(directory(iLoop)); % read the path name
        loadData % load the DICOMs
        
        % load the previously defined annual band coordinates
        load0 = load(fullfile(fileOpen{:},saveName), 'userBands');
        userBands= load0.userBands;
        userBands = flipdim(userBands,3);
        
        userBands = layers-userBands;
        userBands(userBands==layers) = 0;
        
        % set up slab positions: sample core for general location within image
        samp = round(layers/2);
        filteredXcore(:,:) = imfilter(X(:,:,samp), h2, 'replicate');
        level = 0;
        if thresh == 0
            level = graythresh(filteredXcore);
            coralSamp = im2bw(filteredXcore,level);
        else
            coralSamp = ind2sub(size(filteredXcore(:,:)),filteredXcore(:,:) > coreTol);
        end
        
        [r,c] = find(coralSamp);
        
        [val,loc] = max(r);
        topMost = [c(loc),r(loc)];
        [val,loc] = min(r);
        bottomMost = [c(loc),r(loc)];
        [val,loc] = max(c);
        rightMost = [c(loc),r(loc)];
        [val,loc] = min(c);
        leftMost = [c(loc),r(loc)];
        
        center = [round((rightMost(1)+leftMost(1))/2),round((topMost(2)+bottomMost(2))/2)];
        
        spacing = linspace((numSlabs+1)/2-numSlabs,(numSlabs-1)/2,numSlabs)...
            .*round((mean(topMost(2))-mean(bottomMost(2)))/numSlabs-1);
        if numSlabs == 3
            spacing = round(spacing*0.9); % if only 3 slabs, bring a little closer
        end
        slab1 = round(mean(center)+spacing);
        slab2 = round(mean(center));
        
        totBandsIn = length(userBands(1,1,:));
        
        ldbDraw1 = zeros(row,col,layers);
        
        contra = [-1000 1800];
        thick = 15;
        proj = 'min';
        mid = round((numSlabs+1)/2);
        ords = [mid, numSlabs+1, 1:mid-1 mid+1:numSlabs];
        numLoop  = numSlabs+1;
        iBand = 0;
        setCo = 0;
        iBand = iBand+1;
        i = ords(iBand);
    
        if i < numSlabs+1
            slabDraw3 = zeros(row,1,layers);
            if strcmp(proj,'min')
                slabDraw3(:,:,1:layers)  = min(X(:,slab1(i)-thick:slab1(i)+thick,:),[],2);
            elseif strcmp(proj,'mean')
                slabDraw3(:,:,1:layers)  = mean(X(:,slab1(i)-thick:slab1(i)+thick,:),2);
            elseif strcmp(proj,'max')
                slabDraw3(:,:,1:layers)  = max(X(:,slab1(i)-thick:slab1(i)+thick,:),[],2);
            end
            slabDraw3 = permute(slabDraw3,[3,1,2]);
        else % the last slab is perpendicular to the others
            slabDraw3 = zeros(row,1,layers);
            if strcmp(proj,'min')
                slabDraw3(:,:,1:layers)  = min(X(slab2-thick:slab2+thick,:,:),[],1);
            elseif strcmp(proj,'mean')
                slabDraw3(:,:,1:layers)  = mean(X(slab2-thick:slab2+thick,:,:),1);
            elseif strcmp(proj,'max')
                slabDraw3(:,:,1:layers)  = max(X(slab2-thick:slab2+thick,:,:),[],1);
            end
            slabDraw3 = permute(slabDraw3,[3,1,2]);
        end
    
    
        f(iBand) = figure;
        %imagesc(flipud(slabDraw3))
        pcolor(slabDraw3)
        shading interp
        colormap(bone)
        set(f(iBand),'Units','inches');
        set(f(iBand), 'Position', [1 1 6 10]);
        %axis equal
        set(gca,'CLim',contra)
    
        topC = layers;

        int = get(gca,'YTick');
        sp = int(1)-int(2);
        set(gca,'YTick',layers-topC:-sp:layers)
        set(gca,'YTickLabel',(0:-sp:layers)*pxS/10)
        ylabel('Distance Downcore (cm)','FontWeight','bold','FontSize',16)
        set(gca,'FontWeight','bold','FontSize',16)
    
        button = 1; % reset button
        j = 0; % reset band number
    
    
        while setCo == 0;
            prompt = {'Enter min HU:','Enter max HU:','"max", "min", or "mean"',...
                'slab thickness (mm)','slab position (mm)','Enter "1" when done'};
            dlg_title = 'Input';
            num_lines = 1;
            if i ~= numSlabs+1
                slabPos = (slab1(i))*hpxS;
            else
                slabPos = (slab2)*hpxS;
            end
            def = {num2str(contra(1)),num2str(contra(2)),proj,...
                num2str((thick+1)*2*hpxS),num2str(slabPos),'0'};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            contra(1) = str2num(answer{1});
            contra(2) = str2num(answer{2});
            setCo = str2num(answer{end});
            thick = round((str2num(answer{4})/2/hpxS)-1);
            if i ~= numSlabs+1
                slab1(i) = round(str2num(answer{5})/hpxS);
            else
                slab2 = round(str2num(answer{5})/hpxS);
            end
            proj = answer{3};
        
            if i < numSlabs+1
                slabDraw3 = zeros(row,1,layers);
                if strcmp(proj,'min')
                    slabDraw3(:,:,1:layers)  = min(X(:,slab1(i)-thick:slab1(i)+thick,:),[],2);
                elseif strcmp(proj,'mean')
                    slabDraw3(:,:,1:layers)  = mean(X(:,slab1(i)-thick:slab1(i)+thick,:),2);
                elseif strcmp(proj,'max')
                    slabDraw3(:,:,1:layers)  = max(X(:,slab1(i)-thick:slab1(i)+thick,:),[],2);
                end
                slabDraw3 = permute(slabDraw3,[3,1,2]);
            else % the last slab is perpendicular to the others
                slabDraw3 = zeros(row,1,layers);
                if strcmp(proj,'min')
                    slabDraw3(:,:,1:layers)  = min(X(slab2-thick:slab2+thick,:,:),[],1);
                elseif strcmp(proj,'mean')
                    slabDraw3(:,:,1:layers)  = mean(X(slab2-thick:slab2+thick,:,:),1);
                elseif strcmp(proj,'max')
                    slabDraw3(:,:,1:layers)  = max(X(slab2-thick:slab2+thick,:,:),[],1);
                end
                slabDraw3 = permute(slabDraw3,[3,1,2]);
            end
            %imagesc(flipud(slabDraw3))
            pcolor(slabDraw3)
            shading interp

            topC = layers;

            int = get(gca,'YTick');
            sp = int(1)-int(2);
            set(gca,'YTick',layers-topC:-sp:layers)
            set(gca,'YTickLabel',(0:-sp:layers)*pxS/10)
            ylabel('Distance Downcore (cm)','FontWeight','bold','FontSize',16)
            set(gca,'FontWeight','bold','FontSize',16)
            %axis equal
            set(gca,'CLim',contra)
            drawnow
        end
        
        hashMarks = 1;
        

        updateBands
        
        prompt = {'Delete Band #:','Add Band Above #:','Enter 1 to add bands to bottom'};
            dlg_title = 'Do only one';
            num_lines = 1;
            def = {'0','0','0'};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            bandRemove = str2num(answer{1});
            bandAdd = str2num(answer{2});
            bandAddToBottom = str2num(answer{3});
        
            if bandRemove > 0
                userBands(:,:,bandRemove) = [];
            end
            
            if bandAdd > 0
                userBands(:,:,bandAdd+1:totBandsIn+1) = userBands(:,:,bandAdd:totBandsIn);
                userBandHoldOn = userBands;
                manualBandID4loop
                userBands = layers-userBands;
                userBands(userBands==layers) = 0;
                userBandHoldOn(:,:,bandAdd) = userBands(:,:,1);
                userBands = userBandHoldOn;
            end
            
            if bandAddToBottom > 0
                userBandHoldOn = userBands;
                manualBandID4loop
                userBands = layers-userBands;
                userBands(userBands==layers) = 0;
                userBandHoldOn(:,:,totBandsIn+1:totBandsIn+1+length(userBands(1,1,:))-1) = userBands;
                userBands = userBandHoldOn;
            end
            
            userBands = layers-userBands;
            userBands(userBands==layers) = 0;
        
            userBands = flipdim(userBands,3);
            
                    
    bands3D4loop % interpolate bands in 3D from user coordinates
    
    save(fullfile(fileOpen{:},saveName), 'userBands');
    
    % create slab
    slabDraw1 = zeros(row,1,layers);
    slabDraw1(:,:,1:layers)  = min(X(row/2:row/2+25,:,:));
    slabDraw1 = permute(slabDraw1,[3,1,2]);
    ylab = (floor(layers/labInt)*labInt:-labInt:0)*(pxS/10);
    ylabS = {num2str(ylab(1))};
    for i = 2:length(ylab)
        ylabS{i} = num2str(ylab(i));
    end

    ldbDraw1 = zeros(row,col,layers);
    ldbDraw2 = zeros(1,row,layers);
    for i = 1:length(LDBdata(1,1,:))
        for j = row/2-15:row/2+15
            for k = 1:col
                if LDBdata(j,k,i) > 1
                    ldbDraw1(j,k,LDBdata(j,k,i)) = 3000;
                end
            end
        end
    end
    ldbDraw2(:,:,1:layers)  = max(ldbDraw1(row/2:row/2+25,:,:));
    ldbDraw2 = permute(ldbDraw2,[3,2,1]);


    f3 = figure;
    %imagesc(flipud(slabDraw1+ldbDraw2))
    pcolor((slabDraw1+ldbDraw2))
    shading interp
    colormap(bone)
    set(f3,'Units','inches');
    set(f3, 'Position', [1 1 6 10]);
    axis equal
    set(gca,'CLim',[-1200 1800])
    set(gca,'YTick',rem(layers,labInt):labInt:ceil(layers/labInt)*labInt,'YTickLabel',ylabS)
    set(gca,'XTick',[])
    ylabel('height (cm)   ','FontWeight','bold','FontSize',16)
    set(gca,'FontWeight','bold','FontSize',16);
    ylabel('Distance Downcore (cm)')
    
    extractName
    
    print(gcf,'-dpdf',fullfile(fileOpen{:},strcat('band_IDs_1_',name)),sprintf('-r%d',150)) % save a pdf of the slabs
    close  
    print(gcf,'-dpdf',fullfile(fileOpen{:},strcat('band_IDs_2_',name)),sprintf('-r%d',150)) % save a pdf of the slabs
    userBands = []; % to make sure the band map is reset
    close 
            
    end
end
       


function coralSave
    
    global fileOpen fileSave
    
    % script to save workspace while working with coral CT scans
    
    warning('off','all'); % suppress warnings that large 3D matrices not saved
    
    lastL = fileOpen{:};
    if lastL(end) ~= '/'
        fileOpen{:} = [fileOpen{:},'/'];
    end
    fileSave = fullfile(fileOpen{:},strcat(saveName,'_processed'));
    save(fileSave)
    
    warning('on','all'); % turn warnings back on

end
        


function bands3D4loop
    
    global userBands row col h3
    global LDBdata
    
    % Script to do 2-D interpolation of user-identified annual density bands in
    % coral CT scans
    
    % Last modified 5/19/2013 TMD
    
    blockLDBstore = userBands;
    LDBdata = zeros(row,col,length(blockLDBstore(1,1,:)));
    LDBfilter = zeros(row,col,length(blockLDBstore(1,1,:)));
    
    % create row and col attachments here
    [rowMesh ,colMesh] = meshgrid(1:row,1:col);
    
    % check for and delete and bands ony defined on one slab
    toRem = [];
    for i = 1:length(blockLDBstore(1,1,:))
        
        [r,c] = find(blockLDBstore(:,:,i));
        
        v = [];
        for j = 1:length(r)
            v(j) = blockLDBstore(r(j),c(j),i);
        end
        
        % try 2-d interpolation
        check = 0;
        warning('off','all');
        a = round(griddata(r,c,v,rowMesh,colMesh));
        warning('on','all');
        if ~isempty(a)
            check = 1;
        end
        
        if check == 0
            toRem = [toRem, i];
        end
        
    end
    
    userBands(:,:,toRem) = [];
    blockLDBstore = userBands;
    
    % Loop to do all the surface trending
    for i = 1:length(blockLDBstore(1,1,:))
        
        [r,c] = find(blockLDBstore(:,:,i));
        
        v = [];
        for j = 1:length(r)
            v(j) = blockLDBstore(r(j),c(j),i);
        end
        
        % 2-d interpolation
        LDBdata(:,:,i) = round(griddata(r,c,v,rowMesh,colMesh));
        LDBfilter(:,:,i) = round(imfilter(LDBdata(:,:,i), h3));
    end
    LDBdata = permute(LDBfilter,[2,1,3]);

end

    

function coralStandardCurve
    
    global Xstan layStan densityStd processNumber titleName stdCurve 
    global consisStdDensityMeasured
    
    % Script to create the density standard curve for coral CT analysis

    processNumber = processNumber+1;
    titleName = 'Standard Curve';

    % nominal standard values
    if inputStandards==1
        nominal = [269 0.8095; 270 0.9927; 225 1.0887; 167 1.1655; 237 1.3550; ...
            220 1.2794; 214 1.3862; 221 1.5374; 235 1.3221];
    else
        nominal = nominalStd;
    end
    
    stdIds = dir(standardOpen);
    stdCurve = []; % column 1 is accepted density, column 2 is measured HU
    stdCount = 0;
    
    for iStd = 1:length(stdIds)
        lastL = standardOpen;
        if lastL(end) ~= '/'
            standardOpen = [standardOpen,'/'];
        end
        
        thisStd = [];
        for iNom = 1:length(nominal)
            if str2num(stdIds(iStd).name) == nominal(iNom,1)
                thisStd = iNom;
                [Xstan,metadataStan,slice] = read_dcm3([standardOpen,stdIds(iStd).name],2,voxRz);
                [rowStan,colStan,layStan] = size(Xstan);
                hpxS = metadataStan.PixelSpacing(1);
                stdCount = stdCount + 1;
                volumeDensityStd;
                stdCurve(stdCount,:) = [nominal(thisStd,2),densityStd];
                
            end
            
            
        end
        
        progressUpdate(iStd/length(stdIds))
        
    end
    
    HU2dens = polyfit(stdCurve(:,1),stdCurve(:,2),1); % slope, intercept
    r2std = corrcoef(stdCurve(:,1),stdCurve(:,2));
    r2std = r2std(1,2);
    
    if ~isempty(cStd)
        lastL = cStd;
        if lastL(end) ~= '/'
            cStd = [cStd,'/'];
        end
        [Xstan,metadataStan,slice] = read_dcm3(cStd,2,voxRz);
        [rowStan,colStan,layStan] = size(Xstan);
        hpxS = metadataStan.PixelSpacing(1);
        volumeDensityStd;
        consisStdDensityMeasured = densityStd;
        densityStd = (densityStd-HU2dens(2))/HU2dens(1);
        densAdjust = cStdDens-densityStd;
        huAdjust = densAdjust*HU2dens(1);
        HU2dens(2) = HU2dens(2) - huAdjust;
    end

end



function volumeDensityStd
    global Xstan layStan pxS hpxS volumeStd densityStd
    % Compute total volume and mean density of core
    
    % first, build the rest of the core not already built:
    coralStan = zeros(size(Xstan));
    
    hStd = fspecial('gaussian',12, 4);
    
    % sample the core to see if density is very low
    filteredXcore(:,:) = imfilter(Xstan(:,:,round(layStan/2)), hStd);
    level = graythresh(filteredXcore);
    coralSample = im2bw(filteredXcore,level);
    Xsamp = Xstan(:,:,round(layStan/2));
    densitySample = mean(Xsamp(coralSample==1));
    
    for i = 1:layStan
        filteredXcore(:,:) = imfilter(Xstan(:,:,i), hStd);
        level = graythresh(filteredXcore);
        coralStan(:,:,i) = im2bw(filteredXcore,level);
        if densitySample < 650 % too low for Otsu to work properly
            coralStan(:,:,i) = ind2sub(size(filteredXcore(:,:)),filteredXcore(:,:) > -200);
        end
    end
    
    % now calculate volume and density
    
    pxTstd = sum(sum(sum(coralStan))); % total voxels in core
    
    volumeStd = pxTstd*((hpxS/10)*(hpxS/10)*(pxS/10)); % volume in cm^3
    
    densityStd = mean(Xstan(coralStan==1));

end



function buildCore
    
    global fileOpen coreMovieName X layers row col h2 processNumber titleName corA center
    global coral cracks crackLayers
    
    % Script to build a map of where a core exists in a 3D coral CT scan
    
    processNumber = processNumber + 1;
    titleName = 'Mapping Core';
    
    warning('off','all');
    
    % if user input to display and save a movie
    if inpMov == 2
        writerObj = VideoWriter(strcat(fileOpen{:},coreMovieName{:})); % movie object
        writerObj.FrameRate = 10;
        open(writerObj);
        f1 = figure;
        set(f1,'Units','inches');
        set(f1, 'Position', [1 1 12 5]);
    end
    
    % if just displaying movie
    if inpMov == 1
        f1 = figure;
        set(f1,'Units','inches');
        set(f1, 'Position', [1 1 12 5]);
    end
    
    % scan the entire CT image
    bottomCore = 1;
    topCore = layers;
    
    % initialize filtered image
    filteredXcore = zeros(row,col);
    
    % initialize matrix to store indices of where coral, cracks exist
    coral = zeros(row,col,topCore);
    cracks = zeros(row,col,topCore);
    
    % keep track of size of coral in voxels and if an approximate cylinder
    ellA = zeros(length(bottomCore:topCore),1);
    corA = zeros(length(bottomCore:topCore),1);
    
    % store coordinates of extremes of outer edges of location of core
    center = zeros(topCore-bottomCore,2); % to store center of core in each layer
    leftMost = zeros(topCore-bottomCore,2);
    rightMost = zeros(topCore-bottomCore,2);
    topMost = zeros(topCore-bottomCore,2);
    bottomMost = zeros(topCore-bottomCore,2);
    
    % initialize crack variables and storage
    crackA = zeros(row,col);
    crackCheck = 0;
    crackExist = 0;
    crackLayers = zeros(1,layers);
    crackStart = 1;
    
    for i = bottomCore:topCore % loop through DICOM images
        
        % filter DICOM image
        filteredXcore(:,:) = imfilter(X(:,:,i), h2, 'replicate');
        
        % Use thresholding to identify core region. If 'thresh' == 0, uses
        % Otsu's method, otherwise user can set defined threshold value (in HU)
        level = 0;
        if thresh == 0
            level = graythresh(((filteredXcore-min(min(filteredXcore)))./max(max(filteredXcore-min(min(filteredXcore))))));%.*max(max(filteredXcore-min(min(filteredXcore))))+min(min(filteredXcore));%graythresh(filteredXcore);
            if (max(max(filteredXcore))-min(min(filteredXcore)))*level+min(min(filteredXcore)) > -800
                coral(:,:,i) = im2bw(((filteredXcore-min(min(filteredXcore)))./max(max(filteredXcore-min(min(filteredXcore))))),level);
            end
        else
            coral(:,:,i) = ind2sub(size(filteredXcore(:,:)),filteredXcore(:,:) > coreTol);
        end
        
        % set borders to 0
        coral(:,1,i) = 0;
        coral(:,row,i) = 0;
        coral(1,:,i) = 0;
        coral(col,:,i) = 0;
        
        if (max(max(filteredXcore))-min(min(filteredXcore)))*level+min(min(filteredXcore)) > -800 || sum(sum(coral(:,:,i))) > 500
             
            % coordinates of core
            [r,c] = find(coral(:,:,i));
            
            [val,loc] = max(r);
            topMost(i,1:2) = [c(loc),r(loc)];
            [val,loc] = min(r);
            bottomMost(i,1:2) = [c(loc),r(loc)];
            [val,loc] = max(c);
            rightMost(i,1:2) = [c(loc),r(loc)];
            [val,loc] = min(c);
            leftMost(i,1:2) = [c(loc),r(loc)];
            
            % determine the upperright, upperleft, lowerright, lowerleft
            [val,loc1] = max(r.*c); % note that 'val' does not matter for our purposes
            [val,loc2] = max((row-r).*(col-c));
            [val,loc3] = max((row-r).*c);
            [val,loc4] = max(r.*(col-c));
            
            % define the center.
            center(i,:) = [round((rightMost(i,1)+leftMost(i,1))/2),round((topMost(i,2)+bottomMost(i,2))/2)];
            
            % check if an ellipse can be defined
            check = 0;
            try fitellipse([leftMost(i,1:2);topMost(i,1:2);rightMost(i,1:2);bottomMost(i,1:2);...
                    c(loc1),r(loc1);c(loc2),r(loc2);c(loc3),r(loc3);c(loc4),r(loc4)]);
                check = 1;
            catch
            end
            
            if check == 1
                
                [z, a, b, alpha] = fitellipse([leftMost(i,1:2);topMost(i,1:2);rightMost(i,1:2);bottomMost(i,1:2);...
                    c(loc1),r(loc1);c(loc2),r(loc2);c(loc3),r(loc3);c(loc4),r(loc4)]);
                
                ellA(i) = pi*a*b; % calculate area of ellipse
                corA(i) = sum(sum(coral(:,:,i))); % area of coral
                
                if i > 1 && (ellA(i)-corA(i))/ellA(i)*100 > crackTol % is percent crack > tolerance% ?
                    if crackCheck == 0 % not previous in a crack
                        crackCheck = 2; % start potential crack
                        crackStart = i; % note which layer the crack started in
                        crackA = coral(:,:,crackStart-1)-coral(:,:,i); % define crack area
                        cracks(:,:,i) = coral(:,:,crackStart-1)-coral(:,:,i); % store crack
                    else % previousy in a crack
                        crackA = crackA+coral(:,:,crackStart-1)-coral(:,:,i); % add on to crack area
                        crackA(crackA>0) = 1; % reset crack area to binary
                        crackA(crackA<0) = 0;
                        cracks(:,:,i) = coral(:,:,crackStart-1)-coral(:,:,i); % store crack
                    end
                elseif crackCheck > 0 % potentially exited a crack
                    crackCheck = crackCheck - 1; % allow one layer to miss tolerance before reseting
                    if crackCheck == 0 % exited a crack
                        crackA = zeros(row,col); % reset crack area
                        if crackExist == 0 % if it was not a crack
                            cracks(:,:,crackStart:i) = 0; % delete the potential crack
                        end
                        crackExist = 0; % reset if crack was exited
                    end
                end
                
                if sum(crackA(coral(:,:,crackStart)==1)) > 0.75*sum(sum(sum(coral(:,:,i))))
                    crackExist = 1; % crack found
                    crackLayers(crackStart-1:i+1) = 1; % note which layers are crack layers
                end
                
                % plot
                if inpMov == 1 || inpMov == 2
                    figure(f1);
                    subplot(1,2,1)
                    imagesc(coral(:,:,i));
                    colorbar
                    hold on
                    try
                        plotellipse(z,a,b,alpha);
                    catch
                    end
                    text(10,10,strcat('layer  ', num2str(i), ' of  ', num2str(layers)),'Color','white')
                    scatter([center(i,1),leftMost(i,1),topMost(i,1),rightMost(i,1),bottomMost(i,1),...
                        c(loc1),c(loc2),c(loc3),c(loc4)],...
                        [center(i,2),leftMost(i,2),topMost(i,2),rightMost(i,2),bottomMost(i,2),...
                        r(loc1),r(loc2),r(loc3),r(loc4)])
                    hold off
                    subplot(1,2,2)
                    imagesc(X(:,:,i))
                    colormap bone
                    colorbar
                    set(gca,'CLim',[-500 1800])
                    drawnow
                    if inpMov == 2
                        frame = getframe(f1);
                        writeVideo(writerObj,frame)
                    end
                    
                end
                
            end
            
        end
        
        progressUpdate((i-bottomCore+1)/((topCore-bottomCore)))
        
    end
    
    X(cracks==1) = NaN; % delete the crack regions in X
    
    warning('on','all'); % turn warnings back on

end



function buildPolyps
    
    global Xfull X hpxS h1 crackLayers cracks coral row col layers
    global topCore bottomCore distanceTol processNumber titleName
    global minimaInd polypsX polypsY polypsXshort polypsYshort polypHouns
    global fileOpen coreMovieName
        
    % Script to trace corallite tracks through a coral core
    
    % Last modified TMD 5/19/2013
    
    processNumber = processNumber + 1;
    titleName = 'Tracing corallites';
    
    % if user input to display and save a movie
    if inpMov == 2
        writerObj = VideoWriter(strcat(fileOpen{:},coreMovieName{:},'_corallites')); % movie object
        writerObj.FrameRate = 10;
        open(writerObj);
        f1 = figure;
        set(f1,'Units','inches');
        %set(f1, 'Position', [1 1 12 5]);
    end
    
    % if just displaying movie
    if inpMov == 1
        f1 = figure;
        set(f1,'Units','inches');
        %set(f1, 'Position', [1 1 12 5]);
    end
    
    % initialize filtered image
    filteredX = zeros(row,col);
    
    % initialize matrix to store minima indices
    minimaInd = zeros(1,2,topCore-bottomCore);
    
    % initialize polyp matrices
    if hpxS < 0.1
        polypsX = NaN((topCore-bottomCore)*100,topCore); % x-coordinates
        polypsY = NaN((topCore-bottomCore)*100,topCore); % y-coordinates
    else
        polypsX = NaN((topCore-bottomCore)*500,topCore); % x-coordinates
        polypsY = NaN((topCore-bottomCore)*500,topCore); % y-coordinates
    end
    polypHouns = zeros(size(polypsX)); % stores density
    
    [squareX,squareY] = meshgrid(1:row,1:col); % put indices in a grid
    
    % keep track of where the corallites are actually growing within the polyp
    % matrices. Initialize with estimates of max number of corallites in a
    % single DICOM image, for speed
    if hpxS < 0.1
        stillGrowingMatrix = zeros(800,layers);
    else
        stillGrowingMatrix = zeros(2500,layers);
    end
    
    % keep track of number of new corallites
    births = 0;
    
    jumpStore2 = maxJump;
    
    for i = bottomCore:topCore % loop through DICOM images
        
        % filter the X data
        filteredX(:,:) = imfilter(Xfull(:,:,i), h1);
        
        % if favia, use double filtering technique
        if strcmp(gen,'Favia')
            h4 = fspecial('gaussian',20,6);
            filteredX(:,:) = imfilter(filteredX,h4);
        end
        
        % if diploastrea, use double filtering technique
        if strcmp(gen,'Diploastrea')
            h4 = fspecial('gaussian',110,15);
            filteredX(:,:) = imfilter(filteredX,h4);
        end
        
        % Find local watersheds (based on density)
        allMins = imregionalmin(filteredX);
        
        % clip watersheds outside of core
        min2 = coral(:,:,i) + allMins;
        min2(min2==1) = 0;
        
        % assign watershed centers to storage matrix
        [ro,co] = find(min2);
        count = length(ro);
        minimaInd(:,:,i) = 0;
        minimaInd(1:count,1,i) = ro;
        minimaInd(1:count,2,i) = co;
        
        
        % if there is a crack, continue building corallites. Corallites continue
        % growing directly upwards through the crack.
        if crackLayers(i) == 1
            for j = 2:row-1 % no mimima on image edges
                for k = 2:col-1 % no mimima on image edges
                    if max(cracks(j,k,i-1:i+5)) == 1
                        if min(abs(minimaInd(:,1,i-1)-j)) == 0
                            locs = find(minimaInd(:,1,i-1)==j);
                            if min(abs(minimaInd(locs,2,i-1)-k)) == 0
                                if count == 0 % no polyps assigned in this layer yet
                                    count = count+1;
                                    minimaInd(count,:,i) = [j, k]; % 1st column stores row, 2nd column stores col
                                elseif (min(abs(minimaInd(:,1,i)-j))+min(abs(minimaInd(:,2,i)-k))) < maxJump/2
                                    count = count+1;
                                    minimaInd(count,:,i) = [j, k]; % 1st column stores row, 2nd column stores col
                                end
                            end
                        end
                    end
                end
            end
        end
        
        % double the allowable jump during a crack
        if crackLayers(i) == 1
            maxJump = jumpStore2*2;
        else
            maxJump = jumpStore2;
        end
        
        noDuple2 = zeros(row,col); % used to check for new corallites
        stillGrowingInd = find(polypsX(:,i-1)>0); % keep track of where corallites exist
        
        % build a matrix to keep track of where corallites exist in matrices
        if isempty(stillGrowingInd) < 1
            stillGrowingMatrix(1:length(stillGrowingInd),i) = stillGrowingInd;
        end
        
        if length(stillGrowingInd) > 1
            
            % find the index within the minimaInd matrix that has x and y
            % coordinates most similar to those of the previous layer
            [nextInd, jumps] = knnsearch([minimaInd(:,1,i) minimaInd(:,2,i)],...
                [polypsX(stillGrowingInd,i-1) polypsY(stillGrowingInd,i-1)]);
            
            % polyps can only move a set distance horizontally
            nextInd(jumps>maxJump) = [];
            stillGrowingInd(jumps>maxJump) = [];
            
            % make sure nextInd is unique.
            while mode(nextInd)>min(nextInd)
                r = find(nextInd==mode(nextInd));
                % kill the shorter corallite
                age = [];
                for j = 1:length(r)
                    age(j) = sum((isnan(polypsX(r(j),find(polypsX(r(j),:))))~=1));
                end
                [val,loc] = min(age);
                stillGrowingInd(r(loc)) = [];
                nextInd(r(loc)) = [];
            end
            r = find(nextInd==1);
            while length(r) > 1
                age = [];
                for j = 1:length(r)
                    age(j) = sum((isnan(polypsX(r(j),find(polypsX(r(j),:))))~=1));
                end
                [val,loc] = min(age);
                stillGrowingInd(r(loc)) = [];
                nextInd(r(loc)) = [];
                r = find(nextInd==1);
            end
            
            % assign indices of corallite assignments
            if length(nextInd) > 1
                polypsX(stillGrowingInd,i) = minimaInd(nextInd,1,i);
                polypsY(stillGrowingInd,i) = minimaInd(nextInd,2,i);
                
                % keep track of which corallites were just assigned
                for j = 1:length(nextInd)
                    
                    noDuple2(minimaInd(nextInd(j),1,i),minimaInd(nextInd(j),2,i)) = 1;
                end
            end
        end
        
        % Corallite is 'born' if there exists a corallite within a layer that was not
        % connected to an existing corallite track, and if not right on the edge
        for k = 1:length(minimaInd(:,1,i)) % loop for number of corallite centers in this layer
            if minimaInd(k,1,i) > 0 % loop through corallites in this layer
                if coral(minimaInd(k,1,i),minimaInd(k,2,i),i) == 1
                    
                    % check that corallite has not already been assigned
                    if noDuple2(minimaInd(k,1,i),minimaInd(k,2,i)) == 0
                        % assign 'new' corallite to corallite vectors
                        births = births+1;
                        polypsX(births,i) = minimaInd(k,1,i);
                        polypsY(births,i) = minimaInd(k,2,i);
                    end
                end
            end
        end
        
        % update stillGrowingInd afters births
        stillGrowingInd = find(polypsX(:,i)>0); % keep track of where polyps exist
        
        % assign each voxel in this layer that is part of the core to a polyp
        [polypInd, jumps] = knnsearch(...
            [polypsY(stillGrowingInd,i),polypsX(stillGrowingInd,i)],...
            [squareX((coral(:,:,i)+cracks(:,:,i)==1)), ...
            squareY((coral(:,:,i)+cracks(:,:,i)==1))]);
        
        % assign X data of this layer to a new variable X1
        X1 = X(:,:,i);
        
        % 'unroll' X1 and only use data that is inside the core
        X2 = X1(((coral(:,:,i)+cracks(:,:,i))==1));
        
        % assign mean density for each corallite
        if ~isempty(polypInd)
            maxVox = length(polypInd(polypInd==mode(polypInd)));
            voxStore = nan(length(stillGrowingInd),maxVox);
            voxCount = ones(length(stillGrowingInd),1);
            if length(X2) > 2
                for k = 1:length(polypInd)
                    voxStore(polypInd(k),voxCount(polypInd(k))) = X2(k);
                    voxCount(polypInd(k)) = voxCount(polypInd(k))+1;
                end
                
                for k = 1:length(stillGrowingInd)
                    polypHouns(stillGrowingInd(k),i) = nanmean(voxStore(k,:));
                end
            end
        end
        
        % plot
        if inpMov == 1 || inpMov == 2
            figure(f1);
            hold off
            imagesc(X(:,:,i))
            colorbar
            axis equal
            set(gca,'CLim',[-500 1800])
            existband = (isnan(polypsY(:,i))~=1);
            hold on
            voronoi([1;row;1;row;polypsY(existband,i)],[1;col;col;1;polypsX(existband,i)]);
            colormap(bone)
            xlabel('x-coordinate (voxels)    ','FontWeight','bold','FontSize',16)
            ylabel('y-coordinate (voxels)    ','FontWeight','bold','FontSize',16)
            TT = colorbar('peer',gca);
            set(get(TT,'ylabel'),'String', 'Density (HU)','FontWeight','bold','FontSize',16);
            set(gca,'FontWeight','bold','FontSize',14)
            set(gca,'XTick',0:20:500)
            set(gca,'YTick',0:20:500)
            drawnow
            if inpMov == 2
                frame = getframe;
                writeVideo(writerObj,frame)
            end
        end
        
        progressUpdate((i-bottomCore+1)/((topCore-bottomCore)))
        
    end
    
    % Remove corallites that are very short
    
    % Remove short corallites
    nullTest = sort(abs(polypHouns),2,'descend'); % sort polypHouns by columns
    removeCount = 0;
    distanceTol2 = round(distanceTol);
    toRemove = [];
    for i = 1:length(polypHouns(:,1))
        if min(nullTest(i,1:distanceTol2)) == 0 % if there is a 0 (not connected
            % to a polyp center) within the tolerance distance, delete this
            % corallite
            removeCount = removeCount+1; % keep track of how many were removed
            toRemove(removeCount) = i; % indices of which corallites to remove
        end
    end
    polypHouns(toRemove,:) = []; % remove the short polyps from the matrix
    
    % keep tabs on which original corallites were retained:
    polypsXshort = polypsX;
    polypsYshort = polypsY;
    polypsXshort(toRemove,:) = [];
    polypsYshort(toRemove,:) = [];

end



function bands3D4loopNoFilt
    
    global userBands row col
    global LDBdata blockLDBstore
    
    % Script to do 2-D interpolation of user-identified annual density bands in
    % coral CT scans
    
    blockLDBstore = userBands;
    LDBdata = zeros(row,col,length(blockLDBstore(1,1,:)));
    LDBfilter = zeros(row,col,length(blockLDBstore(1,1,:)));
    % create row and col attachments here
    [rowMesh ,colMesh] = meshgrid(1:row,1:col);
    
    % check for and delete and bands ony defined on one slab
    toRem = [];
    for i = 1:length(blockLDBstore(1,1,:))
        
        [r,c] = find(blockLDBstore(:,:,i));
        
        v = [];
        for j = 1:length(r)
            v(j) = blockLDBstore(r(j),c(j),i);
        end
        
        % try 2-d interpolation
        check = 0;
        warning('off','all');
        a = round(griddata(r,c,v,rowMesh,colMesh));
        warning('on','all');
        if ~isempty(a)
            check = 1;
        end
        
        if check == 0
            toRem = [toRem, i];
        end
        
    end
    
    userBands(:,:,toRem) = [];
    blockLDBstore = userBands;
    
    % Loop to do all the surface trending
    for i = 1:length(blockLDBstore(1,1,:))
        
        [r,c] = find(blockLDBstore(:,:,i));
        
        v = [];
        for j = 1:length(r)
            v(j) = blockLDBstore(r(j),c(j),i);
        end
        
        % 2-d interpolation
        LDBdata(:,:,i) = round(griddata(r,c,v,rowMesh,colMesh));
        LDBfilter(:,:,i) = LDBdata(:,:,i);
    end
    LDBdata = permute(LDBfilter,[2,1,3]);

end



function firstBands
    
    global polypsXshort polypsYshort blockLDBstore bandTol hpxS LDBdata
    global processNumber titleName polypHouns
    global densityMatrix layT adjT layB adjB coralliteRef

    processNumber = processNumber + 1;
    titleName = 'Adjusting bands';
    
    polypsXshort(isnan(polypsXshort))=0;
    polypsYshort(isnan(polypsYshort))=0;
    
    % set the maximum number of corallites - this is done to initialize the
    % storage matrix of all corallites per band in order to improve speed.
    if hpxS < 0.1
        nP = 800;
    else
        nP = 2000;
    end
    polypIDs = zeros(length(blockLDBstore(1,1,:))-1,nP);
    
    % initialize matrix to store individual corallite density tracks
    densityMatrix = 0;
    coralliteRef = 0;
    
    % initialize variables to store extension, density, calcification per band
    calcifPolyp = 0;
    travelPolyp = 0;
    densityPolyp = 0;
    
    % loop through all bands
    for i = 1:length(blockLDBstore(1,1,:))
        
        count = 1; % counter for corallites between annual bands
        
        % loop through all corallites
        for j = 1:length(polypsXshort(:,1))
            
            % identify where the corallite exists
            thisBandExist = find(polypsXshort(j,1:end));
            
            % check if this corallite exists in this annual band
            if min(thisBandExist) > max(max(blockLDBstore(:,:,i))) || ...
                    max(thisBandExist) < min(min(blockLDBstore(:,:,i)))
                
            else % the corallite exists somewhere in this annual band
                
                % now check if the corallite starts below and ends above this band.
                % Since blockLDBstore is only sparsely filled with data
                % points, locate the data point in blockLDBstore nearest to
                % each corallite track. If the corallite track passes through the
                % band, it will count towards adjusting the bands. Since the
                % corallite tracks may not pass through the exact points of
                % blockLDBstore, use a cost function. If the polyp track passes
                % within 1 voxel of the band, considered that it crossed the band.
                costFunBot = nan(1,max(thisBandExist));
                for m = min(thisBandExist):max(thisBandExist)
                    costFunBot(m) = abs(min(min(m-LDBdata(polypsXshort(j,m),polypsYshort(j,m),i))));
                end
                
                % check if corallite cross band
                if min(costFunBot(min(thisBandExist):max(thisBandExist))) <= 1 ...
                        
                % This corallite will count towards band adjustment because it
                % extends across a band. Now, find the layer indices
                % where this corallite is nearest to the annual band.
                % First, tag which corallite this is
                polypIDs(i,count) = j;
                % Then find correct region to consider by identifying
                % where the corallite track is closest to the nearby
                % surface of an annual band.
                [b,layerBottom]  = min(costFunBot(min(thisBandExist):max(thisBandExist)));
                layerBottom = layerBottom+min(thisBandExist)-1; % need to add this because of how min returns index
                
                % compute distance traveled (extension rate)
                count2 = 0; % counter for steps per corallite
                
                % density counter (reset for each corallite counted)
                density = 0;
                adjustB = 0;
                if polypsXshort(j,layerBottom-bandTol) > 0
                    layerBottom = layerBottom-bandTol;
                    adjustB = bandTol;
                else
                    [val loc] = min(find(polypsXshort(j,layerBottom-bandTol:layerBottom)));
                    if ~isempty(val)
                        layerBottom = layerBottom-bandTol+val-1;
                        adjustB = bandTol-val+1;
                    end
                end
                
                layerTop = layerBottom+adjustB;
                adjustT = 0;
                if length(polypsXshort(1,:)) < layerBottom+adjustB+bandTol
                    layerTop = length(polypsXshort(1,:));
                    adjustT = length(polypsXshort(1,:))-(layerBottom+adjustB);
                else
                    if polypsXshort(j,layerBottom+adjustB+bandTol) > 0
                        layerTop = layerBottom+adjustB+bandTol;
                        adjustT = bandTol;
                    else
                        [val loc] = min(find(polypsXshort(j,layerBottom+adjustB:layerBottom+adjustB+bandTol)));
                        if ~isempty(val)
                            layerTop = layerTop+val;
                            adjustT = val-1;
                        end

                    end
                end
                
                for n = layerBottom:layerBottom+adjustB+adjustT % how many DICOM layers the corallite traveled
                    
                    count2 = count2+1; % counter for steps per corallite
                    
                    % corallite density
                    density(count2) = polypHouns(j,n); % HU
                    
                end
                
                % update storage matrix
                densityMatrix(count,1:length(density),i) = density; %(HU)
                coralliteRef(count,i) = j;
                adjB(count,i) = adjustB;
                layB(count,i) = layerBottom;
                adjT(count,i) = adjustT;
                layT(count,i) = layerTop;
                count = count+1; % counter of corallites per band
                end
            end
        end
        
        
        progressUpdate((i/length(blockLDBstore(1,1,:))))
        
    end

end



function adjustBands
    
    global userBands densityMatrix blockLDBstore adjB adjT layB
    global polypsXshort polypsYshort coralliteRef
    global adjustedBands
    
    % Adjust band location to match true density minima
    
    % Last modified 5/29/2013 TMD
    
    originalBands = userBands;
    adjustedBands = zeros(size(userBands));
    
    for i = 1:length(blockLDBstore(1,1,:)) % loop through bands
        
        if length(densityMatrix(1,1,:)) > i-1
            
            inds = find(densityMatrix(:,1,i));
            
            for k = 1:length(inds) % loop through corallites
                j = inds(k);
                
                % location of density minima at bottom band
                if adjB(j,i)==0 && adjT(j,i)==0
                    loc = 0;
                elseif bandType==1 % low-density
                    [val loc] = min(densityMatrix(j,1:adjB(j,i)+adjT(j,i),i));
                elseif bandType==2
                    [val loc] = max(densityMatrix(j,1:adjB(j,i)+adjT(j,i),i));
                end
                adjustedBands(polypsXshort(coralliteRef(j,i),layB(j,i)+loc),...
                    polypsYshort(coralliteRef(j,i),layB(j,i)+loc),i) = ...
                    layB(j,i)+loc;
            end
        end
    end
    
    userBands = adjustedBands;

end



function densityBetweenBands
    
    global X LDBdata processNumber titleName coral
    global densityPerYear densityTotal
    
    % Script to compute density between annual density bands in coral CT scans
    
    % Last modified TMD 5/19/2013
    
    processNumber = processNumber + 1;
    titleName = 'Annual density';
    
    % initialize storage
    densityTotal = zeros(length(LDBdata(1,1,:))-1,1);
    
    % loop through denisty bands less one because can only compute density
    % between consecutive bands
    for i = 1:length(LDBdata(1,1,:))-1
        
        % matrix of zeros same size as CT scan
        Xcor = zeros(size(X));
        
        % difference between consectutive bands
        betweenBands1 = LDBdata(:,:,i+1)-LDBdata(:,:,i);
        
        % bottom band
        betweenBands2 = LDBdata(:,:,i);
        
        % define number of DICOM images above bottom band
        betweenBands1 = LDBdata(:,:,i)+betweenBands1;
        
        % set to 0 any coordinate that is not covered by both bands
        betweenBands2(isnan(betweenBands1)==1) = 0;
        betweenBands1(isnan(betweenBands1)==1) = 0;
        
        % find coordinates covered by both bands
        [ro,co] = find(betweenBands1);
        
        % set region between bands to 1 in Xcor, only if inside core
        for j = 1:length(ro)
            Xcor(ro(j),co(j),betweenBands2(ro(j),co(j)):betweenBands1(ro(j),co(j))) ...
                = coral(ro(j),co(j),betweenBands2(ro(j),co(j)):betweenBands1(ro(j),co(j)));
        end
        
        % unroll region between bands
        Xcor2 = X(Xcor==1);
        
        % take mean of region between bands (nanmean excludes any cracks)
        mDens = nanmean(Xcor2);
        
        % assign mean density of this band
        densityTotal(i) = (nanmean(mDens)-HU2dens(2))/HU2dens(1); % mean density
        
        progressUpdate((i/(length(LDBdata(1,1,:))-1)))
        
    end
    
    densityPerYear = flipud(densityTotal); % reverse order

end



function volumeDensity
    
    global bottomCore topCore layers X h2 coral pxS hpxS
    global volume densityWholeCore SA row col
    
    % Script to compute volume and mean density of a coral CT scan. Requires
    % that 'buildCore' and 'coralStandardCurve' have already been executed.
    
    % check for where the core has not been built yet
    unBuilt = [1:bottomCore-1, topCore+1:layers];
    
    % add on to existing core where it is not already built
    for k = 1:length(unBuilt)
        i = unBuilt(k);
        filteredXcore(:,:) = imfilter(X(:,:,i), h2);
        if thresh == 0
            level = graythresh(((filteredXcore-min(min(filteredXcore)))./max(max(filteredXcore-min(min(filteredXcore))))));%.*max(max(filteredXcore-min(min(filteredXcore))))+min(min(filteredXcore));%graythresh(filteredXcore);
            if (max(max(filteredXcore))-min(min(filteredXcore)))*level+min(min(filteredXcore)) > -800
                coral(:,:,i) = im2bw(((filteredXcore-min(min(filteredXcore)))./max(max(filteredXcore-min(min(filteredXcore))))),level);
            end
        else
            coral(:,:,i) = ind2sub(size(filteredXcore(:,:)),filteredXcore(:,:) > coreTol);
        end
    end
    
    % now calculate volume and density
    
    pxT = sum(sum(sum(coral))); % total voxels in core
    
    volume = pxT*((hpxS/10)*(hpxS/10)*(pxS/10)); % volume in cm^3
    
    % mean whole-core density
    densityWholeCore = (mean(X(coral==1))-HU2dens(2))/HU2dens(1);

    %surface area (using imSurface by David Legland)
    [xi,yi,zi] = meshgrid(1:row,1:col,linspace(1,layers,layers*pxS/hpxS));
    eqPx = round(interp3(coral,xi,yi,zi));
    SA = imSurface(eqPx,13,[hpxS,hpxS,hpxS])/100;
    
end





function calcification
    
    global blockLDBstore row col layers polypsXshort polypsYshort polypHouns
    global densityMatrix processNumber titleName hpxS pxS LDBdata
    global polyps4plot tracks4plot travelPolyp densityPolyp calcifPolyp
    global travelTotal densityTotal calcifTotal polypIDs

    processNumber = processNumber + 1;
    titleName = 'Calcification';
    
    % initialize some storage matrices here, see their usage below
    calcifTotal = zeros(length(blockLDBstore(1,1,:))-1,1);
    travelTotal = zeros(length(blockLDBstore(1,1,:))-1,1);
    densityTotal = zeros(length(blockLDBstore(1,1,:))-1,1);
    polyps4plot = zeros(row,col,layers);
    tracks4plot = zeros(row,col,layers);
    polypsXshort(isnan(polypsXshort))=0;
    polypsYshort(isnan(polypsYshort))=0;
    
    % set the maximum number of corallites - this is done to initialize the
    % storage matrix of all corallites per band in order to improve speed.
    if hpxS < 0.1
        nP = 800;
    else
        nP = 2000;
    end
    polypIDs = zeros(length(blockLDBstore(1,1,:))-1,nP);
    
    % initialize matrix to store individual corallite density tracks
    densityMatrix = 0;
    
    % initialize variables to store extension, density, calcification per band
    calcifPolyp = 0;
    travelPolyp = 0;
    densityPolyp = 0;
    
    % loop through all bands less one because can only measure rates between
    % consecutive bands
    for i = 1:length(blockLDBstore(1,1,:))-1
        
        % initialize for speed
        travel = 0;
        calcif = 0;
        count = 1; % counter for corallites between annual bands
        
        % loop through all corallites
        for j = 1:length(polypsXshort(:,1))
            
            % identify where the corallite exists
            thisBandExist = find(polypsXshort(j,1:end));
            
            % check if this corallite exists in this annual band
            if min(thisBandExist) > max(max(blockLDBstore(:,:,i))) || ...
                    max(thisBandExist) < min(min(blockLDBstore(:,:,i+1)))
                
            else % the corallite exists somewhere in this annual band
                
                % now check if the corallite starts below and ends above this band.
                % Since blockLDBstore is only sparsely filled with data
                % points, locate the data point in blockLDBstore nearest to
                % each corallite track. If the corallite track passes through the
                % consecutive bands, it will count towards calcification. Since the
                % corallite tracks may not pass through the exact points of
                % blockLDBstore, use a cost function. If the polyp track passes
                % within 1 voxel of the band, considered that it crossed the band.
                costFunTop = nan(1,max(thisBandExist));
                costFunBot = nan(1,max(thisBandExist));
                for m = min(thisBandExist):max(thisBandExist)
                    costFunBot(m) = abs(min(min(m-LDBdata(polypsXshort(j,m),polypsYshort(j,m),i))));
                    costFunTop(m) = abs(min(min(m-LDBdata(polypsXshort(j,m),polypsYshort(j,m),i+1))));
                end
                
                % check if corallite cross both bands
                if min(costFunBot(min(thisBandExist):max(thisBandExist))) <= 1 ...
                        && min(costFunTop(min(thisBandExist):max(thisBandExist))) <= 1
                    
                    % This corallite will count towards calcification because it
                    % extends across consecutive LDBs. Now, find the layer indices
                    % where this corallite is nearest to the annual bands.
                    % First, tag which corallite this is
                    polypIDs(i,count) = j;
                    % Then find correct region to consider by identifying
                    % where the corallite track is closest to the nearby
                    % surface of an annual band.
                    [b,layerBottom]  = min(costFunBot(min(thisBandExist):max(thisBandExist)));
                    layerBottom = layerBottom+min(thisBandExist)-1; % need to add this because of how min returns index
                    [c,layerTop]  = min(costFunTop(min(thisBandExist):max(thisBandExist)));
                    layerTop = layerTop+min(thisBandExist)-1; % need to add this because of how min returns index
                    
                    % store location for plotting diagnostics
                    polyps4plot(polypsXshort(j,layerBottom),polypsYshort(j,layerBottom),layerBottom-4:layerBottom+4) = 3000;
                    polyps4plot(polypsXshort(j,layerTop),polypsYshort(j,layerTop),layerTop-4:layerTop+4) = 3000;
                    
                    % compute distance traveled (extension rate)
                    count2 = 0; % counter for steps per corallite
                    
                    % reset travel, calcif, and density counters (these
                    % are reset for each corallite counted)
                    travel = zeros(length(layerBottom:layerTop-1),1);
                    calcif = 0;
                    density = 0;
                    lastMeasure = layerBottom; % keep track of where measured
                    for n = layerBottom:layerTop-1 % how many DICOM layers the corallite traveled
                        
                        % store locations for tracing corallite tracks
                        tracks4plot(polypsXshort(j,n),polypsYshort(j,n),n) = 3000;
                        
                        count2 = count2+1; % counter for steps per corallite
                        
                        % measure extension every 10 layers, or when the end is
                        % reached to determine the distance since last
                        % measurement
                        if (mod(n,10) == 0 && strcmp(gen,'Porites')) || n == layerTop-1
                            dist = n-lastMeasure(end); % distance in layers
                            
                            if dist > 0
                                
                                % measure distance as square root of sum of
                                % squares, which is the Euclidian distance
                                
                                % perform SVD to fit line to data over previous
                                % 10 layers
                                
                                % define points in 'ext'
                                ext = [polypsXshort(j,n+1-dist:n+1);...
                                    polypsYshort(j,n+1-dist:n+1);n+1-dist:n+1]';
                                
                                x0 = mean(ext)'; % centroid of points
                                
                                % translate points
                                A = [(ext(:, 1) - x0(1)) (ext(:, 2) - x0(2))...
                                    (ext(:, 3) - x0(3))];
                                
                                % perform singular value decomposition
                                [U, S, V] = svd(A, 0);
                                
                                % maximum value of S
                                [s, ind] = max(diag(S));
                                
                                % corresponding vector
                                a = V(:, ind);
                                
                                % define value of t, as in f(t), as 0 (initial)
                                t = 0;
                                
                                % define coordinates at this t
                                p0 = a.*t;
                                
                                % define value of t, as in f(t), that corresponds
                                % to the top DICOM layer considered
                                t = dist/a(3);
                                
                                % define coordinates at this t
                                p1 = a.*t;
                                
                                % measure distance as square root of sum of
                                % squares, which is the Euclidian distance
                                travel(count2) = sqrt(sum(((abs(p0(1)...
                                    -p1(1)))*hpxS)^2+((abs(p0(2)-...
                                    p1(2)))*hpxS)^2+((abs(p0(3)-...
                                    p1(3)))*pxS)^2));
                                
                            end
                            
                            lastMeasure(end+1) = n; % reset
                        end
                        
                        % corallite calcification and density
                        calcif(count2) = polypHouns(j,n+1)*travel(count2);
                        density(count2) = polypHouns(j,n+1); % HU
                        
                    end
                    
                    % update storage matrices
                    travelPolyp(i,count) = sum(travel); % (mm)
                    densityPolyp(i,count) = nanmean(density(density>0)); % (HU)
                    calcifPolyp(i,count) = travelPolyp(count)*densityPolyp(count);
                    densityMatrix(count,1:length(density),i) = density; %(HU)
                    count = count+1; % counter of corallites per band
                end
            end
        end
        
        if count > 1
            % update annual extension, density, and calcification
            travelTotal(i) = nanmean(travelPolyp(i,find(travelPolyp(i,:)))); % mean polyp extension per band (mm)
            densityTotal(i) = (nanmean(((densityPolyp(i,find(densityPolyp(i,:))))))-HU2dens(2))./HU2dens(1); % mean density per band (Hounsfeld units)
            calcifTotal(i) = travelTotal(i)*densityTotal(i);
        end        
        
        progressUpdate((i/length(blockLDBstore(1,1,:))))
        
    end

end



function progressUpdate(fraction)
        
    global processNumber processTotal titleName pFig iLoop
        
    if processNumber == 0
        pFig = figure;
        set(pFig,'position',[1128 338+(468-468/7*processTotal) 286 468/7*processTotal])
    else
        figure(pFig)
        subplot(processTotal,1,processNumber)
        title([titleName,': ',num2str(round(fraction*100)),'%'])
        axis([0 1 0 1])
        if mod(iLoop,2)==0
            patch([0 0 fraction fraction],[0 1 1 0],[0 0 1])
        else
            patch([0 0 fraction fraction],[0 1 1 0],[1 0 0])
        end
        set(gca,'Xtick',[],'YTick',[])
        drawnow
    end        
end
 


function displaySlabBands4loop 
    

    global contra LDBdata X row col layers pxS slab1 slab2 thick proj mid

    slabDraw1 = zeros(row,1,layers);
    slabDraw1(:,:,1:layers)  = min(X(row/2:row/2+25,:,:));
    slabDraw1 = permute(slabDraw1,[3,1,2]);
 
    LDBdata2 = LDBdata;
    LDBdata2(isnan(LDBdata2)) = 0;

    ylab = (floor(layers/labInt)*labInt:-labInt:0)*(pxS/10);
    ylabS = {num2str(ylab(1))};
    for i = 2:length(ylab)
        ylabS{i} = num2str(ylab(i));
    end

    ldbDraw1 = zeros(row,col,layers);
    ldbDraw2 = zeros(1,row,layers);
    for i = 1:length(LDBdata(1,1,:))
        for j = slab2-thick:slab2+thick
            for k = 1:col
                if LDBdata(j,k,i) > 1
                    ldbDraw1(j,k,LDBdata(j,k,i)) = 3000;
                end
            end
        end
    end
    ldbDraw2(:,:,1:layers)  = max(ldbDraw1(slab2-thick:slab2+thick,:,:),[],1);
    ldbDraw2 = permute(ldbDraw2,[3,2,1]);
    

    slabDraw3 = zeros(row,1,layers);
    if strcmp(proj,'min')
        slabDraw3(:,:,1:layers)  = min(X(slab2-thick:slab2+thick,:,:),[],1);
    elseif strcmp(proj,'mean')
        slabDraw3(:,:,1:layers)  = mean(X(slab2-thick:slab2+thick,:,:),1);
    elseif strcmp(proj,'max')
        slabDraw3(:,:,1:layers)  = max(X(slab2-thick:slab2+thick,:,:),[],1);
    end
    slabDraw3 = permute(slabDraw3,[3,1,2]);

    f3 = figure;
    %imagesc(flipud(slabDraw3+ldbDraw2))
    pcolor((slabDraw3+ldbDraw2))
    shading interp
    colormap(bone)
    set(f3,'Units','inches');
    set(f3, 'Position', [1 1 6 10]);
    axis equal
    set(gca,'CLim',contra)
    
    topC = layers;

    int = get(gca,'YTick');
    sp = int(1)-int(2);
    set(gca,'YTick',layers-topC:-sp:layers)
    set(gca,'YTickLabel',(0:-sp:layers)*pxS/10)
    ylabel('Distance Downcore (cm)','FontWeight','bold','FontSize',16)
    set(gca,'FontWeight','bold','FontSize',16)   
    
    slabDraw3 = zeros(row,1,layers);
    if strcmp(proj,'min')
        slabDraw3(:,:,1:layers)  = min(X(:,slab1(mid)-thick:slab1(mid)+thick,:),[],2);
    elseif strcmp(proj,'mean')
        slabDraw3(:,:,1:layers)  = mean(X(:,slab1(mid)-thick:slab1(mid)+thick,:),2);
    elseif strcmp(proj,'max')
        slabDraw3(:,:,1:layers)  = max(X(:,slab1(mid)-thick:slab1(mid)+thick,:),[],2);
    end
    slabDraw3 = permute(slabDraw3,[3,1,2]);        
    
    ldbDraw1 = zeros(row,col,layers);
    ldbDraw2 = zeros(1,row,layers);
    for i = 1:length(LDBdata(1,1,:))
        for j = 1:row
            for k = slab1(mid)-thick:slab1(mid)+thick
                if LDBdata(j,k,i) > 1
                    ldbDraw1(j,k,LDBdata(j,k,i)) = 3000;
                end
            end
        end
    end
    ldbDraw2(:,:,1:layers)  = max(ldbDraw1(:,slab1(mid)-thick:slab1(mid)+thick,:),[],2);
    ldbDraw2 = permute(ldbDraw2,[3,2,1]);
            
    
    f3 = figure;
    %imagesc(flipud(slabDraw3+ldbDraw2))
    pcolor((slabDraw3+ldbDraw2))
    shading interp
    colormap(bone)
    set(f3,'Units','inches');
    set(f3, 'Position', [1 1 6 10]);
    axis equal
    set(gca,'CLim',contra)
    
    topC = layers;

    int = get(gca,'YTick');
    sp = int(1)-int(2);
    set(gca,'YTick',layers-topC:-sp:layers)
    set(gca,'YTickLabel',(0:-sp:layers)*pxS/10)
    ylabel('Distance Downcore (cm)','FontWeight','bold','FontSize',16)
    set(gca,'FontWeight','bold','FontSize',16)

end

        

function displayPolypSlab
    
    global layers pxS row col LDBdata tracks4plot slabDraw1 corA

    ylab = (floor(layers/labInt)*labInt:-labInt:0)*(pxS/10);
    ylabS = {num2str(ylab(1))};
    for i = 2:length(ylab)
        ylabS{i} = num2str(ylab(i));
    end
    
    ldbDraw1 = zeros(row,col,layers);
    ldbDraw2 = zeros(1,row,layers);
    for i = 1:length(LDBdata(1,1,:))
        for j = row/2-15:row/2+15
            for k = 1:col
                if LDBdata(j,k,i) > 1
                    ldbDraw1(j,k,LDBdata(j,k,i)) = 3000;
                end
            end
        end
    end
    ldbDraw2(:,:,1:layers)  = max(ldbDraw1(row/2:row/2+25,:,:));
    ldbDraw2 = permute(ldbDraw2,[3,2,1]);
    
    trackDraw1 = zeros(row,1,layers);
    trackDraw1(:,:,1:layers)  = max(tracks4plot(row/2:row/2+25,:,:));
    trackDraw1 = permute(trackDraw1,[3,1,2]);
    
    f3 = figure;
    %imagesc(flipud(slabDraw1+trackDraw1+ldbDraw2))
    pcolor((slabDraw1+trackDraw1+ldbDraw2))
    shading interp
    colormap(bone)
    set(f3,'Units','inches');
    set(f3, 'Position', [1 1 6 10]);
    axis equal
    TT = colorbar('peer',gca);
    set(gca,'CLim',[-1500 2000])
    set(get(TT,'ylabel'),'String', 'HU','FontWeight','bold','FontSize',16);
    set(gca,'YTick',rem(layers,labInt):labInt:ceil(layers/labInt)*labInt,'YTickLabel',ylabS)
    set(gca,'XTick',[])
    ylabel('height (cm)   ','FontWeight','bold','FontSize',16)
    set(gca,'FontWeight','bold','FontSize',16);

    topLim = 30000;
    topC = layers;
    for i2 = layers:-1:1
        if corA(i2) > topLim
            topC = i2;
            break
        end
    end
    int = get(gca,'YTick');
    if length(int) > 1
        sp = int(1)-int(2);
        set(gca,'YTick',layers-topC:-sp:layers)
        set(gca,'YTickLabel',(0:-sp:layers)*pxS/10)
        ylabel('Distance Downcore (cm)')
    end

end


function coralLoopBoringPlot
    
    global coral layers X h2 row col fileOpen slabDraw1 slab1 slab2 i name
        
    % User defined borings locations

    bottomCore = 1;

    i = round(layers/3);
    filteredXcore(:,:) = imfilter(X(:,:,i), h2);
    coral(:,:,i) = ind2sub(size(filteredXcore(:,:)),filteredXcore(:,:) > coreTol);

    coral(:,1,i) = 0;
    coral(:,row,i) = 0;
    coral(1,:,i) = 0;
    coral(col,:,i) = 0;
   
    [r,c] = find(coral(:,:,i));
    
    [val,loc] = max(r);
    topMost = [c(loc),r(loc)];
    [val,loc] = min(r);
    bottomMost = [c(loc),r(loc)];
    [val,loc] = max(c);
    rightMost = [c(loc),r(loc)];
    [val,loc] = min(c);
    leftMost = [c(loc),r(loc)];
    
    center = [round((rightMost(1)+leftMost(1))/2),round((topMost(2)+bottomMost(2))/2)];
    
    % set up slab positions
    spacing = linspace((numSlabs+1)/2-numSlabs,(numSlabs-1)/2,numSlabs)...
        .*round((topMost(2)-bottomMost(2))/(numSlabs-1));
    slab1 = round(center(2)+spacing);
    slab2 = round(center(2));


    for i = 1:numSlabs+1
    
    
        if i < numSlabs+1
        
            if slab1(i)+12 > row
                slab1(i) = row-12;
            end
            if slab1(i)-12 < 1
                slab1(i) = 13;
            end
        
            slabDraw1 = zeros(row,1,layers);
            slabDraw1(:,:,1:layers)  = min(X(:,slab1(i)-12:slab1(i)+12,:),[],2);
            slabDraw1 = permute(slabDraw1,[3,1,2]);
        else % the last slab is perpendicular to the others
            slabDraw1 = zeros(row,1,layers);
            slabDraw1(:,:,1:layers)  = min(X(slab2-12:slab2+12,:,:),[],1);
            slabDraw1 = permute(slabDraw1,[3,1,2]);
        end

        displayBorings
        
        print(gcf,'-dpdf',fullfile(fileOpen{:},strcat('borings_',name,'_',num2str(i)))) % save a pdf of the slabs
    end

end



function displayBorings

    global pxS row layers slab1 slab2 i borings slabDraw1 hpxS

    ylab = (floor(layers/labInt)*labInt:-labInt:0)*(pxS/10);
    ylabS = {num2str(ylab(1))};
    for j = 2:length(ylab)
        ylabS{j} = num2str(ylab(j));
    end

    borings(:,:,end+1:layers) = 0; % make same size as X

    if i < numSlabs+1
        borDraw1 = zeros(row,1,layers);
        borDraw1(:,:,1:layers)  = 3000.*max(borings(:,slab1(i)-12:slab1(i)+12,:),[],2);
        borDraw1 = permute(borDraw1,[3,1,2]);
    else
        borDraw1 = zeros(row,1,layers);
        borDraw1(:,:,1:layers)  = 3000.*max(borings(slab2-12:slab2+12,:,:),[],1);
        borDraw1 = permute(borDraw1,[3,1,2]);
    end    

    f5 = figure;
    subplot(121)
    imagesc(pxS,hpxS,flipud(slabDraw1))
    colormap(bone)
    set(f5,'Units','inches');
    set(f5, 'Position', [1 1 6 10]);
    set(gca,'CLim',[-1500 2000])
    set(gca,'YTick',rem(layers,labInt):labInt:ceil(layers/labInt)*labInt,'YTickLabel',ylabS)
    set(gca,'XTick',[])
    ylabel('Height (cm)   ','FontWeight','bold','FontSize',16)
    set(gca,'FontWeight','bold','FontSize',16);
    subplot(122)
    imagesc(pxS,hpxS,flipud(slabDraw1+borDraw1))
    colormap(bone)
    set(f5,'Units','inches');
    set(f5, 'Position', [1 1 6 10]);
    set(gca,'CLim',[-1500 2000])
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    set(gca,'FontWeight','bold','FontSize',16);
end



function displayBandsPolyps
    
    global  pxS row col LDBdata hpxS polypsXshort polypsYshort

    fRecon = figure;
    
    for i = 1:length(LDBdata(1,1,:))
        surf((1:row).*hpxS/10,(1:col).*hpxS/10,LDBdata(:,:,i).*pxS/10)
        hold on
        shading interp
        colormap(bone)
        if min(min(LDBdata(:,:,i))).*pxS/10 < max(max(LDBdata(:,:,i)))
            set(gca,'CLim',[min(min(LDBdata(:,:,i))).*pxS/10, max(max(LDBdata(:,:,i))).*pxS/10])
        end
        
        freezeColors
        
    end
    axis equal
    set(gca,'FontWeight','bold','FontSize',14)
    set(gca,'XTick',[],'YTick',[])
    zlabel('Distance Downcore (cm)   ','FontWeight','bold','FontSize',16)
    zs = get(gca,'ztick');
    set(gca,'zticklabel',fliplr(zs))
    set(gca,'FontName','calibri')
    
    polypsXshort(polypsXshort==0)=NaN;
    polypsYshort(polypsYshort==0)=NaN;
    
    for i=1:5:length(polypsXshort(:,1)) % only draw 20% of the polyp tracks
        plot3(polypsYshort(i,find(isnan(polypsYshort(i,:))~=1)).*hpxS/10,polypsXshort(i,find(isnan(polypsXshort(i,:))~=1)).*hpxS/10,find(isnan(polypsXshort(i,:))~=1).*pxS/10,'k')
        hold on
    end
    grid off

end



function resultsPlot
        
    global blockLDBstore polypIDs densityTotal calcifTotal travelTotal
    global polypsPerBand

    % Results
    
    % How many polyps were used in measuring calcification in each layer
    polypsPerBand = zeros(length(blockLDBstore(1,1,:))-1,1);
    for i = 1:length(blockLDBstore(1,1,:))-1
        polypsPerBand(i) = length(find(polypIDs(i,:))); % non-zero entries are the polyp IDs, so just count them
    end
    
    
    densityTotal(travelTotal==0) = NaN;
    calcifTotal(travelTotal==0) = NaN;
    travelTotal(travelTotal==0) = NaN;
    
    densCalc = corrcoef(densityTotal(find(isnan(travelTotal)~=1)),calcifTotal(find(isnan(travelTotal)~=1)));
    extCalc = corrcoef(travelTotal(find(isnan(travelTotal)~=1)),calcifTotal(find(isnan(travelTotal)~=1)));
    
    if length(densCalc)>1
        densCalc = densCalc(1,2)^2;
        extCalc = extCalc(1,2)^2;
    end
    
    fSum = figure;
    set(fSum,'Units','inches');
    set(fSum, 'Position', [1 1 8 12]);
    
    subplot(321)
    plot(calcifTotal/10,'LineWidth',2)
    xlabel('band number (from bottom)','FontWeight','bold','FontSize',12)
    ylabel('Calcification (g cm^-^2 yr^-^1)  ','FontWeight','bold','FontSize',12)
    
    subplot(322)
    plot(travelTotal/10,'LineWidth',2)
    xlabel('band number (from bottom)','FontWeight','bold','FontSize',12)
    ylabel('extension rate (cm yr^-^1) ','FontWeight','bold','FontSize',12)
    
    subplot(323)
    plot(densityTotal,'LineWidth',2)
    xlabel('band number (from bottom)','FontWeight','bold','FontSize',12)
    ylabel('density (g cm^-^3)  ','FontWeight','bold','FontSize',12)
    
    subplot(324)
    bar(polypsPerBand)
    xlabel('band # (between LDBs)','FontWeight','bold','FontSize',12)
    ylabel('# of corallites','FontWeight','bold','FontSize',12)
    
    subplot(325)
    scatter(travelTotal(travelTotal>0)/10,calcifTotal(travelTotal>0)/10,'.')
    xlabel('extension rate (cm yr^-^1)  ','FontWeight','bold','FontSize',12)
    ylabel('Calcification (g cm^-^2 yr^-^1)  ','FontWeight','bold','FontSize',12)
    title(strcat('r^2 = ',num2str(extCalc)),'FontWeight','bold','FontSize',12)
    
    subplot(326)
    scatter(densityTotal(travelTotal>0),calcifTotal(travelTotal>0)/10,'.')
    xlabel('density (g cm^-^3)  ','FontWeight','bold','FontSize',12)
    ylabel('Calcification (g cm^-^2 yr^-^1)  ','FontWeight','bold','FontSize',12)
    title(strcat('r^2 = ',num2str(densCalc)),'FontWeight','bold','FontSize',12)

end



function resultsPlotDensity
        
    global densityTotal

    figure;

    plot(densityTotal,'LineWidth',2)
    xlabel('band number (from bottom)','FontWeight','bold','FontSize',16)
    ylabel('density (g cm^-^3)  ','FontWeight','bold','FontSize',16)

end



function processBoringsID
    
    global layers borAllow volume densityWholeCore borPer hpxS boringCount
    global name processTotal processNumber iLoop fileOpen borings SA pxS
    
    processTotal = 3;
    processNumber = 0;
    progressUpdate

    % main loop

    fid = fopen('coralSummary.csv','w');
    fprintf(fid,'%s\n',['name, ','mean HU, ','boring percent, ',...
    'volume, ','Surface area (cm^2),','layers, ','horizontal pixel spacing (mm), ','vertical pixel spacing (mm), ']);


    warning('off','all'); % suppress warnings that large 3D matrices not saved


    for iLoop = dirSet(1):dirSet(2)
        
        processNumber = 0;
        
        directory = coralDirectory;
    
        fileOpen = strcat(directory(iLoop));
    
        loadData
    
        load0 = load(fullfile(fileOpen{:},saveName),'borAllow');
        borAllow = load0.borAllow;
        
        lastL = fileOpen{:};
        if lastL(end) ~= '/'
            fileOpen{:} = [fileOpen{:},'/'];
        end
        name = [];
        slashCount = 0;
        a = fileOpen{:};
        for iN = length(fileOpen{:}):-1:1
            if a(iN) == '/'
                slashCount = slashCount+1;
                if slashCount == 2
                    name = a(iN+1:end-1);
                    break
                end
            end
        end
        
        chooseCoreFilter
    
        buildCoreBorings4loop
        
        filterUserBorings
    
        coralLoopBoringPlot
        
        if ~densEqCheck
            coralStandardCurve
        end
    
        volumeDensity
    
        fprintf(fid, '%s\n',  strcat(directory{iLoop},', ',num2str(densityWholeCore),...
            ', ',num2str(borPer),', ',num2str(volume),', ',num2str(SA),', ',num2str(layers),', ',num2str(hpxS),', ',num2str(pxS)));
        
        % save output in coral folder
        fid2 = fopen(fullfile(fileOpen{:},'bioerosion_output.csv'),'w');
        fprintf(fid2,'%s\n',['name, ','mean HU, ','boring percent, ',...
        'volume, ','Surface area (cm^2),','layers, ','horizontal pixel spacing (mm), ','vertical pixel spacing (mm), ']);
        fprintf(fid2, '%s\n',  strcat(directory{iLoop},', ',num2str(densityWholeCore),...
            ', ',num2str(borPer),', ',num2str(volume),', ',num2str(SA),', ',num2str(layers),', ',num2str(hpxS),', ',num2str(pxS)));
    
        close all  

        coralSave
    
    end
end



function processBoringsAll
    
    global layers fileOpen iLoop volume densityWholeCore borPer hpxS boringCount
    global processNumber processTotal name SA pxS
  
    processTotal = 3;
    processNumber = 0;
    progressUpdate
    
    % main loop

    fid = fopen('coralSummary.csv','w');
    fprintf(fid,'%s\n',['name, ','mean HU, ','boring percent, ',...
    'volume, ','Surface area (cm^2),','layers, ','horizontal pixel spacing (mm), ','vertical pixel spacing (mm), ']);


    warning('off','all'); % suppress warnings that large 3D matrices not saved


    for iLoop = dirSet(1):dirSet(2)
        
        processNumber = 0;
        
        directory = coralDirectory;
    
        fileOpen = strcat(directory(iLoop));
    
        loadData
        
        lastL = fileOpen{:};
        if lastL(end) ~= '/'
            fileOpen{:} = [fileOpen{:},'/'];
        end
        name = [];
        slashCount = 0;
        a = fileOpen{:};
        for iN = length(fileOpen{:}):-1:1
            if a(iN) == '/'
                slashCount = slashCount+1;
                if slashCount == 2
                    name = a(iN+1:end-1);
                    break
                end
            end
        end
        
        chooseCoreFilter
    
        buildCoreBorings4loop
        
        allBorings
    
        coralLoopBoringPlot
        
        if ~densEqCheck
            coralStandardCurve
        end
    
        volumeDensity
    
        fprintf(fid, '%s\n',  strcat(directory{iLoop},', ',num2str(densityWholeCore),...
            ', ',num2str(borPer),', ',num2str(volume),', ',num2str(SA),', ',num2str(layers),', ',num2str(hpxS),', ',num2str(pxS)));
        
        % save output in coral folder
        fid2 = fopen(fullfile(fileOpen{:},'bioerosion_output.csv'),'w');
        fprintf(fid2,'%s\n',['name, ','mean HU, ','boring percent, ',...
        'volume, ','Surface area (cm^2),','layers, ','horizontal pixel spacing (mm), ','vertical pixel spacing (mm), ']);
        fprintf(fid2, '%s\n',  strcat(directory{iLoop},', ',num2str(densityWholeCore),...
            ', ',num2str(borPer),', ',num2str(volume),', ',num2str(SA),', ',num2str(layers),', ',num2str(hpxS),', ',num2str(pxS)));
    
        close all  

        coralSave
    
    end
end



function buildCoreBorings4loop
        
    global row col layers fileOpen coreMovieName X h2
    global borings borBottom corA borA coral processNumber titleName
        
    processNumber = processNumber + 1;
    titleName = 'Mapping Borings';
        
        % locate the core and check for cracks and borings

    if inpMov == 2
        writerObj = VideoWriter(strcat(fileOpen,coreMovieName{:})); % movie object
        writerObj.FrameRate = 10;
        open(writerObj);
        f1 = figure;
        set(f1,'Units','inches');
        set(f1, 'Position', [1 1 12 5]);
    elseif inpMov ~= 1
    end

    bottomCore = 1;
    topCore = layers;

    filteredXcore = zeros(row,col); 

    % initialize matrix to store indices of where coral, cracks, borings exists
    coral = zeros(row,col,topCore);
    cracks = zeros(row,col,topCore);
    borings = zeros(row,col,topCore);

    % keep track of size of coral, cylinder, and borings
    ellA = zeros(length(bottomCore:topCore),1);
    corA = zeros(length(bottomCore:topCore),1);
    borA = zeros(length(bottomCore:topCore),1);

    % store center location of core
    center = zeros(topCore-bottomCore,2); % to store center of core in each layer
    leftMost = zeros(topCore-bottomCore,2);
    rightMost = zeros(topCore-bottomCore,2);
    topMost = zeros(topCore-bottomCore,2);
    bottomMost = zeros(topCore-bottomCore,2);

    % initialize crack variables
    crackA = zeros(row,col);
    crackCheck = 0;
    crackExist = 0;
    crackLayers = zeros(1,layers);
    crackStart = 1;

    borV = 0; % keep track of total volume of borings

    % flag for identifying once the full core has been reached - no borings
    % allowed until this occurs
    borBottom = zeros(size(borings));
    borLayer = zeros(row,col);

    for i = bottomCore:topCore
        
        % filter DICOM image
        filteredXcore(:,:) = imfilter(X(:,:,i), h2, 'replicate');
    
        % Use thresholding to identify core region. If 'thresh' == 0, uses
        % Otsu's method, otherwise user can set defined threshold value (in HU)
        level = 0;
        if thresh == 0
            level = graythresh(((filteredXcore-min(min(filteredXcore)))./max(max(filteredXcore-min(min(filteredXcore))))));%.*max(max(filteredXcore-min(min(filteredXcore))))+min(min(filteredXcore));%graythresh(filteredXcore);
            if (max(max(filteredXcore))-min(min(filteredXcore)))*level+min(min(filteredXcore)) > -800
                coral(:,:,i) = im2bw(((filteredXcore-min(min(filteredXcore)))./max(max(filteredXcore-min(min(filteredXcore))))),level);
            end
        else
            coral(:,:,i) = ind2sub(size(filteredXcore(:,:)),filteredXcore(:,:) > coreTol);
        end
    
        borLayer((coral(:,:,i)+borLayer)>0) = 1;

        borBottom(:,:,i) = borLayer;
    
        if (max(max(filteredXcore))-min(min(filteredXcore)))*level+min(min(filteredXcore)) > -800 % coral exists in this layer
   
            [r,c] = find(coral(:,:,i));
    
            [val,loc] = max(r);
            topMost(i,1:2) = [c(loc),r(loc)];
            [val,loc] = min(r);
            bottomMost(i,1:2) = [c(loc),r(loc)];
            [val,loc] = max(c);
            rightMost(i,1:2) = [c(loc),r(loc)];
            [val,loc] = min(c);
            leftMost(i,1:2) = [c(loc),r(loc)];
    
            % calculate the upperright, upperleft, lowerright, lowerleft
            [val,loc1] = max(r.*c); % note that 'val' does not matter for our purposes
            [val,loc2] = max((row-r).*(col-c));
            [val,loc3] = max((row-r).*c);
            [val,loc4] = max(r.*(col-c));
    
            % define the center. This currently is not used in this code, but could
            % consider using this center instead of 'z'
            center(i,:) = [round((rightMost(i,1)+leftMost(i,1))/2),round((topMost(i,2)+bottomMost(i,2))/2)];
    
            try 
            [z, a, b, alpha] = fitellipse([leftMost(i,1:2)+1;topMost(i,1:2)-1;rightMost(i,1:2)-1;bottomMost(i,1:2)+1;...
            c(loc1)-1,r(loc1)-1;c(loc2)+1,r(loc2)+1;c(loc3)-1,r(loc3)+1;c(loc4)+1,r(loc4)-1]);
    
            ellA(i) = pi*a*b; % calculate area of ellipse
            corA(i) = sum(sum(coral(:,:,i))); % area of coral    
        
            [ro,co] = find(coral(:,:,i)-1); % find where coral does not exist
        
            % determine if voxels are inside ellipse by translating and
            % rotating cords and then aligning with ellipse
            inEll = ((ro-z(2))*cos(alpha)+(co-z(1))*sin(alpha)).^2./a.^2 + ...
                (-(ro-z(2))*sin(alpha)+(co-z(1))*cos(alpha)).^2./b.^2;
        
            borInd = find(inEll<1);
        
            % assign to borings matrix
            for j = 1:length(borInd)
                borings(ro(borInd(j)),co(borInd(j)),i) = 1;
            end
        
            indBor = bwlabel(borings(:,:,i));
        
            for j = 1:max(max(indBor))
                if (sum(sum(indBor==j))/ellA(i))*100 < borTol
                    indBor(indBor==j) = 0;
                end
            end
        
        
            % reassign borings
        
            borings(:,:,i) = 0;
        
            [ro,co] = find(indBor>0);
        
            % assign to borings matrix
            for j = 1:length(ro)
                borings(ro(j),co(j),i) = 1;
            end
        
            borA(i) = sum(sum(borings(:,:,i))); % determine area of borings
            borV = borV + borA(i); % add on to total borings volume
        
            % if not borings, reset borLayer
            if sum(sum(borings(:,:,i))) == 0
                borLayer = coral(:,:,i);
            end
    
        if (ellA(i)-corA(i))/ellA(i)*100 > crackTol % is percent crack > tolerance% ?
            if crackCheck == 0 % not previous in a crack
                crackCheck = 2; % start potential crack
                crackStart = i; % note which layer the crack started in
                crackA = coral(:,:,crackStart-1)-coral(:,:,i); % define crack area
                cracks(:,:,i) = coral(:,:,crackStart-1)-coral(:,:,i); % store crack
            else % previousy in a crack
                crackA = crackA+coral(:,:,crackStart-1)-coral(:,:,i); % add on to crack area
                crackA(crackA>0) = 1; % reset crack area to binary
                crackA(crackA<0) = 0;
                cracks(:,:,i) = coral(:,:,crackStart-1)-coral(:,:,i); % store crack
            end
        elseif crackCheck > 0 % potentially exited a crack
            crackCheck = crackCheck - 1; % allow one layer to miss tolerance before reseting
            if crackCheck == 0 % exited a crack
                crackA = zeros(row,col); % reset crack area
                if crackExist == 0 % if it was not a crack
                    cracks(:,:,crackStart:i) = 0; % delete the potential crack
                end
                crackExist = 0; % reset if crack was exited
            end
        end
    
        if sum(crackA(coral(:,:,crackStart)==1)) > 0.5*sum(sum(sum(coral(:,:,i))))
            crackExist = 1; % crack found
            crackLayers(crackStart-1:i+1) = 1; % note which layers are crack layers
        end
    
        % plot
        if inpMov == 1 || inpMov == 2
            subplot(1,2,1)
            imagesc(coral(:,:,i));
            colorbar
            hold on
            plotellipse(z,a,b,alpha);
            text(10,10,strcat('layer  ', num2str(i), ' of  ', num2str(layers)),'Color','white')
            scatter([center(i,1),leftMost(i,1),topMost(i,1),rightMost(i,1),bottomMost(i,1),...
                c(loc1(i)),c(loc2(i)),c(loc3(i)),c(loc4(i))],...
                [center(i,2),leftMost(i,2),topMost(i,2),rightMost(i,2),bottomMost(i,2),...
                r(loc1(i)),r(loc2(i)),r(loc3(i)),r(loc4(i))])
            hold off
            subplot(1,2,2)
            imagesc(X(:,:,i))
            colormap bone
            colorbar
            set(gca,'CLim',[-500 1800])
            drawnow
            if inpMov == 2
                frame = getframe;
                writeVideo(writerObj,frame)
            end
        end
    
            catch
            end
    
        end
        
        progressUpdate((i-bottomCore+1)/((topCore-bottomCore)))
        
    end
    %close
end



function filterUserBorings
        
    global borings borBottom borAllow boringCount borA borT corT borPer layers corA

    borings(borBottom==0) = 0;

    L = bwlabeln(borings); % separate borings

    userAllow = [];
    [ro co] = find(borAllow);
    for i = 1:length(ro)
        borUser = L(ro(i),co(i),borAllow(ro(i),co(i)));
        userAllow = [userAllow borUser];
    end

    userAllow = unique(userAllow);
    userAllow(userAllow==0) = [];

    boringCount = length(userAllow);

    borings = zeros(size(borings));

    for i = 1:length(userAllow)
        borings(L==userAllow(i)) = 1;
    end

    for i = 1:layers
        borA(i) = sum(sum(borings(:,:,i)));
    end


    borT = sum(sum(sum(borings)));
    corT = sum(corA);

    borPer = 100*(borT/(borT+corT));

end



function allBorings

    global layers borings borBottom boringCount borA borT corT borPer corA

    borings(borBottom==0) = 0;

    L = bwlabeln(borings); % separate borings
    boringCount = max(L);

    for i = 1:layers
        borA(i) = sum(sum(borings(:,:,i)));
    end

    borT = sum(sum(sum(borings)));
    corT = sum(corA);

    borPer = 100*(borT/(borT+corT));

end



function selectBorings
    
    global X layers row col h2 borAllow fileOpen processTotal processNumber iBand f

    processTotal = 1;
    processNumber = 0;
    progressUpdate
    
    % User defined borings locations
    
    for iLoop = dirSet(1):dirSet(2)
        
        processNumber = 0;
        
        fileOpen = strcat(directory(iLoop)); % read the path name
        
        % make sure the filename ends with a '/'
        lastL = fileOpen{:};
        if lastL(end) ~= '/'
            fileOpen{:} = [fileOpen{:},'/'];
        end
        
        loadData
    
        chooseCoreFilter
        
        i = round(layers/3);
        % filter DICOM image
        filteredXcore(:,:) = imfilter(X(:,:,i), h2, 'replicate');
    
        % Use thresholding to identify core region. If 'thresh' == 0, uses
        % Otsu's method, otherwise user can set defined threshold value (in HU)
        level = 0;
        if thresh == 0
            level = graythresh(((filteredXcore-min(min(filteredXcore)))./max(max(filteredXcore-min(min(filteredXcore))))));%.*max(max(filteredXcore-min(min(filteredXcore))))+min(min(filteredXcore));%graythresh(filteredXcore);
            if (max(max(filteredXcore))-min(min(filteredXcore)))*level+min(min(filteredXcore)) > -800
                coral(:,:,i) = im2bw(((filteredXcore-min(min(filteredXcore)))./max(max(filteredXcore-min(min(filteredXcore))))),level);
            end
        else
            coral(:,:,i) = ind2sub(size(filteredXcore(:,:)),filteredXcore(:,:) > coreTol);
        end
        coral(:,1,i) = 0;
        coral(:,row,i) = 0;
        coral(1,:,i) = 0;
        coral(col,:,i) = 0;
   
        [r,c] = find(coral(:,:,i));
    
        [val,loc] = max(r);
        topMost = [c(loc),r(loc)];
        [val,loc] = min(r);
        bottomMost = [c(loc),r(loc)];
        [val,loc] = max(c);
        rightMost = [c(loc),r(loc)];
        [val,loc] = min(c);
        leftMost = [c(loc),r(loc)];
    
        center = [round((rightMost(1)+leftMost(1))/2),round((topMost(2)+bottomMost(2))/2)];
    
        % set up slab positions
        spacing = linspace((numSlabs+1)/2-numSlabs,(numSlabs-1)/2,numSlabs)...
        .*round((topMost(2)-bottomMost(2))/(numSlabs-1));
        slab1 = round(center(2)+spacing);
        slab2 = round(center(2));

        borAllow = zeros(row,col);

        for i = 1:numSlabs+1
    
    
    
            if i < numSlabs+1
        
                if slab1(i)+12 > row
                    slab1(i) = row-12;
                end
                if slab1(i)-12 < 1
                    slab1(i) = 13;
                end
                if slab2+12 > row
                    slab2 = row-12;
                end
                if slab2-12 < 1
                    slab2 = 13;
                end
        
                slabDraw3 = zeros(row,1,layers);
                slabDraw3(:,:,1:layers)  = min(X(:,slab1(i)-12:slab1(i)+12,:),[],2);
                slabDraw3 = permute(slabDraw3,[3,1,2]);
            else % the last slab is perpendicular to the others
                slabDraw3 = zeros(row,1,layers);
                slabDraw3(:,:,1:layers)  = min(X(slab2-12:slab2+12,:,:),[],1);
                slabDraw3 = permute(slabDraw3,[3,1,2]);
            end

            f(iBand) = figure;
            imagesc(flipud(slabDraw3))
            colormap(bone)
            set(f(iBand),'Units','inches');
            set(f(iBand), 'Position', [1 1 6 10]);
            axis equal
            set(gca,'CLim',[-500 2000])
    
            button = 1; % reset button
            j = 0; % reset band number
        
            uistack(f(iBand),'top')
            figure(f(iBand))
        
            [x,y,b] = ginput;
            
            % delete clicks off image
            inds = find(round(x)<1 | round(x)>row);
            y(inds) = [];
            x(inds) = [];
            inds2 = find((layers-round(y))<1 | round(y)<1);
            y(inds2) = [];
            x(inds2) = [];

            for k = 1:length(x)
                if i < numSlabs+1 % for first set of slabs (all but last one)
                    borAllow(round(x(k)),slab1(i)-12:slab1(i)+12) = layers-round(y(k));
                else % reverse x and y assignment order
                    borAllow(slab2-12:slab2+12,round(x(k))) = layers-round(y(k));
                end
            end
        
        end

        close all
        drawnow
    
        save(fullfile(fileOpen{:},saveName), 'borAllow');
        
    end
        
end



function chooseCoreFilter
        
    global h2 hpxS
        
    % filter for core identification
    % first number is size (in voxels), second number is standard deviation
    % (in voxels). Enter 'hpxS' in workspace after loading a core to
    % see the length of one voxel (in mm). Enter 'surf(h2)' in workspace to
    % visualize what the filter looks like.
    if strcmp(gen,'Porites')
        h2 = fspecial('gaussian',round(30/hpxS*0.1), 10);
    elseif strcmp(gen,'Siderastrea')
        h2 = fspecial('gaussian',round(30/hpxS*0.1), 10);
    elseif strcmp(gen,'Acropora')
        h2 = fspecial('gaussian',round(7/hpxS*0.1), 10);
    elseif strcmp(gen,'Montastrea')
        h2 = fspecial('gaussian',round(12/hpxS*0.1), 5);
        thresh = 1;
        coreTol = -300;
    elseif strcmp(gen,'Diploastrea')
        h2 = fspecial('gaussian',round(20/hpxS*0.1), 10);
        thresh = 1;
        coreTol = -400;
    elseif strcmp(gen,'Diploria')
        h2 = fspecial('gaussian',round(15/hpxS*0.1), 7);
        thresh = 1;
        coreTol = -400;
    elseif strcmp(gen,'Favia')
        h2 = fspecial('gaussian',round(40/hpxS*0.1), 10);
        thresh = 1;
        coreTol = -400;
    end

end



function extractName
    
    global fileOpen name
    
    lastL = fileOpen{:};
    if lastL(end) ~= '/'
        fileOpen{:} = [fileOpen{:},'/'];
    end
    name = [];
    slashCount = 0;
    a = fileOpen{:};
    for iN = length(fileOpen{:}):-1:1
        if a(iN) == '/'
            slashCount = slashCount+1;
            if slashCount == 2
                name = a(iN+1:end-1);
                break
            end
        end
    end
end



function [z, a, b, alpha] = fitellipse(x, varargin)
    %FITELLIPSE   least squares fit of ellipse to 2D data
    %
    %   [Z, A, B, ALPHA] = FITELLIPSE(X)
    %       Fit an ellipse to the 2D points in the 2xN array X. The ellipse is
    %       returned in parametric form such that the equation of the ellipse
    %       parameterised by 0 <= theta < 2*pi is:
    %           X = Z + Q(ALPHA) * [A * cos(theta); B * sin(theta)]
    %       where Q(ALPHA) is the rotation matrix
    %           Q(ALPHA) = [cos(ALPHA), -sin(ALPHA);
    %                       sin(ALPHA), cos(ALPHA)]
    %
    %       Fitting is performed by nonlinear least squares, optimising the
    %       squared sum of orthogonal distances from the points to the fitted
    %       ellipse. The initial guess is calculated by a linear least squares
    %       routine, by default using the Bookstein constraint (see below)
    %
    %   [...]            = FITELLIPSE(X, 'linear')
    %       Fit an ellipse using linear least squares. The conic to be fitted
    %       is of the form
    %           x'Ax + b'x + c = 0
    %       and the algebraic error is minimised by least squares with the
    %       Bookstein constraint (lambda_1^2 + lambda_2^2 = 1, where
    %       lambda_i are the eigenvalues of A)
    %
    %   [...]            = FITELLIPSE(..., 'Property', 'value', ...)
    %       Specify property/value pairs to change problem parameters
    %          Property                  Values
    %          =================================
    %          'constraint'              {|'bookstein'|, 'trace'}
    %                                    For the linear fit, the following
    %                                    quadratic form is considered
    %                                    x'Ax + b'x + c = 0. Different
    %                                    constraints on the parameters yield
    %                                    different fits. Both 'bookstein' and
    %                                    'trace' are Euclidean-invariant
    %                                    constraints on the eigenvalues of A,
    %                                    meaning the fit will be invariant
    %                                    under Euclidean transformations
    %                                    'bookstein': lambda1^2 + lambda2^2 = 1
    %                                    'trace'    : lambda1 + lambda2     = 1
    %
    %           Nonlinear Fit Property   Values
    %           ===============================
    %           'maxits'                 positive integer, default 200
    %                                    Maximum number of iterations for the
    %                                    Gauss Newton step
    %
    %           'tol'                    positive real, default 1e-5
    %                                    Relative step size tolerance
    %   Example:
    %       % A set of points
    %       x = [1 2 5 7 9 6 3 8;
    %            7 6 8 7 5 7 2 4];
    %
    %       % Fit an ellipse using the Bookstein constraint
    %       [zb, ab, bb, alphab] = fitellipse(x, 'linear');
    %
    %       % Find the least squares geometric estimate
    %       [zg, ag, bg, alphag] = fitellipse(x);
    %
    %       % Plot the results
    %       plot(x(1,:), x(2,:), 'ro')
    %       hold on
    %       % plotellipse(zb, ab, bb, alphab, 'b--')
    %       % plotellipse(zg, ag, bg, alphag, 'k')
    %
    %   See also PLOTELLIPSE
    
    % Copyright Richard Brown, this code can be freely used and modified so
    % long as this line is retained
    error(nargchk(1, 5, nargin, 'struct'))
    
    % Default parameters
    params.fNonlinear = true;
    params.constraint = 'bookstein';
    params.maxits     = 200;
    params.tol        = 1e-5;
    
    % Parse inputs
    [x, params] = parseinputs(x, params, varargin{:});
    
    % Constraints are Euclidean-invariant, so improve conditioning by removing
    % centroid
    centroid = mean(x, 2);
    x        = x - repmat(centroid, 1, size(x, 2));
    
    % Obtain a linear estimate
    switch params.constraint
        % Bookstein constraint : lambda_1^2 + lambda_2^2 = 1
        case 'bookstein'
            [z, a, b, alpha] = fitbookstein(x);
            
            % 'trace' constraint, lambda1 + lambda2 = trace(A) = 1
        case 'trace'
            [z, a, b, alpha] = fitggk(x);
    end % switch
    
    % Minimise geometric error using nonlinear least squares if required
    if params.fNonlinear
        % Initial conditions
        z0     = z;
        a0     = a;
        b0     = b;
        alpha0 = alpha;
        
        % Apply the fit
        [z, a, b, alpha, fConverged] = ...
            fitnonlinear(x, z0, a0, b0, alpha0, params);
        
        % Return linear estimate if GN doesn't converge
        if ~fConverged
            warning('fitellipse:FailureToConverge', ...'
                'Gauss-Newton did not converge, returning linear estimate');
            z = z0;
            a = a0;
            b = b0;
            alpha = alpha0;
        end
    end
    
    % Add the centroid back on
    z = z + centroid;

end % fitellipse



function [z, a, b, alpha] = fitbookstein(x)
    %FITBOOKSTEIN   Linear ellipse fit using bookstein constraint
    %   lambda_1^2 + lambda_2^2 = 1, where lambda_i are the eigenvalues of A
    
    % Convenience variables
    m  = size(x, 2);
    x1 = x(1, :)';
    x2 = x(2, :)';
    
    % Define the coefficient matrix B, such that we solve the system
    % B *[v; w] = 0, with the constraint norm(w) == 1
    B = [x1, x2, ones(m, 1), x1.^2, sqrt(2) * x1 .* x2, x2.^2];
    
    % To enforce the constraint, we need to take the QR decomposition
    [Q, R] = qr(B);
    
    % Decompose R into blocks
    R11 = R(1:3, 1:3);
    R12 = R(1:3, 4:6);
    R22 = R(4:6, 4:6);
    
    % Solve R22 * w = 0 subject to norm(w) == 1
    [U, S, V] = svd(R22);
    w = V(:, 3);
    
    % Solve for the remaining variables
    v = -R11 \ R12 * w;
    
    % Fill in the quadratic form
    A        = zeros(2);
    A(1)     = w(1);
    A([2 3]) = 1 / sqrt(2) * w(2);
    A(4)     = w(3);
    bv       = v(1:2);
    c        = v(3);
    
    % Find the parameters
    [z, a, b, alpha] = conic2parametric(A, bv, c);

end % fitellipse



function [z, a, b, alpha] = fitggk(x)
    % Linear least squares with the Euclidean-invariant constraint Trace(A) = 1
    
    % Convenience variables
    m  = size(x, 2);
    x1 = x(1, :)';
    x2 = x(2, :)';
    
    % Coefficient matrix
    B = [2 * x1 .* x2, x2.^2 - x1.^2, x1, x2, ones(m, 1)];
    
    v = B \ -x1.^2;
    
    % For clarity, fill in the quadratic form variables
    A        = zeros(2);
    A(1,1)   = 1 - v(2);
    A([2 3]) = v(1);
    A(2,2)   = v(2);
    bv       = v(3:4);
    c        = v(5);
    
    % find parameters
    [z, a, b, alpha] = conic2parametric(A, bv, c);

end



function [z, a, b, alpha, fConverged] = fitnonlinear(x, z0, a0, b0, alpha0, params)
    % Gauss-Newton least squares ellipse fit minimising geometric distance
    
    % Get initial rotation matrix
    Q0 = [cos(alpha0), -sin(alpha0); sin(alpha0) cos(alpha0)];
    m = size(x, 2);
    
    % Get initial phase estimates
    phi0 = angle( [1 i] * Q0' * (x - repmat(z0, 1, m)) )';
    u = [phi0; alpha0; a0; b0; z0];
    
    % Iterate using Gauss Newton
    fConverged = false;
    for nIts = 1:params.maxits
        % Find the function and Jacobian
        [f, J] = sys(u);
        
        % Solve for the step and update u
        h = -J \ f;
        u = u + h;
        
        % Check for convergence
        delta = norm(h, inf) / norm(u, inf);
        if delta < params.tol
            fConverged = true;
            break
        end
    end
    
    alpha = u(end-4);
    a     = u(end-3);
    b     = u(end-2);
    z     = u(end-1:end);
        

    function [f, J] = sys(u)
        % SYS : Define the system of nonlinear equations and Jacobian. Nested
        % function accesses X (but changeth it not)
        % from the FITELLIPSE workspace

        % Tolerance for whether it is a circle
        circTol = 1e-5;
        
        % Unpack parameters from u
        phi   = u(1:end-5);
        alpha = u(end-4);
        a     = u(end-3);
        b     = u(end-2);
        z     = u(end-1:end);
        
        % If it is a circle, the Jacobian will be singular, and the
        % Gauss-Newton step won't work. 
        %TODO: This can be fixed by switching to a Levenberg-Marquardt
        %solver
        if abs(a - b) / (a + b) < circTol
            warning('fitellipse:CircleFound', ...
                'Ellipse is near-circular - nonlinear fit may not succeed')
        end
        
        % Convenience trig variables
        c = cos(phi);
        s = sin(phi);
        ca = cos(alpha);
        sa = sin(alpha);
        
        % Rotation matrices
        Q    = [ca, -sa; sa, ca];
        Qdot = [-sa, -ca; ca, -sa];
        
        % Preallocate function and Jacobian variables
        f = zeros(2 * m, 1);
        J = zeros(2 * m, m + 5);
        for i = 1:m
            rows = (2*i-1):(2*i);
            % Equation system - vector difference between point on ellipse
            % and data point
            f((2*i-1):(2*i)) = x(:, i) - z - Q * [a * cos(phi(i)); b * sin(phi(i))];
            
            % Jacobian
            J(rows, i) = -Q * [-a * s(i); b * c(i)];
            J(rows, (end-4:end)) = ...
                [-Qdot*[a*c(i); b*s(i)], -Q*[c(i); 0], -Q*[0; s(i)], [-1 0; 0 -1]];
        end
    end
end % fitnonlinear



function [z, a, b, alpha] = conic2parametric(A, bv, c)
    % Diagonalise A - find Q, D such at A = Q' * D * Q
    [Q, D] = eig(A);
    Q = Q';
    
    % If the determinant < 0, it's not an ellipse
    if prod(diag(D)) <= 0
        error('fitellipse:NotEllipse', 'Linear fit did not produce an ellipse');
    end
    
    % We have b_h' = 2 * t' * A + b'
    t = -0.5 * (A \ bv);
    
    c_h = t' * A * t + bv' * t + c;
    
    z = t;
    a = sqrt(-c_h / D(1,1));
    b = sqrt(-c_h / D(2,2));
    alpha = atan2(Q(1,2), Q(1,1));
end % conic2parametric



function [x, params] = parseinputs(x, params, varargin)
    % PARSEINPUTS put x in the correct form, and parse user parameters
    
    % CHECK x
    % Make sure x is 2xN where N > 3
    if size(x, 2) == 2
        x = x';
    end
    if size(x, 1) ~= 2
        error('fitellipse:InvalidDimension', ...
            'Input matrix must be two dimensional')
    end
    if size(x, 2) < 6
        error('fitellipse:InsufficientPoints', ...
            'At least 6 points required to compute fit')
    end
    
    
    % Determine whether we are solving for geometric (nonlinear) or algebraic
    % (linear) distance
    if ~isempty(varargin) && strncmpi(varargin{1}, 'linear', length(varargin{1}))
        params.fNonlinear = false;
        varargin(1)       = [];
    else
        params.fNonlinear = true;
    end
    
    % Parse property/value pairs
    if rem(length(varargin), 2) ~= 0
        error('fitellipse:InvalidInputArguments', ...
            'Additional arguments must take the form of Property/Value pairs')
    end
    
    % Cell array of valid property names
    properties = {'constraint', 'maxits', 'tol'};
    
    while length(varargin) ~= 0
        % Pop pair off varargin
        property      = varargin{1};
        value         = varargin{2};
        varargin(1:2) = [];
        
        % If the property has been supplied in a shortened form, lengthen it
        iProperty = find(strncmpi(property, properties, length(property)));
        if isempty(iProperty)
            error('fitellipse:UnknownProperty', 'Unknown Property');
        elseif length(iProperty) > 1
            error('fitellipse:AmbiguousProperty', ...
                'Supplied shortened property name is ambiguous');
        end
        
        % Expand property to its full name
        property = properties{iProperty};
        
        % Check for irrelevant property
        if ~params.fNonlinear && ismember(property, {'maxits', 'tol'})
            warning('fitellipse:IrrelevantProperty', ...
                'Supplied property has no effect on linear estimate, ignoring');
            continue
        end
        
        % Check supplied property value
        switch property
            case 'maxits'
                if ~isnumeric(value) || value <= 0
                    error('fitcircle:InvalidMaxits', ...
                        'maxits must be an integer greater than 0')
                end
                params.maxits = value;
            case 'tol'
                if ~isnumeric(value) || value <= 0
                    error('fitcircle:InvalidTol', ...
                        'tol must be a positive real number')
                end
                params.tol = value;
            case 'constraint'
                switch lower(value)
                    case 'bookstein'
                        params.constraint = 'bookstein';
                    case 'trace'
                        params.constraint = 'trace';
                    otherwise
                        error('fitellipse:InvalidConstraint', ...
                            'Invalid constraint specified')
                end
        end % switch property
    end % while

end % parseinputs



function varargout = plotellipse(varargin)
    %PLOTELLIPSE   Plot parametrically specified ellipse
    %
    %   PLOTELLIPSE(Z, A, B, ALPHA) Plots the ellipse specified by Z, A, B,
    %       ALPHA (as returned by FITELLIPSE)
    %
    %       A, B are positive scalars, Z a 2x1 column vector, and
    %       ALPHA a rotation angle, such that the equation of the ellipse is:
    %           X = Z + Q(ALPHA) * [A * cos(phi); B * sin(phi)]
    %       where Q(ALPHA) is the rotation matrix
    %           Q(ALPHA) = [cos(ALPHA) -sin(ALPHA);
    %                       sin(AlPHA) cos(ALPHA)]
    %
    %   PLOTELLIPSE(..., LineSpec) passes the LineSpec string to the plot
    %       command (e.g. 'r--')
    %
    %   PLOTELLIPSE(Hax, ...) plots into the axes specified by the axes handle
    %       Hax
    %
    %   H = PLOTELLIPSE(...) returns a handle to the created lineseries object
    %       created by the plot command
    %
    %   Example:
    %       % Ellipse centred at 10,10, with semiaxes 5 and 3, rotated by pi/4
    %       a = 5;
    %       b = 3;
    %       z = [10; 10]
    %       alpha = pi/4;
    %       plotellipse(z, a, b, alpha)
    %
    %   See also FITELLIPSE
    
    % Copyright Richard Brown. This code can be freely reused and modified so
    % long as it retains this copyright clause
    
    error(nargchk(4, 6, nargin, 'struct'));
    error(nargchk(0, 1, nargout, 'struct'));
    
    % Parse and check inputs
    if ishandle(varargin{1})
        hAx = varargin{1};
        varargin(1) = [];
    else
        hAx = gca();
    end
    
    % Ellipse centre
    z = varargin{1};
    z = z(:);
    if length(z) ~= 2
        error('plotellipse:InvalidCentre', ...
            'Ellipse center must be a 2 element column vector');
    end
    
    a = varargin{2};
    b = varargin{3};
    if ~isscalar(a) || ~isscalar(b) || a < 0 || b < 0
        error('plotellipse:InvalidAxes', ...
            'A, B must be real, positive scalars');
    end
    
    alpha = varargin{4};
    if ~isscalar(alpha)
        error('plotellipse:InvalidConstant', ...
            'Rotation angle alpha must be a real scalar, in radians');
    end
    varargin(1:4) = [];
    
    % See if a linespec is supplied
    if ~isempty(varargin)
        linespec = varargin{1};
    else
        linespec = '';
    end
    
    % form the parameter vector
    npts = 100;
    t = linspace(0, 2*pi, npts);
    
    % Rotation matrix
    Q = [cos(alpha), -sin(alpha); sin(alpha) cos(alpha)];
    % Ellipse points
    X = Q * [a * cos(t); b * sin(t)] + repmat(z, 1, npts);
    
    % The actual plotting one-liner
    h = plot(hAx, X(1,:), X(2,:), linespec);
    
    % Return the handle if asked for
    if nargout == 1
        varargout = {h};
    end

end



function freezeColors(varargin)
    % freezeColors  Lock colors of plot, enabling multiple colormaps per figure. (v2.3)
    %
    %   Problem: There is only one colormap per figure. This function provides
    %       an easy solution when plots using different colomaps are desired
    %       in the same figure.
    %
    %   freezeColors freezes the colors of graphics objects in the current axis so
    %       that subsequent changes to the colormap (or caxis) will not change the
    %       colors of these objects. freezeColors works on any graphics object
    %       with CData in indexed-color mode: surfaces, images, scattergroups,
    %       bargroups, patches, etc. It works by converting CData to true-color rgb
    %       based on the colormap active at the time freezeColors is called.
    %
    %   The original indexed color data is saved, and can be restored using
    %       unfreezeColors, making the plot once again subject to the colormap and
    %       caxis.
    %
    %
    %   Usage:
    %       freezeColors        applies to all objects in current axis (gca),
    %       freezeColors(axh)   same, but works on axis axh.
    %
    %   Example:
    %       subplot(2,1,1); imagesc(X); colormap hot; freezeColors
    %       subplot(2,1,2); imagesc(Y); colormap hsv; freezeColors etc...
    %
    %       Note: colorbars must also be frozen. Due to Matlab 'improvements' this can
    %				no longer be done with freezeColors. Instead, please
    %				use the function CBFREEZE by Carlos Adrian Vargas Aguilera
    %				that can be downloaded from the MATLAB File Exchange
    %				(http://www.mathworks.com/matlabcentral/fileexchange/24371)
    %
    %       h=colorbar; cbfreeze(h), or simply cbfreeze(colorbar)
    %
    %       For additional examples, see test/test_main.m
    %
    %   Side effect on render mode: freezeColors does not work with the painters
    %       renderer, because Matlab doesn't support rgb color data in
    %       painters mode. If the current renderer is painters, freezeColors
    %       changes it to zbuffer. This may have unexpected effects on other aspects
    %	      of your plots.
    %
    %       See also unfreezeColors, freezeColors_pub.html, cbfreeze.
    %
    %
    %   John Iversen (iversen@nsi.edu) 3/23/05
    %
    
    %   Changes:
    %   JRI (iversen@nsi.edu) 4/19/06   Correctly handles scaled integer cdata
    %   JRI 9/1/06   should now handle all objects with cdata: images, surfaces,
    %                scatterplots. (v 2.1)
    %   JRI 11/11/06 Preserves NaN colors. Hidden option (v 2.2, not uploaded)
    %   JRI 3/17/07  Preserve caxis after freezing--maintains colorbar scale (v 2.3)
    %   JRI 4/12/07  Check for painters mode as Matlab doesn't support rgb in it.
    %   JRI 4/9/08   Fix preserving caxis for objects within hggroups (e.g. contourf)
    %   JRI 4/7/10   Change documentation for colorbars
    
    % Hidden option for NaN colors:
    %   Missing data are often represented by NaN in the indexed color
    %   data, which renders transparently. This transparency will be preserved
    %   when freezing colors. If instead you wish such gaps to be filled with
    %   a real color, add 'nancolor',[r g b] to the end of the arguments. E.g.
    %   freezeColors('nancolor',[r g b]) or freezeColors(axh,'nancolor',[r g b]),
    %   where [r g b] is a color vector. This works on images & pcolor, but not on
    %   surfaces.
    %   Thanks to Fabiano Busdraghi and Jody Klymak for the suggestions. Bugfixes
    %   attributed in the code.
    
    % Free for all uses, but please retain the following:
    %   Original Author:
    %   John Iversen, 2005-10
    %   john_iversen@post.harvard.edu
    
    appdatacode = 'JRI__freezeColorsData';
    
    [h, nancolor] = checkArgs(varargin);
    
    %gather all children with scaled or indexed CData
    cdatah = getCDataHandles(h);
    
    %current colormap
    cmap = colormap;
    nColors = size(cmap,1);
    cax = caxis;
    
    % convert object color indexes into colormap to true-color data using
    %  current colormap
    for hh = cdatah',
        g = get(hh);
        
        %preserve parent axis clim
        parentAx = getParentAxes(hh);
        originalClim = get(parentAx, 'clim');
        
        %   Note: Special handling of patches: For some reason, setting
        %   cdata on patches created by bar() yields an error,
        %   so instead we'll set facevertexcdata instead for patches.
        if ~strcmp(g.Type,'patch'),
            cdata = g.CData;
        else
            cdata = g.FaceVertexCData;
        end
        
        %get cdata mapping (most objects (except scattergroup) have it)
        if isfield(g,'CDataMapping'),
            scalemode = g.CDataMapping;
        else
            scalemode = 'scaled';
        end
        
        %save original indexed data for use with unfreezeColors
        siz = size(cdata);
        setappdata(hh, appdatacode, {cdata scalemode});
        
        %convert cdata to indexes into colormap
        if strcmp(scalemode,'scaled'),
            %4/19/06 JRI, Accommodate scaled display of integer cdata:
            %       in MATLAB, uint * double = uint, so must coerce cdata to double
            %       Thanks to O Yamashita for pointing this need out
            idx = ceil( (double(cdata) - cax(1)) / (cax(2)-cax(1)) * nColors);
        else %direct mapping
            idx = cdata;
            %10/8/09 in case direct data is non-int (e.g. image;freezeColors)
            % (Floor mimics how matlab converts data into colormap index.)
            % Thanks to D Armyr for the catch
            idx = floor(idx);
        end
        
        %clamp to [1, nColors]
        idx(idx<1) = 1;
        idx(idx>nColors) = nColors;
        
        %handle nans in idx
        nanmask = isnan(idx);
        idx(nanmask)=1; %temporarily replace w/ a valid colormap index
        
        %make true-color data--using current colormap
        realcolor = zeros(siz);
        for i = 1:3,
            c = cmap(idx,i);
            c = reshape(c,siz);
            c(nanmask) = nancolor(i); %restore Nan (or nancolor if specified)
            realcolor(:,:,i) = c;
        end
        
        %apply new true-color color data
        
        %true-color is not supported in painters renderer, so switch out of that
        if strcmp(get(gcf,'renderer'), 'painters'),
            set(gcf,'renderer','zbuffer');
        end
        
        %replace original CData with true-color data
        if ~strcmp(g.Type,'patch'),
            set(hh,'CData',realcolor);
        else
            set(hh,'faceVertexCData',permute(realcolor,[1 3 2]))
        end
        
        %restore clim (so colorbar will show correct limits)
        if ~isempty(parentAx),
            set(parentAx,'clim',originalClim)
        end
        
    end %loop on indexed-color objects

end


% ============================================================================ %
% Local functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getCDataHandles -- get handles of all descendents with indexed CData
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hout = getCDataHandles(h)
    % getCDataHandles  Find all objects with indexed CData
    
    %recursively descend object tree, finding objects with indexed CData
    % An exception: don't include children of objects that themselves have CData:
    %   for example, scattergroups are non-standard hggroups, with CData. Changing
    %   such a group's CData automatically changes the CData of its children,
    %   (as well as the children's handles), so there's no need to act on them.
    
    error(nargchk(1,1,nargin,'struct'))
    
    hout = [];
    if isempty(h),return;end
    
    ch = get(h,'children');
    for hh = ch'
        g = get(hh);
        if isfield(g,'CData'),     %does object have CData?
            %is it indexed/scaled?
            if ~isempty(g.CData) && isnumeric(g.CData) && size(g.CData,3)==1,
                hout = [hout; hh]; %#ok<AGROW> %yes, add to list
            end
        else %no CData, see if object has any interesting children
            hout = [hout; getCDataHandles(hh)]; %#ok<AGROW>
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getParentAxes -- return handle of axes object to which a given object belongs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hAx = getParentAxes(h)
    % getParentAxes  Return enclosing axes of a given object (could be self)
    
    error(nargchk(1,1,nargin,'struct'))
    %object itself may be an axis
    if strcmp(get(h,'type'),'axes'),
        hAx = h;
        return
    end
    
    parent = get(h,'parent');
    if (strcmp(get(parent,'type'), 'axes')),
        hAx = parent;
    else
        hAx = getParentAxes(parent);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% checkArgs -- Validate input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h, nancolor] = checkArgs(args)
    % checkArgs  Validate input arguments to freezeColors
    
    nargs = length(args);
    error(nargchk(0,3,nargs,'struct'))
    
    %grab handle from first argument if we have an odd number of arguments
    if mod(nargs,2),
        h = args{1};
        if ~ishandle(h),
            error('JRI:freezeColors:checkArgs:invalidHandle',...
                'The first argument must be a valid graphics handle (to an axis)')
        end
        % 4/2010 check if object to be frozen is a colorbar
        if strcmp(get(h,'Tag'),'Colorbar'),
            if ~exist('cbfreeze.m'),
                warning('JRI:freezeColors:checkArgs:cannotFreezeColorbar',...
                    ['You seem to be attempting to freeze a colorbar. This no longer'...
                    'works. Please read the help for freezeColors for the solution.'])
            else
                cbfreeze(h);
                return
            end
        end
        args{1} = [];
        nargs = nargs-1;
    else
        h = gca;
    end
    
    %set nancolor if that option was specified
    nancolor = [nan nan nan];
    if nargs == 2,
        if strcmpi(args{end-1},'nancolor'),
            nancolor = args{end};
            if ~all(size(nancolor)==[1 3]),
                error('JRI:freezeColors:checkArgs:badColorArgument',...
                    'nancolor must be [r g b] vector');
            end
            nancolor(nancolor>1) = 1; nancolor(nancolor<0) = 0;
        else
            error('JRI:freezeColors:checkArgs:unrecognizedOption',...
                'Unrecognized option (%s). Only ''nancolor'' is valid.',args{end-1})
        end
    end

end




function [surf labels] = imSurface(img, varargin)
%IMSURFACE Surface area of a 3D binary structure
%
%   S = imSurface(IMG)
%   Estimates the surface area of the 3D binary structure represented by
%   IMG.
%
%   S = imSurface(IMG, NDIRS)
%   Specifies the number of directions used for estimating surface area.
%   NDIRS can be either 3 or 13, default is 3.
%
%   S = imSurface(..., RESOL)
%   Specifies image resolution. RESOL is a 1-by-3 row vector containing
%   resolution in the X, Y and Z direction (in that order).
%
%   S = imSurface(LBL)
%   [S L] = imSurface(LBL)
%   When LBL is a label image, returns the surface area of each label in
%   the 3D array, and eventually returns the indices of processed labels.
%
%
%   Example
%     % Create a binary image of a ball
%     [x y z] = meshgrid(1:100, 1:100, 1:100);
%     img = sqrt( (x-50.12).^2 + (y-50.23).^2 + (z-50.34).^2) < 40;
%     % compute surface area of the ball
%     S = imSurface(img)
%     S =
%         20108
%     % compare with theoretical value
%     Sth = 4*pi*40^2;
%     100 * (S - Sth) / Sth
%     ans = 
%         0.0090
%
%   See also
%     imVolume, imMeanBreadth, imSurfaceDensity, imJointSurface
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2010-07-26,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


% check image dimension
if ndims(img) ~= 3
    error('first argument should be a 3D image');
end

% in case of a label image, return a vector with a set of results
if ~islogical(img)
    % extract labels (considers 0 as background)
    labels = unique(img);
    labels(labels == 0) = [];
    
    % allocate result array
    nLabels = length(labels);
    surf = zeros(nLabels, 1);

    props = regionprops(img, 'BoundingBox');
    
    % Compute surface area of each label considered as binary image
    % The computation is performed on a subset of the image for reducing
    % memory footprint.
    for i = 1:nLabels
        label = labels(i);
        box = props(label).BoundingBox;
        % convert bounding box to image extent, in x, y and z directions
        i0 = ceil(box([2 1 3]));
        i1 = i0 + box([5 4 6]) - 1;
        % crop image of current label
        bin = img(i0(1):i1(1), i0(2):i1(2), i0(3):i1(3)) == label;
        surf(i) = imSurface(bin, varargin{:});
    end
    
    return;
end


%% Process input arguments

% in case of binary image, compute only one label...
labels = 1;

% default number of directions
ndir = 3;

% default image resolution
delta = [1 1 1];

% Process user input arguments
while ~isempty(varargin)
    var = varargin{1};
    if ~isnumeric(var)
        error('option should be numeric');
    end
    
    % option is either connectivity or resolution
    if isscalar(var)
        ndir = var;
    else
        delta = var;
    end
    varargin(1) = [];
end


%% Initialisations

% distances between a pixel and its neighbours.
d1  = delta(1);
d2  = delta(2);
d3  = delta(3);

% volume of a voxel (used for computing line densities)
vol = d1*d2*d3;


%% Main processing for 3 directions

% number of voxels
nv = sum(img(:));

% number of connected components along the 3 main directions
% (Use Graph-based formula: chi = nVertices - nEdges)
n1 = nv - sum(sum(sum(img(:,1:end-1,:) & img(:,2:end,:))));
n2 = nv - sum(sum(sum(img(1:end-1,:,:) & img(2:end,:,:))));
n3 = nv - sum(sum(sum(img(:,:,1:end-1) & img(:,:,2:end))));

if ndir == 3
    % compute surface area by averaging over the 3 main directions
    surf = 4/3 * (n1/d1 + n2/d2 + n3/d3) * vol;
    return;
end


%% Additional processing for 13 directions

% Number of connected components along diagonals contained in the three
% main planes
n4 = nv - sum(sum(sum(img(2:end,1:end-1,:)   & img(1:end-1,2:end,:))));
n5 = nv - sum(sum(sum(img(1:end-1,1:end-1,:) & img(2:end,2:end,:))));
n6 = nv - sum(sum(sum(img(:,2:end,1:end-1)   & img(:,1:end-1,2:end))));
n7 = nv - sum(sum(sum(img(:,1:end-1,1:end-1) & img(:,2:end,2:end))));
n8 = nv - sum(sum(sum(img(2:end,:,1:end-1)   & img(1:end-1,:,2:end))));
n9 = nv - sum(sum(sum(img(1:end-1,:,1:end-1) & img(2:end,:,2:end))));

%TODO: add the case of 9 directions ?

% Number of connected components along lines corresponding to diagonals of
% the unit cube
n10 = nv - sum(sum(sum(img(1:end-1,1:end-1,1:end-1) & img(2:end,2:end,2:end))));
n11 = nv - sum(sum(sum(img(2:end,1:end-1,1:end-1) & img(1:end-1,2:end,2:end))));
n12 = nv - sum(sum(sum(img(1:end-1,2:end,1:end-1) & img(2:end,1:end-1,2:end))));
n13 = nv - sum(sum(sum(img(2:end,2:end,1:end-1) & img(1:end-1,1:end-1,2:end))));

% space between 2 voxels in each direction
d12  = hypot(d1, d2);
d13  = hypot(d1, d3);
d23  = hypot(d2, d3);
d123 = sqrt(d1^2 + d2^2 + d3^2);

% Compute weights corresponding to surface fraction of spherical caps
% For isotropic case, weights correspond to:
% c1 = 0.04577789120476 * 2;  % Ox
% c2 = 0.04577789120476 * 2;  % Oy
% c3 = 0.04577789120476 * 2;  % Oz
% c4 = 0.03698062787608 * 2;  % Oxy
% c6 = 0.03698062787608 * 2;  % Oxz
% c8 = 0.03698062787608 * 2;  % Oyz
% c10 = 0.03519563978232 * 2;  % Oxyz
c = computeDirectionWeights3d13(delta);

% compute the weighted sum of each direction
% intersection count * direction weight / line density
surf = 4 * vol * (...
    n1*c(1)/d1 + n2*c(2)/d2 + n3*c(3)/d3 + ...
    (n4+n5)*c(4)/d12 + (n6+n7)*c(6)/d13 + (n8+n9)*c(8)/d23 + ...
    (n10 + n11 + n12 + n13)*c(10)/d123 );

end


function res = computeDirectionWeights3d13(varargin)
%COMPUTEDIRECTIONWEIGHTS3D13 Direction weights for 13 directions in 3D
%
%   C = computeDirectionWeights3d13
%   Returns an array of 13-by-1 values, corresponding to directions:
%   C(1)  = [+1  0  0]
%   C(2)  = [ 0 +1  0]
%   C(3)  = [ 0  0 +1]
%   C(4)  = [+1 +1  0]
%   C(5)  = [-1 +1  0]
%   C(6)  = [+1  0 +1]
%   C(7)  = [-1  0 +1]
%   C(8)  = [ 0 +1 +1]
%   C(9)  = [ 0 -1 +1]
%   C(10) = [+1 +1 +1]
%   C(11) = [-1 +1 +1]
%   C(12) = [+1 -1 +1]
%   C(13) = [-1 -1 +1]
%   The sum of the weights in C equals 1.
%   Some values are equal whatever the resolution:
%   C(4)==C(5);
%   C(6)==C(7);
%   C(8)==C(9);
%   C(10)==C(11)==C(12)==C(13);
%
%   C = computeDirectionWeights3d13(DELTA)
%   With DELTA = [DX DY DZ], specifies the resolution of the grid.
%
%   Example
%   c = computeDirectionWeights3d13;
%   sum(c)
%   ans =
%       1.0000
%
%   c = computeDirectionWeights3d13([2.5 2.5 7.5]);
%   sum(c)
%   ans =
%       1.0000
%
%
%   See also
%
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2010-10-18,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


%% Initializations

% grid resolution
delta = [1 1 1];
if ~isempty(varargin)
    delta = varargin{1};
end

% If resolution is [1 1 1], return the pre-computed set of weights
if all(delta == [1 1 1])
    area1 = 0.04577789120476 * 2;
    area2 = 0.03698062787608 * 2;
    area3 = 0.03519563978232 * 2;
    res = [...
        area1; area1; area1; ...
        area2; area2; area2; area2; area2; area2;...
        area3; area3; area3; area3 ];
    return;
end

% Define points of interest in the 26 discrete directions
% format is pt[Xpos][Ypos][Zpos], with [X], [Y] or [Z] being one of 
% 'N' (for negative), 'P' (for Positive) or 'Z' (for Zero)

% points below the OXY plane
ptPNN = normalizeVector([+1 -1 -1].*delta);
ptPZN = normalizeVector([+1  0 -1].*delta);
ptNPN = normalizeVector([-1 +1 -1].*delta);
ptZPN = normalizeVector([ 0 +1 -1].*delta);
ptPPN = normalizeVector([+1 +1 -1].*delta);

% points belonging to the OXY plane
ptPNZ = normalizeVector([+1 -1  0].*delta);
ptPZZ = normalizeVector([+1  0  0].*delta);
ptNPZ = normalizeVector([-1 +1  0].*delta);
ptZPZ = normalizeVector([ 0 +1  0].*delta);
ptPPZ = normalizeVector([+1 +1  0].*delta);

% points above the OXY plane
ptNNP = normalizeVector([-1 -1 +1].*delta);
ptZNP = normalizeVector([ 0 -1 +1].*delta);
ptPNP = normalizeVector([+1 -1 +1].*delta);
ptNZP = normalizeVector([-1  0 +1].*delta);
ptZZP = normalizeVector([ 0  0 +1].*delta);
ptPZP = normalizeVector([+1  0 +1].*delta);
ptNPP = normalizeVector([-1 +1 +1].*delta);
ptZPP = normalizeVector([ 0 +1 +1].*delta);
ptPPP = normalizeVector([+1 +1 +1].*delta);


%% Spherical cap type 1, direction [1 0 0]

% Compute area of voronoi cell for a point on the Ox axis, i.e. a point
% in the 6-neighborhood of the center.
refPoint = ptPZZ;

% neighbours of chosen point, sorted by CCW angle
neighbors = [ptPNN; ptPNZ; ptPNP; ptPZP; ptPPP; ptPPZ; ptPPN; ptPZN];

% compute area of spherical polygon
area1 = sphericalVoronoiDomainArea(refPoint, neighbors);


%% Spherical cap type 1, direction [0 1 0]

% Compute area of voronoi cell for a point on the Oy axis, i.e. a point
% in the 6-neighborhood of the center.
refPoint    = ptZPZ;

% neighbours of chosen point, sorted by angle
neighbors   = [ptPPZ; ptPPP; ptZPP; ptNPP; ptNPZ; ptNPN; ptZPN; ptPPN];

% compute area of spherical polygon
area2 = sphericalVoronoiDomainArea(refPoint, neighbors);


%% Spherical cap type 1, direction [0 0 1]

% Compute area of voronoi cell for a point on the Oz axis, i.e. a point
% in the 6-neighborhood of the center.
refPoint = ptZZP;

% neighbours of chosen point, sorted by angle
neighbors = [ptPZP; ptPPP; ptZPP; ptNPP; ptNZP; ptNNP; ptZNP; ptPNP];

% compute area of spherical polygon
area3 = sphericalVoronoiDomainArea(refPoint, neighbors);


%% Spherical cap type 2, direction [1 1 0]

% Compute area of voronoi cell for a point on the Oxy plane, i.e. a point
% in the 18-neighborhood
refPoint = ptPPZ;

% neighbours of chosen point, sorted by angle
neighbors = [ptPZZ; ptPPP; ptZPZ; ptPPN];

% compute area of spherical polygon
area4 = sphericalVoronoiDomainArea(refPoint, neighbors);


%% Spherical cap type 2, direction [1 0 1]

% Compute area of voronoi cell for a point on the Oxz plane, i.e. a point
% in the 18-neighborhood
refPoint = ptPZP;
% neighbours of chosen point, sorted by angle
neighbors = [ptPZZ; ptPPP; ptZZP; ptPNP];

% compute area of spherical polygon
area5 = sphericalVoronoiDomainArea(refPoint, neighbors);


%% Spherical cap type 2, direction [0 1 1]

% Compute area of voronoi cell for a point on the Oxy plane, i.e. a point
% in the 18-neighborhood
refPoint = ptZPP;
% neighbours of chosen point, sorted by angle
neighbors = [ptZPZ; ptNPP; ptZZP; ptPPP];

% compute area of spherical polygon
area6 = sphericalVoronoiDomainArea(refPoint, neighbors);


%% Spherical cap type 3 (all cubic diagonals)

% Compute area of voronoi cell for a point on the Oxyz diagonal, i.e. a
% point in the 26 neighborhood only
refPoint = ptPPP;
% neighbours of chosen point, sorted by angle
neighbors = [ptPZP; ptZZP; ptZPP; ptZPZ; ptPPZ; ptPZZ];

% compute area of spherical polygon
area7 = sphericalVoronoiDomainArea(refPoint, neighbors);


%% Concatenate results

% return computed areas, formatted as fraction of sphere surface
res = [...
    area1 area2 area3 ...
    area4 area4 area5 area5 area6 area6...
    area7 area7 area7 area7...
    ]/(2*pi);
end


function vn = normalizeVector(v)
%NORMALIZEVECTOR normalize a vector
%
%   V2 = normalizeVector(V);
%   Returns the normalization of vector V, such that ||V|| = 1. V can be
%   either a row or a column vector.
%
%   When V is a MxN array, normalization is performed for each row of the
%   array.
%
%   Example:
%   normalizeVector([3 3])
%       [ 0.7071   0.7071 ]
%
%   See Also:
%   vectors2d, vecnorm
%
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 29/11/2004.
%

%   HISTORY
%   14/01/2005 correct bug
%   22/05/2009 rename as normalizeVector


dim = size(v);

if dim(1)==1 || dim(2)==1
    vn = v/sqrt(sum(v.*v));
else
    vn = v./repmat(sqrt(sum(v.*v, 2)), [1 dim(2)]);
end

end



function area = sphericalVoronoiDomainArea(refPoint, neighbors)
%SPHERICALVORONOIDOMAINAREA Compute area of a spherical voronoi domain
%
%   AREA = sphericalVoronoiDomainArea(GERM, NEIGHBORS)
%   GERM is a 1-by-3 row vector representing cartesian coordinates of a
%   point on the sphere (in X, Y Z order)
%   NEIGHBORS is a N-by-3 array representing cartesian coordinates of the
%   germ neighbors. It is expected that NEIGHBORS contains only neighbors
%   that effectively contribute to the voronoi domain.
%
%   Example
%   sphericalVoronoiDomainArea
%
%   See also
%
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2010-11-17,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

% reference sphere
sphere = [0 0 0 1];

% number of neigbors, and number of sides of the domain
nbSides = size(neighbors, 1);

% compute planes containing separating circles
planes = zeros(nbSides, 9);
for i=1:nbSides
    planes(i,1:9) = normalizePlane(medianPlane(refPoint, neighbors(i,:)));
end

% allocate memory
lines       = zeros(nbSides, 6);
intersects  = zeros(2*nbSides, 3);

% compute circle-circle intersections
for i=1:nbSides
    lines(i,1:6) = intersectPlanes(planes(i,:), ...
        planes(mod(i,nbSides)+1,:));
    intersects(2*i-1:2*i,1:3) = intersectLineSphere(lines(i,:), sphere);
end

% keep only points in the same direction than refPoint
ind = dot(intersects, repmat(refPoint, [2*nbSides 1]), 2)>0;
intersects = intersects(ind,:);
nbSides = size(intersects, 1);

% compute spherical area of each triangle [center  pt[i+1]%4   pt[i] ]
angles = zeros(nbSides, 1);
for i=1:nbSides
    pt1 = intersects(i, :);
    pt2 = intersects(mod(i  , nbSides)+1, :);
    pt3 = intersects(mod(i+1, nbSides)+1, :);
    
    angles(i) = sphericalAngle(pt1, pt2, pt3);
    angles(i) = min(angles(i), 2*pi-angles(i));
end

% compute area of spherical polygon
area = sum(angles) - pi*(nbSides-2);

end



function plane = medianPlane(p1, p2)
%MEDIANPLANE create a plane in the middle of 2 points
%
%   PLANE = medianPlane(P1, P2)
%   Creates a plane in the middle of 2 points.
%   PLANE is perpendicular to line (P1 P2) and contains the midpoint of P1
%   and P2.
%   The direction of the normal of PLANE is the same as the vector from P1
%   to P2.
%
%   See also:
%   planes3d, createPlane
%
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 18/02/2005.
%

%   HISTORY
%   28/06/2007: add doc, and manage multiple inputs

% unify data dimension
if size(p1, 1)==1
    p1 = repmat(p1, [size(p2, 1) 1]);
elseif size(p2, 1)==1
    p2 = repmat(p2, [size(p1, 1) 1]);
elseif size(p1, 1)~=size(p2, 1)    
    error('data should have same length, or one data should have length 1');
end

% middle point
p0  = (p1 + p2)/2;

% normal to plane
n   = p2-p1;

% create plane from point and normal
plane = createPlane(p0, n);

end



function plane = createPlane(varargin)
%CREATEPLANE create a plane in parametrized form
%
%   Creates a plane in the following format: 
%   PLANE = [X0 Y0 Z0  DX1 DY1 DZ1  DX2 DY2 DZ2], where :
%   - (X0, Y0, Z0) is a point belonging to the plane
%   - (DX1, DY1, DZ1) is a first direction vector
%   - (DX2, DY2, DZ2) is a second direction vector
%   The 2 direction vectors are normalized and orthogonal.
%
%
%   PLANE = createPlane(P1, P2, P3) 
%   create a plane containing the 3 points
%
%   PLANE = createPlane(PTS) 
%   The 3 points are packed into a single 3x3 array.
%
%   PLANE = createPlane(P0, N);
%   create a plane from a point and from a normal to the plane. parameter N
%   is given as [THETA PHI], where THETA is the colatitute (angle with the
%   vertical axis) and PHI is angle with Ox axis, counted
%   counter-clockwise. Both are given in radians.
%   
%   See also:
%   planes3d, medianPlane
%   
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 18/02/2005.
%

%   HISTORY :
%   24/11/2005 : add possibility to pack points for plane creation
%   21/08/2006 : return normalized planes
%   06/11/2006 : update doc for planes created from normal

if length(varargin)==1
    var = varargin{1};
    
    if iscell(var)
        plane = zeros([0 9]);
        for i=1:length(var)
            plane = [plane; createPlane(var{i})];
        end
    elseif size(var, 1)==3
        % 3 points in a single array
        p1 = var(1,:);
        p2 = var(2,:);
        p3 = var(3,:);
        
        % create direction vectors
        v1 = p2-p1;
        v2 = p3-p1;

        % create plane
        plane = normalizePlane([p1 v1 v2]);
        return;
    end
    
elseif length(varargin)==2
    % plane origin
    p0 = varargin{1};
    
    % second parameter is either a 3D vector or a 3D angle (2 params)
    var = varargin{2};
    if size(var, 2)==2
        n = sph2cart2([var ones(size(var, 1))]);
    elseif size(var, 2)==3
        n = normalizeVector3d(var);
    else
        error ('wrong number of parameters in createPlane');
    end
    
    % ensure same dimension for parameters
    if size(p0, 1)==1
        p0 = repmat(p0, [size(n, 1) 1]);
    end
    if size(n, 1)==1
        n = repmat(n, [size(p0, 1) 1]);
    end

    % find a vector not colinear to the normal
    v0 = repmat([1 0 0], [size(p0, 1) 1]);
    inds = vectorNorm(cross(n, v0, 2))<1e-14;
    v0(inds, :) = repmat([0 1 0], [sum(inds) 1]);
%     if abs(cross(n, v0, 2))<1e-14
%         v0 = repmat([0 1 0], [size(p0, 1) 1]);
%     end
    
    % create direction vectors
    v1 = normalizeVector3d(cross(n, v0, 2));
    v2 = -normalizeVector3d(cross(v1, n, 2));

    % concatenate result in the array representing the plane
    plane = [p0 v1 v2];
    return;
    
elseif length(varargin)==3
    p1 = varargin{1};    
    p2 = varargin{2};
    p3 = varargin{3};
    
    % create direction vectors
    v1 = p2-p1;
    v2 = p3-p1;
   
    plane = normalizePlane([p1 v1 v2]);
    return;
  
else
    error('wrong number of arguments in "createPlane".');
end


end



function vn = normalizeVector3d(v)
%NORMALIZEVECTOR3D normalize a 3D vector
%
%   V2 = normalizeVector3d(V);
%   Returns the normalization of vector V, such that ||V|| = 1. Vector V is
%   given as a row vector.
%
%   When V is a Nx3 array, normalization is performed for each row of the
%   array.
%
%   See also:
%   vectors3d, vectorNorm3d
%
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 29/11/2004.
%

%   HISTORY
%   30/11/2005 correct a bug
%   19/06/2009 rename as normalizeVector3d

n   = sqrt(v(:,1).*v(:,1) + v(:,2).*v(:,2) + v(:,3).*v(:,3));
vn  = v./[n n n];

end


function n = vectorNorm(v, varargin)
%VECTORNORM compute norm of vector or of set of vectors
%
%   N = vectorNorm(V);
%   Returns the euclidean norm of vector V.
%
%   N = vectorNorm(V, N);
%   Specifies the norm to use. N can be any value greater than 0. 
%   N=1 -> city lock norm
%   N=2 -> euclidean norm
%   N=inf -> compute max coord.
%
%   When V is a MxN array, compute norm for each vector of the array.
%   Vector are given as rows. Result is then a [M*1] array.
%
%   See Also:
%   vectors2d, vectorAngle
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 21/02/2005.
%

%   HISTORY
%   02/05/2006 manage several norms
%   18/09/2007 use 'isempty'
%   15/10/2008 add comments
%   22/05/2009 rename as vectorNorm
%   01/03/2010 fix bug for inf norm


% size of vector
dim = size(v);

% extract the type of norm to compute
d = 2;
if ~isempty(varargin)
    d = varargin{1};
end

if d==2
    % euclidean norm: sum of squared coordinates, and take square root
    if dim(1)==1 || dim(2)==1
        n = sqrt(sum(v.*v));
    else
        n = sqrt(sum(v.*v, 2));
    end
elseif d==1 
    % absolute norm: sum of absolute coordinates
    if dim(1)==1 || dim(2)==1
        n = sum(abs(v));
    else
        n = sum(abs(v), 2);
    end
elseif d==inf
    % infinite norm: uses the maximal corodinate
    if dim(1)==1 || dim(2)==1
        n = max(v);
    else
        n = max(v, [], 2);
    end
else
    % Other norms, use explicit but slower expression  
    if dim(1)==1 || dim(2)==1
        n = power(sum(power(v, d)), 1/d);
    else
        n = power(sum(power(v, d), 2), 1/d);
    end
end

end


function plane2 = normalizePlane(plane1)
%NORMALIZEPLANE normalize parametric form of a plane
%
%   PLANE2 = normalizePlane(PLANE1);
%   Transforms the plane PLANE1 in the following format :
%   [X0 Y0 Z0  DX1 DY1 DZ1  DX2 DY2 DZ2], where :
%   - (X0, Y0, Z0) is a point belonging to the plane
%   - (DX1, DY1, DZ1) is a first direction vector
%   - (DX2, DY2, DZ2) is a second direction vector
%   into another plane, with the same format, but with :
%   - (x0 y0 z0) is the closest point of plane to origin
%   - (DX1 DY1 DZ1) has norm equal to 1
%   - (DX2 DY2 DZ2) has norm equal to 1 and is orthogonal to (DX1 DY1 DZ1)
%   
%   See also:
%   planes3d, createPlane
%
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 21/02/2005.
%

%   HISTORY
%   21/08/2009 compute origin after computation of vectors (more precise)
%       and add support for several planes.

% compute first direction vector
d1  = normalizeVector3d(plane1(:,4:6));

% compute second direction vector
n   = normalizeVector3d(planeNormal(plane1));
d2  = -normalizeVector3d(cross(d1, n));

% compute origin point of the plane
origins = repmat([0 0 0], [size(plane1, 1) 1]);
p0 = projPointOnPlane(origins, [plane1(:,1:3) d1 d2]);

% create the resulting plane
plane2 = [p0 d1 d2];


end



function n = planeNormal(plane)
%PLANENORMAL compute the normal to a plane
%
%   N = planeNormal(PLANE) 
%   compute the normal of the given plane
%   PLANE : [x0 y0 z0 dx1 dy1 dz1 dx2 dy2 dz2]
%   N : [dx dy dz]
%   
%   See also
%   planes3d, createPlane
%
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 17/02/2005.
%

%   HISTORY


% plane normal
n = cross(plane(:,4:6), plane(:, 7:9), 2);

end



function point = projPointOnPlane(point, plane)
%PROJPOINTONPLANE return the projection of a point on a plane
%
%   PT2 = projPointOnPlane(PT1, PLANE);
%   Compute the (orthogonal) projection of point PT1 onto the line PLANE.
%   
%   Function works also for multiple points and planes. In this case, it
%   returns multiple points.
%   Point PT1 is a [N*3] array, and PLANE is a [N*9] array (see createPlane
%   for details). Result PT2 is a [N*3] array, containing coordinates of
%   orthogonal projections of PT1 onto planes PLANE.
%
%   See also:
%   planes3d, points3d, planePosition, intersectLinePlane
%
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 18/02/2005.
%

%   HISTORY
%   21/08/2006 : debug support for multiple points or planes

if size(point, 1)==1
    point = repmat(point, [size(plane, 1) 1]);
elseif size(plane, 1)==1
    plane = repmat(plane, [size(point, 1) 1]);
elseif size(plane, 1) ~= size(point, 1)
    error('projPointOnPlane : size of inputs differ');
end

n = planeNormal(plane);
line = [point n];
point = intersectLinePlane(line, plane);


end



function point = intersectLinePlane(line, plane)
%INTERSECTLINEPLANE return intersection between a plane and a line
%
%   PT = intersectLinePlane(LINE, PLANE)
%   Returns the intersection point of the given line and the given plane.
%   PLANE : [x0 y0 z0 dx1 dy1 dz1 dx2 dy2 dz2]
%   LINE :  [x0 y0 z0 dx dy dz]
%   PT :    [xi yi zi]
%   If LINE and PLANE are parallel, return [NaN NaN NaN].
%   If LINE (or PLANE) is a matrix with 6 (or 9) columns and N rows, result
%   is an array of points with N rows and 3 columns.
%   
%   See also:
%   lines3d, planes3d, points3d, intersectEdgePlane
%
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 17/02/2005.
%

%   HISTORY
%   24/11/2005 add support for multiple input
%   23/06/2006 correction from Songbai Ji
%   14/12/2006 correction for parallel lines and plane normals
%   05/01/2007 fixup for parallel lines and plane normals
%   24/04/2007 rename as 'intersectLinePlane'

%  Songbai Ji (6/23/2006). Bug fixed; also allow one plane, many lines; 
% many planes one line; or N planes and N lines configuration in the input.

% unify sizes of data
if size(line,1) == 1;   % one line and many planes
    line = repmat(line, size(plane,1), 1);
elseif size(plane, 1) == 1;     % one plane possible many lines
    plane = repmat(plane, size(line,1), 1);
elseif (size(plane,1) ~= size(line,1)) ; % N planes and M lines, not allowed for now.
    error('input size not correct, either one/many plane and many/one line, or same # of planes and lines!');
end

% initialize empty array
point = zeros(size(plane, 1), 3);

% plane normal
n = cross(plane(:,4:6), plane(:,7:9), 2);

% get indices of line and plane which are parallel
par = abs(dot(n, line(:,4:6), 2))<1e-14;
point(par,:) = NaN;

% difference between origins of plane and line
dp = plane(:, 1:3) - line(:, 1:3);

% relative position of intersection on line
% Should be Array multiply, original file had a bug. (songbai ji
% 6/23/2006).
% Divide only for non parallel vectors (DL)
t = dot(n(~par,:), dp(~par,:), 2)./dot(n(~par,:), line(~par,4:6), 2);

% compute coord of intersection point
point(~par, :) = line(~par,1:3) + repmat(t,1,3).*line(~par,4:6);


end



function line = intersectPlanes(plane1, plane2)
%INTERSECTPLANES return intersection between 2 planes in space
%
%   LINE = intersectPlanes(PLANE1, PLANE2)
%   Returns the straight line belonging to both planes
%   PLANE : [x0 y0 z0 dx1 dy1 dz1 dx2 dy2 dz2]
%   LINE :  [x0 y0 z0 dx dy dz]
%
%   See also:
%   planes3d, lines3d, intersectLinePlane
%
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 17/02/2005.
%

%   HISTORY


% plane normal
n1 = normalizeVector3d(cross(plane1(:,4:6), plane1(:, 7:9), 2));
n2 = normalizeVector3d(cross(plane2(:,4:6), plane2(:, 7:9), 2));

% test if planes are parallel
if abs(cross(n1, n2, 2))<1e-14
    line = [NaN NaN NaN NaN NaN NaN];
    return;
end

% Uses Hessian form, ie : N.p = d
% I this case, d can be found as : -N.p0, when N is normalized
d1 = dot(n1, plane1(:,1:3), 2);
d2 = dot(n2, plane2(:,1:3), 2);

% compute dot products
dot1 = dot(n1, n1, 2);
dot2 = dot(n2, n2, 2);
dot12 = dot(n1, n2, 2);


det = dot1*dot2 - dot12*dot12;
c1  = (d1*dot2 - d2*dot12)./det;
c2  = (d2*dot1 - d1*dot12)./det;

p0  = c1*n1 + c2*n2;
dp  = cross(n1, n2, 2);

line = [p0 dp];

end



function point = intersectLineSphere(line, sphere)
%INTERSECTLINESPHERE return intersection between a line and a sphere
%
%   GC = intersectLineSphere(LINE, SPHERE);
%   Returns the two points which are the intersection of the given line and
%   sphere. 
%   LINE   : [x0 y0 z0  dx dy dz]
%   SPHERE : [xc yc zc  R]
%   GC     : [x1 y1 z1 ; x2 y2 z2]
%   
%   See also
%   spheres, circles3d
%
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 18/02/2005.
%

%   HISTORY

% difference between centers
dc = line(1:3)-sphere(1:3);

a = sum(line(:, 4:6).*line(:, 4:6), 2);
b = 2*sum(dc.*line(4:6), 2);
c = sum(dc.*dc, 2) - sphere(:,4).*sphere(:,4);

delta = b.*b -4*a.*c;

if delta>1e-14
    % find two roots of second order equation
    u1 = (-b -sqrt(delta))/2/a;
    u2 = (-b +sqrt(delta))/2/a;
    
    % convert into 3D coordinate
    point = [line(1:3)+u1*line(4:6) ; line(1:3)+u2*line(4:6)];

elseif abs(delta) > 1e-14
    % find unique root, and convert to 3D coord.
    u = -b/2./a;    
    point = line(1:3) + u*line(4:6);
    
else
    point = zeros(0, 3);
    return;
end

end


function alpha = sphericalAngle(p1, p2, p3)
%SPHERICALANGLE compute angle on the sphere
%
%   ALPHA = sphericalAngle(P1, P2, P3)
%   compute angle (P1, P2, P2), in radians, between 0 and 2*PI.
%
%   Points are given either as [x y z] (there will be normalized to lie on
%   the unit sphere), or as [phi theta], with phi being the longitude in [0
%   2*PI] and theta being the elevation on horizontal [-pi/2 pi/2].
%
%
%   NOTE: 
%   this is an 'oriented' version of the angle computation, that is, the
%   result of sphericalAngle(P1, P2, P3) equals
%   2*pi-sphericalAngle(P3,P2,P1). To have the more classical relation
%   (with results given betwen 0 and PI), it suffices to take the minimum
%   of angle and 2*pi-angle.
%   
%   See also:
%   angles3d, spheres
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 21/02/2005.
%

%   HISTORY
%   23/05/2006 fix bug for points with angle from center > pi/2

% test if points are given as matlab spherical coordinate
if size(p1, 2) ==2
    [x y z] = sph2cart(p1(:,1), p1(:,2));
    p1 = [x y z];
    [x y z] = sph2cart(p2(:,1), p2(:,2));
    p2 = [x y z];
    [x y z] = sph2cart(p3(:,1), p3(:,2));
    p3 = [x y z];
end

% normalize points
p1  = normalizeVector3d(p1);
p2  = normalizeVector3d(p2);
p3  = normalizeVector3d(p3);

% create the plane tangent to the unit sphere and containing central point
plane = createPlane(p2, p2);

% project the two other points on the plane
pi1 = planePosition(projPointOnPlane(p1, plane), plane);
pi3 = planePosition(projPointOnPlane(p3, plane), plane);

% compute angle on the tangent plane
alpha = angle3Points(pi1, [0 0], pi3);

end






function pos = planePosition(point, plane)
%PLANEPOSITION compute position of a point on a plane
%
%   PT2 = planePosition(POINT, PLANE)
%   POINT has format [X Y Z], and plane has format
%   [X0 Y0 Z0  DX1 DY1 DZ1  DX2 DY2 DZ2], where :
%   - (X0, Y0, Z0) is a point belonging to the plane
%   - (DX1, DY1, DZ1) is a first direction vector
%   - (DX2, DY2, DZ2) is a second direction vector
%
%   Result PT2 has the form [XP YP], with [XP YP] coordinate of the point
%   in the coordinate system of the plane.
%
%   
%   CAUTION :
%   WORKS ONLY FOR PLANES WITH ORTHOGONAL DIRECTION VECTORS
%
%   See also:
%   planes3d, points3d, planePoint
%
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 21/02/2005.
%

%   HISTORY :
%   24/11/2005 add support for multiple input

% unify size of data
if size(point, 1)~=size(plane, 1)
    if size(point, 1)==1
        point = repmat(point, [size(plane, 1) 1]);
    elseif size(plane, 1)==1
        plane = repmat(plane, [size(point, 1) 1]);
    else
        error('point and plane do not have the same dimension');
    end
end


p0 = plane(:,1:3);
d1 = plane(:,4:6);
d2 = plane(:,7:9);

s = dot(point-p0, d1, 2)./vectorNorm3d(d1);
t = dot(point-p0, d2, 2)./vectorNorm3d(d2);

pos = [s t];

end



function n = vectorNorm3d(v)
%VECTORNORM3D compute norm of vector or of set of 3D vectors
%
%   N = vectorNorm3d(V);
%   Returns norm of vector V.
%
%   When V is a Nx3 array, compute norm for each vector of the array.
%   Vector are given as rows. Result is then a [N*1] array.
%
%   NOTE: compute only euclidean norm.
%
%   See Also
%   vectors3d, normalizeVector3d
%
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 21/02/2005.

%   HISTORY
%   19/06/2009 rename as vectorNorm3d

n = sqrt(sum(v.*v, 2));

end

function theta = angle3Points(varargin)
%ANGLE3POINTS return oriented angle made by 3 points
%
%   ALPHA = ANGLE3POINTS(P1, P2, P3).
%   Pi are either [1*2] arrays, or [N*2] arrays, in this case ALPHA is a 
%   [N*1] array. The angle computed is the directed angle between line 
%   (P2P1) and line (P2P3).
%   Result is always given in radians, between 0 and 2*pi.
%
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 23/02/2004.
%

%   HISTORY :
%   25/09/2005 : enable single parameter

if length(varargin)==3
    p1 = varargin{1};
    p2 = varargin{2};
    p3 = varargin{3};
elseif length(varargin)==1
    var = varargin{1};
    p1 = var(1,:);
    p2 = var(2,:);
    p3 = var(3,:);
end    

% angle line (P2 P1)
theta = lineAngle(createLine(p2, p1), createLine(p2, p3));


function theta = lineAngle(varargin)
%LINEANGLE return angle between lines
%
%   a = LINEANGLE(line) return the angle between horizontal, right-axis 
%   and the given line. Angle is fiven in radians, between 0 and 2*pi,
%   in counter-clockwise direction.
%
%   a = LINEANGLE(line1, line2) return the directed angle between the
%   two lines. Angle is given in radians between 0 and 2*pi.
%
%   see createLine for more details on line representation.
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 31/10/2003.
%

%   HISTORY :
%   19/02/2004 : added support for multiple lines.

nargs = length(varargin);
if nargs == 1
    % one line
    line = varargin{1};
    theta = mod(atan2(line(:,4), line(:,3)) + 2*pi, 2*pi);
elseif nargs==2
    % two lines
    theta1 = lineAngle(varargin{1});
    theta2 = lineAngle(varargin{2});
    theta = mod(theta2-theta1+2*pi, 2*pi);
end


end


function line = createLine(varargin)
%CREATELINE create a line with various inputs.
%
%   Line is represented in a parametric form : [x0 y0 dx dy]
%   x = x0 + t*dx
%   y = y0 + t*dy;
%
%
%   l = CREATELINE(p1, p2) return the line going through the two given
%   points.
%   
%   l = CREATELINE(x0, y0, dx, dy) the line going through point (x0, y0)
%   and with direction vector(dx, dy).
%
%   l = CREATELINE(param) where param is an array of 4 values, create the
%   line oing through the point (param(1) param(2)), and with direction
%   vector (param(3) param(4)).
%   
%   l = CREATELINE(theta) create a line originated at (0,0) and
%   with angle theta.
%
%   l = CREATELINE(rho, theta) create a line with normal theta, and with
%   min distance to origin equal to rho. rho can be negative, in this case,
%   the line is the same as with CREATELINE(-rho, theta+pi), but the
%   orientation is different.
%
%
%   Note : in all cases, parameters can be vertical arrays of the same
%   dimension. The result is then an array of lines, of dimensions [N*4].
%
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 31/10/2003.
%

%   HISTORY :
%   18/02/2004 : add more possibilities to create lines (4 parameters,
%      all param in a single tab, and point + dx + dy.
%      Also add support for creation of arrays of lines.


if length(varargin)==1
    % Only one input parameter. It can be :
    % - line angle
    % - array of four parameters
    var = varargin{1};
    
    if size(var, 2)==4
        % 4 parameters of the line in a single array.
        line = var;
    elseif size(var, 2)==1
        % 1 parameter : angle of the line, going through origin.
        line = [zeros(size(var), 1) zeros(size(var, 1)) cos(var) sin(var)];
    else
        error('wrong number of dimension for arg1 : can be 1 or 4');
    end
    
elseif length(varargin)==2    
    % 2 input parameters. They can be :
    % - line angle and signed distance to origin.
    % - 2 points, then 2 arrays of 1*2 double.
    v1 = varargin{1};
    v2 = varargin{2};
    if size(v1, 2)==1
        % first param is signed distance to origin, and second param is
        % angle of line.
        line = [v1.*cos(v2) v1.*sin(v2) -sin(v2) cos(v2)];
    else
        % first input parameter is first point, and second input is the
        % second point.
        line = [v1(:,1), v1(:,2), v2(:,1)-v1(:,1), v2(:,2)-v1(:,2)];    
    end
    
elseif length(varargin)==3
    % 3 input parameters :
    % first one is a point belonging to the line,
    % second and third ones are direction vector of the line (dx and dy).
    p = varargin{1};
    line = [p(:,1) p(:,2) varargin{2} varargin{3}];
   
elseif length(varargin)==4
    % 4 input parameters :
    % they are x0, y0 (point belongng to line) and dx, dy (direction vector
    % of the line).
    % All parameters should have the same size.
    line = [varargin{1} varargin{2} varargin{3} varargin{4}];
else
    error('Wrong number of arguments in ''createLine'' ');
end

end

end


function alpha = anglePoints3d(varargin)
%ANGLEPOINTS3D Compute angle between three 3D points
%
%   ALPHA = anglePoints3d(P1, P2)
%   Computes angle (P1, O, P2), in radians, between 0 and PI.
%
%   ALPHA = anglePoints3d(P1, P2, P3)
%   Computes angle (P1, P2, P3), in radians, between 0 and PI.
%
%   ALPHA = anglePoints3d(PTS)
%   PTS is a 3x3 or 2x3 array containing coordinate of points.
%
%   See also
%   points3d, angles3d
%
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 21/02/2005.
%

%   HISTORY
%   20/09/2005 : add case of single argument for all points
%   04/01/2007: check typo

p2 = [0 0 0];
if length(varargin)==1
    pts = varargin{1};
    if size(pts, 1)==2
        p1 = pts(1,:);
        p0 = [0 0 0];
        p2 = pts(2,:);
    else
        p1 = pts(1,:);
        p0 = pts(2,:);
        p2 = pts(3,:);
    end
elseif length(varargin)==2
    p1 = varargin{1};
    p0 = [0 0 0];
    p2 = varargin{2};
elseif length(varargin)==3
    p1 = varargin{1};
    p0 = varargin{2};
    p2 = varargin{3};
end

% ensure all data have same size
n1 = size(p1, 1);
n2 = size(p2, 1);
n0 = size(p0, 1);
if n1~=n2
    if n1==1
        p1 = repmat(p1, [n2 1]);
    elseif n2==1
        p2 = repmat(p2, [n1 1]);
    else
        error('Arguments P1 and P2 must have the same size');
    end
end
if n1~=n0
    if n1==1
        p1 = repmat(p1, [n0 1]);
    elseif n0==1
        p0 = repmat(p0, [n1 1]);
    else
        error('Arguments P1 and P0 must have the same size');
    end
end

% normalized vectors
p1 = normalizeVector3d(p1-p0);
p2 = normalizeVector3d(p2-p0);

% compute angle
alpha = acos(dot(p1, p2, 2));

end



end
