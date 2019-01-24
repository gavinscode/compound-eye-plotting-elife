%Written by Gavin Taylor, 2017. MIT License

%This is the computational analysis code for Taylor et al 2019
%it is not well formated or commented, and is set to use the same parameters as in the analysis in the paper
%please contact gavin.taylor.01@gmail.com of assistance if required.  


ConfigForBeeCEComp %load structure of analysis paramters
close all

saveDataFolder = 'Add folder here'; %set folder to save data
saveDataAs = 'Add file name here'; %'Add file name here' %set file name for saved data

fullBeeData = fullBeeData([]); %limit analysis to subset of bees in config file, empty is all
%% settings area
%larger is faster, tends to increase density slightly
measurmentIntervals = 25; %radius microns

%starts with this value, will step around as required.
surfaceFitRngTests = 100;  %in microns radius, for use when getting normals

%to ensure there is some height variation in plot surfaces
%program may decided to adjsut if desired rngs cannot be achieved
baseHeightStep = [5/sqrt(2) 5*sqrt(2)];  
    useLevelSets = 1; %aims to prevent errors from flat bits in xy, yz, or xz plans
                      %however, flat bits in non-orthogonal planes can still cause problems  
    minReqIO = 0.85; %will raise fit rng locally until this value is reached
                        %may not be solving underlying cause

%set which plots to show                        
testPlots = 1; 
    plotGroupHexs = 1;
    plotIndHexs = 0; %one plot for each hexagon that was labelled
    plotLevelSet = 0;
    plotTestAngleCalc = 0;
    plotHexagonInterp = 1;
    plotHexagonInterpInd = 0; %not this will produce thousands of windows, best to test with break points
saveFigs = 0; %best to save with this turned off otherwise files will be huge  

%doesn't choose sample points within this distance of the border
stepInFromBorders = 1;
    stepInInterval = 25; 

autoAlignEyes = 1; %refines registration to align primary axis of eyes vertically, and corrects yaw and roll from eye centers.    
       
loadHeads = 0;
    removeHeadAfterPlot = 1; %otherwise the heads will probably take all the memory

%when selecting hexagons

useDirectHexAxes = 1;
%Either of these can be selected if neither tessolation nor direct axes used, if neither selected, best of either vertical or horizontal is chosen
	%but then interpolation will almost certainly do a bad job, so best to pick one
    forceVertical = 0; %forces facets to be chosen as points vertical
    forceHorizontal = 0; %forces facets to be chosen as points horiztonal
    
%for calculating thickness, otherwise consecutive closest points used
useReverseNormal = 1; 
    CCT_Thresh = 100; %if greater than thresh, closest pt picked, equiv to 0 above
    lensT_Thresh = 70;
    limRetinaDist_Thresh = 1.5; %ratio of closest distance to vector length

IDWThresh = 0.1; %interp values can be plus/minus this ratio of max/min before triggering warning...
powerP = 3; %for IDW on eye, larger value reduces influence of distant points?
%increase in powerP tends to bias error to be positive 
    IDWMinRangeFacet = 100; %only applied to facet interp, both have to be in range (voxels)

nearestOnWorld = 0; %otherwise use IDW, but tends to mess up borders and leave pecks
bordersFromNearest = 1; %include border points so values don't sag to mean there
    powerPonSphere = 4; %for IDW on sphere, max distance between sphere sides is pi, so set to 4 to reduce influence of distant points 
        %increase in this tends to smooth projections onto sphere

sphereRad = 10^9; %large to avoid perspective error thingy
useDensityFromIntegral = 1; %when calculating test and interpolaitnt to sphere
    %not implemented for h/v
trimSphere = 1; %limit found points on sphere to match calculated fov  
    contractBordersToTrim = 1; %can either contract from border layer or remove distant points to eye intersects

%used if density is from paralellogram of lattice points, not from integral
limLargeDensities = 1; %can end up really high otherwise
    densityLim = 0.5; %this seems reasonable

aCoef = 0.0067; %in um^-1, for sensitvity equation
    
%color ranges for plots
curvRng = [0 2000]; curvLabs = {sprintf('%i', curvRng(1)), sprintf('%i', round((curvRng(2)-curvRng(1))/2+curvRng(1))), sprintf('%i', curvRng(2))};
diamRng = [10 30]; diamLabs = {sprintf('%i', diamRng(1)), sprintf('%i', round((diamRng(2)-diamRng(1))/2+diamRng(1))), sprintf('%i', diamRng(2))};
intORng = [0 5]; intOLabs = {sprintf('%i', intORng(1)), sprintf('%i', round((intORng(2)-intORng(1))/2+intORng(1))), sprintf('%i', intORng(2))};
axDRng = [0 5000]; axDLabs = {sprintf('%i', axDRng(1)), sprintf('%i', round((axDRng(2)-axDRng(1))/2+axDRng(1))), sprintf('%i', axDRng(2))};
facARng = [0 1000]; facALabs = {sprintf('%i', facARng(1)), sprintf('%i', round((facARng(2)-facARng(1))/2+facARng(1))), sprintf('%i', facARng(2))};
projARng = [0 5*10^6]; projALabs = {sprintf('%i', projARng(1)), sprintf('%i', round((projARng(2)-projARng(1))/2+projARng(1))), sprintf('%i', projARng(2))};
sensRng = [0 0.5]; sensLabs = {sprintf('%i', sensRng(1)), sprintf('%i', round((sensRng(2)-sensRng(1))/2+sensRng(1))), sprintf('%i', sensRng(2))};

lensTRng = [0 100]; lensTLabs = {sprintf('%i', lensTRng(1)), sprintf('%i', round((lensTRng(2)-lensTRng(1))/2+lensTRng(1))), sprintf('%i', lensTRng(2))};
CCTRng = [0 100]; CCTLabs = {sprintf('%i', CCTRng(1)), sprintf('%i', round((CCTRng(2)-CCTRng(1))/2+CCTRng(1))), sprintf('%i', CCTRng(2))};
RetTRng = [0 500]; RetTLabs = {sprintf('%i', RetTRng(1)), sprintf('%i', round((RetTRng(2)-RetTRng(1))/2+RetTRng(1))), sprintf('%i', RetTRng(2))};

%second param controls num points: 3 - 642, 4 - 2562, 5 - 10242
refSphere = SubdivideSphericalMesh(IcosahedronMesh,5); 
    
numBees = size(fullBeeData,2);

for i = 1:numBees
    fullBeeData(i).OrigHexF = NaN; fullBeeData(i).headEyeF = NaN; fullBeeData(i).SamplingF = NaN;      
    fullBeeData(i).histFig = NaN; fullBeeData(i).hexInterpF = NaN; fullBeeData(i).AreaTestF = NaN;    
    fullBeeData(i).testAngF = NaN; fullBeeData(i).testOnEyeF = NaN; fullBeeData(i).testOnEyeFExtra = NaN;     
    fullBeeData(i).testFOV = NaN; fullBeeData(i).testOnWorldF = NaN; fullBeeData(i).testOnWorldFExtra = NaN;
    fullBeeData(i).LSetF = NaN;
end

minDst = 0;
for i = round(1:size(refSphere.X,1)/50:size(refSphere.X,1));
    toTest = 1:size(refSphere.X,1); toTest(i) = [];
    minDst = minDst + min(sqrt((refSphere.X(i,1) - refSphere.X(toTest,1)).^2 + (refSphere.X(i,2) - refSphere.X(toTest,2)).^2 + (refSphere.X(i,3) - refSphere.X(toTest,3)).^2));
end
minDst = minDst/length(1:100:size(refSphere.X,1));

%outer loop per bee
for i = 1:numBees
    %ok so here we have to test the interpolation range for each bee
    itFit = 1;
    accetableFit = 0;
    surfaceFitRng = surfaceFitRngTests(itFit);
    fitResults = [];
    surfacesTested = [];
    htStep = baseHeightStep;
    fprintf('processing %s\n', fullBeeData(i).BeeID);
    forceAccept = 0;
    
    while ~accetableFit
        fprintf('using range %.1f\n', surfaceFitRng);
        %% load, preprocesss, register and do test plot
        %get eye data
        fullBeeData(i).LabelsStack = loadTIFStack( fullBeeData(i).StackFolder, 1);
        fullBeeData(i).eyeVolSize = size(fullBeeData(i).LabelsStack);

        %calculate volume of eye
        fullBeeData(i).EyeVolume = length(getIndicesOfMultipleIDs(fullBeeData(i).LabelsStack,[fullBeeData(i).StackLabels.Lens, fullBeeData(i).StackLabels.Retina, ...
                fullBeeData(i).StackLabels.Cones, fullBeeData(i).StackLabels.LaminaC, fullBeeData(i).StackLabels.LensOuterS, fullBeeData(i).StackLabels.LensInnerS, ...
                fullBeeData(i).StackLabels.RetinaOuterS]))*fullBeeData(i).VoxelSize^3;
        fullBeeData(i).lensVolume = length(getIndicesOfMultipleIDs(fullBeeData(i).LabelsStack,[fullBeeData(i).StackLabels.Lens, fullBeeData(i).StackLabels.LensOuterS, ...
                fullBeeData(i).StackLabels.LensInnerS]))*fullBeeData(i).VoxelSize^3;
        fullBeeData(i).CCVolume = length(getIndicesOfMultipleIDs(fullBeeData(i).LabelsStack,[fullBeeData(i).StackLabels.Cones]))*fullBeeData(i).VoxelSize^3;
        fullBeeData(i).retinaVolume = length(getIndicesOfMultipleIDs(fullBeeData(i).LabelsStack,[fullBeeData(i).StackLabels.Retina, fullBeeData(i).StackLabels.LaminaC, ...
                fullBeeData(i).StackLabels.RetinaOuterS]))*fullBeeData(i).VoxelSize^3;

        if fullBeeData(i).LeftEyeOriginal
            origEyeVolTrans = fullBeeData(i).LeftEyeTrans;
            mirrorEyeVolTrans = fullBeeData(i).RightEyeTrans;
            mirrorEyVolTransOrig = fullBeeData(i).RightEyeTrans; %preserve for later
        else
            origEyeVolTrans = fullBeeData(i).RightEyeTrans;
            mirrorEyeVolTrans = fullBeeData(i).LeftEyeTrans;
            mirrorEyVolTransOrig = fullBeeData(i).LeftEyeTrans; %preserve for later
        end    

        %get cornea subs
        fullBeeData(i).EyeFrontSurfInds = getIndicesOfMultipleIDs(fullBeeData(i).LabelsStack, fullBeeData(i).StackLabels.LensOuterS);
        [LOx,LOy,LOz] = ind2sub(fullBeeData(i).eyeVolSize, fullBeeData(i).EyeFrontSurfInds);

        %new approach, find flat and closed sets in each axis for avoiding eye descritizaion
        if useLevelSets
            if plotLevelSet
                if ishandle(fullBeeData(i).LSetF)
                   close(fullBeeData(i).LSetF);  
                end
                fullBeeData(i).LSetF = figure;
                subplot(1,4,1); 
                hold on; axis equal
                plot3(LOx,LOy,LOz,'k.');
            end     
            xFlatSets = zeros(length(LOx),1);
            yFlatSets = zeros(length(LOy),1);
            zFlatSets = zeros(length(LOz),1);

            convexAreaReq = 0.75;
            %its ok if there is quite a few bits along the bordes, but main regions need to be caught correctly

            %scan through x planes
            xSetsCount = 0;
            for j = 1:fullBeeData(i).eyeVolSize(1)
                if sum(sum(fullBeeData(i).LabelsStack(j,:,:) == fullBeeData(i).StackLabels.LensOuterS)) > 0
                    testIm = permute(fullBeeData(i).LabelsStack(j,:,:), [2 3 1]);
                    testIm(testIm ~=  fullBeeData(i).StackLabels.LensOuterS) = 0;
                    cc = bwconncomp(testIm); 
                    %note bw comp pixels are flipped :/
                    stats = regionprops(cc,'FilledArea','Area','ConvexArea','PixelList');
                    for k = 1:length(stats)
                        %count a level set if its closed (so filled = original areas)
                        %and area is similar to convex area 
                        %%% not ideal test, bit of a fudge
                        if stats(k).FilledArea == stats(k).Area & stats(k).Area > stats(k).ConvexArea*convexAreaReq
                            xSetsCount = xSetsCount + 1;
                            inds = sub2ind(fullBeeData(i).eyeVolSize, ones(stats(k).Area,1)*j, stats(k).PixelList(:,2), stats(k).PixelList(:,1));
                            [~, intInds] = intersect(fullBeeData(i).EyeFrontSurfInds ,inds);
                            xFlatSets(intInds) = xSetsCount;
                            if plotLevelSet
                                  plot3(LOx(intInds), LOy(intInds), LOz(intInds), 'mo');
                            end
                        end
                    end
                end
            end

            %scan through y planes
            ySetsCount = 0;
            for j = 1:fullBeeData(i).eyeVolSize(2)
                if sum(sum(fullBeeData(i).LabelsStack(:,j,:) == fullBeeData(i).StackLabels.LensOuterS)) > 0
                    testIm = permute(fullBeeData(i).LabelsStack(:,j,:), [1 3 2]);
                    testIm(testIm ~=  fullBeeData(i).StackLabels.LensOuterS) = 0;
                    cc = bwconncomp(testIm); 
                    stats = regionprops(cc,'FilledArea','Area','ConvexArea','PixelList');
                    for k = 1:length(stats)
                        if stats(k).FilledArea == stats(k).Area & stats(k).Area > stats(k).ConvexArea*convexAreaReq
                            ySetsCount = ySetsCount + 1;
                            inds = sub2ind(fullBeeData(i).eyeVolSize, stats(k).PixelList(:,2), ones(stats(k).Area,1)*j, stats(k).PixelList(:,1));
                            [~, intInds] = intersect(fullBeeData(i).EyeFrontSurfInds ,inds);
                            yFlatSets(intInds) = ySetsCount;
                            if plotLevelSet
                                plot3(LOx(intInds), LOy(intInds), LOz(intInds), 'co');
                            end
                        end
                    end
                end
            end

            %scan through z planes
            zSetsCount = 0;
            for j = 1:fullBeeData(i).eyeVolSize(3)
                if sum(sum(fullBeeData(i).LabelsStack(:,:,j) == fullBeeData(i).StackLabels.LensOuterS)) > 0
                    testIm = fullBeeData(i).LabelsStack(:,:,j);
                    testIm(testIm ~=  fullBeeData(i).StackLabels.LensOuterS) = 0;
                    cc = bwconncomp(testIm); 
                    stats = regionprops(cc,'FilledArea','Area','ConvexArea','PixelList');
                    for k = 1:length(stats)
                        if stats(k).FilledArea == stats(k).Area & stats(k).Area > stats(k).ConvexArea*convexAreaReq
                            zSetsCount = zSetsCount + 1;
                            inds = sub2ind(fullBeeData(i).eyeVolSize, stats(k).PixelList(:,2), stats(k).PixelList(:,1), ones(stats(k).Area,1)*j);
                            [~, intInds] = intersect(fullBeeData(i).EyeFrontSurfInds ,inds);
                            zFlatSets(intInds) = zSetsCount;
                            if plotLevelSet
                                plot3(LOx(intInds), LOy(intInds), LOz(intInds), 'go');
                            end
                        end
                    end
                end
            end
            if plotLevelSet
                %basis of old level set approach, useful for visualization
                XUnique = unique(LOx); 
                YUnique = unique(LOy);
                ZUnique = unique(LOz);

                subplot(1,4,2); hold on; axis equal
                for j = 1:length(XUnique)
                    inds = find(LOx == XUnique(j));
                    plot3(LOx(inds), LOy(inds), LOz(inds),'.','markersize',5);
                end
                subplot(1,4,3); hold on; axis equal
                for j = 1:length(YUnique)
                    inds = find(LOy == YUnique(j));
                    plot3(LOx(inds), LOy(inds), LOz(inds),'.','markersize',5);
                end
                subplot(1,4,4); hold on; axis equal
                for j = 1:length(ZUnique)
                    inds = find(LOz == ZUnique(j));
                    plot3(LOx(inds), LOy(inds), LOz(inds),'.','markersize',5);
                end
            end
        else
           xFlatSets = [];
           yFlatSets = [];
           zFlatSets = []; 
        end
        
        %get eye border voxels
        %firstly grow lens vol
        tempVol1 = zeros(fullBeeData(i).eyeVolSize);
        tempVol1(getIndicesOfMultipleIDs(fullBeeData(i).LabelsStack, fullBeeData(i).StackLabels.Lens)) = 1;
        tempVol1 = imdilate(tempVol1,strel('sphere',1));
        %then remove bit already in lens and front
        tempVol1(getIndicesOfMultipleIDs(fullBeeData(i).LabelsStack, [[fullBeeData(i).StackLabels.Lens, fullBeeData(i).StackLabels.LensOuterS, ...
                fullBeeData(i).StackLabels.LensInnerS]])) = 0;
        %grow this and front then find intersection with front, is border
        tempVol1 = imdilate(tempVol1,strel('sphere',2));
        fullBeeData(i).borderInds = find(tempVol1(fullBeeData(i).EyeFrontSurfInds));
        
        %remove flat sets that are only present on borders
        if useLevelSets
            for j = 1:xSetsCount
                inds = find(xFlatSets == j);
                if length(intersect(fullBeeData(i).EyeFrontSurfInds(fullBeeData(i).borderInds), fullBeeData(i).EyeFrontSurfInds(inds))) == length(inds)
                    xFlatSets(inds) = 0;
                    if plotLevelSet
                        plot3(LOx(inds), LOy(inds), LOz(inds), 'rx');
                    end
                end
            end
            for j = 1:ySetsCount
                inds = find(yFlatSets == j);
                if length(intersect(fullBeeData(i).EyeFrontSurfInds(fullBeeData(i).borderInds), fullBeeData(i).EyeFrontSurfInds(inds))) == length(inds)
                    yFlatSets(inds) = 0;
                    if plotLevelSet
                        plot3(LOx(inds), LOy(inds), LOz(inds), 'rx');
                    end
                end
            end
            for j = 1:zSetsCount
                inds = find(zFlatSets == j);
                if length(intersect(fullBeeData(i).EyeFrontSurfInds(fullBeeData(i).borderInds), fullBeeData(i).EyeFrontSurfInds(inds))) == length(inds)
                    zFlatSets(inds) = 0;
                    if plotLevelSet
                        plot3(LOx(inds), LOy(inds), LOz(inds), 'rx');
                    end
                end
            end
        end
       
        if ~isnan(fullBeeData(i).AddEyeScale )
            %previously eye was increased in size to match head, remove this from transfrom
            origEyeVolTrans = make_transformation_matrix([0 0 0], [0 0 0], [1 1 1]/fullBeeData(i).AddEyeScale)*origEyeVolTrans;
        end
        %note, -1 -> so subscript origin is always zero
        fullBeeData(i).EyeFrontSurfSubsTrans = [[LOx,LOy,LOz]*fullBeeData(i).VoxelSize-1 ones(size(LOx,1),1)]*origEyeVolTrans';
        fullBeeData(i).EyeFrontSurfSubsTrans(:,4) = [];

        %y is vertical direction
        fullBeeData(i).EyeLength = max(fullBeeData(i).EyeFrontSurfSubsTrans(:,2)) - min(fullBeeData(i).EyeFrontSurfSubsTrans(:,2));

        %get parameters of mirror cords
        newVolBounds = (fullBeeData(i).eyeVolSize-1)*fullBeeData(i).VoxelSize;
        newVolBounds = [ 0 0 0 1;...
          0 0 newVolBounds(3) 1;...
          0 newVolBounds(2) 0 1;...
          0 newVolBounds(2) newVolBounds(3) 1;...
          newVolBounds(1) 0 0 1;...
          newVolBounds(1) newVolBounds(2) 0 1;...
          newVolBounds(1) 0 newVolBounds(3) 1;...
          newVolBounds 1]*origEyeVolTrans';
        volOffset = min(newVolBounds(:,1:3));
        newVolSize = ceil((max(newVolBounds(:,1:3)) - min(newVolBounds(:,1:3)))/fullBeeData(i).VoxelSize);

        %calculate right eye subs using Amiras procedure
        fullBeeData(i).mirrorEyeSubs = fullBeeData(i).EyeFrontSurfSubsTrans-volOffset;
        fullBeeData(i).mirrorEyeSubs(:,1) = -fullBeeData(i).mirrorEyeSubs(:,1)+(newVolSize(1)-1)*fullBeeData(i).VoxelSize;
        fullBeeData(i).mirrorEyeSubs = fullBeeData(i).mirrorEyeSubs+volOffset;

        if ~isnan(fullBeeData(i).AddEyeScale)
            %for mirror eye we need to remove effect of scale only from translations
            %so firstly scale up, then transfrom, then scale down.
            %note matrices applied from right to left
            %great success!
            mirrorEyeVolTrans = make_transformation_matrix([0 0 0], [0 0 0], [1 1 1]/fullBeeData(i).AddEyeScale)*mirrorEyeVolTrans*...
                    make_transformation_matrix([0 0 0], [0 0 0], [1 1 1]*fullBeeData(i).AddEyeScale);
        end
        fullBeeData(i).mirrorEyeSubs = [fullBeeData(i).mirrorEyeSubs ones(size(fullBeeData(i).mirrorEyeSubs,1),1)]*mirrorEyeVolTrans';
        fullBeeData(i).mirrorEyeSubs(:,4) = [];

        %get head center
        fullBeeData(i).HeadCenter = mean([fullBeeData(i).EyeFrontSurfSubsTrans' fullBeeData(i).mirrorEyeSubs']');
        fullBeeData(i).EyeFrontSurfSubsTrans = fullBeeData(i).EyeFrontSurfSubsTrans-fullBeeData(i).HeadCenter;
        fullBeeData(i).mirrorEyeSubs = fullBeeData(i).mirrorEyeSubs-fullBeeData(i).HeadCenter;

        %get eye axis of both eyes
        pcaAx = pca([fullBeeData(i).EyeFrontSurfSubsTrans' fullBeeData(i).mirrorEyeSubs']');
        if autoAlignEyes
            meanPosEye = mean(fullBeeData(i).EyeFrontSurfSubsTrans);
            meanPosMirror = mean(fullBeeData(i).mirrorEyeSubs);
            
            %check 2nd axis is close to y axis
            if sqrt((pcaAx(1,2))^2+(pcaAx(2,2)-1)^2+pcaAx(3,2)^2) < 0.1
                %algin 2nd axis to y to get correct pitch
                extraRotAd = vectorRotationFromAtoB(pcaAx(:,2),[0 1 0]);
                pcaAx = pcaAx*extraRotAd;

                 meanPosEye = meanPosEye*extraRotAd;
                 meanPosMirror = meanPosMirror*extraRotAd;
            else 
                error('axis may not be aligned');
            end

            %then correct rotation in x/y plane (roll)
            eyeFrontalAng = atan2(meanPosEye(2),meanPosEye(1));
            mirrorFrontalAng = atan2(meanPosMirror(2),meanPosMirror(1));
            %one axis should go close to zero, the other close to pi
            %correct such that ccw rotation is positive
            if eyeFrontalAng > pi/2
                eyeFrontalAng = eyeFrontalAng-pi;
            elseif eyeFrontalAng < -pi/2
                eyeFrontalAng = pi+eyeFrontalAng;
            end
            if mirrorFrontalAng > pi/2
                mirrorFrontalAng = mirrorFrontalAng-pi;
            elseif mirrorFrontalAng < -pi/2
                mirrorFrontalAng = pi+mirrorFrontalAng;
            end
            extraRotAd2 = vrrotvec2mat([0 0 1 (mirrorFrontalAng+eyeFrontalAng)/2]);
            meanPosEye = meanPosEye*extraRotAd2;
            meanPosMirror = meanPosMirror*extraRotAd2;

            %now update in xz plane (yaw)
            eyeFrontalAng = atan2(meanPosEye(3),meanPosEye(1));
            mirrorFrontalAng = atan2(meanPosMirror(3),meanPosMirror(1));
             if eyeFrontalAng > pi/2
                eyeFrontalAng = eyeFrontalAng-pi;
            elseif eyeFrontalAng < -pi/2
                eyeFrontalAng = pi+eyeFrontalAng;
            end
            if mirrorFrontalAng > pi/2
                mirrorFrontalAng = mirrorFrontalAng-pi;
            elseif mirrorFrontalAng < -pi/2
                mirrorFrontalAng = pi+mirrorFrontalAng;
            end
            extraRotAd3 = vrrotvec2mat([0 1 0 -(mirrorFrontalAng+eyeFrontalAng)/2]);

            extraRotAd = extraRotAd*extraRotAd2*extraRotAd3;
    %       update subscripts using 3x3 transform - just for plotting
            eyeFrontSurfSubsTransTemp = fullBeeData(i).EyeFrontSurfSubsTrans*extraRotAd;
            mirrorEyeSubsTemp = fullBeeData(i).mirrorEyeSubs*extraRotAd;

    %       make 4x4 transform to update other transofrms, note 4x4 transforms go forwards
            extraRotAd(4,4) = 1;
            %%% note that this is only releveant for calcs to align in world coord frame, so just applied at projection to sphere...

            %take again just for plotting
            pcaAx = pca([fullBeeData(i).EyeFrontSurfSubsTrans' fullBeeData(i).mirrorEyeSubs']');
        else
            extraRotAd = eye(4); %otherwise do nothing.
            eyeFrontSurfSubsTransTemp = fullBeeData(i).EyeFrontSurfSubsTrans;
            mirrorEyeSubsTemp = fullBeeData(i).mirrorEyeSubs;
        end
        if loadHeads
            %load head data
            fullBeeData(i).headStack = loadTIFStack( fullBeeData(i).headStackFolder, 1);
            fullBeeData(i).headVolSize = size(fullBeeData(i).headStack);

            if ~isnan(fullBeeData(i).AddEyeScale )
                %instead, reduce size of head
               headTransToUse = make_transformation_matrix([0 0 0], [0 0 0], [1 1 1]/fullBeeData(i).AddEyeScale)*fullBeeData(i).headTrans;
            else
               headTransToUse = fullBeeData(i).headTrans;
            end
            %get head subs
            [tX,tY,tZ] = ind2sub(fullBeeData(i).headVolSize, find( fullBeeData(i).headStack));
            fullBeeData(i).headStack = [];
            fullBeeData(i).HeadSurfSubsTrans = [[tX, tY, tZ]*fullBeeData(i).headVoxelSize-1 ones(size(tX,1),1)]*headTransToUse';
            fullBeeData(i).HeadSurfSubsTrans(:,4) = [];
            fullBeeData(i).HeadSurfSubsTrans = fullBeeData(i).HeadSurfSubsTrans-fullBeeData(i).HeadCenter;
        end

        %get facet D measurments
        origFacetInfo = load(fullBeeData(i).FacetSizeFile); %loads as size, x1, y1, z1, x2, y2, z2

        %sort to link neighbors in voxel coords
        counted = zeros(size(origFacetInfo,1),1);
        neihbors = zeros(size(origFacetInfo,1)/3,3);
        avgMeas = zeros(size(origFacetInfo,1)/3,1); %x, y
        avgCords = zeros(size(origFacetInfo,1)/3,3);

        if rem(size(origFacetInfo,1),3) ~= 0
             error(sprintf('incorrect number of faceter markers %s',fullBeeData(i).BeeID));
        end

        tempCoords = [mean(origFacetInfo(:,[2 5]),2) mean(origFacetInfo(:,[3 6]),2) mean(origFacetInfo(:,[4 7]),2)]/fullBeeData(i).VoxelSize;
        it = 1;
        while sum(counted == 0) 
            toGetInd = find(counted == 0);
            toGetInd = toGetInd(1);
            [sDists, sDistsI] = sort(sqrt((tempCoords(:,1) - tempCoords(toGetInd,1)).^2 + (tempCoords(:,2) - tempCoords(toGetInd,2)).^2 +...
                (tempCoords(:,3) - tempCoords(toGetInd,3)).^2));

            %1st one will be self, and two neighbors
            neihbors(it, :) = sDistsI(1:3);
            avgCords(it,:) = [mean(tempCoords(sDistsI(1:3),1)) mean(tempCoords(sDistsI(1:3),2)) mean(tempCoords(sDistsI(1:3),3))];

            avgMeas(it) = mean(origFacetInfo(sDistsI(1:3),1))/3; %measured across three hexagons originally
            %mark which were counted
            counted(sDistsI(1:3)) = 1;
            it = it + 1;
        end

        fullBeeData(i).facetSizes = avgMeas; %sizes should be correct scale

        [tX, tY, tZ] = ind2sub(fullBeeData(i).eyeVolSize, fullBeeData(i).EyeFrontSurfInds);

        %snap data onto eye surface in volume
        for j = 1:size(avgCords,1)
            %subtract 1 from from vol subscripts as avgCoords origin is zero from amira
            [~,closestI] = min(sqrt((avgCords(j,1) - tX).^2 + (avgCords(j,2) - tY).^2 + (avgCords(j,3) - tZ).^2));
            avgCords(j,:) = [tX(closestI), tY(closestI), tZ(closestI)];
        end

%       plot3(avgCords(:,1),avgCords(:,2),avgCords(:,3),'ro');

        fullBeeData(i).facetInds = sub2ind(fullBeeData(i).eyeVolSize, avgCords(:,1), avgCords(:,2), avgCords(:,3));
        fullBeeData(i).facetLocsTrans = [avgCords*fullBeeData(i).VoxelSize ones(size(avgCords,1),1)]*origEyeVolTrans';
        fullBeeData(i).facetLocsTrans(:,4) = [];
        fullBeeData(i).facetLocsTrans =fullBeeData(i).facetLocsTrans-fullBeeData(i).HeadCenter;

        %%%
        %try to decompose facet orientations into horizontal and world components
        if plotGroupHexs
            if ishandle(fullBeeData(i).OrigHexF)
               close(fullBeeData(i).OrigHexF);  
            end
            fullBeeData(i).OrigHexF = figure;
            
            subplot(1,2,1); hold on
            xlabel('x'); ylabel('y'); zlabel('z'); axis equal 
            ax = gca; ax.Clipping = 'off'; 

            subplot(1,2,2); hold on
            xlabel('x'); ylabel('y'); zlabel('z'); axis equal
            ax = gca; ax.Clipping = 'off'; 
        end

        lensSurfCenter = mean(fullBeeData(i).EyeFrontSurfSubsTrans);

        colsF = jet(100);
        minF = min(origFacetInfo(:,1)); maxF = max(origFacetInfo(:,1));
        facetMeasCol = (origFacetInfo(:,1)-minF);
        facetMeasCol = round(facetMeasCol/(maxF-minF)*99+1);

        if ~useDirectHexAxes
            fullBeeData(i).distAndAng2Point = zeros(length(fullBeeData(i).facetInds),2);
            fullBeeData(i).distAndAng2Para = zeros(length(fullBeeData(i).facetInds),2);
            fullBeeData(i).angBetweenPerpAndPara = zeros(length(fullBeeData(i).facetInds),2);
        else
            fullBeeData(i).HexagonCoords = zeros(length(fullBeeData(i).facetInds),6,2);
            fullBeeData(i).HexLinkOrder = zeros(length(fullBeeData(i).facetInds),6);
        end
        fullBeeData(i).hexArea = zeros(length(fullBeeData(i).facetInds),1);
        toDelete = [];
        
        for j = 1:size(avgCords,1)
           % calculate normals to align in world
           [tempNorm] = normalToSurface_withCurvature( fullBeeData(i).facetLocsTrans(j,:),...
                    fullBeeData(i).EyeFrontSurfSubsTrans, lensSurfCenter, 1, surfaceFitRng,[xFlatSets, yFlatSets, zFlatSets],fullBeeData(i).borderInds, htStep);
           el = acos(tempNorm(2));
           az = atan2(tempNorm(3),tempNorm(1));

           %rotate lines to world coordinate frames
           rotLine1 = [origFacetInfo(neihbors(j, 1),[2 5])', origFacetInfo(neihbors(j, 1),[3 6])', origFacetInfo(neihbors(j, 1),[4 7])', ones(2,1)]*origEyeVolTrans'-[fullBeeData(i).HeadCenter 0];
           rotLine2 = [origFacetInfo(neihbors(j, 2),[2 5])', origFacetInfo(neihbors(j, 2),[3 6])', origFacetInfo(neihbors(j, 2),[4 7])', ones(2,1)]*origEyeVolTrans'-[fullBeeData(i).HeadCenter 0];
           rotLine3 = [origFacetInfo(neihbors(j, 3),[2 5])', origFacetInfo(neihbors(j, 3),[3 6])', origFacetInfo(neihbors(j, 3),[4 7])', ones(2,1)]*origEyeVolTrans'-[fullBeeData(i).HeadCenter 0];
           combLines = [rotLine1', rotLine2', rotLine3']';

           if plotGroupHexs
               figure(fullBeeData(i).OrigHexF); %plot main
               subplot(1,2,1);
               text(fullBeeData(i).facetLocsTrans(j,1),fullBeeData(i).facetLocsTrans(j,2),fullBeeData(i).facetLocsTrans(j,3),sprintf('%i',j));
               line(combLines(1:2,1), combLines(1:2,2), combLines(1:2,3),'linewidth',1,'color',colsF(facetMeasCol(neihbors(j, 1)),:),'linestyle',':');
               line(combLines(3:4,1), combLines(3:4,2), combLines(3:4,3),'linewidth',1,'color',colsF(facetMeasCol(neihbors(j, 2)),:),'linestyle',':');
               line(combLines(5:6,1), combLines(5:6,2), combLines(5:6,3),'linewidth',1,'color',colsF(facetMeasCol(neihbors(j, 3)),:),'linestyle',':');
           end
           combLines = (combLines(:,1:3)-fullBeeData(i).facetLocsTrans(j,:))*vrrotvec2mat([0 1 0 -az])*vrrotvec2mat([0 0 1 -el+pi/2])+fullBeeData(i).facetLocsTrans(j,:); %

           if plotGroupHexs
               %plot rotated
               subplot(1,2,2);
               text(fullBeeData(i).facetLocsTrans(j,1),fullBeeData(i).facetLocsTrans(j,2),fullBeeData(i).facetLocsTrans(j,3),sprintf('%i',j));
               line(combLines(1:2,1), combLines(1:2,2), combLines(1:2,3),'linewidth',1,'color',colsF(facetMeasCol(neihbors(j, 1)),:),'linestyle',':');
               line(combLines(3:4,1), combLines(3:4,2), combLines(3:4,3),'linewidth',1,'color',colsF(facetMeasCol(neihbors(j, 2)),:),'linestyle',':');
               line(combLines(5:6,1), combLines(5:6,2), combLines(5:6,3),'linewidth',1,'color',colsF(facetMeasCol(neihbors(j, 3)),:),'linestyle',':');
           end

           %do indivdually shrink indvidually for each line before
           m1 = mean(combLines(1:2,:));
           combLines(1:2,:) = (combLines(1:2,:)-m1)*2/3+m1;
           m2 = mean(combLines(3:4,:));
           combLines(3:4,:) = (combLines(3:4,:)-m2)*2/3+m2;
           m3 = mean(combLines(5:6,:));
           combLines(5:6,:) = (combLines(5:6,:)-m3)*2/3+m3;

           meanLines = mean(combLines);
           combLines = combLines-meanLines;

           if plotGroupHexs
                plot3(meanLines(1), meanLines(2), meanLines(3),'k.');
           end  

           %get voroni diagram from lattice pts
           [vx,vy] = voronoi([combLines(:,2)' 0],[combLines(:,3)' 0]);
           vx = vx(:); vy = vy(:);
           %actual verts are included three times
           %not elegant, but first remove one set from unique get far tail points
           [~,IToD] = unique([vx,vy], 'rows');
           vx(IToD) = []; vy(IToD) = [];
           %then take unique again to get one set of verts
           [facetVerts,IToD] = unique([vx,vy], 'rows');

           %here x is el(y), y is az(z)

           if plotIndHexs
               figure;  hold on; axis equal %note these are not deleted per iteration
               voronoi([combLines(:,2)' 0],[combLines(:,3)' 0])
               plot([combLines(:,2)' 0],[combLines(:,3)' 0],'r.');
               plot(facetVerts(:,1),facetVerts(:,2),'rx');
           end

          if size(facetVerts,1) ~= 6
              warning(sprintf('degenerate facet with %i sides ignored %i', size(facetVerts,1)));
              toDelete = [toDelete j]; %need to delete from everything!
          else
               %get angular order of hexagon points - ccw
               angsVInW = zeros(6,1);
               for k = 1:6
                   angsVInW(k) = atan2(facetVerts(k,2),facetVerts(k,1));
               end
               [~, linkOrderV] = sort(angsVInW);

               %this is mainly for tesselation as useful for plotting and getting linkages
               %now find corner angles and side lengths 
               %using schema from http://mathstat.slu.edu/escher/index.php/Tessellations_by_Polygons#cite_note-2 
               angsV = zeros(6,1);
               distsV = zeros(6,1);
               angL = {'A','B','C','D','E','F'};
               disL = {'a','b','c','d','e','f'};
               latLink = zeros(6,1); 
               for k = 1:6
                   ptSet1 = [k-1, k];%in CCW direction
                   ptSet2 = [k+1, k];%in CW direction
                   ptSet1(ptSet1 < 1) = ptSet1(ptSet1 < 1) + 6;
                   ptSet2(ptSet2 > 6) = ptSet2(ptSet2 > 6) - 6;
                   relPt1 = [facetVerts(linkOrderV(ptSet1(1)),1)-facetVerts(linkOrderV(ptSet1(2)),1), facetVerts(linkOrderV(ptSet1(1)),2)-facetVerts(linkOrderV(ptSet1(2)),2)];
                   relPt2 = [facetVerts(linkOrderV(ptSet2(1)),1)-facetVerts(linkOrderV(ptSet2(2)),1), facetVerts(linkOrderV(ptSet2(1)),2)-facetVerts(linkOrderV(ptSet2(2)),2)];

                   %https://math.stackexchange.com/questions/59/calculating-an-angle-from-2-points-in-space
                   angsV(k) = acos((relPt1(1)*relPt2(1)+relPt1(2)*relPt2(2))/sqrt((relPt1(1)^2+relPt1(2)^2)*(relPt2(1)^2+relPt2(2)^2))); 

                   distsV(k) = sqrt((facetVerts(linkOrderV(ptSet2(1)),1)-facetVerts(linkOrderV(ptSet2(2)),1))^2+(facetVerts(linkOrderV(ptSet2(1)),2)-facetVerts(linkOrderV(ptSet2(2)),2))^2);

                   %get linkage to lattice
                   sideCent = [mean(facetVerts(linkOrderV(ptSet2),1)), mean(facetVerts(linkOrderV(ptSet2),2))];
                   angDiff = zeros(6,1);
                   for m = 1:6
                      angDiff(m) = acos((combLines(m,2)*sideCent(1)+combLines(m,3)*sideCent(2))/sqrt((combLines(m,2)^2+combLines(m,3)^2)*(sideCent(1)^2+sideCent(2)^2)));
                   end
                   [~,latLink(k)] = min(angDiff);

                   %find area of this triangle from heronds formula
                   lenA = sqrt((0-facetVerts(linkOrderV(ptSet2(1)),1))^2+(0-facetVerts(linkOrderV(ptSet2(1)),2))^2);
                   lenB = sqrt((0-facetVerts(linkOrderV(ptSet2(2)),1))^2+(0-facetVerts(linkOrderV(ptSet2(2)),2))^2);
                   lenC = sqrt((facetVerts(linkOrderV(ptSet2(1)),1)-facetVerts(linkOrderV(ptSet2(2)),1))^2+(facetVerts(linkOrderV(ptSet2(1)),2)-facetVerts(linkOrderV(ptSet2(2)),2))^2);
                   p = (lenA+lenB+lenC)/2;
                   triArea = sqrt(p*(p-lenA)*(p-lenB)*(p-lenC));

                   fullBeeData(i).hexArea(j) = fullBeeData(i).hexArea(j) + triArea;

                   if plotIndHexs
                       text(facetVerts(linkOrderV(k),1), facetVerts(linkOrderV(k),2), sprintf('%s : %.1f', angL{k}, angsV(k)/pi*180));
                       text(mean(facetVerts(linkOrderV(ptSet2),1)), mean(facetVerts(linkOrderV(ptSet2),2)), sprintf('%s : %.1f', disL{k}, distsV(k)));
                       text(combLines(latLink(k),2), combLines(latLink(k),3), sprintf('%s', disL{k}));

                       meanC1 = (0 + facetVerts(linkOrderV(ptSet2(1)),1) + facetVerts(linkOrderV(ptSet2(2)),1))/3;
                       meanC2 = (0 + facetVerts(linkOrderV(ptSet2(1)),2) + facetVerts(linkOrderV(ptSet2(2)),2))/3;
                       text(meanC1, meanC2, sprintf('%.1f', triArea));
    %                      text(facetVerts(linkOrderV(k),1), facetVerts(linkOrderV(k),2), sprintf('%s : %.1f : %.1f', angL{k}, angsVInW(linkOrderV(k))/pi*180, angsV(k)/pi*180));
                   end
               end

               if useDirectHexAxes
                   fullBeeData(i).HexagonCoords(j,:,1) = facetVerts(:,1);
                   fullBeeData(i).HexagonCoords(j,:,2) = facetVerts(:,2);
                   fullBeeData(i).HexLinkOrder(j,:) = linkOrderV;
               else
                  %calculations from parallel/perp lines are the same for all types

                   %note explicit divide by 2 to go from lattice spacing to d values
                   paraLine = paraLine/2;
                   perpLine = perpLine/2; %note length of perp/point line extends further than hexagon, to next lattice line

                   %flip as require to get correct line up
                   %this ensures that perp points up (pos z) or right (pos y)
                        %hence para points left (neg y) or up (pos z) 
                   if ppMatch(minAngI)
                       %if para chosen and points vert, perp should point right
                       if axisMatch(minAngI) & perpLine(1) < 0
                           pointTestLine = -pointTestLine;
                           paraLine = -paraLine;
                           perpLine = -perpLine;
                           %or if it points horiz, perp should point up
                       elseif ~axisMatch(minAngI) & perpLine(2) < 0
                           pointTestLine = -pointTestLine;
                           paraLine = -paraLine;
                           perpLine = -perpLine;
                       end
                   else
                       %if point chosen, point up if on vert, or right if on horiz
                       if axisMatch(minAngI) & perpLine(2) < 0
                           pointTestLine = -pointTestLine;
                           paraLine = -paraLine;
                           perpLine = -perpLine;
                       elseif ~axisMatch(minAngI) & perpLine(1) < 0
                           pointTestLine = -pointTestLine;
                           paraLine = -paraLine;
                           perpLine = -perpLine;
                       end
                   end

                   %angle between vectors
                   axisAng = acos((paraLine(1)*perpLine(1)+paraLine(2)*perpLine(2))/sqrt((paraLine(1)^2+paraLine(2)^2)*(perpLine(1)^2+perpLine(2)^2)));
                   %we always want the acute angle
                   if axisAng > pi/2
                       axisAng = pi-axisAng;
                   end

                   %note angles here aren't rotated... problem?
                   paraAng = atan2(paraLine(2), paraLine(1));
                   perpAng = atan2(perpLine(2), perpLine(1));
                   paraDist = sqrt(paraLine(1)^2+paraLine(2)^2);
                   perpDist = sqrt(perpLine(1)^2+perpLine(2)^2);
                   fullLen = sqrt((paraLine(1)+perpLine(1))^2+(paraLine(2)+perpLine(2))^2);

                   %%%not 100% sure this is the way to do things as both end up being quite increased if not aligned to world.
                   %is rotating d to world and then calculating r in world the same as:
                   %calculating r as rotated and then calcuting io in world?

                   %make paralellogram to get world coordinates - each is max across pt range
                   parrelogram = [0 0; paraLine; perpLine; perpLine+paraLine];
                   horizLen = max(parrelogram(:,1))-min(parrelogram(:,1));
                   vertLen = max(parrelogram(:,2))-min(parrelogram(:,2));
                   worldLen = sqrt(horizLen^2+vertLen^2); %this can be unrepresentative...

                   %store things
                   fullBeeData(i).distAndAng2Point(j,:) = [perpDist, perpAng];
                   fullBeeData(i).distAndAng2Para(j,:) = [paraDist, paraAng];

                   if plotIndHexs
                       %rotate so vector sum always pointing in positive x
                       if paraLine(1)+perpLine(1) < 0
                           paraLine = -paraLine;
                           perpLine = -perpLine;
                       end

                       %plot lines and write params 
                       line([0 paraLine(1)], [0 paraLine(2)], 'color', 'g');
                       line([0 perpLine(1)], [0 perpLine(2)], 'color', 'm');
                       line([0 (perpLine(1)+paraLine(1))], [0 (perpLine(2)+paraLine(2))], 'color', 'c');
                       line([0 horizLen], [0 0], 'color','k');
                       line([0 0], [0 vertLen], 'color','k');

                       text((perpLine(1)+paraLine(1)), (perpLine(2)+paraLine(2)), sprintf('V %.1f : M %.1f : W%.1f',fullLen, fullBeeData(i).facetSizes(j), worldLen));
                       text(paraLine(1), paraLine(2), sprintf('%.1f',paraDist));
                       text(perpLine(1), perpLine(2), sprintf('%.1f',perpDist));
                       text(0, 0, sprintf('A %.1f : R %.1f',axisAng/pi*180, paraDist/perpDist));
                   end
               end

               if plotGroupHexs
                   figure(fullBeeData(i).OrigHexF); %plot main

                   %now we should plot hexagon on both % so we need to rotate it back
                   %and then plot major axis used...
                   linkOrderVT = [linkOrderV' linkOrderV(1)]';

                   %plot unrotated
                   subplot(1,2,2);
                   for k = 1:6
                        line([0 0]+meanLines(1), [facetVerts(linkOrderVT(k),1) facetVerts(linkOrderVT(k+1),1)]+meanLines(2),...
                            [facetVerts(linkOrderVT(k),2) facetVerts(linkOrderVT(k+1),2)]+meanLines(3),'linewidth',1,'color','k','linestyle','--');
                   end

                   %extra bits 
                   if ~useDirectHexAxes 
                       line([0 0]+meanLines(1),[0 pointTestLine(1)]+meanLines(2), [0 pointTestLine(2)]+meanLines(3), 'color', 'r','linewidth',2);
                       line([0 0]+meanLines(1),[0 paraLine(1)]+meanLines(2), [0 paraLine(2)]+meanLines(3), 'color', 'g');
                       if ppMatch(minAngI)
                           plot3(meanLines(1),paraLine(1)+meanLines(2), paraLine(2)+meanLines(3), 'g*');
                       end
                       line([0 0]+meanLines(1),[0 perpLine(1)]+meanLines(2), [0 perpLine(2)]+meanLines(3), 'color', 'm');
                       if ~ppMatch(minAngI)
                           plot3(meanLines(1),perpLine(1)+meanLines(2), perpLine(2)+meanLines(3), 'm*');
                       end
                       line([0 0]+meanLines(1),[0 (perpLine(1)+paraLine(1))]+meanLines(2), [0 (perpLine(2)+paraLine(2))]+meanLines(3), 'color', 'c');
                       line([0 0]+meanLines(1),[0 horizLen]+meanLines(2), [0 0]+meanLines(3), 'color','k');
                       line([0 0]+meanLines(1),[0 0]+meanLines(2), [0 vertLen]+meanLines(3), 'color','k');    

                       text(meanLines(1), paraLine(1)+meanLines(2), paraLine(2)+meanLines(3), sprintf('%.0f',paraAng/pi*180));
                       text(meanLines(1), perpLine(1)+meanLines(2), perpLine(2)+meanLines(3), sprintf('%.0f',perpAng/pi*180));
                     end

                   %plot rotated (back)
                   facetVerts = [zeros(6,1), facetVerts]+meanLines;
                   facetVerts = (facetVerts-fullBeeData(i).facetLocsTrans(j,:))*vrrotvec2mat([0 0 1 -(-el+pi/2)])*vrrotvec2mat([0 1 0 az])+fullBeeData(i).facetLocsTrans(j,:);

                   subplot(1,2,1);
                   title(sprintf('%s',fullBeeData(i).BeeID));

                   for k = 1:6
                        line([facetVerts(linkOrderVT(k),1) facetVerts(linkOrderVT(k+1),1)], [facetVerts(linkOrderVT(k),2) facetVerts(linkOrderVT(k+1),2)],...
                            [facetVerts(linkOrderVT(k),3) facetVerts(linkOrderVT(k+1),3)],'linewidth',1,'color','k','linestyle','--');
                   end

                   %extra bits rotated
                   if ~useDirectHexAxes 
                       tempLines = [zeros(1,3); 0, pointTestLine]+meanLines;
                       tempLines = (tempLines-fullBeeData(i).facetLocsTrans(j,:))*vrrotvec2mat([0 0 1 -(-el+pi/2)])*vrrotvec2mat([0 1 0 az])+fullBeeData(i).facetLocsTrans(j,:);
                       line(tempLines(:,1), tempLines(:,2),tempLines(:,3), 'color','r','linewidth',2);
                       tempLines = [zeros(1,3); 0, paraLine]+meanLines;
                       tempLines = (tempLines-fullBeeData(i).facetLocsTrans(j,:))*vrrotvec2mat([0 0 1 -(-el+pi/2)])*vrrotvec2mat([0 1 0 az])+fullBeeData(i).facetLocsTrans(j,:);
                       line(tempLines(:,1), tempLines(:,2),tempLines(:,3), 'color','g');
                       if ppMatch(minAngI)
                           plot3(tempLines(2,1), tempLines(2,2),tempLines(2,3), 'g*');
                       end
                       tempLines = [zeros(1,3); 0, perpLine]+meanLines;
                       tempLines = (tempLines-fullBeeData(i).facetLocsTrans(j,:))*vrrotvec2mat([0 0 1 -(-el+pi/2)])*vrrotvec2mat([0 1 0 az])+fullBeeData(i).facetLocsTrans(j,:);
                       line(tempLines(:,1), tempLines(:,2),tempLines(:,3), 'color','m');
                       if ~ppMatch(minAngI)
                           plot3(tempLines(2,1), tempLines(2,2),tempLines(2,3), 'm*');
                       end
                       tempLines = [zeros(1,3); 0, (perpLine+paraLine)]+meanLines;
                       tempLines = (tempLines-fullBeeData(i).facetLocsTrans(j,:))*vrrotvec2mat([0 0 1 -(-el+pi/2)])*vrrotvec2mat([0 1 0 az])+fullBeeData(i).facetLocsTrans(j,:);
                       line(tempLines(:,1), tempLines(:,2),tempLines(:,3), 'color','c');
                       tempLines = [zeros(1,3); 0, horizLen, 0]+meanLines;
                       tempLines = (tempLines-fullBeeData(i).facetLocsTrans(j,:))*vrrotvec2mat([0 0 1 -(-el+pi/2)])*vrrotvec2mat([0 1 0 az])+fullBeeData(i).facetLocsTrans(j,:);
                       line(tempLines(:,1), tempLines(:,2),tempLines(:,3), 'color','k');
                       tempLines = [zeros(1,3); 0, 0, vertLen]+meanLines;
                       tempLines = (tempLines-fullBeeData(i).facetLocsTrans(j,:))*vrrotvec2mat([0 0 1 -(-el+pi/2)])*vrrotvec2mat([0 1 0 az])+fullBeeData(i).facetLocsTrans(j,:);
                       line(tempLines(:,1), tempLines(:,2),tempLines(:,3), 'color','k');
                 end
               end
            end
        end

        %remove degenerate facets
        if ~isempty(toDelete)
           fullBeeData(i).facetSizes(toDelete) = [];
           fullBeeData(i).facetInds(toDelete) = [];
           fullBeeData(i).facetLocsTrans(toDelete,:) = [];
           fullBeeData(i).hexArea(toDelete) = [];
            if ~useDirectHexAxes
                fullBeeData(i).distAndAng2Point(toDelete,:) = [];
                fullBeeData(i).distAndAng2Para(toDelete,:) = [];
                fullBeeData(i).angBetweenPerpAndPara(toDelete,:) = [];
            else
                fullBeeData(i).HexagonCoords(toDelete,:,:) = [];
                fullBeeData(i).HexLinkOrder(toDelete,:,:) = [];
            end  
        end
        
        if testPlots
            if ishandle(fullBeeData(i).headEyeF)
                close(fullBeeData(i).headEyeF);
            end
            fullBeeData(i).headEyeF = figure;
            
            subplot(2,1,1); axis equal; hold on; view([0, -90])
            plot3(eyeFrontSurfSubsTransTemp(:,1), eyeFrontSurfSubsTransTemp(:,2), eyeFrontSurfSubsTransTemp(:,3),'r.');
            plot3(mirrorEyeSubsTemp(:,1), mirrorEyeSubsTemp(:,2), mirrorEyeSubsTemp(:,3),'b.');

            %%%head is a bit slow because of many points to plot
            if loadHeads
                toPInds = 1:100:size(fullBeeData(i).HeadSurfSubsTrans,1);
                plot3(fullBeeData(i).HeadSurfSubsTrans(toPInds,1), fullBeeData(i).HeadSurfSubsTrans(toPInds,2), fullBeeData(i).HeadSurfSubsTrans(toPInds,3),'k.');
                if removeHeadAfterPlot
                    fullBeeData(i).HeadSurfSubsTrans = [];
                end
            end

            title(sprintf('%s : check eyes are corretly aligned (on head)',fullBeeData(i).BeeID));

            line([0 pcaAx(1,1)]*1000, [0 pcaAx(2,1)]*1000, [0 pcaAx(3,1)]*1000, 'color','g');
            line([0 pcaAx(1,2)]*1000, [0 pcaAx(2,2)]*1000, [0 pcaAx(3,2)]*1000, 'color','b');
            line([0 pcaAx(1,3)]*1000, [0 pcaAx(2,3)]*1000, [0 pcaAx(3,3)]*1000, 'color','r');
            
            subplot(2,1,2); axis equal; hold on; view([0, -90])
            tempCoords = [mean(origFacetInfo(:,[2 5]),2) mean(origFacetInfo(:,[3 6]),2) mean(origFacetInfo(:,[4 7]),2)]/fullBeeData(i).VoxelSize;
            plot3(tempCoords(:,1), tempCoords(:,2),tempCoords(:,3),'go','markersize',5,'MarkerFaceColor','g');
            [tX, tY, tZ] = ind2sub(fullBeeData(i).eyeVolSize, fullBeeData(i).EyeFrontSurfInds);
            plot3(tX, tY, tZ,'r.');
            title(sprintf('facet measurments on original volume : check they are close : num %i',size(origFacetInfo,1)/3));
        end

        %% get sample points and build conencted neighborhood
        %then get interpolation points on eye

        %calculate by picking a point and then masking out surronding points
        fullBeeData(i).interpPoints = zeros(length(fullBeeData(i).EyeFrontSurfInds),1);
        if stepInFromBorders
            %firstly mask out borders, by stepping in
            for j = 1:length(fullBeeData(i).borderInds)
                indsToCut = find(sqrt((fullBeeData(i).EyeFrontSurfSubsTrans(:,1)-fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).borderInds(j),1)).^2+...
                            (fullBeeData(i).EyeFrontSurfSubsTrans(:,2)-fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).borderInds(j),2)).^2+...
                            (fullBeeData(i).EyeFrontSurfSubsTrans(:,3)-fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).borderInds(j),3)).^2) < stepInInterval);
                fullBeeData(i).interpPoints(indsToCut) = -1;
            end
        end
        while sum(fullBeeData(i).interpPoints == 0)
            indsToTest = find(fullBeeData(i).interpPoints == 0);
            %work down from top.... which if from negative Y in usual coord frame
            [~,maxI] = min(fullBeeData(i).EyeFrontSurfSubsTrans(indsToTest,2));

            %find points within range
            indsToCut = find(sqrt((fullBeeData(i).EyeFrontSurfSubsTrans(indsToTest,1)-fullBeeData(i).EyeFrontSurfSubsTrans(indsToTest(maxI),1)).^2+...
                        (fullBeeData(i).EyeFrontSurfSubsTrans(indsToTest,2)-fullBeeData(i).EyeFrontSurfSubsTrans(indsToTest(maxI),2)).^2+...
                        (fullBeeData(i).EyeFrontSurfSubsTrans(indsToTest,3)-fullBeeData(i).EyeFrontSurfSubsTrans(indsToTest(maxI),3)).^2) < measurmentIntervals);
            fullBeeData(i).interpPoints(indsToTest(indsToCut)) = -1;
            fullBeeData(i).interpPoints(indsToTest(maxI)) = 1;
        end    
        fullBeeData(i).interpPoints = find(fullBeeData(i).interpPoints == 1);    
        
        if testPlots
            if ishandle(fullBeeData(i).SamplingF)
                close(fullBeeData(i).SamplingF)
            end
            fullBeeData(i).SamplingF = figure; 
            
            subplot(1,2,1); axis equal; hold on; view([0, -90]); title(sprintf('%s',fullBeeData(i).BeeID));
            plot3(fullBeeData(i).EyeFrontSurfSubsTrans(:,1), fullBeeData(i).EyeFrontSurfSubsTrans(:,2), fullBeeData(i).EyeFrontSurfSubsTrans(:,3),'r.');
            plot3(fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,1), fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,2),...
                    fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,3), 'go');
                
            subplot(1,2,2); axis equal; hold on; view([0, -90]);
            plot3(LOx(fullBeeData(i).interpPoints), LOy(fullBeeData(i).interpPoints),LOz(fullBeeData(i).interpPoints), 'go');
        end

        %connectedness section - this is not really used now that normals are calcualted in following section
        %note that nearest point networks are often arranged in hexagons, but this has no relation to the hexagonal facet lens pattern

        %determine min dists
        minDist = zeros(length(fullBeeData(i).interpPoints),1);
        for j = 1:length(fullBeeData(i).interpPoints)
            dists = sort(sqrt((fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),1) - fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,1)).^2 + ...
                (fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),2) - fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,2)).^2 +...
                (fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),3) - fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,3)).^2));
            minDist(j) = dists(2);
        end
        minDistL = min(minDist)*1.5;

        %calculate connectedness of each point - from amber gnat/scallop
        fullBeeData(i).samplePointConnections = cell(length(fullBeeData(i).interpPoints),1);
        for j = 1:length(fullBeeData(i).interpPoints)
            %get sorted list of dists to other sample points
            distInd = find(sqrt((fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),1) - fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,1)).^2 + ...
                (fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),2) - fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,2)).^2 +...
                (fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),3) - fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,3)).^2)< minDistL);
            distInd(distInd == j) = []; %remove self
            fullBeeData(i).samplePointConnections{j} = distInd;
        end

        %check that you are in all your connected connections
        for j = 1:length(fullBeeData(i).interpPoints)
            tempConns = fullBeeData(i).samplePointConnections{j};
            for k = 1:length(tempConns)
                fullBeeData(i).samplePointConnections{tempConns(k)} = unique([fullBeeData(i).samplePointConnections{tempConns(k)}' j]');
            end
        end

        %find the borders...
        fullBeeData(i).numConnections = zeros(length(fullBeeData(i).interpPoints),1);
        for j = 1:length(fullBeeData(i).interpPoints)
            fullBeeData(i).numConnections(j) = size(fullBeeData(i).samplePointConnections{j},1);
            %remove furthest if more than 6 conns
            if fullBeeData(i).numConnections(j) > 6
                tempConns = fullBeeData(i).samplePointConnections{j};
                [~,SInd] = sort(sqrt((fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),1) - fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(fullBeeData(i).samplePointConnections{j}),1)).^2 + ...
                        (fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),2) - fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(fullBeeData(i).samplePointConnections{j}),2)).^2 +...
                        (fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),3) - fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(fullBeeData(i).samplePointConnections{j}),3)).^2));

                %self should never be in here.    
                fullBeeData(i).samplePointConnections{j} = tempConns(SInd(1:6));
                fullBeeData(i).numConnections(j) = 6;

                %delete from former friends
                for k = 7:length(SInd)
                    testCons = fullBeeData(i).samplePointConnections{tempConns(SInd(k))};
                    testCons(testCons == j) = [];
                    fullBeeData(i).samplePointConnections{tempConns(SInd(k))} = testCons;
                    fullBeeData(i).numConnections(tempConns(SInd(k))) = length(testCons);
                end
            end
        end

        %get normals
        fullBeeData(i).normals = zeros(length(fullBeeData(i).interpPoints),3);
        lensSurfCenter = mean(fullBeeData(i).EyeFrontSurfSubsTrans);

        %for calcuating intersects
        fullBeeData(i).lensThickness = zeros(length(fullBeeData(i).interpPoints),1)*NaN;
        fullBeeData(i).CCThickness = zeros(length(fullBeeData(i).interpPoints),1)*NaN;
        fullBeeData(i).RetThickness = zeros(length(fullBeeData(i).interpPoints),1)*NaN;
        
        if useReverseNormal
            grid3D.nx = fullBeeData(i).eyeVolSize(1); grid3D.ny = fullBeeData(i).eyeVolSize(2); grid3D.nz = fullBeeData(i).eyeVolSize(3);
            grid3D.minBound = [1 1 1]';
            grid3D.maxBound = [fullBeeData(i).eyeVolSize(1), fullBeeData(i).eyeVolSize(2), fullBeeData(i).eyeVolSize(3)]';
            jumpForward = 3;
            
            lensInnerInds = find(fullBeeData(i).LabelsStack == fullBeeData(i).StackLabels.LensInnerS);
            retinaOuterInds = find(fullBeeData(i).LabelsStack == fullBeeData(i).StackLabels.RetinaOuterS);
            laminaOuterInds = find(fullBeeData(i).LabelsStack == fullBeeData(i).StackLabels.LaminaC);
        end
        [LIx,LIy,LIz] = ind2sub(fullBeeData(i).eyeVolSize, find(fullBeeData(i).LabelsStack == fullBeeData(i).StackLabels.LensInnerS));
        [ROx,ROy,ROz] = ind2sub(fullBeeData(i).eyeVolSize, find(fullBeeData(i).LabelsStack == fullBeeData(i).StackLabels.RetinaOuterS));
        [Lax,Lay,Laz] = ind2sub(fullBeeData(i).eyeVolSize, find(fullBeeData(i).LabelsStack == fullBeeData(i).StackLabels.LaminaC));
        
        for j = 1:length(fullBeeData(i).interpPoints)
            %%%note that calculating curvature with this method is not implemented correctly
            fullBeeData(i).normals(j,:) = normalToSurface_withCurvature( fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),:),...
                    fullBeeData(i).EyeFrontSurfSubsTrans, lensSurfCenter, 1, surfaceFitRng, [xFlatSets, yFlatSets, zFlatSets], fullBeeData(i).borderInds, htStep);

            lens_CP = 0;
            CC_CP = 0;
            ret_CP = 0;
            
            if useReverseNormal
                %draw line back in volume - transform it back first
                transNormal = (origEyeVolTrans^-1*[fullBeeData(i).normals(j,:) 0]')';

                [xCords, yCords, zCords] = amanatidesWooAlgorithm_preAlloc([LOx(fullBeeData(i).interpPoints(j)) LOy(fullBeeData(i).interpPoints(j)) LOz(fullBeeData(i).interpPoints(j))],...
                        -1*transNormal(1:3), grid3D, 0);
                inds = sub2ind(fullBeeData(i).eyeVolSize, xCords, yCords, zCords);
                backLensInteresects = find(fullBeeData(i).LabelsStack(inds) == fullBeeData(i).StackLabels.LensInnerS);
                backLensInteresects(backLensInteresects < jumpForward) = [];
                
                %%%another option to get missing inds would be to make distmap using vector as seed points, 
                %%%and then pick closest surface ind to that
                
                if ~isempty(backLensInteresects)
                   %take top point normal goes through
                   lensT = sqrt((LOx(fullBeeData(i).interpPoints(j))-xCords(backLensInteresects(1))).^2+...
                        (LOy(fullBeeData(i).interpPoints(j))-yCords(backLensInteresects(1))).^2+(LOz(fullBeeData(i).interpPoints(j))-zCords(backLensInteresects(1))).^2);
                   backLensInd = find(lensInnerInds == inds(backLensInteresects(1)));
                end
                if isempty(backLensInteresects) | lensT > lensT_Thresh/fullBeeData(i).VoxelSize
                   %instead, find closest point on entire surface
                   [lensT, backLensInd] = min(sqrt((LOx(fullBeeData(i).interpPoints(j))-LIx).^2+(LOy(fullBeeData(i).interpPoints(j))-LIy).^2+(LOz(fullBeeData(i).interpPoints(j))-LIz).^2));
                   %then get new normal based on this vector
                   transNormal = [LIx(backLensInd)-LOx(fullBeeData(i).interpPoints(j)), LIy(backLensInd)-LOy(fullBeeData(i).interpPoints(j)), LIz(backLensInd)-LOz(fullBeeData(i).interpPoints(j))];
                   [xCords, yCords, zCords] = amanatidesWooAlgorithm_preAlloc([LOx(fullBeeData(i).interpPoints(j)) LOy(fullBeeData(i).interpPoints(j)) LOz(fullBeeData(i).interpPoints(j))],...
                            transNormal(1:3)/norm(transNormal(1:3)), grid3D, 0);
                   inds = sub2ind(fullBeeData(i).eyeVolSize, xCords, yCords, zCords);
                   lens_CP = 1;
                end
                fullBeeData(i).lensThickness(j) = lensT*fullBeeData(i).VoxelSize;
                
                frontRetInteresects = find(fullBeeData(i).LabelsStack(inds) == fullBeeData(i).StackLabels.RetinaOuterS);
                frontRetInteresects(frontRetInteresects < jumpForward) = [];
                
                if ~isempty(frontRetInteresects)
                   %take top point normal goes through
                   CCT = sqrt((LIx(backLensInd)-xCords(frontRetInteresects(1))).^2+...
                        (LIy(backLensInd)-yCords(frontRetInteresects(1))).^2+(LIz(backLensInd)-zCords(frontRetInteresects(1))).^2);
                   backCCInd = find(retinaOuterInds == inds(frontRetInteresects(1)));
                end
                if isempty(frontRetInteresects) | CCT > CCT_Thresh/fullBeeData(i).VoxelSize
                   %instead, find closest point on entire surface'
                   [CCT, backCCInd] = min(sqrt((LIx(backLensInd)-ROx).^2+(LIy(backLensInd)-ROy).^2+(LIz(backLensInd)-ROz).^2));
                   %keep using same normal, if recalculated, should then be used to find lens intersect prior to this as well
                   CC_CP = 1;
                end
                fullBeeData(i).CCThickness(j) = CCT*fullBeeData(i).VoxelSize;
                
                %lamina ints
                laminaInteresects = find(fullBeeData(i).LabelsStack(inds) == fullBeeData(i).StackLabels.LaminaC);
                laminaInteresects(laminaInteresects < jumpForward) = [];
                
                if ~isempty(laminaInteresects)
                   %take top point normal goes through
                   retT2 = sqrt((ROx(backCCInd)-xCords(laminaInteresects(1))).^2+...
                        (ROy(backCCInd)-yCords(laminaInteresects(1))).^2+(ROz(backCCInd)-zCords(laminaInteresects(1))).^2);
                   backRetInd2 = find(laminaOuterInds == inds(laminaInteresects(1)));
                end
                [retT, backRetInd] = min(sqrt((ROx(backCCInd)-Lax).^2+(ROy(backCCInd)-Lay).^2+(ROz(backCCInd)-Laz).^2));
                ret_CP = 1; 
                if ~isempty(laminaInteresects)
                    %if normal intersection is not much further than closest point keep it
                    %not another good threshold for this :/
                    if retT2 < retT*limRetinaDist_Thresh
                        retT = retT2;
                        backRetInd = backRetInd2;
                        ret_CP = 0; 
                    end
                end
                fullBeeData(i).RetThickness(j) = retT*fullBeeData(i).VoxelSize;
            else
                %just take closest point from rear surf to inner lens
                [lensT, backLensInd] = min(sqrt((LOx(fullBeeData(i).interpPoints(j))-LIx).^2+(LOy(fullBeeData(i).interpPoints(j))-LIy).^2+(LOz(fullBeeData(i).interpPoints(j))-LIz).^2));
                fullBeeData(i).lensThickness(j) = lensT*fullBeeData(i).VoxelSize;
                
                %and from there to inner cc
                [CCT, backCCInd] = min(sqrt((LIx(backLensInd)-ROx).^2+(LIy(backLensInd)-ROy).^2+(LIz(backLensInd)-ROz).^2));
                fullBeeData(i).CCThickness(j) = CCT*fullBeeData(i).VoxelSize;
                
                %now from inner cc to lamina interface - always done from closest point on rear of cc
                [retT, backRetInd] = min(sqrt((ROx(backCCInd)-Lax).^2+(ROy(backCCInd)-Lay).^2+(ROz(backCCInd)-Laz).^2));
                fullBeeData(i).RetThickness(j) = retT*fullBeeData(i).VoxelSize;    
            end

            if testPlots 
                subplot(1,2,1);
                line([0 fullBeeData(i).normals(j,1)*1000]+fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),1),...
                        [0 fullBeeData(i).normals(j,2)*1000]+fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),2),...
                        [0 fullBeeData(i).normals(j,3)*1000]+fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),3)); 
                    
                subplot(1,2,2);
                line([LOx(fullBeeData(i).interpPoints(j)) LIx(backLensInd)],[LOy(fullBeeData(i).interpPoints(j)) LIy(backLensInd)],...
                        [LOz(fullBeeData(i).interpPoints(j)) LIz(backLensInd)], 'color', 'r'); 
                line([ROx(backCCInd) LIx(backLensInd)],[ROy(backCCInd) LIy(backLensInd)],...
                        [ROz(backCCInd) LIz(backLensInd)], 'color', 'k');  
                line([ROx(backCCInd) Lax(backRetInd)],[ROy(backCCInd) Lay(backRetInd)],...
                        [ROz(backCCInd) Laz(backRetInd)], 'color', 'b'); 
                if lens_CP
                    plot3(LIx(backLensInd), LIy(backLensInd), LIz(backLensInd), 'mx');
                end
                if CC_CP
                    plot3(ROx(backCCInd), ROy(backCCInd), ROz(backCCInd), 'cx');
                end
                if ret_CP
                    plot3(Lax(backRetInd), Lay(backRetInd), Laz(backRetInd), 'gx');
                end
            end
        end
       
        %show some warnings
        if sum(fullBeeData(i).lensThickness > lensT_Thresh) > length(fullBeeData(i).lensThickness)*0.9
           warning('more than 10% of lens thickness over length limit'); 
        end
        
        if sum(fullBeeData(i).CCThickness > CCT_Thresh) > length(fullBeeData(i).CCThickness)*0.9
           warning('more than 10% of CC thickness over length limit'); 
        end
        %%%lens thickness if frequnely bi-modal, maybe some weird normals
        %should limit length and calc shortest distance if crazy long
        %calcualte normals for border voxels
        %just keep ones that are say, greater than some voxels appart
        [tX,tY,tZ] = ind2sub(fullBeeData(i).eyeVolSize, fullBeeData(i).EyeFrontSurfInds(fullBeeData(i).borderInds));
        for j = 1:length(fullBeeData(i).borderInds)
            [minD, minI] = sort(sqrt((tX-tX(j)).^2+(tY-tY(j)).^2+(tZ-tZ(j)).^2));
            %may need to adjust three depending on eye size. more doesn't hurt
            indsToShift = find(minD > 0 & minD < 3);
            tX(minI(indsToShift)) = -100;
            tY(minI(indsToShift)) = -100;
            tZ(minI(indsToShift)) = -100;
        end

        fullBeeData(i).borderInds(find(tX < 0)) = [];
        fullBeeData(i).borderNormals = zeros(length(fullBeeData(i).borderInds),3);

         for j = 1:length(fullBeeData(i).borderInds)
            %%%note that calculating curvature with this method is not implemented correctly
            fullBeeData(i).borderNormals(j,:) = normalToSurface_withCurvature( fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).borderInds(j),:),...
                    fullBeeData(i).EyeFrontSurfSubsTrans, lensSurfCenter, 1, surfaceFitRng, [xFlatSets, yFlatSets, zFlatSets],fullBeeData(i).borderInds, htStep);

            if testPlots 
                subplot(1,2,1);
                line([0 fullBeeData(i).borderNormals(j,1)*1000]+fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).borderInds(j),1),...
                        [0 fullBeeData(i).borderNormals(j,2)*1000]+fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).borderInds(j),2),...
                        [0 fullBeeData(i).borderNormals(j,3)*1000]+fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).borderInds(j),3), 'color','r'); 
            end
         end
         
        if testPlots 
            if ishandle(fullBeeData(i).histFig)
                close(fullBeeData(i).histFig)
            end
            fullBeeData(i).histFig = figure; 
            
            subplot(1,3,1); hist(fullBeeData(i).lensThickness,20); 
            title(sprintf('%s Lens thickness dist - check no clustered at limt - %i', fullBeeData(i).BeeID, lensT_Thresh));
            subplot(1,3,2); hist(fullBeeData(i).CCThickness,20)
            title(sprintf('CC thickness dist - check no clustered at limt - %i', CCT_Thresh));
            subplot(1,3,3); hist(fullBeeData(i).RetThickness,20)
            title('ret thickness dist');
        end
      
        %% build interpolant for facet size points using IDW
        % then build IDW and do curvature calculations
        %want an array of interp points x num facets
        interpolationDistances = zeros(length(fullBeeData(i).interpPoints),length(fullBeeData(i).facetSizes));
        facetSeperationDistances = zeros(length(fullBeeData(i).facetSizes),length(fullBeeData(i).facetSizes))*NaN; %only interpolate up on thise

        %calc dmap from each facet given geodesic distance path along surface from given facet measurment point
        %then save distance at each interp point

        surfaceMaskVol = zeros(fullBeeData(i).eyeVolSize,'logical');
        surfaceMaskVol(fullBeeData(i).EyeFrontSurfInds) = 1;

        for j = 1:length(fullBeeData(i).facetInds)
            distMap = bwdistgeodesic(surfaceMaskVol,fullBeeData(i).facetInds(j));

            %get distance to sample points
            if sum(isinf(distMap(fullBeeData(i).EyeFrontSurfInds(fullBeeData(i).interpPoints))))
               warning(sprintf('%i infs in surface map for IDW, interpolating to disconected point: %s : bee %i', sum(isinf(distMap(fullBeeData(i).EyeFrontSurfInds(fullBeeData(i).interpPoints)))),...
                   fullBeeData(i).BeeID, i)); 
               indsAsInfs = find(isinf(distMap(fullBeeData(i).EyeFrontSurfInds(fullBeeData(i).interpPoints))));
               indsAreGood = find(~isinf(distMap(fullBeeData(i).EyeFrontSurfInds)));

               %quick fix, get closest with normal dist
               for k = 1:length(indsAsInfs)
                    [~, minI] = min(sqrt((fullBeeData(i).EyeFrontSurfSubsTrans(indsAreGood,1)-fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(indsAsInfs(k)),1)).^2+...
                        (fullBeeData(i).EyeFrontSurfSubsTrans(indsAreGood,2)-fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(indsAsInfs(k)),2)).^2+...
                        (fullBeeData(i).EyeFrontSurfSubsTrans(indsAreGood,3)-fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(indsAsInfs(k)),3)).^2));
                    distMap(fullBeeData(i).EyeFrontSurfInds(fullBeeData(i).interpPoints(indsAsInfs(k)))) = distMap(fullBeeData(i).EyeFrontSurfInds(indsAreGood(minI)));
               end
            end      
            interpolationDistances(:,j) = distMap(fullBeeData(i).EyeFrontSurfInds(fullBeeData(i).interpPoints));

            %get distance between hexagon measrument points
            if sum(isinf(distMap(fullBeeData(i).facetInds)))
               warning(sprintf('%i infs in surface map for IDW, interpolating to disconected point: %s : bee %i', sum(isinf(distMap(fullBeeData(i).EyeFrontSurfInds(fullBeeData(i).facetInds)))),...
                   fullBeeData(i).BeeID, i));
               indsAsInfs = find(isinf(distMap(fullBeeData(i).facetInds)));
               indsAreGood = find(~isinf(distMap(fullBeeData(i).EyeFrontSurfInds)));

               %quick fix, get closest with normal dist
               for k = 1:length(indsAsInfs)
                    [~, minI] = min(sqrt((fullBeeData(i).EyeFrontSurfSubsTrans(indsAreGood,1)-fullBeeData(i).facetLocsTrans(indsAsInfs(k),1)).^2+...
                        (fullBeeData(i).EyeFrontSurfSubsTrans(indsAreGood,2)-fullBeeData(i).facetLocsTrans(indsAsInfs(k),2)).^2+...
                        (fullBeeData(i).EyeFrontSurfSubsTrans(indsAreGood,3)-fullBeeData(i).facetLocsTrans(indsAsInfs(k),3)).^2));
                    distMap(fullBeeData(i).facetInds(indsAsInfs(k))) = distMap(fullBeeData(i).EyeFrontSurfInds(indsAreGood(minI)));
               end
            end      
            facetSeperationDistances(:,j) = distMap(fullBeeData(i).facetInds);
        end

        %now interp facet measurments to interp points using IDW from geodesic distance across lens surface
        if ~useDirectHexAxes
            %original
            interpData = zeros(length(fullBeeData(i).interpPoints),2,2);
        else
            %need to get derivative for each point on initial to target
            
            %slightly inefficient as half will be nan padded
            hexagonPointDerivative = zeros(length(fullBeeData(i).facetInds), length(fullBeeData(i).facetInds), 6, 2);
            for j = 1:length(fullBeeData(i).facetInds)-1
                for k = j+1:length(fullBeeData(i).facetInds)
                    basedCoords = permute(fullBeeData(i).HexagonCoords(j,:,:), [2 3 1]);
                    targetCoords = permute(fullBeeData(i).HexagonCoords(k,:,:), [2 3 1]);

                    %tried aligning shapes with procrustes before hand, but tended to produce weird results as some were rotated by 180

                    %calc diff between each set of coords
                    diffMatrix = zeros(6,6);
                    for l = 1:6
                        for m = 1:6
                            diffMatrix(l,m) = sqrt((basedCoords(fullBeeData(i).HexLinkOrder(j,l),1)-targetCoords(fullBeeData(i).HexLinkOrder(k,m),1))^2+...
                                (basedCoords(fullBeeData(i).HexLinkOrder(j,l),2)-targetCoords(fullBeeData(i).HexLinkOrder(k,m),2))^2);
                        end
                    end

                    %get closest link between both sets
                    [~, linkage] = min(diffMatrix(:));
                    [linkBase,linkTarget] = ind2sub([6, 6], linkage);

                    if plotHexagonInterpInd
                        figure; hold on %not this will not be deleted on each iteration
                        subplot(2,1,1); axis equal; 
                        title(sprintf('%i',facetSeperationDistances(j,k)));
                        for l = 1:6
                            setHere = [l l+1];
                            setHere(setHere > 6) = setHere(setHere > 6) - 6;

                            line([basedCoords(fullBeeData(i).HexLinkOrder(j,setHere(1)),1) basedCoords(fullBeeData(i).HexLinkOrder(j,setHere(2)),1)], ...
                                    [basedCoords(fullBeeData(i).HexLinkOrder(j,setHere(1)),2) basedCoords(fullBeeData(i).HexLinkOrder(j,setHere(2)),2)]);
                            line([basedCoords(fullBeeData(i).HexLinkOrder(j,setHere(1)),1) basedCoords(fullBeeData(i).HexLinkOrder(j,setHere(2)),1)], ...
                                    [basedCoords(fullBeeData(i).HexLinkOrder(j,setHere(1)),2) basedCoords(fullBeeData(i).HexLinkOrder(j,setHere(2)),2)],'color','g');
                            line([targetCoords(fullBeeData(i).HexLinkOrder(k,setHere(1)),1) targetCoords(fullBeeData(i).HexLinkOrder(k,setHere(2)),1)], ...
                                    [targetCoords(fullBeeData(i).HexLinkOrder(k,setHere(1)),2) targetCoords(fullBeeData(i).HexLinkOrder(k,setHere(2)),2)],'color','r');
                        end
                    end

                    %step around and get derivatives
                    for l = 1:6
                        indBase = linkBase+l-1; if indBase > 6; indBase = indBase-6; end
                        indTarget = linkTarget+l-1; if indTarget > 6; indTarget = indTarget-6; end

                        %correct in case link + l greater than 
                        hexagonPointDerivative(j,k,fullBeeData(i).HexLinkOrder(j,indBase),1) = (targetCoords(fullBeeData(i).HexLinkOrder(k,indTarget),1)-basedCoords(fullBeeData(i).HexLinkOrder(j,indBase),1))/...
                                facetSeperationDistances(j,k);
                        hexagonPointDerivative(j,k,fullBeeData(i).HexLinkOrder(j,indBase),2) = (targetCoords(fullBeeData(i).HexLinkOrder(k,indTarget),2)-basedCoords(fullBeeData(i).HexLinkOrder(j,indBase),2))/...
                                facetSeperationDistances(j,k);
                    end

                    if plotHexagonInterpInd
                        subplot(2,1,2); axis equal; hold on
                        for l = 1:6
                            setHere = [l l+1];
                            setHere(setHere > 6) = setHere(setHere > 6) - 6;

                            line([basedCoords(fullBeeData(i).HexLinkOrder(j,setHere(1)),1) basedCoords(fullBeeData(i).HexLinkOrder(j,setHere(2)),1)], ...
                                    [basedCoords(fullBeeData(i).HexLinkOrder(j,setHere(1)),2) basedCoords(fullBeeData(i).HexLinkOrder(j,setHere(2)),2)]);
                            line([targetCoords(fullBeeData(i).HexLinkOrder(k,setHere(1)),1) targetCoords(fullBeeData(i).HexLinkOrder(k,setHere(2)),1)]+facetSeperationDistances(j,k), ...
                                    [targetCoords(fullBeeData(i).HexLinkOrder(k,setHere(1)),2) targetCoords(fullBeeData(i).HexLinkOrder(k,setHere(2)),2)],'color','r','linewidth',3);
                        end
                        for m = (1:10)/10
                            basedCoordsNew = basedCoords;
                            for l = 1:6
                                basedCoordsNew(l,1) =  basedCoordsNew(l,1) + hexagonPointDerivative(j,k,l,1)*facetSeperationDistances(j,k)*m;
                                basedCoordsNew(l,2) =  basedCoordsNew(l,2) + hexagonPointDerivative(j,k,l,2)*facetSeperationDistances(j,k)*m;
                            end

                            for l = 1:6
                                setHere = [l l+1];
                                setHere(setHere > 6) = setHere(setHere > 6) - 6;

                                line([basedCoordsNew(fullBeeData(i).HexLinkOrder(j,setHere(1)),1) basedCoordsNew(fullBeeData(i).HexLinkOrder(j,setHere(2)),1)]+facetSeperationDistances(j,k)*m, ...
                                        [basedCoordsNew(fullBeeData(i).HexLinkOrder(j,setHere(1)),2) basedCoordsNew(fullBeeData(i).HexLinkOrder(j,setHere(2)),2)],'color','g');
                            end
                        end
                    end
                end
            end

            interpData = zeros(length(fullBeeData(i).interpPoints),6,2);
            estimatedArea = zeros(length(fullBeeData(i).interpPoints),1);
        end

        fullBeeData(i).interpFacetArea = zeros(length(fullBeeData(i).interpPoints),1);

        if plotHexagonInterp
            if ishandle(fullBeeData(i).hexInterpF)
                close(fullBeeData(i).hexInterpF)
            end
            fullBeeData(i).hexInterpF = figure; 
            
            hold on; axis equal; ax = gca; ax.Clipping = 'off';
            title(sprintf('%s',fullBeeData(i).BeeID));

            for j = 1:length(fullBeeData(i).facetInds)
                basedCoords = permute(fullBeeData(i).HexagonCoords(j,:,:), [2 3 1]);
                for k = 1:6
                    setHere = [k k+1]; setHere(setHere > 6) = setHere(setHere > 6) - 6;

                    line([0 0]+fullBeeData(i).facetLocsTrans(j,1),...
                            [basedCoords(fullBeeData(i).HexLinkOrder(j,setHere(1)),1) basedCoords(fullBeeData(i).HexLinkOrder(j,setHere(2)),1)]+fullBeeData(i).facetLocsTrans(j,2), ...
                            [basedCoords(fullBeeData(i).HexLinkOrder(j,setHere(1)),2) basedCoords(fullBeeData(i).HexLinkOrder(j,setHere(2)),2)]+fullBeeData(i).facetLocsTrans(j,3));
                end
            end
        end

        for j = 1:length(fullBeeData(i).interpPoints)
            %test if point at loc
            ind = find(interpolationDistances(j,:) == 0);
            if ~isempty(ind)
                %just copy
                if ~useDirectHexAxes
                    interpData(j,1,:) = fullBeeData(i).distAndAng2Point(ind,:);
                    interpData(j,2,:) = fullBeeData(i).distAndAng2Para(ind,:);
                else
                    estimatedArea(j) = fullBeeData(i).hexArea(ind);
                    
                    %calcaulte parameters from hexagon as we haven't done it elsewhere
                    basedCoords = permute(fullBeeData(i).HexagonCoords(ind,:,:), [2 3 1]);
                    
                   %sort as in original section
                   angsVInW = zeros(6,1);
                   for k = 1:6
                       angsVInW(k) = atan2(basedCoords(k,2),basedCoords(k,1));
                   end
                   [~, linkOrderV] = sort(angsVInW);

                   %get line lens as in original section
                   lineSetTemp = zeros(6,2);
                   for k = 1:6
                       sideSet1 = [k, k+1]; %adjacent point e.g. a,b

                       sideSet1(sideSet1 > 6) = sideSet1(sideSet1 > 6) - 6;
                       sideCent1 = [mean(basedCoords(linkOrderV(sideSet1),1)), mean(basedCoords(linkOrderV(sideSet1),2))];

                       lineSetTemp(k,1) = sqrt(sideCent1(1)^2+sideCent1(2)^2)*2;
                       lineSetTemp(k,2) = atan2(sideCent1(2), sideCent1(1));
                   end
                   %sort based on angles
                   [~, sInds] = sort(lineSetTemp(:,2));  
                   lineSetTemp = lineSetTemp(sInds,:); 

                   %copy to interpolator!
                   interpData(j,:,:) = lineSetTemp;
                end

                fullBeeData(i).interpFacetArea(j) = fullBeeData(i).hexArea(ind);
            else
                weights = 1./(interpolationDistances(j,:).^powerP);
                
                if ~useDirectHexAxes
                    interpData(j,1,1) = sum(fullBeeData(i).distAndAng2Point(:,1).*weights(:)')/sum(weights(:));
                    interpData(j,1,2) = sum(fullBeeData(i).distAndAng2Point(:,2).*weights(:)')/sum(weights(:));
                    interpData(j,2,1) = sum(fullBeeData(i).distAndAng2Para(:,1).*weights(:)')/sum(weights(:));
                    interpData(j,2,2) = sum(fullBeeData(i).distAndAng2Para(:,2).*weights(:)')/sum(weights(:));
                else
                    %interpolate hexagon shape - party times
                    predCoords = zeros(length(fullBeeData(i).facetInds)^2,6,2)*NaN;
                    coordWeights = zeros(length(fullBeeData(i).facetInds)^2,6)*NaN;
                    goodIts = 1;

                    %get interpolation estimates from all pair combinations
                    for k = 1:length(fullBeeData(i).facetInds)-1
                        for l = k+1:length(fullBeeData(i).facetInds)
                              %including this test is probably good 
                              %but then need to 
                                if interpolationDistances(j,k) < IDWMinRangeFacet & interpolationDistances(j,l) < IDWMinRangeFacet
                                        stepBetween = facetSeperationDistances(k,l)*interpolationDistances(j,k)/(interpolationDistances(j,k)+interpolationDistances(j,l));

                                        basedCoords = permute(fullBeeData(i).HexagonCoords(k,:,:), [2 3 1]);
                                        for m = 1:6
                                            basedCoords(m,1) =  basedCoords(m,1) + hexagonPointDerivative(k,l,m,1)*stepBetween;
                                            basedCoords(m,2) =  basedCoords(m,2) + hexagonPointDerivative(k,l,m,2)*stepBetween;
                                        end
                                        predCoords(goodIts,:,:) = basedCoords;
                                        coordWeights(goodIts,:) = (weights(k)+weights(l))/2;
                                        goodIts = goodIts + 1;
                                end
                        end
                    end
                    nInds = find(isnan(coordWeights(:,1)));
                    predCoords(nInds,:,:) = [];
                    coordWeights(nInds,:) = [];

                    %now use k means to break up estimated point could into hexagon verts
                    %start at highest weighted point
                    [~, startInd] = max(coordWeights(:,1));

                    tempCoords = reshape(predCoords, [size(predCoords,1)*6 2]);
                    tempStartCoords = reshape(predCoords(startInd,:,:), [6 2]);

                    %cosine takes everything down to unit circle, so looses radial structure
                    %'Distance', 'cosine' - can still get group back in weighted means though
                    [label,centroid] = kmeans(tempCoords, 6, 'Start', tempStartCoords);

                    %need to take weighted average of center of each group
                    wMCoords = zeros(6,2);
                    for k = 1:6
                       inds = find(label == k);
                       wMCoords(k,:) = [wmean(tempCoords(inds,1),coordWeights(inds)) wmean(tempCoords(inds,2),coordWeights(inds))];
                    end

                   %sort as in original section
                   angsVInW = zeros(6,1);
                   for k = 1:6
                       angsVInW(k) = atan2(wMCoords(k,2),wMCoords(k,1));
                   end
                   [~, linkOrderV] = sort(angsVInW);

                   linkOrderV = [linkOrderV' linkOrderV(1)]';
                   for k = 1:6
                        ptSet = [k+1, k];%in CW direction
                        ptSet(ptSet > 6) = ptSet(ptSet > 6) - 6;

                        if plotHexagonInterp 
                            line([0 0]+fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),1), ...
                                [wMCoords(linkOrderV(k),1) wMCoords(linkOrderV(k+1),1)]+fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),2),...
                                [wMCoords(linkOrderV(k),2) wMCoords(linkOrderV(k+1),2)]+fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),3),'linewidth',0.1,'color','k');
                        end
                        %calcualte triangle area
                        lenA = sqrt((0-wMCoords(linkOrderV(ptSet(1)),1))^2+(0-wMCoords(linkOrderV(ptSet(1)),2))^2);
                        lenB = sqrt((0-wMCoords(linkOrderV(ptSet(2)),1))^2+(0-wMCoords(linkOrderV(ptSet(2)),2))^2);
                        lenC = sqrt((wMCoords(linkOrderV(ptSet(1)),1)-wMCoords(linkOrderV(ptSet(2)),1))^2 + (wMCoords(linkOrderV(ptSet(1)),2)-wMCoords(linkOrderV(ptSet(2)),2))^2);
                        p = (lenA+lenB+lenC)/2;
                        triArea = sqrt(p*(p-lenA)*(p-lenB)*(p-lenC));
                        estimatedArea(j) = estimatedArea(j) + triArea;
                   end

                   %get line lens as in original section
                   lineSetTemp = zeros(6,2);
                   for k = 1:6
                       sideSet1 = [k, k+1]; %adjacent point e.g. a,b

                       sideSet1(sideSet1 > 6) = sideSet1(sideSet1 > 6) - 6;
                       sideCent1 = [mean(wMCoords(linkOrderV(sideSet1),1)), mean(wMCoords(linkOrderV(sideSet1),2))];

                       %note explicit multiply by two to extend from side border to next lattice point
                       lineSetTemp(k,1) = sqrt(sideCent1(1)^2+sideCent1(2)^2)*2;
                       lineSetTemp(k,2) = atan2(sideCent1(2), sideCent1(1));
                       
                       if plotHexagonInterp 
                            line([0 0]+fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),1), ...
                                [sideCent1(1) 0]+fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),2),...
                                [sideCent1(2) 0]+fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),3),'linewidth',0.1,'color','r');
                       end
                   end
                   %sort based on angles
                   [~, sInds] = sort(lineSetTemp(:,2));  
                   lineSetTemp = lineSetTemp(sInds,:); 

                   %copy to interpolator!
                   interpData(j,:,:) = lineSetTemp;
                end

                fullBeeData(i).interpFacetArea(j) = sum(fullBeeData(i).hexArea(:).*weights(:))/sum(weights(:));
            end
        end

        %tests for IDW validity - just if range exceeded, note if too small
        if ~useDirectHexAxes
            if sum(interpData(:,1,1) > max(fullBeeData(i).distAndAng2Point(:,1))*(1+IDWThresh)) | sum(interpData(:,1,1) < min(fullBeeData(i).distAndAng2Point(:,1))*(1-IDWThresh))
                [max(interpData(:,1,1)) min(interpData(:,1,1)) max(fullBeeData(i).distAndAng2Point(:,1)) min(fullBeeData(i).distAndAng2Point(:,1))]
                error('interpData1(:,1) exceeds threshold on bee %i', i)
            end
            %correct range test in case one is negative
            if min(fullBeeData(i).distAndAng2Point(:,2)) < 0; negFlip = -1; else; negFlip = 1; end
            if max(fullBeeData(i).distAndAng2Point(:,2)) < 0; posFlip = -1; else; posFlip = 1; end
            if sum(interpData(:,1,2) > max(fullBeeData(i).distAndAng2Point(:,2))*(1+IDWThresh*posFlip)) | sum(interpData(:,1,2) < min(fullBeeData(i).distAndAng2Point(:,2))*(1-IDWThresh*negFlip))
                [max(interpData(:,1,2)) min(interpData(:,1,2)) max(fullBeeData(i).distAndAng2Point(:,2)) min(fullBeeData(i).distAndAng2Point(:,2))]
                error('interpData1(:,2) exceeds threshold on bee %i', i)
            end

            if sum(interpData(:,2,1) > max(fullBeeData(i).distAndAng2Para(:,1))*(1+IDWThresh)) | sum(interpData(:,2,1) < min(fullBeeData(i).distAndAng2Para(:,1))*(1-IDWThresh))
               [max(interpData(:,2,1)) min(interpData(:,2,1)) max(fullBeeData(i).distAndAng2Para(:,1)) min(fullBeeData(i).distAndAng2Para(:,1))]
               error('interpData2(:,1) exceeds threshold on bee %i', i)
            end
            if min(fullBeeData(i).distAndAng2Para(:,2)) < 0; negFlip = -1; else; negFlip = 1; end
            if max(fullBeeData(i).distAndAng2Para(:,2)) < 0; posFlip = -1; else; posFlip = 1; end
            if sum(interpData(:,2,2) > max(fullBeeData(i).distAndAng2Para(:,2))*(1+IDWThresh*posFlip)) |  sum(interpData(:,2,2) < min(fullBeeData(i).distAndAng2Para(:,2))*(1-IDWThresh*negFlip))
                [max(interpData(:,2,2)) min(interpData(:,2,2)) max(fullBeeData(i).distAndAng2Para(:,2)) min(fullBeeData(i).distAndAng2Para(:,2))]
                error('interpData2(:,2) exceeds threshold on bee %i', i)
            end
        end

        if sum(fullBeeData(i).interpFacetArea > max(fullBeeData(i).hexArea)*(1+IDWThresh)) | sum(fullBeeData(i).interpFacetArea < min(fullBeeData(i).hexArea)*(1-IDWThresh))
           [max(fullBeeData(i).interpFacetArea) min(fullBeeData(i).interpFacetArea) max(fullBeeData(i).hexArea) min(fullBeeData(i).hexArea)]
           error('fullBeeData(i).interpFacetArea exceeds threshold on bee %i', i)
        end

        %for testing facet interpolation
        if testPlots & useDirectHexAxes

            if ishandle(fullBeeData(i).AreaTestF)
                close(fullBeeData(i).AreaTestF)
            end
            fullBeeData(i).AreaTestF = figure; 

            bP = mean(fullBeeData(i).facetLocsTrans);
            areaRng = [min(fullBeeData(i).hexArea), max(fullBeeData(i).hexArea)];
            
            subplot(1,4,1); title(sprintf('area interp %s',fullBeeData(i).BeeID));
            tempVals = (fullBeeData(i).hexArea-areaRng(1))./(areaRng(2)-areaRng(1))*100; 
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            fscatter3([fullBeeData(i).facetLocsTrans(:,1)' bP(1) bP(1)],...
                    [fullBeeData(i).facetLocsTrans(:,2)' bP(2) bP(2)],...
                    [fullBeeData(i).facetLocsTrans(:,3)' bP(3) bP(3)],[tempVals' 0 100],jet(100),30);
            tempVals = (fullBeeData(i).interpFacetArea-areaRng(1))./(areaRng(2)-areaRng(1))*100; 
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            [~,hCB] = fscatter3([fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,1)' bP(1) bP(1)],...
                        [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,2)' bP(2) bP(2)],...
                        [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,3)' bP(3) bP(3)],[tempVals' 0 100],jet(100),20);
            for j = 1:length(fullBeeData(i).hexArea)
                text(fullBeeData(i).facetLocsTrans(j,1), fullBeeData(i).facetLocsTrans(j,2), fullBeeData(i).facetLocsTrans(j,3), sprintf('%i',j));
            end
            view(0, -90);
            
            subplot(1,4,2); title('area from facet interp');
            tempVals = (fullBeeData(i).hexArea-areaRng(1))./(areaRng(2)-areaRng(1))*100; 
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            fscatter3([fullBeeData(i).facetLocsTrans(:,1)' bP(1) bP(1)],...
                    [fullBeeData(i).facetLocsTrans(:,2)' bP(2) bP(2)],...
                    [fullBeeData(i).facetLocsTrans(:,3)' bP(3) bP(3)],[tempVals' 0 100],jet(100),30);
            tempVals = (estimatedArea-areaRng(1))./(areaRng(2)-areaRng(1))*100; 
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            [~,hCB] = fscatter3([fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,1)' bP(1) bP(1)],...
                        [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,2)' bP(2) bP(2)],...
                        [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,3)' bP(3) bP(3)],[tempVals' 0 100],jet(100),20);
            for j = 1:length(fullBeeData(i).hexArea)
                text(fullBeeData(i).facetLocsTrans(j,1), fullBeeData(i).facetLocsTrans(j,2), fullBeeData(i).facetLocsTrans(j,3), sprintf('%i',j));
            end
            view(0, -90);
            
            subplot(1,4,3); title('diff % fairly random cols - higher than mid is pos, less is neg');
            tempVals = (estimatedArea./fullBeeData(i).interpFacetArea*100-100)*5+50; 
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            [~,hCB] = fscatter3([fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,1)' bP(1) bP(1)],...
                        [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,2)' bP(2) bP(2)],...
                        [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,3)' bP(3) bP(3)],[tempVals' 0 100],jet(100),20);
            for j = 1:length(fullBeeData(i).hexArea)
                text(fullBeeData(i).facetLocsTrans(j,1), fullBeeData(i).facetLocsTrans(j,2), fullBeeData(i).facetLocsTrans(j,3), sprintf('%i',j));
            end
            view(0, -90);
            
            subplot(1,4,4); hist(estimatedArea./fullBeeData(i).interpFacetArea,20);
            title('error diff percentage');
        end
        
        if useDirectHexAxes
            fullBeeData(i).areaRatio = estimatedArea./fullBeeData(i).interpFacetArea;

            %tends to be biased towards underestimate
            fprintf('Area error min %.2f, max %.2f, mean %.2f percent\n', min(fullBeeData(i).areaRatio), max(fullBeeData(i).areaRatio), mean(fullBeeData(i).areaRatio));

            %firstly check total hexagon range is similar - slightly more relaxed than other
            if sum(fullBeeData(i).areaRatio > (1+IDWThresh*2.5)) | sum(fullBeeData(i).areaRatio < (1-IDWThresh*2.5))
                [max(estimatedArea) min(estimatedArea) max(fullBeeData(i).interpFacetArea) min(fullBeeData(i).interpFacetArea)]
                error('estimatedArea exceeds fullBeeData(i).interpFacetArea range on bee %i', i)
            end
            
            for j = 1:length(fullBeeData(i).interpPoints)
                interpData(j,:,:) = interpData(j,:,:)./sqrt(fullBeeData(i).areaRatio(j));
            end
        end

        %now calculate normals on front of lattice       
        lensSurfCenter = mean(fullBeeData(i).EyeFrontSurfSubsTrans);

        if ~useDirectHexAxes
            tilVal = 2;
            cols = lines(2);
        else
            tilVal = 6;
            cols = lines(6);
        end
        normalFromSpacing = zeros(length(fullBeeData(i).interpPoints),tilVal)*NaN;
        fullBeeData(i).projectedAngle = zeros(length(fullBeeData(i).interpPoints),1);
        
        rUsed = zeros(length(fullBeeData(i).interpPoints),tilVal)*NaN;
        
        if plotTestAngleCalc
            if ishandle(fullBeeData(i).testAngF)
                close(fullBeeData(i).testAngF)
            end
            fullBeeData(i).testAngF = figure; 
            
            hold on; axis equal
            plot3(fullBeeData(i).EyeFrontSurfSubsTrans(:,1),fullBeeData(i).EyeFrontSurfSubsTrans(:,2),fullBeeData(i).EyeFrontSurfSubsTrans(:,3),'k.');
            axis off; ax = gca; ax.Clipping = 'off';
        end

        for j = 1:length(fullBeeData(i).interpPoints)
            surfFit2Use = surfaceFitRng;
            
            %new test, can enlarge if IO too low on 1st test
            while (nanmean(normalFromSpacing(j,:)) < minReqIO/180*pi | surfFit2Use == surfaceFitRng) & ~isempty(minReqIO)
                 if plotTestAngleCalc %& ...
                        line([0 fullBeeData(i).normals(j,1)]*100+fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),1), [0 fullBeeData(i).normals(j,2)]*100+fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),2),...
                         [0 fullBeeData(i).normals(j,3)]*100+fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),3),'color','m','linewidth',3)

                    text(fullBeeData(i).normals(j,1)*100+fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),1), fullBeeData(i).normals(j,2)*100+fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),2),...
                         fullBeeData(i).normals(j,3)*100+fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),3), sprintf('%i', j));
                 end

                %store some temp variables for later analysis
                normalResults = zeros(6,3); %the calculated normal vector
                normalPoints = zeros(6,3); %the start of the normal vector
                areaPoints = zeros(6,3); %the coords of the start points in the rotated coord frame (so essentially flat in x plane)

                %adjust rotation of central point such that its normal is aligned to the equator
                elOfCentralPoint = acos(fullBeeData(i).normals(j,2));
                azOfCentralPoint = atan2(fullBeeData(i).normals(j,3),fullBeeData(i).normals(j,1));
                centralNormalCorrected = fullBeeData(i).normals(j,:)*vrrotvec2mat([0 1 0 -azOfCentralPoint])*vrrotvec2mat([0 0 1 -elOfCentralPoint+pi/2]);
                %central point remains the same as rotation about it
                %now rotate entire surface into correct cord frame about current test point
                allPointsCorrected = (fullBeeData(i).EyeFrontSurfSubsTrans-fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),:))*...
                            vrrotvec2mat([0 1 0 -azOfCentralPoint])*vrrotvec2mat([0 0 1 -elOfCentralPoint+pi/2]);

                for k = 1:tilVal
                    %get points in each direction along lattice axis
                    pVal = [interpData(j,k,1)*cos(interpData(j,k,2)), interpData(j,k,1)*sin(interpData(j,k,2))];
                    [~,intIndPos] = min(sqrt((allPointsCorrected(:,1)-0).^2+(allPointsCorrected(:,2)-pVal(1)).^2+(allPointsCorrected(:,3)-pVal(2)).^2));
                     %calc distance based on 2D projections, otherwise slopping areas will get dropped
                    intDistPos = sqrt((allPointsCorrected(intIndPos,2)-pVal(1)).^2+(allPointsCorrected(intIndPos,3)-pVal(2)).^2);

                    %get angles between vectors as direct calculation of normal differences
                    %check each is close to desired point
                    if fullBeeData(i).VoxelSize*2 > intDistPos
        %              http://geomalgorithms.com/a07-_distance.html
                         %need to take normal at test points, 

                        %%% currently this is called 7 times per measurment loc, 
                        %%% would be more efficient to save result for facet center and query derivative at hex points from same surface. 
                        %%%should match for nearby points and save comp time
                        [normalIntPos,~,~,~,rUsed(j,k)] = normalToSurface_withCurvature( fullBeeData(i).EyeFrontSurfSubsTrans(intIndPos,:),fullBeeData(i).EyeFrontSurfSubsTrans, lensSurfCenter, 1,...
                                surfFit2Use, [xFlatSets, yFlatSets, zFlatSets],fullBeeData(i).borderInds, htStep);

                        %store temp variables
                        normalResults(k,:) = normalIntPos;
                        normalPoints(k,:) = fullBeeData(i).EyeFrontSurfSubsTrans(intIndPos,:);
                        areaPoints(k,:) = allPointsCorrected(intIndPos,:);

                        %get angle between divided by measurment distance then mult by actual distance of facet
                        distPos = sqrt(allPointsCorrected(intIndPos,2).^2+allPointsCorrected(intIndPos,3).^2);
                        angeBetween_IntPos = 2*asin(sqrt((normalIntPos(1)-fullBeeData(i).normals(j,1)).^2+(normalIntPos(2)-fullBeeData(i).normals(j,2)).^2+...
                              (normalIntPos(3)-fullBeeData(i).normals(j,3)).^2)/2)*interpData(j,k,1)/distPos;
                    else
                        %%% could also be useful to check difference is not zero, which would indicate an interpolation failure...
                        normalResults(k,:) = NaN;
                        normalPoints(k,:) = NaN;
                        areaPoints(k,:) = NaN;

                        normalIntPos = [NaN NaN NaN];
                        angeBetween_IntPos = NaN;
                    end
                    if plotTestAngleCalc  %& ...
                         line([0 normalIntPos(1)]*100+fullBeeData(i).EyeFrontSurfSubsTrans(intIndPos,1), [0 normalIntPos(2)]*100+fullBeeData(i).EyeFrontSurfSubsTrans(intIndPos,2),...
                             [0 normalIntPos(3)]*100+fullBeeData(i).EyeFrontSurfSubsTrans(intIndPos,3),'color',cols(k,:),'linewidth',3)
                    end

                    %do opposite, as not explicitly encoded
                    if ~useDirectHexAxes
                        pVal = -[interpData(j,k,1)*cos(interpData(j,k,2)), interpData(j,k,1)*sin(interpData(j,k,2))];
                        [~,intIndNeg] = min(sqrt((allPointsCorrected(:,1)-0).^2+(allPointsCorrected(:,2)-pVal(1)).^2+(allPointsCorrected(:,3)-pVal(2)).^2));
                        intDistNeg = sqrt((allPointsCorrected(intIndNeg,2)-pVal(1)).^2+(allPointsCorrected(intIndNeg,3)-pVal(2)).^2);

                        if fullBeeData(i).VoxelSize*2 > intDistNeg
                            [normalIntNeg] = normalToSurface_withCurvature( fullBeeData(i).EyeFrontSurfSubsTrans(intIndNeg,:),fullBeeData(i).EyeFrontSurfSubsTrans, lensSurfCenter, 1, ...
                                    surfFit2Use, [xFlatSets, yFlatSets, zFlatSets], fullBeeData(i).borderInds, htStep);    
                            distPos = sqrt(allPointsCorrected(intDistNeg,2).^2+allPointsCorrected(intDistNeg,3).^2);
                            angeBetween_IntNeg = 2*asin(sqrt((normalIntNeg(1)-fullBeeData(i).normals(j,1)).^2+(normalIntNeg(2)-fullBeeData(i).normals(j,2)).^2+...
                                  (normalIntNeg(3)-fullBeeData(i).normals(j,3)).^2)/2)*interpData(j,k,1)/distPos;
                        else
                            normalIntNeg = [NaN NaN NaN];
                            angeBetween_IntNeg = NaN;
                        end
                        if plotTestAngleCalc
                             line([0 normalIntNeg(1)]*100+fullBeeData(i).EyeFrontSurfSubsTrans(intIndNeg,1), [0 normalIntNeg(2)]*100+fullBeeData(i).EyeFrontSurfSubsTrans(intIndNeg,2),...
                                 [0 normalIntNeg(3)]*100+fullBeeData(i).EyeFrontSurfSubsTrans(intIndNeg,3),'color',cols(k,:),'linewidth',3)
                        end
                        normalFromSpacing(j,k) = nanmean([angeBetween_IntPos angeBetween_IntNeg]);
                    else
                        normalFromSpacing(j,k) = angeBetween_IntPos;
                    end
                end
                surfFit2Use = surfFit2Use+50;
                if surfFit2Use > 500
                    error('probably error with regulation surf size to minReqIO parameter');
                end
            end
            
            if useDirectHexAxes
                %with all normals found, figure out their intersections on spehre
                toIgnore = find(isnan(normalResults(:,1)));
                normalResults(toIgnore,:) = [];
                normalPoints(toIgnore,:) = [];
                normalInts = normalResults*0;

                for k = 1:size(normalResults,1)
                    origin = normalPoints(k,:);
                    normal = normalResults(k,:);
                    p = (2*normal(1)* origin(1)+2*normal(2)*origin(2)+2*normal(3)*origin(3))/(normal(1)^2+normal(2)^2+normal(3)^2);
                    q = (origin(1)^2+origin(2)^2+origin(3)^2-sphereRad^2)/(normal(1)^2+normal(2)^2+normal(3)^2);
                    val = -(p/2)+sqrt((p/2)^2-q);
                    normalInts(k,:) = [normal(1)*val+origin(1), normal(2)*val+origin(2), normal(3)*val+origin(3)]/sphereRad;
                end
                %ok, this is already done later - oh well, its fast
                origin = fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),:);
                normal = fullBeeData(i).normals(j,:);
                p = (2*normal(1)* origin(1)+2*normal(2)*origin(2)+2*normal(3)*origin(3))/(normal(1)^2+normal(2)^2+normal(3)^2);
                q = (origin(1)^2+origin(2)^2+origin(3)^2-sphereRad^2)/(normal(1)^2+normal(2)^2+normal(3)^2);
                val = -(p/2)+sqrt((p/2)^2-q);
                centInt= [normal(1)*val+origin(1), normal(2)*val+origin(2), normal(3)*val+origin(3)]/sphereRad;

                %if one is removed, may break sorting...
                if size(normalResults,1) ~= 6
                   angsVInW = zeros(size(normalResults,1),1);
                   for k = 1:size(normalResults,1)
                       angsVInW(k) = atan2(areaPoints(k,3),areaPoints(k,2));
                   end
                   [~, linkOrderV] = sort(angsVInW);
                   areaPoints = areaPoints(linkOrderV,:);
                   normalInts = normalInts(linkOrderV,:);
                end

                %sum up angular areas subtended by lattice points
                area = 0;
                projectedArea = 0;
                for k = 1:size(normalResults,1)
                    ptSet = [k+1, k];%in CW direction
                    ptSet(ptSet > size(normalResults,1)) = ptSet(ptSet > size(normalResults,1)) - size(normalResults,1);

                    %calculate 2D area of triangle
                    lenA = sqrt((0-areaPoints(ptSet(1),2))^2+(0-areaPoints(ptSet(1),3))^2);
                    lenB = sqrt((0-areaPoints(ptSet(2),2))^2+(0-areaPoints(ptSet(2),3))^2);
                    lenC = sqrt((areaPoints(ptSet(1),2)-areaPoints(ptSet(2),2))^2+(areaPoints(ptSet(1),3)-areaPoints(ptSet(2),3))^2);
                    p = (lenA+lenB+lenC)/2;
                    triArea = sqrt(p*(p-lenA)*(p-lenB)*(p-lenC));
                    area = area + triArea;

                    %calculate angular area
                    p1 = centInt;
                    p2 = normalInts(ptSet(1),:);
                    p3 = normalInts(ptSet(2),:);
                    thetaP1P2 = acos((p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3))/sqrt((p1(1)^2+p1(2)^2+p1(3)^2)*(p2(1)^2+p2(2)^2+p2(3)^2))); 
                    thetaP1P3 = acos((p1(1)*p3(1)+p1(2)*p3(2)+p1(3)*p3(3))/sqrt((p1(1)^2+p1(2)^2+p1(3)^2)*(p3(1)^2+p3(2)^2+p3(3)^2))); 
                    thetaP2P3 = acos((p2(1)*p3(1)+p2(2)*p3(2)+p2(3)*p3(3))/sqrt((p2(1)^2+p2(2)^2+p2(3)^2)*(p3(1)^2+p3(2)^2+p3(3)^2))); 
                    thetaS = (thetaP1P2+thetaP1P3+thetaP2P3)/2;
                    projectedArea = projectedArea + 4*atan(sqrt(tan(thetaS/2)*tan((thetaS-thetaP1P2)/2)*tan((thetaS-thetaP1P3)/2)*tan((thetaS-thetaP2P3)/2)));
                end
                %find out proj to area relationship at this location
                
                %%%note that this is based on projected 2D area when later integrated to compare total area, it probably
                %%%undershoots because surface has 3D component from curvature 
                %%%- probably best to get either based on surface mesh, of from fitted surface 
                
                fullBeeData(i).projectedAngle(j) = (projectedArea/area);
                %obviously invert this to get area for projected angle
            else
               error('area density calcs not implemented for h/v yet') 
            end
        end
        fprintf('mean radius used %.1f - min radius used %.1f - max radius used %.1f\n', nanmean(rUsed(:)), min(rUsed(:)), max(rUsed(:)))
        fullBeeData(i).RadiusUsed = [nanmean(rUsed(:)), min(rUsed(:)), max(rUsed(:))];
        
        %interp missing values with idw - possibly none with step in form side
        if ~useDirectHexAxes
            for j = 1:tilVal
                goodDInds = find(~isnan(normalFromSpacing(:,j)));
                missingValsToInterp = find(isnan(normalFromSpacing(:,j)));
                for k = 1:length(missingValsToInterp)
                    dists = sqrt((fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(goodDInds),1)-fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(missingValsToInterp(k)),1)).^2+...
                        (fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(goodDInds),2)-fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(missingValsToInterp(k)),2)).^2+...
                        (fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(goodDInds),3)-fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(missingValsToInterp(k)),3)).^2);
                    weights = 1./(dists.^powerP);
                    normalFromSpacing(missingValsToInterp(k),j) = sum(normalFromSpacing(goodDInds,j).*weights)/sum(weights);
                end
            end
        else
            for j = 1:length(fullBeeData(i).interpPoints)
               missingValsToInterp = find(isnan(normalFromSpacing(j,:)));
               if ~isempty(missingValsToInterp)
                   normalFromSpacing(j,missingValsToInterp) = nanmean(normalFromSpacing(j,:));
               end
            end
        end
        
        fullBeeData(i).avgDiameter = zeros(length(fullBeeData(i).interpPoints),1);
        for j = 1:length(fullBeeData(i).interpPoints)
            if ~useDirectHexAxes
                %get value of d from lattice parralelogram     
                fullBeeData(i).avgDiameter = sqrt((interpData(1,1,1).*cos(interpData(1,1,2))+interpData(1,2,1).*cos(interpData(1,2,2))).^2+...
                    (interpData(1,1,1).*sin(interpData(1,1,2))+interpData(1,2,1).*sin(interpData(1,2,2))).^2);
                %check if calc in other parallelogram orientaiton is greater, if so use that
                altVals = sqrt((interpData(1,1,1).*cos(interpData(1,1,2))-interpData(1,2,1).*cos(interpData(1,2,2))).^2+...
                    (interpData(1,1,1).*sin(interpData(1,1,2))-interpData(1,2,1).*sin(interpData(1,2,2))).^2);
                indsToSwap = find(fullBeeData(i).avgDiameter < altVals);
                fullBeeData(i).avgDiameter(indsToSwap) = altVals(indsToSwap);   
            else
                fullBeeData(i).avgDiameter(j) = mean(interpData(j,:,1));
            end
        end

        %calcualte surface area from half isosurface area

        %create vol and set front surf to 1
        tempBWVol = zeros(fullBeeData(i).eyeVolSize);
        tInds = getIndicesOfMultipleIDs(fullBeeData(i).LabelsStack, [fullBeeData(i).StackLabels.LensOuterS]);
        tempBWVol(tInds) = 1;
        tempSurf = isosurface(tempBWVol,0.5);

        %https://se.mathworks.com/matlabcentral/answers/93023-is-there-a-matlab-function-that-can-compute-the-area-of-my-patch
        verts = tempSurf.vertices*fullBeeData(i).VoxelSize;
        faces = tempSurf.faces;
        a = verts(faces(:, 2), :) - verts(faces(:, 1), :);
        b = verts(faces(:, 3), :) - verts(faces(:, 1), :);
        c = cross(a, b, 2);
        fullBeeData(i).lensSurfaceArea = 1/2 * sum(sqrt(sum(c.^2, 2)))/2;

        fullBeeData(i).NumFacets = fullBeeData(i).lensSurfaceArea/mean(fullBeeData(i).interpFacetArea);
        fullBeeData(i).densityFromArea = 1./(fullBeeData(i).projectedAngle.*fullBeeData(i).interpFacetArea);
        fullBeeData(i).fovFilledFromIntegral = mean(fullBeeData(i).projectedAngle)*fullBeeData(i).lensSurfaceArea/(4*pi);

        fullBeeData(i).WorldIO = zeros(length(fullBeeData(i).interpPoints),2);
        
        %calculate axis density
        if ~useDirectHexAxes
            for j = 1:length(fullBeeData(i).interpPoints)
                %get coords with angles
                IOPointC = [normalFromSpacing(:,1).*cos(interpData(j,1,2)) normalFromSpacing(:,1).*sin(interpData(j,1,2))];
                IOParaC = [normalFromSpacing(:,2).*cos(interpData(j,2,2)) normalFromSpacing(:,2).*sin(interpData(j,2,2))];

                axisAng = acos((IOParaC(:,1).*IOPointC(:,1)+IOParaC(:,2).*IOPointC(:,2))./sqrt((IOParaC(:,1).^2+IOParaC(:,2).^2).*(IOPointC(:,1).^2+IOPointC(:,2).^2)));
                %we always want the acute angle
                axisAng(axisAng > pi/2) = pi-axisAng(axisAng > pi/2);

                %get long axis of parallelogram for interommatidial angle
                fullBeeData(i).calInterO = sqrt((IOPointC(:,1)+IOParaC(:,1)).^2+(IOPointC(:,2)+IOParaC(:,2)).^2);
                %check if other long axis is longer, and swap if so.
                altVals = sqrt((IOPointC(:,1)-IOParaC(:,1)).^2+(IOPointC(:,2)-IOParaC(:,2)).^2);
                indsToSwap = find(fullBeeData(i).calInterO < altVals);
                fullBeeData(i).calInterO(indsToSwap) = altVals(indsToSwap);

                %build full parallelogram for each point to get world el/az resolution

                for k = 1:length(fullBeeData(i).interpPoints)
                    parrelogram = [0 0; IOParaC(k,:); IOPointC(k,:); IOPointC(k,:)+IOParaC(k,:)];

                    %1st is el (vert) but horiz in this coord frame, 2nd is az
                    fullBeeData(i).WorldIO(k,:) = [ max(parrelogram(:,1))-min(parrelogram(:,1)), max(parrelogram(:,2))-min(parrelogram(:,2))];
                end

                %this generally predicts the correct number of facets

                %note formula in Beershema 1977 paper appears to omit factor of 2
                %this assumes not overly irregular grid? I think it assumes one axis of symmetry
                fullBeeData(i).AxisDensity = 1./(2*normalFromSpacing(:,1).*normalFromSpacing(:,2).*sin(axisAng));
            end
        else
            %with direct Hex, its simpler (kind of)
            fullBeeData(i).calInterO = mean(normalFromSpacing')';
            fullBeeData(i).AxisDensity = zeros(length(fullBeeData(i).interpPoints),1);
            
            warning('world IO is probably wrong');

            for j = 1:length(fullBeeData(i).interpPoints)
                %calcualte axis of each to get world coordinates
                IOC = zeros(6,2);
                for k = 1:6
                    IOC(k,:) =  [normalFromSpacing(j,k).*cos(interpData(j,k,2)) normalFromSpacing(j,k).*sin(interpData(j,k,2))];
                end

                %this finds equivelenet to h and v values?
                sortedEl = sort([0 (IOC(1,1)-IOC(4,1))/2 (IOC(2,1)-IOC(5,1))/2 (IOC(3,1)-IOC(6,1))/2]);
                sortedAz = sort([0 (IOC(1,2)-IOC(4,2))/2 (IOC(2,2)-IOC(5,2))/2 (IOC(3,2)-IOC(6,2))/2]);
                fullBeeData(i).WorldIO(j,:) = [max(diff(sortedEl)) max(diff(sortedAz))];

                %this appears to be best method of tested to for getting axis density
                %but to be honest, its not clear exactly what physical parameter we need to calculate here
                angsVInW = zeros(6,1);
                for k = 1:6
                    angsVInW(k) = atan2(IOC(k,2),IOC(k,1));
                end
                [~, sAng] = sort(angsVInW);

                %%%this many not be correct way to calculate axis densitites
                tempArea = zeros(6,1);
                for k = 1:6
                    ptSet1 = [k, k+1, k+2];%in CW direction
                    ptSet1(ptSet1 > 6) = ptSet1(ptSet1 > 6) - 6;
                    ptSet2 = [k, k+2, k+5];%in CW direction
                    ptSet2(ptSet2 > 6) = ptSet2(ptSet2 > 6) - 6;

                    %calcualte triangle area
                    lenA = sqrt((IOC(sAng(ptSet1(3)),1)*0.5-IOC(sAng(ptSet1(1)),1))^2+(IOC(sAng(ptSet1(3)),2)*0.5-IOC(sAng(ptSet1(1)),2))^2);
                    lenB = sqrt((IOC(sAng(ptSet1(3)),1)*0.5-IOC(sAng(ptSet1(2)),1))^2+(IOC(sAng(ptSet1(3)),2)*0.5-IOC(sAng(ptSet1(2)),2))^2);
                    lenC = sqrt((IOC(sAng(ptSet1(1)),1)-IOC(sAng(ptSet1(2)),1))^2 + (IOC(sAng(ptSet1(1)),2)-IOC(sAng(ptSet1(2)),2))^2);
                    p = (lenA+lenB+lenC)/2;
                    triArea1 = sqrt(p*(p-lenA)*(p-lenB)*(p-lenC));

                    lenA = sqrt((IOC(sAng(ptSet2(3)),1)*0.5-IOC(sAng(ptSet2(1)),1))^2+(IOC(sAng(ptSet2(3)),2)*0.5-IOC(sAng(ptSet2(1)),2))^2);
                    lenB = sqrt((IOC(sAng(ptSet2(3)),1)*0.5-IOC(sAng(ptSet2(2)),1)*0.5)^2+(IOC(sAng(ptSet2(3)),2)*0.5-IOC(sAng(ptSet2(2)),2)*0.5)^2);
                    lenC = sqrt((IOC(sAng(ptSet2(1)),1)-IOC(sAng(ptSet2(2)),1)*0.5)^2 + (IOC(sAng(ptSet2(1)),2)-IOC(sAng(ptSet2(2)),2)*0.5)^2);
                    p = (lenA+lenB+lenC)/2;
                    triArea2 = sqrt(p*(p-lenA)*(p-lenB)*(p-lenC));

                    tempArea(k) = (triArea1+triArea2);
                end
                sArea = sort(tempArea);

                if limLargeDensities
                    %assumption of symmetrical hexagon always smaller
                    if (1/(mean(sArea(1:2))))/(1/(sqrt(3)*fullBeeData(i).calInterO(j).^2/2)) > 1+densityLim
                       %if it is really larger than threshold, use mean of all, which is smaller
                       fullBeeData(i).AxisDensity(j) = 1/(mean(sArea));
                       
                       %%%note that a better catch may be if very assymetric IO values by axis, which creates weird rectangle here 
                    else
                       fullBeeData(i).AxisDensity(j) = 1/(mean(sArea(1:2))); 
                    end
                else
                    fullBeeData(i).AxisDensity(j) = 1/(mean(sArea(1:2)));
                end
            end
        end
        
        %calc sensitvity, assuming that projected area of facet equals its acceptance angle
        if useDensityFromIntegral
            fullBeeData(i).sensitvityApproximation = fullBeeData(i).interpFacetArea./fullBeeData(i).densityFromArea.*(aCoef*fullBeeData(i).RetThickness./(2.3+aCoef*fullBeeData(i).RetThickness));
        else
            fullBeeData(i).sensitvityApproximation = fullBeeData(i).interpFacetArea./fullBeeData(i).AxisDensity.*(aCoef*fullBeeData(i).RetThickness./(2.3+aCoef*fullBeeData(i).RetThickness));
        end
        
        %%% would be good to have eye wide metric, but it is not as simple as this
        fullBeeData(i).expectedSensInt = mean(fullBeeData(i).sensitvityApproximation)*fullBeeData(i).NumFacets;
        
        fprintf('mean ratio between axis density methods %.2f\n', mean(fullBeeData(i).AxisDensity./fullBeeData(i).densityFromArea));
        
        %calc radius from effective IO, at least this works for sure!
        fullBeeData(i).avgRadius = fullBeeData(i).avgDiameter./fullBeeData(i).calInterO;

        for j = 1:length(fullBeeData(i).interpPoints)
            fullBeeData(i).avg(j) = mean(interpData(j,:,1)./normalFromSpacing(j,:));
        end
        
        %this is azimuth divided by elevation?
        %if less than 1, az greater, if more than 1, el greater
        fullBeeData(i).worldRatio = fullBeeData(i).WorldIO(:,1)./fullBeeData(i).WorldIO(:,2);

       if testPlots
            if ishandle(fullBeeData(i).testOnEyeF)
                close(fullBeeData(i).testOnEyeF)
            end
            fullBeeData(i).testOnEyeF = figure; 
            
            bP = mean(fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,:));
            subplot(2,3,1); hold on; 
            tempVals = (fullBeeData(i).facetSizes-diamRng(1))./(diamRng(2)-diamRng(1))*100; 
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            fscatter3([fullBeeData(i).facetLocsTrans(:,1)' bP(1) bP(1)],...
                    [fullBeeData(i).facetLocsTrans(:,2)' bP(2) bP(2)],...
                    [fullBeeData(i).facetLocsTrans(:,3)' bP(3) bP(3)],[tempVals' 0 100],jet(100),40);
            tempVals = (fullBeeData(i).avgDiameter-diamRng(1))./(diamRng(2)-diamRng(1))*100;
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            [~,hCB] = fscatter3([fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,1)' bP(1) bP(1)],...
                    [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,2)' bP(2) bP(2)],...
                    [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,3)' bP(3) bP(3)],[tempVals' 0 100],jet(100),20);
            plot3(bP(1), bP(2), bP(3), 'w.','markersize',100);
            set(hCB, 'ticks', [0 50 100], 'TickLabels', diamLabs,'YAxisLocation','right','TickDirection','out');
            view([0, -90]); title(sprintf('%s facet D - min %.1f, max %.1f, mean %.1f',fullBeeData(i).BeeID, min(fullBeeData(i).avgDiameter), max(fullBeeData(i).avgDiameter), mean(fullBeeData(i).avgDiameter)));

            subplot(2,3,2); hold on;
            tempVals = (fullBeeData(i).avgRadius-curvRng(1))./(curvRng(2)-curvRng(1))*100;
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            [~,hCB] = fscatter3([fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,1)' bP(1) bP(1)],...
                    [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,2)' bP(2) bP(2)],...
                    [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,3)' bP(3) bP(3)],[tempVals' 0 100],jet(100),20);
            plot3(bP(1), bP(2), bP(3), 'w.','markersize',50);
            set(hCB, 'Ticks', [0 50 100], 'TickLabels', curvLabs,'YAxisLocation','right','TickDirection','out');
            view([0, -90])
            plot3(fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(missingValsToInterp),1), fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(missingValsToInterp),2), fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(missingValsToInterp),3), 'mx');
            title('curvature');
            
            subplot(2,3,3); hold on; 
            if useDensityFromIntegral
                tempVals = (fullBeeData(i).densityFromArea-axDRng(1))./(axDRng(2)-axDRng(1))*100;
            else
                tempVals = (fullBeeData(i).AxisDensity-axDRng(1))./(axDRng(2)-axDRng(1))*100; 
            end
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            [~,hCB] = fscatter3([fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,1)' bP(1) bP(1)],...
                    [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,2)' bP(2) bP(2)],...
                    [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,3)' bP(3) bP(3)],[tempVals' 0 100],jet(100),20);
            plot3(bP(1), bP(2), bP(3), 'w.','markersize',50); view([0, -90]); 
            set(hCB, 'Ticks', [0 50 100], 'TickLabels', axDLabs,'YAxisLocation','right','TickDirection','out');
            if useDensityFromIntegral
                title(sprintf('AxisDensity, min %.1f, max %.1f, mean %.1f', min(fullBeeData(i).densityFromArea), max(fullBeeData(i).densityFromArea), mean(fullBeeData(i).densityFromArea)));
            else
                title(sprintf('AxisDensity, min %.1f, max %.1f, mean %.1f', min(fullBeeData(i).AxisDensity), max(fullBeeData(i).AxisDensity), mean(fullBeeData(i).AxisDensity)));
            end
            plot3(fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(missingValsToInterp),1), fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(missingValsToInterp),2), fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(missingValsToInterp),3), 'mx');

            subplot(2,3,4); hold on;
            tempVals = (fullBeeData(i).WorldIO(:,1)/pi*180-intORng(1))./(intORng(2)-intORng(1))*100;
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            [~,hCB] = fscatter3([fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,1)' bP(1) bP(1)],...
                    [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,2)' bP(2) bP(2)],...
                    [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,3)' bP(3) bP(3)],[tempVals' 0 100],jet(100),20);
            plot3(bP(1), bP(2), bP(3), 'w.','markersize',50);
            set(hCB, 'Ticks', [0 50 100], 'TickLabels', intOLabs,'YAxisLocation','right','TickDirection','out');
            view([0, -90]); title(sprintf('Partial IO, world El- min %.1f, max %.1f', min(fullBeeData(i).WorldIO(:,1)/pi*180), max(fullBeeData(i).WorldIO(:,1)/pi*180)));
            plot3(fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(missingValsToInterp),1), fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(missingValsToInterp),2), fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(missingValsToInterp),3), 'mx');

            subplot(2,3,5); hold on; 
            tempVals = (fullBeeData(i).WorldIO(:,2)/pi*180-intORng(1))./(intORng(2)-intORng(1))*100;
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            [~,hCB] = fscatter3([fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,1)' bP(1) bP(1)],...
                    [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,2)' bP(2) bP(2)],...
                    [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,3)' bP(3) bP(3)],[tempVals' 0 100],jet(100),20);
            plot3(bP(1), bP(2), bP(3), 'w.','markersize',50);
            set(hCB, 'Ticks', [0 50 100], 'TickLabels', intOLabs,'YAxisLocation','right','TickDirection','out');
            view([0, -90]); title(sprintf('Partial IO, world Az- min %.1f, max %.1f', min(fullBeeData(i).WorldIO(:,2)/pi*180), max(fullBeeData(i).WorldIO(:,2)/pi*180)));
            plot3(fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(missingValsToInterp),1), fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(missingValsToInterp),2), fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(missingValsToInterp),3), 'mx');

            subplot(2,3,6); hold on; 
            tempVals = (fullBeeData(i).calInterO/pi*180-intORng(1))./(intORng(2)-intORng(1))*100;
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            [~,hCB] = fscatter3([fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,1)' bP(1) bP(1)],...
                    [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,2)' bP(2) bP(2)],...
                    [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,3)' bP(3) bP(3)],[tempVals' 0 100],jet(100),20);
            plot3(bP(1), bP(2), bP(3), 'w.','markersize',50);
            set(hCB, 'Ticks', [0 50 100], 'TickLabels', intOLabs,'YAxisLocation','right','TickDirection','out');
                    plot3(fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(missingValsToInterp),1), fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(missingValsToInterp),2), fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(missingValsToInterp),3), 'mx');
            [~,minI] = min(fullBeeData(i).calInterO);
             plot3(fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(minI),1), fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(minI),2), fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(minI),3), 'mx');
            view([0, -90]); title(sprintf('Resultant IO- min %.1f, max %.1f, mean %.1f - minInd %i', min(fullBeeData(i).calInterO/pi*180), max(fullBeeData(i).calInterO/pi*180), mean(fullBeeData(i).calInterO/pi*180), minI));

            if ishandle(fullBeeData(i).testOnEyeFExtra)
                close(fullBeeData(i).testOnEyeFExtra)
            end
            fullBeeData(i).testOnEyeFExtra = figure; 
            
            subplot(2,3,1); hold on;
            tempVals = (fullBeeData(i).lensThickness-lensTRng(1))./(lensTRng(2)-lensTRng(1))*100;
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            [~,hCB] = fscatter3([fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,1)' bP(1) bP(1)],...
                    [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,2)' bP(2) bP(2)],...
                    [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,3)' bP(3) bP(3)],[tempVals' 0 100],jet(100),20);
            plot3(bP(1), bP(2), bP(3), 'w.','markersize',50);
            set(hCB, 'Ticks', [0 50 100], 'TickLabels', lensTLabs,'YAxisLocation','right','TickDirection','out');
            view([0, -90]); title(sprintf('%s Lens thickness',fullBeeData(i).BeeID));
            
            subplot(2,3,2); hold on;
            tempVals = (fullBeeData(i).CCThickness-CCTRng(1))./(CCTRng(2)-CCTRng(1))*100;
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            [~,hCB] = fscatter3([fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,1)' bP(1) bP(1)],...
                    [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,2)' bP(2) bP(2)],...
                    [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,3)' bP(3) bP(3)],[tempVals' 0 100],jet(100),20);
            plot3(bP(1), bP(2), bP(3), 'w.','markersize',50);
            set(hCB, 'Ticks', [0 50 100], 'TickLabels', CCTLabs,'YAxisLocation','right','TickDirection','out');
            view([0, -90]); title(sprintf('CC thickness'));
            
            subplot(2,3,3); hold on;
            tempVals = (fullBeeData(i).RetThickness-RetTRng(1))./(RetTRng(2)-RetTRng(1))*100;
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            [~,hCB] = fscatter3([fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,1)' bP(1) bP(1)],...
                    [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,2)' bP(2) bP(2)],...
                    [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,3)' bP(3) bP(3)],[tempVals' 0 100],jet(100),20);
            plot3(bP(1), bP(2), bP(3), 'w.','markersize',50);
            set(hCB, 'Ticks', [0 50 100], 'TickLabels', RetTLabs,'YAxisLocation','right','TickDirection','out');
            view([0, -90]); title(sprintf('Retina thickness'));
            
            subplot(2,3,4); hold on;
            tempVals = (mean(rUsed')-50)./(300-50)*100;
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            [~,hCB] = fscatter3([fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,1)' bP(1) bP(1)],...
                    [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,2)' bP(2) bP(2)],...
                    [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,3)' bP(3) bP(3)],[tempVals 0 100],jet(100),20);
            plot3(bP(1), bP(2), bP(3), 'w.','markersize',50);
            set(hCB, 'Ticks', [0 50 100],'YAxisLocation','right','TickDirection','out');
            view([0, -90]); title(sprintf('Mean radius for normals'));
            
            subplot(2,3,5);
            plot(fullBeeData(i).calInterO, fullBeeData(i).AxisDensity./(1./(sqrt(3)*fullBeeData(i).calInterO.^2/2)),'b.')
            title(sprintf('diff parralelogram density chosen and symmetrical hex density - std of dif %.2f', std(fullBeeData(i).AxisDensity./(1./(sqrt(3)*fullBeeData(i).calInterO.^2/2)))));
 
            subplot(2,3,6); hold on;
            tempVals = (fullBeeData(i).sensitvityApproximation-sensRng(1))./(sensRng(2)-sensRng(1))*100;
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            [~,hCB] = fscatter3([fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,1)' bP(1) bP(1)],...
                    [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,2)' bP(2) bP(2)],...
                    [fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,3)' bP(3) bP(3)],[tempVals' 0 100],jet(100),20);
            plot3(bP(1), bP(2), bP(3), 'w.','markersize',50);
            set(hCB, 'Ticks', [0 50 100], 'TickLabels', sensLabs,'YAxisLocation','right','TickDirection','out');
            view([0, -90]); title(sprintf('Sensitvity approximation'));
       end

        %% plot onto world to look at sampling density
        fullBeeData(i).sphereIntersect = zeros(length(fullBeeData(i).interpPoints),3);
        for j = 1:length(fullBeeData(i).interpPoints) 
            origin = fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints(j),:);
            normal = fullBeeData(i).normals(j,:);

            p = (2*normal(1)* origin(1)+2*normal(2)*origin(2)+2*normal(3)*origin(3))/(normal(1)^2+normal(2)^2+normal(3)^2);
            q = (origin(1)^2+origin(2)^2+origin(3)^2-sphereRad^2)/(normal(1)^2+normal(2)^2+normal(3)^2);
            val = -(p/2)+sqrt((p/2)^2-q);% +-?
            fullBeeData(i).sphereIntersect(j,:) = [normal(1)*val+origin(1), normal(2)*val+origin(2), normal(3)*val+origin(3)];
        end
        fullBeeData(i).sphereIntersect = fullBeeData(i).sphereIntersect/sphereRad;

        fullBeeData(i).sphereIntersectBorders = zeros(length(fullBeeData(i).borderInds),3);
        for j = 1:length(fullBeeData(i).borderInds) 
            origin = fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).borderInds(j),:);
            normal = fullBeeData(i).borderNormals(j,:);

            p = (2*normal(1)* origin(1)+2*normal(2)*origin(2)+2*normal(3)*origin(3))/(normal(1)^2+normal(2)^2+normal(3)^2);
            q = (origin(1)^2+origin(2)^2+origin(3)^2-sphereRad^2)/(normal(1)^2+normal(2)^2+normal(3)^2);
            val = -(p/2)+sqrt((p/2)^2-q);% +-?
            fullBeeData(i).sphereIntersectBorders(j,:) = [normal(1)*val+origin(1), normal(2)*val+origin(2), normal(3)*val+origin(3)];
        end
        fullBeeData(i).sphereIntersectBorders = fullBeeData(i).sphereIntersectBorders/sphereRad;

        %apply auto align if present (note done to subs used to calc normals)
        fullBeeData(i).sphereIntersect = [fullBeeData(i).sphereIntersect zeros(size(fullBeeData(i).sphereIntersect,1),1)]*extraRotAd';
        fullBeeData(i).sphereIntersectBorders = [fullBeeData(i).sphereIntersectBorders  zeros(size(fullBeeData(i).sphereIntersectBorders,1),1)]*extraRotAd';
        fullBeeData(i).sphereIntersect = fullBeeData(i).sphereIntersect(:,1:3);
        fullBeeData(i).sphereIntersectBorders = fullBeeData(i).sphereIntersectBorders(:,1:3);
            
        %need to get expected IO from nearest normal point for the borders
        expectedIO = zeros(length(fullBeeData(i).borderInds),1);
        for j = 1:length(fullBeeData(i).borderInds) 
            [~,minI] = min(sqrt((fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,1)-fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).borderInds(j),1)).^2 +...
                    (fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,2)-fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).borderInds(j),2)).^2 +...
                    (fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).interpPoints,3)-fullBeeData(i).EyeFrontSurfSubsTrans(fullBeeData(i).borderInds(j),3)).^2));
            expectedIO(j) = fullBeeData(i).calInterO(minI);
        end

        %determine points on ref sphere field of view for this eye
        fullBeeData(i).inFOV = zeros(size(refSphere.X,1),1);
        distToTargets = zeros(size(refSphere.X,1),1);
        for j = 1:size(refSphere.X,1)
            %get central angle between this and all sphere intersects
            deltaAng = 2*asin(sqrt((fullBeeData(i).sphereIntersect(:,1)-refSphere.X(j,1)).^2+(fullBeeData(i).sphereIntersect(:,2)-refSphere.X(j,2)).^2+...
                      (fullBeeData(i).sphereIntersect(:,3)-refSphere.X(j,3)).^2)/2);  

            deltaAng2 = 2*asin(sqrt((fullBeeData(i).sphereIntersectBorders(:,1)-refSphere.X(j,1)).^2+(fullBeeData(i).sphereIntersectBorders(:,2)-refSphere.X(j,2)).^2+...
                      (fullBeeData(i).sphereIntersectBorders(:,3)-refSphere.X(j,3)).^2)/2);  
            %wow vectorized asin is sooooo much faster than sequential atan vector calcs   

            %store closest ratio of dist to interO to nearerst main for use in trimming
            %if min angle used directly, removal is biased to points in low res fov
            [minD, minI] = min(deltaAng);
            distToTargets(j) = minD/fullBeeData(i).calInterO(minI);
            
            %test if angular differences are less than IO 2x of given point
            %so essentially, world points within +1 IO will be included
            %2 is a bit of an arbitrary choice, but appears to catch all points between test ones
            if sum(deltaAng < fullBeeData(i).calInterO*2) > 0 | sum(deltaAng2 < expectedIO*2) > 0
                fullBeeData(i).inFOV(j) = 1;
            end
        end

        %this is the option for multi region fov
%       [leftMapResults, numLeftReg] = calculateRegionsInMap( refSphere.X, fullBeeData(i).inFOV, minDst*4);

        %this is only to deal with single region fov - fills holes
        [ regionIds, regLengths, regVals]  = findRegionsOnMesh( refSphere.X, fullBeeData(i).inFOV, minDst*4);
        outsideInds = find(regVals == 0); [~,maxI] = max(regLengths(outsideInds));
        fullBeeData(i).inFOV(:) = 1; fullBeeData(i).inFOV(regionIds == outsideInds(maxI)) = 0;

        %get coords of sphere border
        [~, borderInds] = plotSortedIsoLine(refSphere.X, fullBeeData(i).inFOV, 0.5,0,minDst*4);
        
        if trimSphere
            %as ratio is currently stored, balance point removal between high and low res fov
            
            %go through and remove border points until smaller fov than expected from calc
            while (sum(fullBeeData(i).inFOV)/size(refSphere.X,1)) > fullBeeData(i).fovFilledFromIntegral
                if contractBordersToTrim
                    %remover border ind further from actual point first
                    %this is guaranteed not to leave holes in fov
                    [~,maxI] = max(distToTargets(borderInds));
                    fullBeeData(i).inFOV(borderInds(maxI)) = 0;
                    borderInds(maxI) = [];

                    if isempty(borderInds)
                        %no more border, get next layer
                        [~, borderInds] = plotSortedIsoLine(refSphere.X, fullBeeData(i).inFOV, 0.5,0,minDst*3);
                    end
                else
                    %this could well leave holes in fov
                    %points will probably end up being removed from regions with largest interO first
                    indsT = find(fullBeeData(i).inFOV);
                    [~,maxI] = max(distToTargets(indsT));
                    fullBeeData(i).inFOV(indsT(maxI)) = 0;
                end
                if sum(fullBeeData(i).inFOV) == 0
                   error('entire fov deleted!'); 
                end
            end
            
            if ~contractBordersToTrim
                %fill holes again... will increase fov a bit, but wtf
                [ regionIds, regLengths, regVals] = findRegionsOnMesh( refSphere.X, fullBeeData(i).inFOV, minDst*4);
                outsideInds = find(regVals == 0); [~,maxI] = max(regLengths(outsideInds));
                fullBeeData(i).inFOV(:) = 1; fullBeeData(i).inFOV(regionIds == outsideInds(maxI)) = 0;
            end
        else
            %border defines aproximate limit, but remove it from FOV
            fullBeeData(i).inFOV(borderInds) = 0;
        end

        fullBeeData(i).inFOV = find(fullBeeData(i).inFOV);

        %get opposing fov with transform for opposite eye
        tempSphereInts = refSphere.X(fullBeeData(i).inFOV,:);

       %firstly need to remove extraRotAdj applied to original
        tempSphereInts = (extraRotAd^-1*[tempSphereInts zeros(size(tempSphereInts,1),1)]')';
       %then flip x
        tempSphereInts(:,1) = -tempSphereInts(:,1);
         
        %apply transforms
        %then apply eye mirroring transform
        tempSphereInts = (tempSphereInts*mirrorEyVolTransOrig')*extraRotAd';
        tempSphereInts = tempSphereInts(:,1:3);
        
        %correct for scale changes
        for j = 1:size(tempSphereInts,1)
            tempSphereInts(j,:) = tempSphereInts(j,:)/norm(tempSphereInts(j,:));
        end
   
        %snap transformed points back onto sphere
        oppositeFOV = zeros(length(fullBeeData(i).inFOV),1);
        for j = 1:length(fullBeeData(i).inFOV)
            [~, oppositeFOV(j)] = min(sqrt((refSphere.X(:,1)-tempSphereInts(j,1)).^2+...
                    (refSphere.X(:,2)-tempSphereInts(j,2)).^2+(refSphere.X(:,3)-tempSphereInts(j,3)).^2));
        end
        
        %not sure this is guaranteed to be closed (some holes could pop up)
        tempFOV = zeros(size(refSphere.X,1),1); 
        tempFOV(oppositeFOV) = 1;
        [ regionIds, regLengths, regVals]  = findRegionsOnMesh( refSphere.X, tempFOV, minDst*4);
        outsideInds = find(regVals == 0); [~,maxI] = max(regLengths(outsideInds));
        tempFOV(:) = 1; tempFOV(regionIds == outsideInds(maxI)) = 0;
        oppositeFOV = find(tempFOV);
        
        if fullBeeData(i).LeftEyeOriginal
           fullBeeData(i).inFOVRight = oppositeFOV;
        else
           %need to swap so inFOV is always left
           temp = fullBeeData(i).inFOV;
           fullBeeData(i).inFOV = oppositeFOV;
           fullBeeData(i).inFOVRight = temp;
           
           %need to adjust sphere intersects here after flip, not ideal :/
           fullBeeData(i).sphereIntersect = (extraRotAd^-1*[fullBeeData(i).sphereIntersect zeros(size(fullBeeData(i).sphereIntersect,1),1)]')';
           fullBeeData(i).sphereIntersect(:,1) = -fullBeeData(i).sphereIntersect(:,1);
           fullBeeData(i).sphereIntersect = (fullBeeData(i).sphereIntersect*mirrorEyVolTransOrig')*extraRotAd';
           fullBeeData(i).sphereIntersect = fullBeeData(i).sphereIntersect(:,1:3);
            
           fullBeeData(i).sphereIntersectBorders = (extraRotAd^-1*[fullBeeData(i).sphereIntersectBorders zeros(size(fullBeeData(i).sphereIntersectBorders,1),1)]')';
           fullBeeData(i).sphereIntersectBorders(:,1) = -fullBeeData(i).sphereIntersectBorders(:,1);
           fullBeeData(i).sphereIntersectBorders = (fullBeeData(i).sphereIntersectBorders*mirrorEyVolTransOrig')*extraRotAd';
           fullBeeData(i).sphereIntersectBorders = fullBeeData(i).sphereIntersectBorders(:,1:3);
        end

        %get final border after swap       
        tempFOV = zeros(size(refSphere.X,1),1); tempFOV(fullBeeData(i).inFOV) = 1;
        fullBeeData(i).inFOVBorderSubs = plotSortedIsoLine(refSphere.X, tempFOV, 0.5,0,minDst*4);
        fullBeeData(i).inFOVBorderSubs = sortPointsWithTSP(fullBeeData(i).inFOVBorderSubs);
        fullBeeData(i).inFOVBorderSubs = vertcat(fullBeeData(i).inFOVBorderSubs,fullBeeData(i).inFOVBorderSubs(1,:)); %close line...       
        
        fullBeeData(i).inFOVBino = intersect(fullBeeData(i).inFOV, fullBeeData(i).inFOVRight);
        
        %get closest point in FOV to straight ahead (neg z)
        [~, fullBeeData(i).frontalPoint] = min(2*asin(sqrt((refSphere.X(fullBeeData(i).inFOV,1)-0).^2+(refSphere.X(fullBeeData(i).inFOV,2)-0).^2+(refSphere.X(fullBeeData(i).inFOV,3)-(-1)).^2)/2));

        fullBeeData(i).fovFilledFromPoints = length(fullBeeData(i).inFOV)/size(refSphere.X,1);
        
        fprintf('fov on points %.2f, fov from calc %.2f\n', fullBeeData(i).fovFilledFromPoints, fullBeeData(i).fovFilledFromIntegral);

        if testPlots
            if ishandle(fullBeeData(i).testFOV)
                close(fullBeeData(i).testFOV)
            end
            fullBeeData(i).testFOV = figure; 
            
            subplot(1,2,1); axis equal; hold on; view([0, -90])
            plot3(fullBeeData(i).sphereIntersect(:,1), fullBeeData(i).sphereIntersect(:,2), fullBeeData(i).sphereIntersect(:,3), 'r.')
            plot3(refSphere.X(fullBeeData(i).inFOV,1), refSphere.X(fullBeeData(i).inFOV,2), refSphere.X(fullBeeData(i).inFOV,3), 'bo')

            plot3(fullBeeData(i).inFOVBorderSubs(:,1), fullBeeData(i).inFOVBorderSubs(:,2), fullBeeData(i).inFOVBorderSubs(:,3), 'g-')

            plot3(refSphere.X(fullBeeData(i).inFOV(fullBeeData(i).frontalPoint),1), refSphere.X(fullBeeData(i).inFOV(fullBeeData(i).frontalPoint),2),...
                    refSphere.X(fullBeeData(i).inFOV(fullBeeData(i).frontalPoint),3), 'g*','markersize',10)
            
            patch('Faces',refSphere.Triangulation,'Vertices',refSphere.X,'FaceAlpha',0.2,'EdgeAlpha',0,'FaceColor',[0.8 0.8 0.8]);
        
            title(sprintf('%s',fullBeeData(i).BeeID));
            
            subplot(1,2,2); hold on; axis equal
            plot3(refSphere.X(fullBeeData(i).inFOV,1), refSphere.X(fullBeeData(i).inFOV,2), refSphere.X(fullBeeData(i).inFOV,3), 'b+')
            plot3(refSphere.X(fullBeeData(i).inFOVRight,1), refSphere.X(fullBeeData(i).inFOVRight,2), refSphere.X(fullBeeData(i).inFOVRight,3), 'm*')
            plot3(refSphere.X(fullBeeData(i).inFOVBino,1), refSphere.X(fullBeeData(i).inFOVBino,2), refSphere.X(fullBeeData(i).inFOVBino,3), 'ro')
            view([0, -90])
            title('points in left right and bino');
        end
        
        %% interpolate all the parameters onto the world
        fullBeeData(i).diameterOnSphere = zeros(size(fullBeeData(i).inFOV,1),1);
        fullBeeData(i).curvatureOnSphere = zeros(size(fullBeeData(i).inFOV,1),1);

        fullBeeData(i).InterOOnSphere = zeros(size(fullBeeData(i).inFOV,1),1);
        fullBeeData(i).AxisDensityOnSphere = zeros(size(fullBeeData(i).inFOV,1),1);
        fullBeeData(i).areaProjectionOnSphere = zeros(size(fullBeeData(i).inFOV,1),1);
        fullBeeData(i).facetAreaOnSphere = zeros(size(fullBeeData(i).inFOV,1),1);
        fullBeeData(i).SensitvityOnSphere = zeros(size(fullBeeData(i).inFOV,1),1);
        
        fullBeeData(i).InterOWorldOnSphere = zeros(size(fullBeeData(i).inFOV,1),2);
        fullBeeData(i).InterORatioOnSphere = zeros(size(fullBeeData(i).inFOV,1),1);

        fullBeeData(i).lensTOnSphere = zeros(size(fullBeeData(i).inFOV,1),1);
        fullBeeData(i).CCTOnSphere = zeros(size(fullBeeData(i).inFOV,1),1);
        fullBeeData(i).retTOnSphere = zeros(size(fullBeeData(i).inFOV,1),1);
        
        %interp paramters onto world using IDW for angular/circle distance

        %option to include border points based on their nearest neighbor
        if bordersFromNearest
            %don't do last as its a duplicate of first
            quickNeighborLink = zeros(length(fullBeeData(i).inFOVBorderSubs)-1,1);
            %get the nearest normal neighbor of each border point
            for j = 1:length(fullBeeData(i).inFOVBorderSubs)-1
               [~,quickNeighborLink(j)] = min(2*asin(sqrt((fullBeeData(i).sphereIntersect(:,1)-fullBeeData(i).inFOVBorderSubs(j,1)).^2+(fullBeeData(i).sphereIntersect(:,2)-fullBeeData(i).inFOVBorderSubs(j,2)).^2+...
                      (fullBeeData(i).sphereIntersect(:,3)-fullBeeData(i).inFOVBorderSubs(j,3)).^2)/2));
            end
        end    
            
        %%% would be good just to make IDW function sometime, this is getting long
        for j = 1:size(fullBeeData(i).inFOV,1)
            %get angle to each sphere point - note for r = 1 this equal the great circle distance
            deltaAng = 2*asin(sqrt((fullBeeData(i).sphereIntersect(:,1)-refSphere.X(fullBeeData(i).inFOV(j),1)).^2+(fullBeeData(i).sphereIntersect(:,2)-refSphere.X(fullBeeData(i).inFOV(j),2)).^2+...
                      (fullBeeData(i).sphereIntersect(:,3)-refSphere.X(fullBeeData(i).inFOV(j),3)).^2)/2);  
            if bordersFromNearest
                 deltaAngBorder = 2*asin(sqrt((fullBeeData(i).inFOVBorderSubs(1:end-1,1)-refSphere.X(fullBeeData(i).inFOV(j),1)).^2+(fullBeeData(i).inFOVBorderSubs(1:end-1,2)-refSphere.X(fullBeeData(i).inFOV(j),2)).^2+...
                      (fullBeeData(i).inFOVBorderSubs(1:end-1,3)-refSphere.X(fullBeeData(i).inFOV(j),3)).^2)/2);  
            end
                  
            %%%improvement for borders could be to use nearest neihbors just for border pix intersects, and then use these in main interp. 
            %otherwise clear drop down to mean on borders
                  
            if ~nearestOnWorld
                %test if point at loc
                ind = find(deltaAng == 0);
                indB = find(deltaAngBorder == 0);
                
                if ~isempty(ind) | ~isempty(indB)
                    if ~isempty(ind)
                        ind2Get = ind;
                    elseif ~isempty(indB)
                        ind2Get = quickNeighborLink(indB); %two inds for last point
                        %could both not be empty???
                    end
                    %border thing will not matter here
                    fullBeeData(i).diameterOnSphere(j) = fullBeeData(i).avgDiameter(ind2Get);
                    fullBeeData(i).curvatureOnSphere(j) = fullBeeData(i).avgRadius(ind2Get);
                    fullBeeData(i).InterOOnSphere(j) = fullBeeData(i).calInterO(ind2Get);

                    fullBeeData(i).InterOWorldOnSphere(j,1) = fullBeeData(i).WorldIO(ind2Get,1);
                    fullBeeData(i).InterOWorldOnSphere(j,2) = fullBeeData(i).WorldIO(ind2Get,2);

                    fullBeeData(i).InterORatioOnSphere(j) = fullBeeData(i).worldRatio(ind2Get);
                    if useDensityFromIntegral
                        fullBeeData(i).AxisDensityOnSphere(j) = fullBeeData(i).densityFromArea(ind2Get);
                    else
                        fullBeeData(i).AxisDensityOnSphere(j) = fullBeeData(i).AxisDensity(ind2Get);
                    end
                    fullBeeData(i).areaProjectionOnSphere(j) = 1./fullBeeData(i).projectedAngle(ind2Get);
                    fullBeeData(i).facetAreaOnSphere(j) = fullBeeData(i).interpFacetArea(ind2Get);
                    fullBeeData(i).SensitvityOnSphere(j) = fullBeeData(i).sensitvityApproximation(ind2Get);
                    
                    fullBeeData(i).lensTOnSphere(j) = fullBeeData(i).lensThickness(ind2Get);
                    fullBeeData(i).CCTOnSphere(j) = fullBeeData(i).CCThickness(ind2Get);
                    fullBeeData(i).retTOnSphere(j) = fullBeeData(i).RetThickness(ind2Get);
                else
                    if bordersFromNearest
                        weightsBorder = 1./(deltaAngBorder.^powerPonSphere); 
                    else
                        weightsBorder = deltaAngBorder*0;
                    end
                    weights = 1./(deltaAng.^powerPonSphere);
                    
                    fullBeeData(i).diameterOnSphere(j) = (sum(fullBeeData(i).avgDiameter.*weights)+sum(fullBeeData(i).avgDiameter(quickNeighborLink).*weightsBorder))/(sum(weights)+sum(weightsBorder));
                    fullBeeData(i).curvatureOnSphere(j) = (sum(fullBeeData(i).avgRadius.*weights)+sum(fullBeeData(i).avgRadius(quickNeighborLink).*weightsBorder))/(sum(weights)+sum(weightsBorder));
                    fullBeeData(i).InterOOnSphere(j) = (sum(fullBeeData(i).calInterO.*weights)+sum(fullBeeData(i).calInterO(quickNeighborLink).*weightsBorder))/(sum(weights)+sum(weightsBorder));

                    fullBeeData(i).InterOWorldOnSphere(j,1) = (sum(fullBeeData(i).WorldIO(:,1).*weights)+sum(fullBeeData(i).WorldIO(quickNeighborLink,1).*weightsBorder))/(sum(weights)+sum(weightsBorder));
                    fullBeeData(i).InterOWorldOnSphere(j,2) = (sum(fullBeeData(i).WorldIO(:,2).*weights)+sum(fullBeeData(i).WorldIO(quickNeighborLink,2).*weightsBorder))/(sum(weights)+sum(weightsBorder));

                    fullBeeData(i).InterORatioOnSphere(j) = (sum(fullBeeData(i).worldRatio.*weights)+sum(fullBeeData(i).worldRatio(quickNeighborLink).*weightsBorder))/(sum(weights)+sum(weightsBorder));
                    if useDensityFromIntegral
                        fullBeeData(i).AxisDensityOnSphere(j) = (sum(fullBeeData(i).densityFromArea.*weights)+sum(fullBeeData(i).densityFromArea(quickNeighborLink).*weightsBorder))/(sum(weights)+sum(weightsBorder));
                    else
                        fullBeeData(i).AxisDensityOnSphere(j) = (sum(fullBeeData(i).AxisDensity.*weights)+sum(fullBeeData(i).AxisDensity(quickNeighborLink).*weightsBorder))/(sum(weights)+sum(weightsBorder));
                    end
                    fullBeeData(i).areaProjectionOnSphere(j) = 1./((sum(fullBeeData(i).projectedAngle.*weights)+sum(fullBeeData(i).projectedAngle(quickNeighborLink).*weightsBorder))/(sum(weights)+sum(weightsBorder)));
                    fullBeeData(i).facetAreaOnSphere(j) = (sum(fullBeeData(i).interpFacetArea.*weights)+sum(fullBeeData(i).interpFacetArea(quickNeighborLink).*weightsBorder))/(sum(weights)+sum(weightsBorder));
                    fullBeeData(i).SensitvityOnSphere(j) = (sum(fullBeeData(i).sensitvityApproximation.*weights)+sum(fullBeeData(i).sensitvityApproximation(quickNeighborLink).*weightsBorder))/(sum(weights)+sum(weightsBorder));
                    
                    fullBeeData(i).lensTOnSphere(j) = (sum(fullBeeData(i).lensThickness.*weights)+sum(fullBeeData(i).lensThickness(quickNeighborLink).*weightsBorder))/(sum(weights)+sum(weightsBorder));
                    fullBeeData(i).CCTOnSphere(j) = (sum(fullBeeData(i).CCThickness.*weights)+sum(fullBeeData(i).CCThickness(quickNeighborLink).*weightsBorder))/(sum(weights)+sum(weightsBorder));
                    fullBeeData(i).retTOnSphere(j) = (sum(fullBeeData(i).RetThickness.*weights)+sum(fullBeeData(i).RetThickness(quickNeighborLink).*weightsBorder))/(sum(weights)+sum(weightsBorder));
                end
            else
                %this just doesn't use borders, as closest point will be the same
                [~, ind] = min(deltaAng);
                fullBeeData(i).diameterOnSphere(j) = fullBeeData(i).avgDiameter(ind);
                fullBeeData(i).curvatureOnSphere(j) = fullBeeData(i).avgRadius(ind);
                fullBeeData(i).InterOOnSphere(j) = fullBeeData(i).calInterO(ind);

                fullBeeData(i).InterOWorldOnSphere(j,1) = fullBeeData(i).WorldIO(ind,1);
                fullBeeData(i).InterOWorldOnSphere(j,2) = fullBeeData(i).WorldIO(ind,2);

                fullBeeData(i).InterORatioOnSphere(j) = fullBeeData(i).worldRatio(ind);
                if useDensityFromIntegral
                    fullBeeData(i).AxisDensityOnSphere(j) = fullBeeData(i).densityFromArea(ind);
                else
                    fullBeeData(i).AxisDensityOnSphere(j) = fullBeeData(i).AxisDensity(ind);
                end
                fullBeeData(i).areaProjectionOnSphere(j) = 1./fullBeeData(i).projectedAngle(ind);
                fullBeeData(i).facetAreaOnSphere(j) = fullBeeData(i).interpFacetArea(ind);
                fullBeeData(i).SensitvityOnSphere(j) = fullBeeData(i).sensitvityApproximation(ind);    
                
                fullBeeData(i).lensTOnSphere(j) = fullBeeData(i).lensThickness(ind);
                fullBeeData(i).CCTOnSphere(j) = fullBeeData(i).CCThickness(ind);
                fullBeeData(i).retTOnSphere(j) = fullBeeData(i).RetThickness(ind);
            end
        end

        %test IDW validity - just if range exceeded, not if too small
        if sum(fullBeeData(i).diameterOnSphere > max(fullBeeData(i).avgDiameter)*(1+IDWThresh)) |  sum(fullBeeData(i).diameterOnSphere < min(fullBeeData(i).avgDiameter)*(1-IDWThresh))
            [max(fullBeeData(i).diameterOnSphere) min(fullBeeData(i).diameterOnSphere) max(fullBeeData(i).avgDiameter) min(fullBeeData(i).avgDiameter)] 
            error('fullBeeData(i).diameterOnSphere exceeds threshold on bee %i', i)
        end
        if sum(fullBeeData(i).curvatureOnSphere > max(fullBeeData(i).avgRadius)*(1+IDWThresh)) |   sum(fullBeeData(i).curvatureOnSphere < min(fullBeeData(i).avgRadius)*(1-IDWThresh))
            [max(fullBeeData(i).curvatureOnSphere) min(fullBeeData(i).curvatureOnSphere) max(fullBeeData(i).avgRadius) min(fullBeeData(i).avgRadius)] 
            error('fullBeeData(i).curvatureOnSphere exceeds threshold on bee %i', i)
        end
         if sum(fullBeeData(i).InterOOnSphere > max(fullBeeData(i).calInterO)*(1+IDWThresh)) |  sum(fullBeeData(i).InterOOnSphere < min(fullBeeData(i).calInterO)*(1-IDWThresh))
             [max(fullBeeData(i).InterOOnSphere) min(fullBeeData(i).InterOOnSphere) max(fullBeeData(i).calInterO) min(fullBeeData(i).calInterO)] 
             error('fullBeeData(i).InterOOnSphere exceeds threshold on bee %i', i)
         end
        if sum(fullBeeData(i).InterOWorldOnSphere(:,1) > max(fullBeeData(i).WorldIO(:,1))*(1+IDWThresh)) |  sum(fullBeeData(i).InterOWorldOnSphere(:,1) < min(fullBeeData(i).WorldIO(:,1))*(1-IDWThresh))
            [max(fullBeeData(i).InterOWorldOnSphere(:,1)) min(fullBeeData(i).InterOWorldOnSphere(:,1)) max(fullBeeData(i).WorldIO(:,1)) min(fullBeeData(i).WorldIO(:,1))] 
            error('fullBeeData(i).InterOWorldOnSphere(:,1) exceeds threshold on bee %i', i)
        end
        if sum(fullBeeData(i).InterOWorldOnSphere(:,2) > max(fullBeeData(i).WorldIO(:,2))*(1+IDWThresh)) |  sum(fullBeeData(i).InterOWorldOnSphere(:,2) < min(fullBeeData(i).WorldIO(:,2))*(1-IDWThresh))
            [max(fullBeeData(i).InterOWorldOnSphere(:,2)) min(fullBeeData(i).InterOWorldOnSphere(:,2)) max(fullBeeData(i).WorldIO(:,2)) min(fullBeeData(i).WorldIO(:,2))] 
            error('fullBeeData(i).InterOWorldOnSphere(:,2) exceeds threshold on bee %i', i)
        end
        if sum(fullBeeData(i).InterORatioOnSphere > max(fullBeeData(i).worldRatio)*(1+IDWThresh)) |  sum(fullBeeData(i).InterORatioOnSphere < min(fullBeeData(i).worldRatio)*(1-IDWThresh))
           [max(fullBeeData(i).InterORatioOnSphere) min(fullBeeData(i).InterORatioOnSphere) max(fullBeeData(i).worldRatio) min(fullBeeData(i).worldRatio)] 
            error('fullBeeData(i).InterORatioOnSphere exceeds threshold on bee %i', i)
        end
        if useDensityFromIntegral
            if sum(fullBeeData(i).AxisDensityOnSphere > max(fullBeeData(i).densityFromArea)*(1+IDWThresh)) |  sum(fullBeeData(i).AxisDensityOnSphere < min(fullBeeData(i).densityFromArea)*(1-IDWThresh))
               [max(fullBeeData(i).AxisDensityOnSphere) min(fullBeeData(i).AxisDensityOnSphere) max(fullBeeData(i).densityFromArea) min(fullBeeData(i).densityFromArea)] 
                error('fullBeeData(i).AxisDensityOnSphere exceeds threshold on bee %i', i)
            end
        else
            if sum(fullBeeData(i).AxisDensityOnSphere > max(fullBeeData(i).AxisDensity)*(1+IDWThresh)) |  sum(fullBeeData(i).AxisDensityOnSphere < min(fullBeeData(i).AxisDensity)*(1-IDWThresh))
               [max(fullBeeData(i).AxisDensityOnSphere) min(fullBeeData(i).AxisDensityOnSphere) max(fullBeeData(i).AxisDensity) min(fullBeeData(i).AxisDensity)] 
                error('fullBeeData(i).AxisDensityOnSphere exceeds threshold on bee %i', i)
            end
        end
        if sum(fullBeeData(i).areaProjectionOnSphere > max(1./fullBeeData(i).projectedAngle)*(1+IDWThresh)) |  sum(fullBeeData(i).areaProjectionOnSphere < min(1./fullBeeData(i).projectedAngle)*(1-IDWThresh))
           [max(fullBeeData(i).areaProjectionOnSphere) min(fullBeeData(i).areaProjectionOnSphere) max(1./fullBeeData(i).projectedAngle) min(1./fullBeeData(i).projectedAngle)] 
            error('fullBeeData(i).areaProjectionOnSphere exceeds threshold on bee %i', i)
        end 
        if sum(fullBeeData(i).facetAreaOnSphere > max(fullBeeData(i).facetAreaOnSphere)*(1+IDWThresh)) |  sum(fullBeeData(i).interpFacetArea < min(fullBeeData(i).interpFacetArea)*(1-IDWThresh))
           [max(fullBeeData(i).facetAreaOnSphere) min(fullBeeData(i).facetAreaOnSphere) max(fullBeeData(i).interpFacetArea) min(fullBeeData(i).interpFacetArea)] 
            error('fullBeeData(i).facetAreaOnSphere exceeds threshold on bee %i', i)
        end 
        if sum(fullBeeData(i).SensitvityOnSphere > max(fullBeeData(i).SensitvityOnSphere)*(1+IDWThresh)) |  sum(fullBeeData(i).sensitvityApproximation < min(fullBeeData(i).sensitvityApproximation)*(1-IDWThresh))
           [max(fullBeeData(i).SensitvityOnSphere) min(fullBeeData(i).SensitvityOnSphere) max(fullBeeData(i).sensitvityApproximation) min(fullBeeData(i).sensitvityApproximation)] 
            error('fullBeeData(i).SensitvityOnSphere exceeds threshold on bee %i', i)
        end 
        if sum(fullBeeData(i).lensTOnSphere > max(fullBeeData(i).lensThickness)*(1+IDWThresh)) |  sum(fullBeeData(i).lensTOnSphere < min(fullBeeData(i).lensThickness)*(1-IDWThresh))
           [max(fullBeeData(i).lensTOnSphere) min(fullBeeData(i).lensTOnSphere) max(fullBeeData(i).lensThickness) min(fullBeeData(i).lensThickness)] 
            error('fullBeeData(i).lensTOnSphere exceeds threshold on bee %i', i)
        end
        if sum(fullBeeData(i).CCTOnSphere > max(fullBeeData(i).CCThickness)*(1+IDWThresh)) |  sum(fullBeeData(i).CCTOnSphere < min(fullBeeData(i).CCThickness)*(1-IDWThresh))
           [max(fullBeeData(i).CCTOnSphere) min(fullBeeData(i).CCTOnSphere) max(fullBeeData(i).CCThickness) min(fullBeeData(i).CCThickness)] 
            error('fullBeeData(i).CCTOnSphere exceeds threshold on bee %i', i)
        end
        if sum(fullBeeData(i).retTOnSphere > max(fullBeeData(i).RetThickness)*(1+IDWThresh)) |  sum(fullBeeData(i).retTOnSphere < min(fullBeeData(i).RetThickness)*(1-IDWThresh))
           [max(fullBeeData(i).retTOnSphere) min(fullBeeData(i).retTOnSphere) max(fullBeeData(i).RetThickness) min(fullBeeData(i).RetThickness)] 
            error('fullBeeData(i).retTOnSphere exceeds threshold on bee %i', i)
        end
        
        if testPlots
            if ishandle(fullBeeData(i).testOnWorldF)
                close(fullBeeData(i).testOnWorldF)
            end
            fullBeeData(i).testOnWorldF = figure; 
            
            subplot(2,3,1);axis equal; hold on;
            bP = mean(fullBeeData(i).inFOVBorderSubs);

            tempVals = (fullBeeData(i).diameterOnSphere-diamRng(1))./(diamRng(2)-diamRng(1))*100; 
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            plot3(fullBeeData(i).inFOVBorderSubs(:,1), fullBeeData(i).inFOVBorderSubs(:,2), fullBeeData(i).inFOVBorderSubs(:,3), 'g-')
            [~,hCB] = fscatter3([refSphere.X(fullBeeData(i).inFOV,1)' bP(1) bP(1)],...
                [refSphere.X(fullBeeData(i).inFOV,2)' bP(2) bP(2)], ...
                [refSphere.X(fullBeeData(i).inFOV,3)' bP(3) bP(3)], [tempVals' 0 100], jet(100),20);
            set(hCB, 'Ticks', [0 50 100], 'TickLabels', diamLabs,'YAxisLocation','right','TickDirection','out');
            plot3(bP(1), bP(2), bP(3), 'w.','markersize',50);
            patch('Faces',refSphere.Triangulation,'Vertices',refSphere.X,'FaceAlpha',0.2,'EdgeAlpha',0,'FaceColor',[0.8 0.8 0.8]);
            view([0, -90]); title(sprintf('%s %i facet D',fullBeeData(i).BeeID, round(fullBeeData(i).EyeVolume/10^3)));

            subplot(2,3,2);axis equal; hold on; view([0, -90]);
            tempVals = (fullBeeData(i).curvatureOnSphere-curvRng(1))./(curvRng(2)-curvRng(1))*100;
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            plot3(fullBeeData(i).inFOVBorderSubs(:,1), fullBeeData(i).inFOVBorderSubs(:,2), fullBeeData(i).inFOVBorderSubs(:,3), 'g-')
            [~,hCB] = fscatter3([refSphere.X(fullBeeData(i).inFOV,1)' bP(1) bP(1)],...
                [refSphere.X(fullBeeData(i).inFOV,2)' bP(2) bP(2)],...
                [refSphere.X(fullBeeData(i).inFOV,3)' bP(3) bP(3)], [tempVals' 0 100], jet(100),20);
            set(hCB, 'Ticks', [0 50 100], 'TickLabels', curvLabs,'YAxisLocation','right','TickDirection','out');
            plot3(bP(1), bP(2), bP(3), 'w.','markersize',50);
            patch('Faces',refSphere.Triangulation,'Vertices',refSphere.X,'FaceAlpha',0.2,'EdgeAlpha',0,'FaceColor',[0.8 0.8 0.8]);
            view([0, -90]); title('curvature');

            subplot(2,3,3);axis equal; hold on; view([0, -90]);
            tempVals = (fullBeeData(i).AxisDensityOnSphere-axDRng(1))./(axDRng(2)-axDRng(1))*100;
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            plot3(fullBeeData(i).inFOVBorderSubs(:,1), fullBeeData(i).inFOVBorderSubs(:,2), fullBeeData(i).inFOVBorderSubs(:,3), 'g-')
            [~,hCB] = fscatter3([refSphere.X(fullBeeData(i).inFOV,1)' bP(1) bP(1)],...
                [refSphere.X(fullBeeData(i).inFOV,2)' bP(2) bP(2)],...
                [refSphere.X(fullBeeData(i).inFOV,3)' bP(3) bP(3)], [tempVals' 0 100], jet(100),20);
            set(hCB, 'Ticks', [0 50 100], 'TickLabels', axDLabs,'YAxisLocation','right','TickDirection','out');
            plot3(bP(1), bP(2), bP(3), 'w.','markersize',50);
            patch('Faces',refSphere.Triangulation,'Vertices',refSphere.X,'FaceAlpha',0.2,'EdgeAlpha',0,'FaceColor',[0.8 0.8 0.8]);
            view([0, -90]);
            title(sprintf('facet density min %.1f - max %.1f', min(fullBeeData(i).AxisDensityOnSphere), max(fullBeeData(i).AxisDensityOnSphere)));

            subplot(2,3,4);axis equal; hold on; 
            plot3(fullBeeData(i).inFOVBorderSubs(:,1), fullBeeData(i).inFOVBorderSubs(:,2), fullBeeData(i).inFOVBorderSubs(:,3), 'g-')
            tempVals = (fullBeeData(i).InterOWorldOnSphere(:,1)/pi*180-intORng(1))./(intORng(2)-intORng(1))*100;
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            [~,hCB] = fscatter3([refSphere.X(fullBeeData(i).inFOV,1)' bP(1) bP(1)],...
                [refSphere.X(fullBeeData(i).inFOV,2)' bP(2) bP(2)],...
                [refSphere.X(fullBeeData(i).inFOV,3)' bP(3) bP(3)], [tempVals' 0 100], jet(100),20);
            set(hCB, 'Ticks', [0 50 100], 'TickLabels', intOLabs,'YAxisLocation','right','TickDirection','out');
            plot3(bP(1), bP(2), bP(3), 'w.','markersize',50);
            patch('Faces',refSphere.Triangulation,'Vertices',refSphere.X,'FaceAlpha',0.2,'EdgeAlpha',0,'FaceColor',[0.8 0.8 0.8]);
            view([0, -90]); title('partio IO world el');

            subplot(2,3,5);axis equal; hold on; 
            plot3(fullBeeData(i).inFOVBorderSubs(:,1), fullBeeData(i).inFOVBorderSubs(:,2), fullBeeData(i).inFOVBorderSubs(:,3), 'g-')
            tempVals = (fullBeeData(i).InterOWorldOnSphere(:,2)/pi*180-intORng(1))./(intORng(2)-intORng(1))*100;
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            [~,hCB] = fscatter3([refSphere.X(fullBeeData(i).inFOV,1)' bP(1) bP(1)],...
                [refSphere.X(fullBeeData(i).inFOV,2)' bP(2) bP(2)],...
                [refSphere.X(fullBeeData(i).inFOV,3)' bP(3) bP(3)], [tempVals' 0 100], jet(100),20);
            set(hCB, 'Ticks', [0 50 100], 'TickLabels', intOLabs,'YAxisLocation','right','TickDirection','out');
            plot3(bP(1), bP(2), bP(3), 'w.','markersize',50);
            patch('Faces',refSphere.Triangulation,'Vertices',refSphere.X,'FaceAlpha',0.2,'EdgeAlpha',0,'FaceColor',[0.8 0.8 0.8]);
            view([0, -90]); title('partial IO world az');

            subplot(2,3,6);axis equal; hold on; 
            plot3(fullBeeData(i).inFOVBorderSubs(:,1), fullBeeData(i).inFOVBorderSubs(:,2), fullBeeData(i).inFOVBorderSubs(:,3), 'g-')
            tempVals = (fullBeeData(i).InterOOnSphere/pi*180-intORng(1))./(intORng(2)-intORng(1))*100;
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            [~,hCB] = fscatter3([refSphere.X(fullBeeData(i).inFOV,1)' bP(1) bP(1)],...
                [refSphere.X(fullBeeData(i).inFOV,2)' bP(2) bP(2)],...
                [refSphere.X(fullBeeData(i).inFOV,3)' bP(3) bP(3)], [tempVals' 0 100], jet(100),20);
            set(hCB, 'Ticks', [0 50 100], 'TickLabels', intOLabs,'YAxisLocation','right','TickDirection','out');
            plot3(bP(1), bP(2), bP(3), 'w.','markersize',50);
            patch('Faces',refSphere.Triangulation,'Vertices',refSphere.X,'FaceAlpha',0.2,'EdgeAlpha',0,'FaceColor',[0.8 0.8 0.8]);
            view([0, -90]); title('IO angle')
            
            if ishandle(fullBeeData(i).testOnWorldFExtra)
                close(fullBeeData(i).testOnWorldFExtra)
            end
            fullBeeData(i).testOnWorldFExtra = figure; 
            
            subplot(2,3,1); axis equal; hold on; 
            plot3(fullBeeData(i).inFOVBorderSubs(:,1), fullBeeData(i).inFOVBorderSubs(:,2), fullBeeData(i).inFOVBorderSubs(:,3), 'g-')
            tempVals = (fullBeeData(i).lensTOnSphere-lensTRng(1))./(lensTRng(2)-lensTRng(1))*100;
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            [~,hCB] = fscatter3([refSphere.X(fullBeeData(i).inFOV,1)' bP(1) bP(1)],...
                [refSphere.X(fullBeeData(i).inFOV,2)' bP(2) bP(2)],...
                [refSphere.X(fullBeeData(i).inFOV,3)' bP(3) bP(3)], [tempVals' 0 100], jet(100),20);
            set(hCB, 'Ticks', [0 50 100], 'TickLabels', lensTLabs,'YAxisLocation','right','TickDirection','out');
            plot3(bP(1), bP(2), bP(3), 'w.','markersize',50);
            patch('Faces',refSphere.Triangulation,'Vertices',refSphere.X,'FaceAlpha',0.2,'EdgeAlpha',0,'FaceColor',[0.8 0.8 0.8]);
            view([0, -90]); title(sprintf('%s lens thickness on world',fullBeeData(i).BeeID));
            
            subplot(2,3,2); axis equal; hold on; 
            plot3(fullBeeData(i).inFOVBorderSubs(:,1), fullBeeData(i).inFOVBorderSubs(:,2), fullBeeData(i).inFOVBorderSubs(:,3), 'g-')
            tempVals = (fullBeeData(i).CCTOnSphere-CCTRng(1))./(CCTRng(2)-CCTRng(1))*100;
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            [~,hCB] = fscatter3([refSphere.X(fullBeeData(i).inFOV,1)' bP(1) bP(1)],...
                [refSphere.X(fullBeeData(i).inFOV,2)' bP(2) bP(2)],...
                [refSphere.X(fullBeeData(i).inFOV,3)' bP(3) bP(3)], [tempVals' 0 100], jet(100),20);
            set(hCB, 'Ticks', [0 50 100], 'TickLabels', CCTLabs,'YAxisLocation','right','TickDirection','out');
            plot3(bP(1), bP(2), bP(3), 'w.','markersize',50);
            patch('Faces',refSphere.Triangulation,'Vertices',refSphere.X,'FaceAlpha',0.2,'EdgeAlpha',0,'FaceColor',[0.8 0.8 0.8]);
            view([0, -90]); title('cc thickness on world');
            
            subplot(2,3,3); axis equal; hold on; 
            plot3(fullBeeData(i).inFOVBorderSubs(:,1), fullBeeData(i).inFOVBorderSubs(:,2), fullBeeData(i).inFOVBorderSubs(:,3), 'g-')
            tempVals = (fullBeeData(i).retTOnSphere-RetTRng(1))./(RetTRng(2)-RetTRng(1))*100;
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            [~,hCB] = fscatter3([refSphere.X(fullBeeData(i).inFOV,1)' bP(1) bP(1)],...
                [refSphere.X(fullBeeData(i).inFOV,2)' bP(2) bP(2)],...
                [refSphere.X(fullBeeData(i).inFOV,3)' bP(3) bP(3)], [tempVals' 0 100], jet(100),20);
            set(hCB, 'Ticks', [0 50 100], 'TickLabels', RetTLabs,'YAxisLocation','right','TickDirection','out');
            plot3(bP(1), bP(2), bP(3), 'w.','markersize',50);
            patch('Faces',refSphere.Triangulation,'Vertices',refSphere.X,'FaceAlpha',0.2,'EdgeAlpha',0,'FaceColor',[0.8 0.8 0.8]);
            view([0, -90]); title('retina thickness on world');
            
            subplot(2,3,4); axis equal; hold on; 
            plot3(fullBeeData(i).inFOVBorderSubs(:,1), fullBeeData(i).inFOVBorderSubs(:,2), fullBeeData(i).inFOVBorderSubs(:,3), 'g-')
            tempVals = (fullBeeData(i).SensitvityOnSphere-sensRng(1))./(sensRng(2)-sensRng(1))*100;
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            [~,hCB] = fscatter3([refSphere.X(fullBeeData(i).inFOV,1)' bP(1) bP(1)],...
                [refSphere.X(fullBeeData(i).inFOV,2)' bP(2) bP(2)],...
                [refSphere.X(fullBeeData(i).inFOV,3)' bP(3) bP(3)], [tempVals' 0 100], jet(100),20);
            set(hCB, 'Ticks', [0 50 100], 'TickLabels', sensLabs,'YAxisLocation','right','TickDirection','out');
            plot3(bP(1), bP(2), bP(3), 'w.','markersize',50);
            patch('Faces',refSphere.Triangulation,'Vertices',refSphere.X,'FaceAlpha',0.2,'EdgeAlpha',0,'FaceColor',[0.8 0.8 0.8]);
            view([0, -90]); title('sensitvity');
            
            subplot(2,3,5); axis equal; hold on; 
            plot3(fullBeeData(i).inFOVBorderSubs(:,1), fullBeeData(i).inFOVBorderSubs(:,2), fullBeeData(i).inFOVBorderSubs(:,3), 'g-')
            tempVals = (fullBeeData(i).areaProjectionOnSphere-projARng(1))./(projARng(2)-projARng(1))*100;
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            [~,hCB] = fscatter3([refSphere.X(fullBeeData(i).inFOV,1)' bP(1) bP(1)],...
                [refSphere.X(fullBeeData(i).inFOV,2)' bP(2) bP(2)],...
                [refSphere.X(fullBeeData(i).inFOV,3)' bP(3) bP(3)], [tempVals' 0 100], jet(100),20);
            set(hCB, 'Ticks', [0 50 100], 'TickLabels', projALabs,'YAxisLocation','right','TickDirection','out');
            plot3(bP(1), bP(2), bP(3), 'w.','markersize',50);
            patch('Faces',refSphere.Triangulation,'Vertices',refSphere.X,'FaceAlpha',0.2,'EdgeAlpha',0,'FaceColor',[0.8 0.8 0.8]);
            view([0, -90]); title('area projection on sphere');
            
            subplot(2,3,6); axis equal; hold on; 
            plot3(fullBeeData(i).inFOVBorderSubs(:,1), fullBeeData(i).inFOVBorderSubs(:,2), fullBeeData(i).inFOVBorderSubs(:,3), 'g-')
            tempVals = (fullBeeData(i).facetAreaOnSphere-facARng(1))./(facARng(2)-facARng(1))*100;
            tempVals( tempVals < 1 ) = 1; tempVals( tempVals > 100 ) = 100;
            [~,hCB] = fscatter3([refSphere.X(fullBeeData(i).inFOV,1)' bP(1) bP(1)],...
                [refSphere.X(fullBeeData(i).inFOV,2)' bP(2) bP(2)],...
                [refSphere.X(fullBeeData(i).inFOV,3)' bP(3) bP(3)], [tempVals' 0 100], jet(100),20);
            set(hCB, 'Ticks', [0 50 100], 'TickLabels', facALabs,'YAxisLocation','right','TickDirection','out');
            plot3(bP(1), bP(2), bP(3), 'w.','markersize',50);
            patch('Faces',refSphere.Triangulation,'Vertices',refSphere.X,'FaceAlpha',0.2,'EdgeAlpha',0,'FaceColor',[0.8 0.8 0.8]);
            view([0, -90]); title('facet area on world');
        end

        %% integrate over sphere, just density for now´
        tempInFOVLinks = zeros(size(refSphere.X,1),1);
        tempInFOVLinks(fullBeeData(i).inFOV) = 1:length(fullBeeData(i).inFOV);

        fullBeeData(i).integratedAreaProj = 0;
        fullBeeData(i).integratedAngleDensity = 0;
        fullBeeData(i).integratedSensitvity = 0;
        
        for j = 1:size(refSphere.Triangulation,1)
            axisDensityOnTri = tempInFOVLinks(refSphere.Triangulation(j,:));
            projAreaOnTri = axisDensityOnTri;
            projSensOnTri = axisDensityOnTri;
            %if not zero it is an index to a value, so get the value
            projAreaOnTri(projAreaOnTri > 0) = fullBeeData(i).areaProjectionOnSphere(projAreaOnTri(projAreaOnTri > 0));
            axisDensityOnTri(axisDensityOnTri > 0) = fullBeeData(i).AxisDensityOnSphere(axisDensityOnTri(axisDensityOnTri > 0));
            projSensOnTri(projSensOnTri > 0) = fullBeeData(i).SensitvityOnSphere(projSensOnTri(projSensOnTri > 0));
            projSensOnTri = projSensOnTri.*axisDensityOnTri;
            
            p1 = refSphere.X(refSphere.Triangulation(j,1),:);
            p2 = refSphere.X(refSphere.Triangulation(j,2),:);
            p3 = refSphere.X(refSphere.Triangulation(j,3),:);
            tC = mean( refSphere.X(refSphere.Triangulation(j,:),:));
            tC = tC/norm(tC);

            %interpolate value at center - angular version of herons formula
            %https://classes.soe.ucsc.edu/cmps160/Fall10/resources/barycentricInterpolation.pdf
            %A1 to match p1
            thetaP2P3 = acos((p2(1)*p3(1)+p2(2)*p3(2)+p2(3)*p3(3))/sqrt((p2(1)^2+p2(2)^2+p2(3)^2)*(p3(1)^2+p3(2)^2+p3(3)^2))); 
            thetaCP2 = acos((tC(1)*p2(1)+tC(2)*p2(2)+tC(3)*p2(3))/sqrt((tC(1)^2+tC(2)^2+tC(3)^2)*(p2(1)^2+p2(2)^2+p2(3)^2)));
            thetaCP3 = acos((tC(1)*p3(1)+tC(2)*p3(2)+tC(3)*p3(3))/sqrt((tC(1)^2+tC(2)^2+tC(3)^2)*(p3(1)^2+p3(2)^2+p3(3)^2)));
            thetaS = (thetaP2P3+thetaCP2+thetaCP3)/2;
            areaCP2P3 = 4*atan(sqrt(tan(thetaS/2)*tan((thetaS-thetaP2P3)/2)*tan((thetaS-thetaCP2)/2)*tan((thetaS-thetaCP3)/2))); %*(180/pi)^2

            %A2 to match p2
            thetaP1P3 = acos((p1(1)*p3(1)+p1(2)*p3(2)+p1(3)*p3(3))/sqrt((p1(1)^2+p1(2)^2+p1(3)^2)*(p3(1)^2+p3(2)^2+p3(3)^2))); 
            thetaCP1 = acos((tC(1)*p1(1)+tC(2)*p1(2)+tC(3)*p1(3))/sqrt((tC(1)^2+tC(2)^2+tC(3)^2)*(p1(1)^2+p1(2)^2+p1(3)^2)));
            thetaS = (thetaP1P3+thetaCP1+thetaCP3)/2;
            areaCP1P3 = 4*atan(sqrt(tan(thetaS/2)*tan((thetaS-thetaP1P3)/2)*tan((thetaS-thetaCP1)/2)*tan((thetaS-thetaCP3)/2))); %*(180/pi)^2

            %A3 to match p3
            thetaP1P2 = acos((p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3))/sqrt((p1(1)^2+p1(2)^2+p1(3)^2)*(p2(1)^2+p2(2)^2+p2(3)^2))); 
            thetaS = (thetaP1P2+thetaCP1+thetaCP2)/2;
            areaCP1P2 = 4*atan(sqrt(tan(thetaS/2)*tan((thetaS-thetaP1P2)/2)*tan((thetaS-thetaCP1)/2)*tan((thetaS-thetaCP2)/2))); %*(180/pi)^2

            %calculate 
            fullBeeData(i).integratedAngleDensity = fullBeeData(i).integratedAngleDensity + (axisDensityOnTri(1)*areaCP2P3 + axisDensityOnTri(2)*areaCP1P3 + axisDensityOnTri(3)*areaCP1P2);
            fullBeeData(i).integratedAreaProj =  fullBeeData(i).integratedAreaProj + (projAreaOnTri(1)*areaCP2P3 + projAreaOnTri(2)*areaCP1P3 + projAreaOnTri(3)*areaCP1P2);
            fullBeeData(i).integratedSensitvity = fullBeeData(i).integratedSensitvity + (projSensOnTri(1)*areaCP2P3 + projSensOnTri(2)*areaCP1P3 + projSensOnTri(3)*areaCP1P2);
        end 
        
        %% test if fit is satisfactory
        %store current results
        fitResults(itFit) = fullBeeData(i).integratedAngleDensity./fullBeeData(i).NumFacets;
        surfacesTested(itFit) = surfaceFitRng;

        fprintf('current fit -for integral %.2f -given r %.1f (area int dif: %.2f, sens int dif: %.2f)\n', fitResults(itFit), surfaceFitRng, fullBeeData(i).integratedAreaProj/fullBeeData(i).lensSurfaceArea, fullBeeData(i).integratedSensitvity/fullBeeData(i).expectedSensInt);
        %facet num measurment tends to produce more, so aim to produce less
        %to neutralize bias
        if fitResults(itFit) < 1 & fitResults(itFit) > 0.95 | forceAccept  
            accetableFit = 1;
            fullBeeData(i).fitUsed = surfaceFitRng;
            fullBeeData(i).facetError = fitResults(itFit);
            fullBeeData(i).htStepUsed = htStep;
        else
            %otherwise choose next parameters to test
            if itFit == 1
                if fitResults(itFit) > 1
                    %should reduce ratio
                    surfaceFitRng = surfaceFitRng+50;
                elseif fitResults(itFit) < 0.95
                    %should increase ratio
                    surfaceFitRng = surfaceFitRng-50;
                end
            elseif itFit == 2
                %interpolate best choice from two closest tests
                ptsLess = find(1 - fitResults > 0);
                ptsGreater = find(1 - fitResults < 0);
                if ~isempty(ptsLess) & ~isempty(ptsGreater)
                   [~, ptLow] = max(fitResults(ptsLess)); ptLow = ptsLess(ptLow);
                   [~, ptHigh] = min(fitResults(ptsGreater)); ptHigh = ptsGreater(ptHigh);
                else if isempty(ptsLess)
                        %in this case take closest two above
                        [~, tempI] = sort(fitResults(ptsGreater));
                        ptLow = ptsGreater(tempI(1));
                        ptHigh = ptsGreater(tempI(2));
                    else if isempty(ptsGreater)
                            %in this case take closest two below
                            [~, tempI] = sort(fitResults(ptsLess));
                            ptLow = ptsLess(tempI(end-1));
                            ptHigh = ptsLess(tempI(end));
                end; end; end

                errorGrad = (fitResults(ptHigh) - fitResults(ptLow))/(surfacesTested(ptHigh) - surfacesTested(ptLow));
                if ~isempty(ptsLess) & ~isempty(ptsGreater) | isempty(ptsLess) 
                    surfaceFitRng = (1-fitResults(ptLow))/errorGrad+surfacesTested(ptLow);
                else
                    surfaceFitRng = (1-fitResults(ptHigh))/errorGrad+surfacesTested(ptHigh);
                end
            elseif itFit == 3
                %just find best so far
                [~, minI] = min(fitResults);
                forceAccept = 1;
                surfaceFitRng = surfacesTested(minI);
                warning(sprintf('!!! fit not achieved for bee %i - %s - will fit %.2f to rng %.1f', i, fullBeeData(i).BeeID, fitResults(minI), surfaceFitRng));
            end

            if surfaceFitRng < 15 | surfaceFitRng > 500
                error('probably error in fit prediction');
            end
            itFit = itFit + 1;
            if itFit > 10
                fitResults
                surfacesTested
                error('probably a problem with finding correct surface fit range');
            end
        end
    end
end

if ~isempty(saveDataAs)
    
    %probably not needed
    if ~saveFigs    
        for i = 1:numBees
            fullBeeData(i).OrigHexF = NaN; fullBeeData(i).headEyeF = NaN; fullBeeData(i).SamplingF = NaN;      
            fullBeeData(i).histFig = NaN; fullBeeData(i).hexInterpF = NaN; fullBeeData(i).AreaTestF = NaN;    
            fullBeeData(i).testAngF = NaN; fullBeeData(i).testOnEyeF = NaN; fullBeeData(i).testOnEyeFExtra = NaN;     
            fullBeeData(i).testFOV = NaN; fullBeeData(i).testOnWorldF = NaN; fullBeeData(i).testOnWorldFExtra = NaN;
            fullBeeData(i).LSetF = NaN;
        end
    end
    fprintf('saved data for %i bees', numBees);
    
    cd(saveDataFolder);
    save(saveDataAs, 'fullBeeData', 'refSphere', 'measurmentIntervals', 'useDirectHexAxes','forceVertical', 'forceHorizontal',...
        'powerP','powerPonSphere','IDWMinRangeFacet','useReverseNormal','CCT_Thresh','lensT_Thresh','limRetinaDist_Thresh',...
        'surfaceFitRngTests','baseHeightStep','stepInFromBorders','stepInInterval','limLargeDensities','densityLim','nearestOnWorld','bordersFromNearest');
end
