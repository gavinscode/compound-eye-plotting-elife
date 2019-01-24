function [ norm, radiusX, radiusY, gaussRadius, rng] = normalToSurface_withCurvature( pointCoords, surface, referenceCenter, flipAway, rng, levelSets, edgeSet, heightReq)
   %Written by Gavin Taylor, 2017. MIT License

    %%% this function needs a massive cleanup!
    
    %%%note that none of the curvature calculations appear to be correct

    %this technique is quite similar to the cross method combined with the coordinate transform methods in
    %CURVATURE MEASUREMENT OF 3D OBJECTS: EVALUATION AND COMPARISON OF THREE METHODS
    %and: SURFACE PARAMETERIZATION AND CURVATURE MEASUREMENT OF ARIBTRARY 3-D OBJECTS: FIVE PRACTICAL METHODS
    %except I used a full surface
    %unfortunately this is considered limited accuracy and slow :(
    
    %maybe a good modfication would be to find the center of the levels
    %after rotation, as done for 2d. discritization may have less of error
    %with relation to flat places, although their total error would be minimized.
    
    %mostly pulled from AngleOfRefractedRay
    %%% ray-tracing scripts need updating based on this
    
    debugPlot = 0;
    
    if isempty(rng)
        rng = 50;    %this may need tunning...
    end
    
    if debugPlot
        figure;
    end
    
    %ok, so this updated to correct for flat areas
    %flat areas often cause solver problems, because we know that it is a
    %concave structure, and solver should produce positive coeefs for x^2,
    %y^2 and xy. But it is expensive to run solver to check these are
    %produced, so this is a heuristic to approxiamte it.
    
    %this currently checks for the mean height of the border in some number
    %of sectors and ensure they are all have some variation
    
    %current solution is just to expand area radially, but could expand in
    %direction that is flat as well. Not sure if this would bias solver?
    
    %it also prevents areas that are flat from discritization in one of the
    %orignal axis from being too present on the border
    
    %might be godo to have these as external funcitons...
    pixSize = 5;
    
    %1st val is for normal areas, second val is for borders
    %these are painful tuning variables :/
    %now allocated dynamically in main loop
        %5/sqrt(2) also worked well here... and 5/2 for 5
    if ~isempty(heightReq)
        pixStepReq = heightReq; %height step must be greater than this, 
    else
        pixStepReq = [5/sqrt(2) 5*sqrt(2)];
    end
    numSectors = 8; %division of circle for height check
    testheights = 1;
    
    percentFlatSetOnBorder = 1/8; %checks that amount of border inds from flat area is less than this
    %and also does this for percentage of total num points, either one above percent fails
    
    origRng = rng;
    referenceCenterOrig = referenceCenter;

    satisfactoryBorders = 0;
    while ~satisfactoryBorders 
        surfInds = find(sqrt((surface(:,1)-pointCoords(1)).^2 + (surface(:,2)-pointCoords(2)).^2 + (surface(:,3)-pointCoords(3)).^2) < rng);
        [~, orignalPointInd] = min(sqrt((surface(surfInds,1)-pointCoords(1)).^2 + (surface(surfInds,2)-pointCoords(2)).^2 + (surface(surfInds,3)-pointCoords(3)).^2));
        
        if ~isempty(levelSets)
            xSets = levelSets(surfInds,1);
            xSets(xSets == 0) = [];
            xSets = unique(xSets);

            ySets = levelSets(surfInds,2);
            ySets(ySets == 0) = [];
            ySets = unique(ySets);
            
            zSets = levelSets(surfInds,3);
            zSets(zSets == 0) = [];
            zSets = unique(zSets);   
        end
        
        if ~isempty(edgeSet)
            ptsOnEdge = intersect(surfInds, edgeSet);
        else
            ptsOnEdge = []; 
        end
      
        if debugPlot
            subplot(1,3,1); hold on; axis equal
            plot3(surface(surfInds,1), surface(surfInds,2), surface(surfInds,3), '.');
           
            if ~isempty(levelSets)
                for i = 1:length(xSets)
                    inds = find(levelSets(surfInds,1) == xSets(i));
                    plot3(surface(surfInds(inds),1), surface(surfInds(inds),2), surface(surfInds(inds),3), 'rx');
                end
                
                for i = 1:length(ySets)
                    inds = find(levelSets(surfInds,2) == ySets(i));
                    plot3(surface(surfInds(inds),1), surface(surfInds(inds),2), surface(surfInds(inds),3), 'rx');
                end
                
                for i = 1:length(zSets)
                    inds = find(levelSets(surfInds,3) == zSets(i));
                    plot3(surface(surfInds(inds),1), surface(surfInds(inds),2), surface(surfInds(inds),3), 'rx');
                end
            end
            
            if ~isempty(ptsOnEdge)
                plot3(surface(ptsOnEdge,1), surface(ptsOnEdge,2), surface(ptsOnEdge,3), 'mx');
            end
        end

        %rotate surface to axis of max variation
%         surfPts = [surface(surfInds,1)-pointCoords(1),surface(surfInds,2)-pointCoords(2),surface(surfInds,3)-pointCoords(3)];
        centM = mean(surface(surfInds,:));
        surfPts = surface(surfInds,:)-centM;
        pcaVecs = pca(surfPts);
    
        rer = vectorRotationFromAtoB([0 0 1]',pcaVecs(:,3)');
        
        %test if ref center below, if so rotate up
        referenceCenter = referenceCenterOrig-centM; %[referenceCenterOrig(1) - pointCoords(1), referenceCenterOrig(2) - pointCoords(2), referenceCenterOrig(3) - pointCoords(3)];
        referenceCenter = referenceCenter*rer;
        if referenceCenter(3) < 0
            rer = vectorRotationFromAtoB([0 0 -1]',pcaVecs(:,3)');
            rer2 = vectorRotationFromAtoB(pcaVecs(:,3)',[0 0 -1]');
            referenceCenter = referenceCenterOrig-centM;  %[referenceCenterOrig(1) - pointCoords(1), referenceCenterOrig(2) - pointCoords(2), referenceCenterOrig(3) - pointCoords(3)];
            referenceCenter = referenceCenter*rer;
        else
            rer2 = vectorRotationFromAtoB(pcaVecs(:,3)',[0 0 1]');  
        end
        
        newPts = surfPts*rer;    
        
        bInds = find(sqrt((newPts(:,1)-newPts(orignalPointInd,1)).^2 + (newPts(:,2)-newPts(orignalPointInd,2)).^2 + (newPts(:,3)-newPts(orignalPointInd,3)).^2) > rng-pixSize*sqrt(2));
          
        %% test the flat sets are not on border
        if ~isempty(levelSets)
            setsFlag = 1;
            
            toRemove = zeros(length(surfInds),1);
            for i = 1:length(xSets)
                inds = find(levelSets(surfInds,1) == xSets(i));
                if ~isempty(intersect(bInds, inds))
                    toRemove(inds) = 1;
                end
            end
            for i = 1:length(ySets)
                inds = find(levelSets(surfInds,2) == ySets(i));
                if ~isempty(intersect(bInds, inds))
                    toRemove(inds) = 1;
                end
            end
            for i = 1:length(zSets)
                inds = find(levelSets(surfInds,3) == zSets(i));
                if ~isempty(intersect(bInds, inds))
                    toRemove(inds) = 1;
                end
            end
            
            numFlat = sum(toRemove);
            [~,toRemove,~] = intersect(bInds,find(toRemove));
            if length(toRemove)/length(bInds) > percentFlatSetOnBorder | ...
                    numFlat/length(surfInds) > percentFlatSetOnBorder
                setsFlag = 0;
            end
        else
            setsFlag = 1;
        end
        
        %% test that there is sufficent height variation along border
        if testheights
            ptHts = newPts(bInds,3)-min(newPts(:,3));
            ptAngs = zeros(length(bInds),1);
            for i = 1:length(bInds)
                ptAngs(i) = atan2(newPts(bInds(i),2),newPts(bInds(i),1));
            end

            %if looking for flat bits - remove from borders...
            %not sure this does much, as flat bit will block anyway
            if ~isempty(levelSets)
                if ~isempty(toRemove)
                    ptHts(toRemove) = [];
                    ptAngs(toRemove) = [];
                end
            end
    %         figure; 
    %         subplot(1,2,1); hold on
    %         plot3(newPts(bInds,1), newPts(bInds,2), newPts(bInds,3), 'cx');
    %         plot3(newPts(:,1), newPts(:,2), newPts(:,3), 'r.','markersize',5);
    %         for i = 1:length(bInds)
    %             text(newPts(bInds(i),1), newPts(bInds(i),2), newPts(bInds(i),3), sprintf('%.1f', ptHts(i)));
    %         end
    %         subplot(1,2,2); plot(ptAngs, ptHts, 'r.');

            %%%method based on sectors
            %-idea is to ensure some texture on all sides
            testRanges = -pi:2*pi/numSectors:pi;
            sectorMeans = zeros(length(testRanges)-1,1)*NaN;
            for i = 1:length(testRanges)-1
                inds = find(ptAngs >= testRanges(i) & ptAngs < testRanges(i+1));
                if ~isempty(inds)
                    sectorMeans(i) = mean(ptHts(inds));
                end
            end

            %edge requires large pix height diff than remainder
            if ~isempty(ptsOnEdge) 
                acceptableSector = sectorMeans > pixStepReq(2);
            else
                acceptableSector = sectorMeans > pixStepReq(1);
            end
            %if sector mean is nan, no values in sector because edge extendsbeyond, so its ok
            acceptableSector(isnan(sectorMeans)) = 1;

            %just check all sectors ok
            if sum(acceptableSector) ~= length(acceptableSector)
                sectorsFlag = 0;
            else
               sectorsFlag = 1; 
            end
            %step around and check that no two adjacent sectors have max below pixel size
%             for i = 1:length(acceptableSector)
%                 if i == length(acceptableSector)
%                     if ~acceptableSector(i) & ~acceptableSector(1)
%                         sectorsFlag = 0;
%                     end
%                 else
%                     if ~acceptableSector(i) & ~acceptableSector(i+1)
%                         sectorsFlag = 0;
%                     end
%                 end
%             end

            if debugPlot
                plot3(surface(surfInds(bInds),1), surface(surfInds(bInds),2), surface(surfInds(bInds),3), 'cx');
            end
        else
            sectorsFlag = 1;
        end
        if  setsFlag & sectorsFlag
            satisfactoryBorders = 1;
%             rng
        end
            
        %shift base to zero as cf tool can't translate horizontally
        %this might not influence the results much except for edges
        [~, minimaInd] = min(newPts(:,3));
        originOffset = newPts(minimaInd,:);
        newPts = newPts - originOffset; 
        
        %approach to find enough unique values as indicate of depth
        if ~satisfactoryBorders
            rng = rng + 10;
        end
        
        if rng > 500 %10*origRng 
           rng
           error('range has increased a lot, probably should check'); 
        end
    end

    if debugPlot
        title(sprintf('%.1f', rng));
         subplot(1,3,2); hold on; axis equal
         plot3(newPts(:,1), newPts(:,2), newPts(:,3), 'r.','markersize',20);
         line([0 referenceCenter(1)], [0 referenceCenter(2)], [0 referenceCenter(3)], 'color', 'r');
%          plot3(newPts(bInds,1), newPts(bInds,2), newPts(bInds,3), 'cx');
%          plot3(newPts(cInds,1), newPts(cInds,2), newPts(cInds,3), 'mx');
    end
    
    %fit the plane
    %force to touch origin
%    
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    %for fixed at zero - can be a large offset
    %not ideal, as point of interest may lie above zero in z
%   opts.Lower = [0 -Inf -Inf -Inf -Inf -Inf]; %
%   opts.Upper = [0 Inf Inf Inf Inf Inf];

%   for even
%   opts.Lower = [-Inf 0 0 -Inf 0 -Inf]; %
%   opts.Upper = [Inf 0 0 Inf 0 Inf];

%%%currently used!
%   for force concave around lens - x*y can be weird if limited area
%   but then need to make sure orientation is correct
    opts.Lower = [-Inf -Inf -Inf 0 0 0]; %
    opts.Upper = [Inf Inf Inf Inf Inf Inf];

    fittedPoly = fit(newPts(:,1:2), newPts(:,3),'poly22', opts);
     
    %fittedPoly = csaps(newPts(:,1:2), newPts(:,3));
    if debugPlot
        zlim([-100 100]); xlim([-500 500]); ylim([-500 500]);
        
        plot(fittedPoly);
    end
    
    %find normal at point that original was shifted to
    closeP = newPts(orignalPointInd,:);
    
    [dX, dY, dXX, dXY, dYY] = differentiate(fittedPoly, closeP(1), closeP(2));% first and second differential

    if isnan(dX)
        dX = 0;
    end
    if isnan(dY)
        dY = 0;
    end
    if isnan(dXX)
       dXX = 0;
    end
    if isnan(dXY)
       dXY = 0;
   end 
    if isnan(dYY)
       dYY = 0;
    end 
    
    %%%while in principle cuvature can just be calcualted on a plane bisecting a surface,
    %I don't think the equation for plane curves are correct
%     Kx= abs(dXX)/((1+dX^2)^1.5);
%     Ky= abs(dYY)/((1+dY^2)^1.5);

    % http://www.mathpages.com/rr/s5-03/5-03.htm
    %%%relativity page puts the tangent parallel to xy plane, so no linear terms contribut to curvature  apparently terms higher than 2nd order also do not contirbute to curvature at origin
    %but this is only correct if tangent is parallel to xy plane, which ours is not
%     Kx = fittedPoly.p20;
%     Ky =  fittedPoly.p02;
%     radiusX=(1/Kx);
%     radiusY=(1/Ky);
%     if isinf(radiusX)
%         radiusX = 10^12;
%     end
%     if isinf(radiusY)
%         radiusY = 10^12;
%     end
%     
    %%%I think this relates to the curvature of a space curve, 
    %at least where f(x,y)
    curvature = abs(dX*dYY-dY*dXX)/((dX^2+dY^2).^(3/2));
    %my derivation matches at least
%     curvature = abs(2*fittedPoly.p10*fittedPoly.p02-2*fittedPoly.p01*fittedPoly.p20)/((fittedPoly.p01.^2+fittedPoly.p10^2)^1.5)
    
    %extrinsic section suggest this equation when surface is rotated to align tangent to xy plane
%     curvature = 4*fittedPoly.p02*fittedPoly.p20-fittedPoly.p11^2

    %from intrinsic section this is the equaiton for a surface
    %with a miscligned tangent
%     curvature = (4*fittedPoly.p02*fittedPoly.p20-fittedPoly.p11^2)/(1+fittedPoly.p01^2+fittedPoly.p10^2);
%     curvature = (dXX*dYY-dXY^2)/(1+dX^2+dY^2);

    %not sure what curvatures are predicted here?
    radius = 1/curvature;
    radiusX = [];
    radiusY = [];
    gaussRadius = [];
    %plot circle
%     if debugPlot
%         th = 0:pi/50:2*pi;
%         xx = radiusX * cos(th);
%         xz = radiusX * sin(th)+radiusX;   
% 
%         xy = radiusY * cos(th);
%         yy = radiusY * sin(th);
% 
%         toPlotInd = find(xz< 20);
%          subplot(1,3,2)
%          plot3(xx(toPlotInd), zeros(length(toPlotInd),1),xz(toPlotInd));
%          
%         yy = radiusY * cos(th);
%         yz = radiusY * sin(th)+radiusY;
%     
%         toPlotInd = find(yz< 20);
%         plot3(zeros(length(toPlotInd),1), yy(toPlotInd), yz(toPlotInd))
%         xlabel('x');ylabel('y');zlabel('z');
%     end

    N = [dX, dY, -1]; N = N/sqrt(N(1)^2+N(2)^2+N(3)^2);
    referenceCenter = referenceCenter/sqrt(referenceCenter(1)^2+referenceCenter(2)^2+referenceCenter(3)^2);
    %flip normal to point towards lens
    if ~flipAway
        if sqrt((referenceCenter(1)-N(1))^2+(referenceCenter(2)-N(2))^2+(referenceCenter(3)-N(3))^2) > sqrt(2)
            N = N*-1;
        end
        
    else
        if sqrt((referenceCenter(1)-N(1))^2+(referenceCenter(2)-N(2))^2+(referenceCenter(3)-N(3))^2) < sqrt(2)
            N = N*-1;
        end
    end
    
    if debugPlot
        line([closeP(1) closeP(1)+N(1)*100],[closeP(2) closeP(2)+N(2)*100] , [0 N(3)*100]);
    end
    
    norm = N*rer2;
    
    %were previously inputs, now removed
%     if ~isempty(angle1) & ~isempty(angle2)
%         %get az and el of relocated normal
%         el = acos(norm(2));
%         az = atan2(norm(3),norm(1));
%         
%         %get lines for angles in flat space
%         %noet that these represent y,z in orignal coord space
%         line1 = [cos(angle1), sin(angle1)];
%         line2 = [cos(angle2), sin(angle2)];
%         
%         %rotate to align with normal vector
%         line1 = [0, line1]*vrrotvec2mat([0 0 1 -(-el+pi/2)])*vrrotvec2mat([0 1 0 az]);
%         line2 = [0, line2]*vrrotvec2mat([0 0 1 -(-el+pi/2)])*vrrotvec2mat([0 1 0 az]);
%         
%         %now rotate to coord frame of pca surf
%         line1 = line1*rer;
%         line2 = line2*rer;
%         
%         %http://www.mathpages.com/rr/s5-03/5-03.htm
%         %now to get ratios between line components (slope)
%         %this is an extrinsic approach
%          
%         ratio1 = line1(2)/line1(1);
%         ratio2 = line2(2)/line2(1);
%         
%         %maybe have to choose this such that it doesn't grow excessively
%         %large? i.e. swap ratio if one tends to infinity...
%         
%         
%         %this is only correct if surface tangent aligned to xy plane
%         %in principle rotations or intrinsic approach will solve fully
%         kurv1 = 2*(fittedPoly.p20+fittedPoly.p11*ratio1+fittedPoly.p02*ratio1^2)/(1+ratio1^2);
%         kurv2 = 2*(fittedPoly.p20+fittedPoly.p11*ratio2+fittedPoly.p02*ratio2^2)/(1+ratio2^2);
%         radius1=(1/kurv1);
%         radius2=(1/kurv2);
%     end
           
    if debugPlot
         subplot(1,3,1);
         line([pointCoords(1) pointCoords(1)+norm(1)*50],[pointCoords(2) pointCoords(2)+norm(2)*50] , [pointCoords(3) pointCoords(3)+norm(3)*50]);
    end
end

% normalToSurface(latOcelli.frontRetinaSurfaceSubs(50,:), latOcelli.frontRetinaSurfaceSubs);
