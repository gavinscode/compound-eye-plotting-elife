%Written by Gavin Taylor, 2017. MIT License

%%%This file calcualtes and saves plots to describe compound eye vision
%%%calcuated by volumetric analysis, as descrbied in Taylor et al. 2019

%%%% Settings %%%%
%set parameters here before running the script

clear; clc; close all
cd('Add here'); %folder with data file - modfy as needed
sDat = load('DataForPlots.mat'); %file to plot from

toPlotMaps = [1 3]; %indices of specific bees to plot maps for
saveAllThePlots = 0; %set to save all plots to the folder below
saveFolder = 'C:\SaveDirectory'; %folder to save to

numHistBins = 20; %number of bins in histograms (will be num+1)
scaleEyeVols = 1; %will convert eye volumes to um and take cube root

imageConvRes = 1; %resolution of plotted maps in degrees
angStep = 10; %resolution of plotted profiles in degrees
remIrregularSpatialBands = 1; 

regressionGroups = {{'AM'},{'BT'}}; %set species groups for performing regression
groupNames = {'Apis','Bombus'}; %set names to plot for species groups
groupLineStyles = {'-','-.'}; %set line type for plotting species groups
groupMarkerStyles = {'o','s'}; %set marker type for plotting species groups
volCols = cool(100); %set color for plotting volumes

plotBinoBorderOnMap = 1; %plot border of binocular field of color maps
plotFOVBorderOnMap = 1; %plot border of FOV on color maps

%parameters for plotting colour maps
t100 = 0:1/101:1;

minRangeIO = 1/180*pi; %set min IO on scale
maxRangeIO = 3/180*pi; %set max IO on scale
IOColMap = flipud([t100'*255 t100'*51 t100'*51]/255); %color map to use

minRangeD = 16; %set min diameter on scale
maxRangeD = 28; %set max diameter on scale
DColMap = [t100'*51 t100'*153 t100'*255]/255;

minRangeRL = 50; %set min retina thickness on scale
maxRangeRL = 350; %set max retina thickness on scale
RLColMap = [t100'*153 t100'*255 t100'*51]/255;

minRangeEyePara = 0.3; %set min P on scale
maxRangeEyePara = 1.5; %set max P on scale
EyeParaColMap = [t100'*255 t100'*153 t100'*51]/255;

minRangeLensT = 20; %set min Lens thickness on scale
maxRangeLensT = 70; %set max lens thickness on scale
LensTColMap = [t100'*51 t100'*255 t100'*255]/255;

minRangeCCT = 25; %set min CC thickness on scale
maxRangeCCT = 75; %set max CC thickness on scale
CCTColMap = [t100'*255 t100'*151 t100'*53]/255;

minRangeCurv = 200; %set min curvature on scale
maxRangeCurv = 1400; %set max curvature on scale
CurvColMap = [t100'*192 t100'*192 t100'*192]/255;

maxRangeScale = 3; %set min scaling coeff on scale
minRangeScale = -1; %set max scaling coeff on scale
ScaleColMap = jet(100);
ScaleColMap(1:25,:) = [];

%%% Note: Other colors and plotting parameters are set in main code

%% end settings, start plotting code
numBees = size(sDat.fullBeeData,2);

if saveAllThePlots
    cd(saveFolder);
end 

eyeVols = [sDat.fullBeeData.EyeVolume];
if scaleEyeVols
    eyeVols = (eyeVols/10^9).^(1/3);
end
volRng = [min(eyeVols) max(eyeVols)];

regGroup = zeros(numBees,1);

%save min and max for histograms
IORange = zeros(numBees,2);
DRange = zeros(numBees,2);
LensLRange = zeros(numBees,2);
CCLRange = zeros(numBees,2);
RetLRange = zeros(numBees,2);
SensRange = zeros(numBees,2);
DensityRange = zeros(numBees,2);
ProjARange = zeros(numBees,2);

%get ranges and single values first
for i = 1:numBees
    IORange(i,:) = [min(sDat.fullBeeData(i).calInterO) max(sDat.fullBeeData(i).calInterO)];
    DRange(i,:) = [min(sDat.fullBeeData(i).avgDiameter) max(sDat.fullBeeData(i).avgDiameter)];
    LensLRange(i,:) = [min(sDat.fullBeeData(i).lensThickness) max(sDat.fullBeeData(i).lensThickness)];
    CCLRange(i,:) = [min(sDat.fullBeeData(i).CCThickness) max(sDat.fullBeeData(i).CCThickness)];
    RetLRange(i,:) = [min(sDat.fullBeeData(i).RetThickness) max(sDat.fullBeeData(i).RetThickness)];
    SensRange(i,:) = [min(sDat.fullBeeData(i).sensitvityApproximation) max(sDat.fullBeeData(i).sensitvityApproximation)];
    DensityRange(i,:) = [min(sDat.fullBeeData(i).AxisDensityOnSphere) max(sDat.fullBeeData(i).AxisDensityOnSphere)];

    %get what group this bee belongs to
    for j = 1:length(regressionGroups)
        tempGroup = regressionGroups{j};
        for k = 1:length(tempGroup)
            if strcmp(sDat.fullBeeData(i).BeeSpecies, tempGroup{k})
                regGroup(i) = j;
            end
        end
    end
end

%get approx area per point on sphere
stPPoint = 4*pi/size(sDat.refSphere.X,1);

%  resample into facet reference frame - req for histograms and means
facetWeightInds = cell(numBees,1);
projWeightInds = cell(numBees,1);
for i = 1:numBees
    facetWeightInds{i} = sDat.fullBeeData(i).NumFacets./sDat.fullBeeData(i).interpFacetArea./sum(1./sDat.fullBeeData(i).interpFacetArea); 
    projWeightInds{i} = sDat.fullBeeData(i).AxisDensityOnSphere*stPPoint;
end

%% summary plots
% plot of itw vs. vol and surface area
volAllomF = figure; 
yyaxis left; hold on
areaCol = [0 0 0]/255;
legPoints = [0,0];
for i = 1:length(regressionGroups)
    tempX = [sDat.fullBeeData(regGroup == i).ITW];
    tempY = ([sDat.fullBeeData(regGroup == i).lensSurfaceArea]/10^6).^(1/2);
    legPoints(i) = plot(tempX, tempY, groupMarkerStyles{i}, 'color', areaCol,'markersize',10,'markerfacecolor','w');
end
ylim([0 3]);
ylabel('Eye Area^0^.^5 (mm)'); 
set(gca, 'YTick', [0 1 2 3]);
set(gca,'ycolor',areaCol)
%do fit on last one, it is bombus group
%log transform fits log(y) = log(a)+b*log(y)
fitArea = polyfit(log10(tempX),log10(tempY),1);
fitArea(2) = 10.^fitArea(2);
plotX = min(tempX):0.001:max(tempX);
legLine1 = plot(plotX, fitArea(2)*plotX.^fitArea(1),'color',areaCol,'linewidth',2);
[R,P] = corrcoef(log10(tempX),log10(tempY));
text(3, 1,sprintf('Y=%.2fX^%.2f R^2=%.2f, p=%.3f',fitArea(2),fitArea(1), R(2)^2, P(2)),'color', areaCol, 'Interpreter', 'none')
uistack(legLine1,'down');

yyaxis right; hold on
volCol = [51 153 255]/255;
for i = 1:length(regressionGroups)
    tempX = [sDat.fullBeeData(regGroup == i).ITW];
    tempY = ([sDat.fullBeeData(regGroup == i).EyeVolume]/10^9).^(1/3); % retinaVolume
    hl = plot(tempX,tempY, groupMarkerStyles{i}, 'color', volCol,'markerfacecolor','w','markersize',10);
end
ylabel('Eye Volume^0^.^3 (mm)'); xlabel('ITW (mm)');

fitVol = polyfit(log10(tempX),log10(tempY),1);
fitVol(2) = 10.^fitVol(2);
plotX = min(tempX):0.001:max(tempX);
legLine2 = plot(plotX, fitVol(2)*plotX.^fitVol(1),'color',volCol, 'linewidth',2);
[R,P] = corrcoef(log10(tempX),log10(tempY));
text(1.5, 0.8,sprintf('Y=%.2fX^%.2f R^2=%.2f, p=%.3f',fitVol(2),fitVol(1), R(2)^2, P(2)),'color', volCol, 'Interpreter', 'none')
%uistack does not work on second y axis...
delete(hl)
plot(tempX,tempY, groupMarkerStyles{i}, 'color', volCol,'markerfacecolor','w','markersize',10)

legend([legPoints], {groupNames{1}, groupNames{2}},'Location','SouthWest');
ylabel('Volume^0^.^3 (mm)');
xlim([1 7]); ylim([0 1]);
set(gca, 'YTick', [0 0.25 0.5 0.75 1]);
set(gca,'XTick',[1 3 5 7],'TickDir','out');
set(gca,'ycolor',volCol)
set(gcf, 'Position', [50 50 500*1.04 500/1.61]);

% plot of vol vs facets
facetAllomF = figure; hold on

facetCol = [0 0 0]/255;
for i = 1:length(regressionGroups)
    tempX = ([sDat.fullBeeData(regGroup == i).EyeVolume]/10^9).^(1/3);
    tempY = ([sDat.fullBeeData(regGroup == i).NumFacets]);
    plot(tempX,tempY, groupMarkerStyles{i}, 'color', facetCol,'markerfacecolor','w','markersize',10);
end
xlabel('Eye Volume^0^.^3 (mm)'); ylabel('# Facets');

fitVol = polyfit(log10(tempX),log10(tempY),1);
fitVol(2) = 10.^fitVol(2);
plotX = min(tempX):0.001:max(tempX);
legLine = plot(plotX, fitVol(2)*plotX.^fitVol(1),'color',facetCol, 'linewidth',2);
[R,P] = corrcoef(log10(tempX),log10(tempY));
text(0.6,3000,sprintf('Y=%.2fX^%.2f R^2=%.2f, p=%.3f',fitVol(2),fitVol(1), R(2)^2, P(2)),'color', facetCol, 'Interpreter', 'none')
uistack(legLine,'down');

for i = 1:length(regressionGroups)
    tempX = ([sDat.fullBeeData(regGroup == i).EyeVolume]/10^9).^(1/3);
    inds = find(regGroup == i);
    tempY = zeros(length(inds),1)';
    for j = 1:length(inds)
        [~,binoRefInds] = intersect(sDat.fullBeeData(inds(j)).inFOV, sDat.fullBeeData(inds(j)).inFOVBino);
        tempY(j) = sum(sDat.fullBeeData(inds(j)).AxisDensityOnSphere(binoRefInds)*stPPoint);
    end
    plot(tempX,tempY, groupMarkerStyles{i}, 'color', facetCol,'markerfacecolor','w','markersize',10);
end

% legend(groupNames);
xlim([0.5 0.9]); ylim([0 8000]);
set(gca, 'YTick', [0 2000 4000 6000 8000]);
set(gca,'XTick',[0.5 0.6 0.7 0.8 0.9],'TickDir','out');
set(gcf, 'Position', [50 50 500 500/1.61]);

% plot sensitivty
sensAllomF = figure; hold on
for i = 1:length(regressionGroups)
    tempX = ([sDat.fullBeeData(regGroup == i).EyeVolume]/10^9).^(1/3);
    inds = find(regGroup == i);
    tempY = zeros(length(inds),1)';
    for j = 1:length(inds)
        tempY(j) = wmean(sDat.fullBeeData(inds(j)).sensitvityApproximation,facetWeightInds{inds(j)});
    end
    legs3 = plot(tempX,tempY, groupMarkerStyles{i}, 'color', facetCol,'markerfacecolor','w','markersize',10);
end
xlabel('Eye Volume^0^.^3 (mm)'); ylabel('Sensitvity (um.str)');

% legend(groupNames);
xlim([0.5 0.9]); ylim([0 0.2]);
set(gca, 'YTick', [0 0.05 0.1 0.15 0.2]);
set(gca,'XTick',[0.5 0.6 0.7 0.8 0.9],'TickDir','out');
set(gcf, 'Position', [50 50 500 500/1.61]);

% plot for area/fov ratio and eye paramter for comp
fovRatioAndEyeP = figure; 
yyaxis left; hold on
xlabel('Eye Volume^0^.^3 (mm)'); ylabel('Eye proportion (um/rad)');
ylim([0 1200]);
set(gca, 'YTick', [0 400 800 1200]);

facetCol = [0 0 0];
for i = 1:length(regressionGroups)
    tempX = ([sDat.fullBeeData(regGroup == i).EyeVolume]/10^9).^(1/3);
    inds = find(regGroup == i);
    tempY = zeros(length(inds),1)';
    for j = 1:length(inds)
        tempY(j) = 1./wmean((sDat.fullBeeData(inds(j)).projectedAngle).^0.5,facetWeightInds{inds(j)});
    end
    legs2 = plot(tempX,tempY, groupMarkerStyles{i}, 'color', facetCol,'markerfacecolor','w','markersize',10);
end

fitVol = polyfit(log10(tempX),log10(tempY),1);
fitVol(2) = 10.^fitVol(2);
plotX = min(tempX):0.001:max(tempX);
legLine = plot(plotX, fitVol(2)*plotX.^fitVol(1),'color',facetCol, 'linewidth',2);
[R,P] = corrcoef(log10(tempX),log10(tempY));
text(0.55,1100,sprintf('Y=%.2fX^%.2f R^2=%.2f, p=%.3f',fitVol(2),fitVol(1), R(2)^2, P(2)),'color', facetCol, 'Interpreter', 'none')
uistack(legLine,'down');
set(gca,'ycolor','k')

facetCol = [102 102 255]/255;
for i = 1:length(regressionGroups)
    tempX = ([sDat.fullBeeData(regGroup == i).EyeVolume]/10^9).^(1/3);
    tempY = (([sDat.fullBeeData(regGroup == i).lensSurfaceArea])./(([sDat.fullBeeData(regGroup == i).fovFilledFromIntegral]*4*pi))).^0.5;
    legs1 = plot(tempX,tempY, groupMarkerStyles{i}, 'color', facetCol,'markerfacecolor',facetCol,'markersize',8);
end

yyaxis right; hold on
facetCol = [255 153 51]/255;
for i = 1:length(regressionGroups)
    tempX = ([sDat.fullBeeData(regGroup == i).EyeVolume]/10^9).^(1/3);
    inds = find(regGroup == i);
    tempY = zeros(length(inds),1)';
    for j = 1:length(inds)
        tempY(j) = wmean(sDat.fullBeeData(inds(j)).calInterO.*sDat.fullBeeData(inds(j)).avgDiameter,facetWeightInds{inds(j)});
    end
    legs3 = plot(tempX,tempY, groupMarkerStyles{i}, 'color', facetCol,'markerfacecolor','w','markersize',10);
end
ylabel('Eye parameter (um.rad)');
ylim([0 1]);
set(gca, 'YTick', [0 0.25 0.5 0.75 1]);

fitVol = polyfit(log10(tempX),log10(tempY),1);
fitVol(2) = 10.^fitVol(2);
plotX = min(tempX):0.001:max(tempX);
legLine = plot(plotX, fitVol(2)*plotX.^fitVol(1),'color',facetCol, 'linewidth',2);
[R,P] = corrcoef(log10(tempX),log10(tempY));
text(0.65,0.5,sprintf('Y=%.2fX^%.2f R^2=%.2f, p=%.3f',fitVol(2),fitVol(1), R(2)^2, P(2)),'color', facetCol, 'Interpreter', 'none')
delete(hl)
plot(tempX,tempY, groupMarkerStyles{i}, 'color', facetCol,'markerfacecolor','w','markersize',10);

legend([legs2, legs1, legs3], {'Global eye ratio','Average eye ratio','Eye parameter'},'Location','SouthWest');
xlim([0.5 0.9]); 
set(gca,'XTick',[0.5 0.6 0.7 0.8 0.9],'TickDir','out');
set(gca,'ycolor',facetCol)
set(gcf, 'Position', [50 50 500 500/1.61]);

%plot for combined and binocular fov
fovAllomF = figure; 
yyaxis left; hold on
monFOVCol = [102 102 255]/255;
legPoints = [0,0];
%%%fov is in square radians, but should not need to normalize with sqrt?
for i = 1:length(regressionGroups)
    tempX = ([sDat.fullBeeData(regGroup == i).EyeVolume]/10^9).^(1/3);
    tempY = (([sDat.fullBeeData(regGroup == i).fovFilledFromIntegral])*4*pi);
    legPoints(i) = plot(tempX, tempY, groupMarkerStyles{i}, 'color', monFOVCol,'markerfacecolor','w','markersize',10);
end

fitArea = polyfit(log10(tempX),log10(tempY),1);
fitArea(2) = 10.^fitArea(2);
plotX = min(tempX):0.001:max(tempX);
legLines1 = plot(plotX, fitArea(2)*plotX.^fitArea(1),'-','color',monFOVCol,'linewidth',2);
[R,P] = corrcoef(log10(tempX),log10(tempY));
text(0.51,4,sprintf('Y=%.2fX^%.2f R^2=%.2f, p=%.3f',fitArea(2),fitArea(1), R(2)^2, P(2)),'color', monFOVCol, 'Interpreter', 'none')
uistack(legLines1,'down');
set(gca,'ycolor',[0 0 0]);

BinoFOVCol = [102 255 102]/255;
for i = 1:length(regressionGroups)
    tempX = ([sDat.fullBeeData(regGroup == i).EyeVolume]/10^9).^(1/3);
    inds = find(regGroup == i);
    tempY = zeros(length(inds),1)';
    for j = 1:length(inds)
        tempY(j) = (length(sDat.fullBeeData(inds(j)).inFOVBino)./size(sDat.refSphere.X,1)*4*pi);
    end
    plot(tempX, tempY, groupMarkerStyles{i}, 'color', BinoFOVCol,'markerfacecolor','w','markersize',10);
end
fitArea = polyfit(log10(tempX),log10(tempY),1);
fitArea(2) = 10.^fitArea(2);
plotX = min(tempX):0.001:max(tempX);
legLines2 = plot(plotX, fitArea(2)*plotX.^fitArea(1),'-','color',BinoFOVCol,'linewidth',2);
[R,P] = corrcoef(log10(tempX),log10(tempY));
text(0.51,1.5,sprintf('Y=%.2fX^%.2f R^2=%.2f, p=%.3f',fitArea(2),fitArea(1), R(2)^2, P(2)),'color', BinoFOVCol, 'Interpreter', 'none')
uistack(legLines2,'down');

CompFOVCol = [255 102 102]/255; 
for i = 1:length(regressionGroups)
    tempX = ([sDat.fullBeeData(regGroup == i).EyeVolume]/10^9).^(1/3);
     inds = find(regGroup == i);
    tempY = zeros(length(inds),1)';
    for j = 1:length(inds)
        tempY(j) = (length(unique([sDat.fullBeeData(inds(j)).inFOV' sDat.fullBeeData(inds(j)).inFOVRight']))./size(sDat.refSphere.X,1)*4*pi);
    end
    plot(tempX,tempY, groupMarkerStyles{i}, 'color', CompFOVCol,'markerfacecolor','w','markersize',10);
end
xlabel('Eye Volume^0^.^3 (mm)'); ylabel('Corneal FOV (str)');

fitVol = polyfit(log10(tempX),log10(tempY),1);
fitVol(2) = 10.^fitVol(2);
plotX = min(tempX):0.001:max(tempX);
legLines3 = plot(plotX, fitVol(2)*plotX.^fitVol(1),'-','color',CompFOVCol,'linewidth',2);
[R,P] = corrcoef(log10(tempX),log10(tempY));
text(0.51,6.5,sprintf('Y=%.2fX^%.2f R^2=%.2f, p=%.3f',fitVol(2),fitVol(1), R(2)^2, P(2)),'color', CompFOVCol, 'Interpreter', 'none')
uistack(legLines3,'down');

legend([legPoints, legLines3, legLines1, legLines2], {groupNames{1}, groupNames{2}, 'Complete', 'Monocular', 'Binocular'});
xlim([0.5 0.9]); ylim([0 7]);
set(gca, 'YTick', [0 1 2 3 4 5 6 7]);
set(gca,'XTick',[0.5 0.6 0.7 0.8 0.9 1],'TickDir','out');

yyaxis right; hold on;
ylim([0 7]/(4*pi)*100);
set(gca, 'YTick', [0 0.1 0.2 0.3 0.4 0.5]*100);
set(gca,'ycolor',[0 0 0]);
ylabel('% of total visual field');
set(gcf, 'Position', [50 50 500*1.02 500/1.61*2*1.06]);

%% plot indvidual of means and violin hists
%for IO
IOHistsF = figure; hold on

smthL = 3; %smooth histogram over this num of bins

HistBins = 0:0.1:5; 
flipChange = [0 0 1 0 0 0 0 0]; %bit of a fudge, indicate which to flip if needed
histCol = [0.5 0.5 0.5];

%histograms are scalled up to hist width, such that hist cal would equal hist width
histW = 0.018; 
histCal = 0.1;
labels = zeros(length(regressionGroups),1);

for i = 1:length(regressionGroups)
    tempX = ([sDat.fullBeeData(regGroup == i).EyeVolume]/10^9).^(1/3);
    inds = find(regGroup == i);
    tempY = zeros(length(inds),1);
    tempY2 = sqrt(23818./[sDat.fullBeeData(regGroup == i).NumFacets]);
    tempY3 = sqrt([sDat.fullBeeData(regGroup == i).fovFilledFromIntegral]*4*pi./[sDat.fullBeeData(regGroup == i).NumFacets])/pi*180;
    tempY4 = zeros(length(inds),1);
    
    for j = 1:length(inds)
        tempY(j) = wmean(sDat.fullBeeData(inds(j)).calInterO,facetWeightInds{inds(j)})/pi*180; 
        tempY4(j) = max(sDat.fullBeeData(inds(j)).calInterO)/pi*180;
        [n,x] = histwc(sDat.fullBeeData(inds(j)).calInterO/pi*180,facetWeightInds{inds(j)},HistBins);
        
        relN = smooth(n/sum(n),smthL);
        %adjust for diff between max n range and intended cal
        xAdj = (max(relN)/histCal);
        
        if ~flipChange(inds(j))
            dirToP = 'left';
            histWToUse = histW*xAdj;
            xToP = tempX(j)-histWToUse/2;
        else
            dirToP = 'right';
            histWToUse = histW*xAdj;
            xToP = tempX(j)+histWToUse/2;
        end
        %note, most options for calculaiton don't seem to work well with weighted histogram?
        distributionPlot({[x',relN]}, 'xValues', xToP,'distWidth', histWToUse, 'showMM',0,'histOri',dirToP,'globalNorm',3,'color',histCol)
    end
    labels(i) = plot(tempX, tempY, groupMarkerStyles{i}, 'color', 'k', 'markersize',10,'markerfacecolor','w');
    label2 = plot(tempX, tempY2, groupMarkerStyles{i}, 'color', [204 102 0]/255, 'markersize',5,'markerfacecolor',[204 102 0]/255);
    label3 = plot(tempX, tempY3, groupMarkerStyles{i}, 'color', [102 102 255]/255, 'markersize',5,'markerfacecolor',[102 102 255]/255);
end

fitVol = polyfit(log10(tempX),log10(tempY)',1);
fitVol(2) = 10.^fitVol(2);
plotX = min(tempX):0.001:max(tempX);
legLines3 = plot(plotX, fitVol(2)*plotX.^fitVol(1),'-','color','k','linewidth',2);
[R,P] = corrcoef(log10(tempX),log10(tempY));
text(0.51,3.5,sprintf('Y=%.2fX^%.2f R^2=%.2f, p=%.3f',fitVol(2),fitVol(1), R(2)^2, P(2)),'color', 'k', 'Interpreter', 'none')
uistack(legLines3,'down',3); %get behind other three markers

xlim([0.5 0.9]); ylim([0 4]);
set(gca, 'YTick', [0 1 2 3 4]);
set(gca,'XTick',[0.5 0.6 0.7 0.8 0.9],'TickDir','out');
xlabel('Eye Volume^0^.^3 (mm)'); ylabel('Corneal interommatidial angle (^o)');
legend([labels(1) labels(2) label2 label3], {groupNames{1}, groupNames{2},'HemiApprox','FOVApprox'});
% title('Magenta are from hemisphere div num eqn - blue are from fov div');
set(gcf, 'Position', [50 50 500 500/1.61]);

% plot summary for facet d
DiamHistsF = figure; hold on

HistBins = 0:0.25:30; 
labels = zeros(length(regressionGroups),1);

for i = 1:length(regressionGroups)
    tempX = ([sDat.fullBeeData(regGroup == i).EyeVolume]/10^9).^(1/3);
    inds = find(regGroup == i);
    tempY = zeros(length(inds),1);
    tempY2 = zeros(length(inds),1);
    for j = 1:length(inds)
        tempY(j) = wmean(sDat.fullBeeData(inds(j)).avgDiameter,facetWeightInds{inds(j)}); 
        tempY2(j) = max(sDat.fullBeeData(inds(j)).avgDiameter);
        [n,x] = histwc(sDat.fullBeeData(inds(j)).avgDiameter,facetWeightInds{inds(j)},HistBins);
        
        relN = smooth(n/sum(n),smthL);
        %adjust for diff between max n range and intended cal
        xAdj = (max(relN)/histCal);
        
        if ~flipChange(inds(j))
            dirToP = 'left';
            histWToUse = histW*xAdj;
            xToP = tempX(j)-histWToUse/2;
        else
            dirToP = 'right';
            histWToUse = histW*xAdj;
            xToP = tempX(j)+histWToUse/2;
        end
        distributionPlot({[x',relN]}, 'xValues', xToP,'distWidth', histWToUse, 'showMM',0,'histOri',dirToP,'globalNorm',3,'color',histCol)
    end
    labels(i) = plot(tempX, tempY, groupMarkerStyles{i}, 'color', 'k', 'markersize',10,'markerfacecolor','w');
end
fitVol = polyfit(log10(tempX),log10(tempY)',1);
fitVol(2) = 10.^fitVol(2);

fitVol2 = polyfit(log10(tempX),log10(tempY2)',1);
fitVol2(2) = 10.^fitVol2(2);
plotX = min(tempX):0.001:max(tempX);
legLines3 = plot(plotX, (fitVol(2)*plotX.^fitVol(1)),'-','color','k','linewidth',2);
[R,P] = corrcoef(log10(tempX),log10(tempY));
text(0.6,10,sprintf('Y=%.2fX^%.2f R^2=%.2f, p=%.3f',fitVol(2),fitVol(1), R(2)^2, P(2)),'color', 'k','Interpreter', 'none')
uistack(legLines3,'down');

xlim([0.5 0.9]); ylim([0 30]);
set(gca, 'YTick', [0 10 20 30]);
set(gca,'XTick',[0.5 0.6 0.7 0.8 0.9],'TickDir','out');
xlabel('Eye Volume^0^.^3 (mm)'); ylabel('Facet diameter (um)');
legend(labels, groupNames,'Location','South');
set(gcf, 'Position', [50 50 500 500/1.61 ]);

% plot summary for rhabdoms thickness
RhabdomsHistsF = figure; hold on

HistBins = 0:5:500; 
labels = zeros(length(regressionGroups),1);

for i = 1:length(regressionGroups)
    tempX = ([sDat.fullBeeData(regGroup == i).EyeVolume]/10^9).^(1/3);
    inds = find(regGroup == i);
    tempY = zeros(length(inds),1);
    for j = 1:length(inds)
        tempY(j) = wmean(sDat.fullBeeData(inds(j)).RetThickness,facetWeightInds{inds(j)}); 
        [n,x] = histwc(sDat.fullBeeData(inds(j)).RetThickness,facetWeightInds{inds(j)},HistBins);
        
        relN = smooth(n/sum(n),smthL);
        %adjust for diff between max n range and intended cal
        xAdj = (max(relN)/histCal);
        
        if ~flipChange(inds(j))
            dirToP = 'left';
            histWToUse = histW*xAdj;
            xToP = tempX(j)-histWToUse/2;
        else
            dirToP = 'right';
            histWToUse = histW*xAdj;
            xToP = tempX(j)+histWToUse/2;
        end
        distributionPlot({[x',relN]}, 'xValues', xToP,'distWidth', histWToUse, 'showMM',0,'histOri',dirToP,'globalNorm',3,'color',histCol)
    end
    labels(i) = plot(tempX, tempY, groupMarkerStyles{i}, 'color', 'k', 'markersize',10,'markerfacecolor','w');
end
fitVol = polyfit(log10(tempX),log10(tempY)',1);
fitVol(2) = 10.^fitVol(2);
plotX = min(tempX):0.001:max(tempX);
legLines3 = plot(plotX, (fitVol(2)*plotX.^fitVol(1)),'-','color','k','linewidth',2);
[R,P] = corrcoef(log10(tempX),log10(tempY));
text(0.6,350,sprintf('Y=%.2fX^%.2f R^2=%.2f, p=%.3f',fitVol(2),fitVol(1), R(2)^2, P(2)),'color', 'k','Interpreter', 'none')
uistack(legLines3,'down');

xlim([0.5 0.9]); ylim([0 500]);
set(gca, 'YTick', [0 100 200 300 400 500]);
set(gca,'XTick',[0.5 0.6 0.7 0.8 0.9],'TickDir','out');
xlabel('Eye Volume^0^.^3 (mm)'); ylabel('Retina Thickness (um)');
legend(labels, groupNames,'Location','NorthWest');
set(gcf, 'Position', [50 50 500 500/1.61 ]);

% plot summary for cornea thickness
corneaHistF = figure; hold on

HistBins = 0:1:100; 
labels = zeros(length(regressionGroups),1);

for i = 1:length(regressionGroups)
    tempX = ([sDat.fullBeeData(regGroup == i).EyeVolume]/10^9).^(1/3);
    inds = find(regGroup == i);
    tempY = zeros(length(inds),1);
    for j = 1:length(inds)
        tempY(j) = wmean(sDat.fullBeeData(inds(j)).lensThickness,facetWeightInds{inds(j)}); 
        [n,x] = histwc(sDat.fullBeeData(inds(j)).lensThickness,facetWeightInds{inds(j)},HistBins);
        
        relN = smooth(n/sum(n),smthL);
        %adjust for diff between max n range and intended cal
        xAdj = (max(relN)/histCal);
        
        if ~flipChange(inds(j))
            dirToP = 'left';
            histWToUse = histW*xAdj;
            xToP = tempX(j)-histWToUse/2;
        else
            dirToP = 'right';
            histWToUse = histW*xAdj;
            xToP = tempX(j)+histWToUse/2;
        end
        distributionPlot({[x',relN]}, 'xValues', xToP,'distWidth', histWToUse, 'showMM',0,'histOri',dirToP,'globalNorm',3,'color',histCol)

    end
    labels(i) = plot(tempX, tempY, groupMarkerStyles{i}, 'color', 'k', 'markersize',10,'markerfacecolor','w');
end

fitVol = polyfit(log10(tempX),log10(tempY)',1);
fitVol(2) = 10.^fitVol(2);
plotX = min(tempX):0.001:max(tempX);
legLines3 = plot(plotX, (fitVol(2)*plotX.^fitVol(1)),'-','color','k','linewidth',2);
[R,P] = corrcoef(log10(tempX),log10(tempY));
text(0.51,90,sprintf('Y=%.2fX^%.2f R^2=%.2f, p=%.3f',fitVol(2),fitVol(1), R(2)^2, P(2)),'color', 'k','Interpreter', 'none')
uistack(legLines3,'down');

xlim([0.5 0.9]); ylim([0 100]);
set(gca, 'YTick', [0 25 50 75 100]);
set(gca,'XTick',[0.5 0.6 0.7 0.8 0.9],'TickDir','out');
xlabel('Eye Volume^0^.^3 (mm)'); ylabel('Lens Thickness (um)');
legend(labels, groupNames);
set(gcf, 'Position', [50 50 500 500/1.61 ]);

% plot summary for CC thickness
CCHistF = figure; hold on

HistBins = 0:1:100; 
labels = zeros(length(regressionGroups),1);

for i = 1:length(regressionGroups)
    tempX = ([sDat.fullBeeData(regGroup == i).EyeVolume]/10^9).^(1/3);
    inds = find(regGroup == i);
    tempY = zeros(length(inds),1);
    for j = 1:length(inds)
        tempY(j) = wmean(sDat.fullBeeData(inds(j)).CCThickness,facetWeightInds{inds(j)}); 
        [n,x] = histwc(sDat.fullBeeData(inds(j)).CCThickness,facetWeightInds{inds(j)},HistBins);
        
        relN = smooth(n/sum(n),smthL);
        %adjust for diff between max n range and intended cal
        xAdj = (max(relN)/histCal);
        
        if ~flipChange(inds(j))
            dirToP = 'left';
            histWToUse = histW*xAdj;
            xToP = tempX(j)-histWToUse/2;
        else
            dirToP = 'right';
            histWToUse = histW*xAdj;
            xToP = tempX(j)+histWToUse/2;
        end
        distributionPlot({[x',relN]}, 'xValues', xToP,'distWidth', histWToUse, 'showMM',0,'histOri',dirToP,'globalNorm',3,'color',histCol)

    end
    labels(i) = plot(tempX, tempY, groupMarkerStyles{i}, 'color', 'k', 'markersize',10,'markerfacecolor','w');
end
fitVol = polyfit(log10(tempX),log10(tempY)',1);
fitVol(2) = 10.^fitVol(2);
plotX = min(tempX):0.001:max(tempX);
legLines3 = plot(plotX, (fitVol(2)*plotX.^fitVol(1)),'-','color','k','linewidth',2);
[R,P] = corrcoef(log10(tempX),log10(tempY));
text(0.51,90,sprintf('Y=%.2fX^%.2f R^2=%.2f, p=%.3f',fitVol(2),fitVol(1), R(2)^2, P(2)),'color', 'k','Interpreter', 'none')
uistack(legLines3,'down');

xlim([0.5 0.9]); ylim([0 100]);
set(gca, 'YTick', [0 25 50 75 100]);
set(gca,'XTick',[0.5 0.6 0.7 0.8 0.9],'TickDir','out');
xlabel('Eye Volume^0^.^3 (mm)'); ylabel('CC Thickness (um)');
legend(labels, groupNames);
set(gcf, 'Position', [50 50 500 500/1.61 ]);

% plot summary for curvature
curvHistF = figure; hold on

HistBins = 0:40:2000; 
labels = zeros(length(regressionGroups),1);

for i = 1:length(regressionGroups)
    tempX = ([sDat.fullBeeData(regGroup == i).EyeVolume]/10^9).^(1/3);
    inds = find(regGroup == i);
    tempY = zeros(length(inds),1);
    for j = 1:length(inds)
        tempY(j) = wmean(sDat.fullBeeData(inds(j)).avgRadius,facetWeightInds{inds(j)}); 
        [n,x] = histwc(sDat.fullBeeData(inds(j)).avgRadius,facetWeightInds{inds(j)},HistBins);
        
        relN = smooth(n/sum(n),smthL);
        %adjust for diff between max n range and intended cal
        xAdj = (max(relN)/histCal);
        
        if ~flipChange(inds(j))
            dirToP = 'left';
            histWToUse = histW*xAdj;
            xToP = tempX(j)-histWToUse/2;
        else
            dirToP = 'right';
            histWToUse = histW*xAdj;
            xToP = tempX(j)+histWToUse/2;
        end
        distributionPlot({[x',relN]}, 'xValues', xToP,'distWidth', histWToUse, 'showMM',0,'histOri',dirToP,'globalNorm',3,'color',histCol)

    end
    labels(i) = plot(tempX, tempY, groupMarkerStyles{i}, 'color', 'k', 'markersize',10,'markerfacecolor','w');
end
fitVol = polyfit(log10(tempX),log10(tempY)',1);
fitVol(2) = 10.^fitVol(2);
plotX = min(tempX):0.001:max(tempX);
legLines3 = plot(plotX, (fitVol(2)*plotX.^fitVol(1)),'-','color','k','linewidth',2);
[R,P] = corrcoef(log10(tempX),log10(tempY));
text(0.51,1700,sprintf('Y=%.2fX^%.2f R^2=%.2f, p=%.3f',fitVol(2),fitVol(1), R(2)^2, P(2)),'color', 'k','Interpreter', 'none')
uistack(legLines3,'down');

xlim([0.5 0.9]); ylim([0 2000]);
set(gca, 'YTick', [0 500 1000 1500 2000]);
set(gca,'XTick',[0.5 0.6 0.7 0.8 0.9],'TickDir','out');
xlabel('Eye Volume^0^.^3 (mm)'); ylabel('Curvature (um)');
legend(labels, groupNames);
set(gcf, 'Position', [50 50 500 500/1.61 ]);

% plot summary for eye parameter
eyeParamHistF = figure; hold on

HistBins = 0:0.025:2; 
labels = zeros(length(regressionGroups),1);

for i = 1:length(regressionGroups)
    tempX = ([sDat.fullBeeData(regGroup == i).EyeVolume]/10^9).^(1/3);
    inds = find(regGroup == i);
    tempY = zeros(length(inds),1);
    for j = 1:length(inds)
        
        tempY(j) = wmean(sDat.fullBeeData(inds(j)).calInterO.*sDat.fullBeeData(inds(j)).avgDiameter,facetWeightInds{inds(j)}); 
        [n,x] = histwc(sDat.fullBeeData(inds(j)).calInterO.*sDat.fullBeeData(inds(j)).avgDiameter,facetWeightInds{inds(j)},HistBins);
        
        relN = smooth(n/sum(n),smthL);
        %adjust for diff between max n range and intended cal
        xAdj = (max(relN)/histCal);
        
        if ~flipChange(inds(j))
            dirToP = 'left';
            histWToUse = histW*xAdj;
            xToP = tempX(j)-histWToUse/2;
        else
            dirToP = 'right';
            histWToUse = histW*xAdj;
            xToP = tempX(j)+histWToUse/2;
        end
        distributionPlot({[x',relN]}, 'xValues', xToP,'distWidth', histWToUse, 'showMM',0,'histOri',dirToP,'globalNorm',3,'color',histCol)

    end
    labels(i) = plot(tempX, tempY, groupMarkerStyles{i}, 'color', 'k', 'markersize',10,'markerfacecolor','w');
end
fitVol = polyfit(log10(tempX),log10(tempY)',1);
fitVol(2) = 10.^fitVol(2);
plotX = min(tempX):0.001:max(tempX);
legLines3 = plot(plotX, (fitVol(2)*plotX.^fitVol(1)),'-','color','k','linewidth',2);
[R,P] = corrcoef(log10(tempX),log10(tempY));
text(0.51,1.75,sprintf('Y=%.2fX^%.2f R^2=%.2f, p=%.3f',fitVol(2),fitVol(1), R(2)^2, P(2)),'color', 'k','Interpreter', 'none')
uistack(legLines3,'down');

xlim([0.5 0.9]); ylim([0 2]);
set(gca, 'YTick', [0 0.5 1 1.5 2]);
set(gca,'XTick',[0.5 0.6 0.7 0.8 0.9],'TickDir','out');
xlabel('Eye Volume^0^.^3 (mm)'); ylabel('Eye parameter (um.rad)');
legend(labels, groupNames);
set(gcf, 'Position', [50 50 500 500/1.61 ]);
%% map fov onto sphere and then onto projections
% make world plots on sphere showing borders of fov for mono, bino and full, and integrate across visual space.

smoothRngFOV = 3;

%map original sphere to az and el - generic
sphereEls = acos(sDat.refSphere.X(:,2))/pi*180-90;
sphereAzs = -atan2(sDat.refSphere.X(:,1),-sDat.refSphere.X(:,3))/pi*180; %note that this is negated to get correct az in image
topInd = find(sphereEls == 90);

azBins = -180:angStep:180;
elBins = -90:angStep:90;

azSumT = zeros(length(azBins)-1,1);
elSumT = zeros(length(elBins)-1,1);
for j = 1:length(azBins)-1
    azSumT(j) = length(find(sphereAzs > azBins(j) & sphereAzs <= azBins(j+1)));
end
for j = 1:length(elBins)-1
    elSumT(j) = length(find(sphereEls > elBins(j) & sphereEls <= elBins(j+1)));
end

% plot full fov on sphere
worldSphereF = figure; 
set(gcf, 'Position', [50 50 850 850]);
 hold on; axis off
patch('Faces',sDat.refSphere.Triangulation,'Vertices',sDat.refSphere.X,'FaceAlpha',0.8,'EdgeAlpha',0,'FaceColor',[0.95 0.95 0.95]);
    view([0, -90]); axis equal
set(gcf,'color','w');
for i = [-180 -150 -120 -90 -30 -60 0 30 60 90 120 150 180]
    azs = i*ones(1,361);
    els = -180:180;
    if i == 180 | i == 0 
        lw = 3;
    else
        lw = 1.5;
    end
    plot3( sin(els/180*pi).*sin(azs/180*pi), cos(els/180*pi),sin(els/180*pi).*cos(azs/180*pi) , '-','color', [0.75 0.75 0.75],'linewidth',lw)  
end
for i = [-90 -60 -30 0 30 60 90]
    azs = -180:180;
    els = i*ones(1,361);
    if i == 0
        lw = 4;
    else
        lw =1.5;
    end
    plot3( sin(els/180*pi).*sin(azs/180*pi),cos(els/180*pi),sin(els/180*pi).*cos(azs/180*pi), '-','color', [0.75 0.75 0.75],'linewidth',lw)  
    %do second half of hemisphere
    if i < 0
        els = (i-90)*ones(1,361);
    else
        els = (i+90)*ones(1,361);
    end
    plot3( sin(els/180*pi).*sin(azs/180*pi),cos(els/180*pi),sin(els/180*pi).*cos(azs/180*pi),'-', 'color', [0.75 0.75 0.75],'linewidth',lw)  
end
for i = 1:8   
    if regGroup(i) == 1
        beeCol = 'k';
    else
        beeCol = volCols(round((eyeVols(i)-volRng(1))./(volRng(2)-volRng(1))*(size(volCols,1)-1))+1,:);
    end 
    %remove last ind (the closing one) before smooth then add after, otherwise smooth opens it up again
    tempSubs = [circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,1)',smoothRngFOV), circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,2)',smoothRngFOV), circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,3)',smoothRngFOV)];
    tempSubs = vertcat(tempSubs,tempSubs(1,:));
    plot3(tempSubs(:,1), tempSubs(:,2), tempSubs(:,3),'-','color',beeCol','linewidth',4);
end

%plot full FOV on proj
projFullFov = figure; 
hold on; xlim([-180 180]); ylim([-90 90]);
xlabel('Azimuth (^o)');
ylabel('Elevation (^o)');
set(gcf, 'Position', [50 50 1000 500]); axis off
set(gca, 'YTick', [],'XTick',[]);
for i = [-180 -150 -120 -90 -30 -60 0 30 60 90 120 150 180]
      plot((i*ones(1,181)).*cos((-90:90)/180*pi), -90:90, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
end
for i = [-90 -60 -30 0 30 60 90]
     plot((-180:180)*cos(i/180*pi), i*ones(361,1), '-','color', [0.75 0.75 0.75],'linewidth',1.5)
end
%add some labels to the meridians
for i = [-180 -120 -60 0 60 120 180]
   text(i,0,sprintf('%i^o',i)); 
end
for i = [-90 -60 -30 30 60 90]
   text(-180*cos(i/180*pi),i,sprintf('%i^o',i)); 
end
for i = 1:8   
    if regGroup(i) == 1
        beeCol = 'k';
    else
        beeCol = volCols(round((eyeVols(i)-volRng(1))./(volRng(2)-volRng(1))*(size(volCols,1)-1))+1,:);
    end
    tempSubs = [circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,1)',smoothRngFOV), circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,2)',smoothRngFOV), circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,3)',smoothRngFOV)];
    tempSubs = vertcat(tempSubs,tempSubs(1,:));
    lineEls = acos(tempSubs(:,2))/pi*180-90;
    lineAzs = -atan2(tempSubs(:,1),-tempSubs(:,3))/pi*180;
    
    %check if line crosses top of pole and insert 90o point if so
    if sum(sDat.fullBeeData(i).inFOV == topInd)
        %if so find point on left of linkage
        tempInds = find(lineAzs < 0);
        [~,leftLink] = max(lineEls(tempInds));
        leftLink = tempInds(leftLink);
        %and then on right
        tempInds = find(lineAzs >= 0);
        [~,rightLink] = max(lineEls(tempInds));
        rightLink = tempInds(rightLink);
        
        if rightLink > leftLink
            lineEls = [lineEls(1:leftLink)' 90 lineEls(rightLink:end)'];    
            lineAzs = [lineAzs(1:leftLink)' 0 lineAzs(rightLink:end)'];
        else
            lineEls = [lineEls(1:rightLink)' 90 lineEls(leftLink:end)'];    
            lineAzs = [lineAzs(1:rightLink)' 0 lineAzs(leftLink:end)'];
        end
    end
  
    plot(lineAzs.*cos(lineEls/180*pi),lineEls,'-','color',beeCol','linewidth',2);
end

minDst = 0;
for i = round(1:size(sDat.refSphere.X,1)/50:size(sDat.refSphere.X,1));
    toTest = 1:size(sDat.refSphere.X,1); toTest(i) = [];
    minDst = minDst + min(sqrt((sDat.refSphere.X(i,1) - sDat.refSphere.X(toTest,1)).^2 + (sDat.refSphere.X(i,2) - sDat.refSphere.X(toTest,2)).^2 + (sDat.refSphere.X(i,3) - sDat.refSphere.X(toTest,3)).^2));
end
minDst = minDst/length(1:100:size(sDat.refSphere.X,1));

% plot bincular FOV on projection
projBinoFOV = figure; 
hold on; xlim([-180 180]); ylim([-90 90]);
xlabel('Azimuth (^o)');
ylabel('Elevation (^o)');
set(gcf, 'Position', [50 50 1000 500]); axis off
set(gca, 'YTick', [],'XTick',[]);
for i = [-180 -150 -120 -90 -30 -60 0 30 60 90 120 150 180]
      plot((i*ones(1,181)).*cos((-90:90)/180*pi), -90:90, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
end
for i = [-90 -60 -30 0 30 60 90]
     plot((-180:180)*cos(i/180*pi), i*ones(361,1), '-','color', [0.75 0.75 0.75],'linewidth',1.5)
end
%add some labels to the meridians
for i = [-180 -120 -60 0 60 120 180]
   text(i,0,sprintf('%i^o',i)); 
end
for i = [-90 -60 -30 30 60 90]
   text(-180*cos(i/180*pi),i,sprintf('%i^o',i)); 
end

mapResultsForBees = cell(8,1);
numRegionsForBees = zeros(8,1);

for i = 1:8   
    if regGroup(i) == 1
        beeCol = 'k';
    else
        beeCol = volCols(round((eyeVols(i)-volRng(1))./(volRng(2)-volRng(1))*(size(volCols,1)-1))+1,:);
    end
    %get iso line
    tempFOV = zeros(size(sDat.refSphere.X,1),1); tempFOV(sDat.fullBeeData(i).inFOVBino) = 1;
    %account for holes
    [leftMapResults, numLeftReg] = calculateRegionsInMap( sDat.refSphere.X, tempFOV, minDst*4 );
    mapResultsForBees{i} = leftMapResults;
    numRegionsForBees(i) = numLeftReg;
    
    for j = 1:numLeftReg
        tempSubs = leftMapResults{j,1};
        tempSubs = [circularSmooth(tempSubs(1:end-1,1)',smoothRngFOV), circularSmooth(tempSubs(1:end-1,2)',smoothRngFOV), circularSmooth(tempSubs(1:end-1,3)',smoothRngFOV)];
        tempSubs = vertcat(tempSubs,tempSubs(1,:));
        lineEls = acos(tempSubs(:,2))/pi*180-90;
        lineAzs = -atan2(tempSubs(:,1),-tempSubs(:,3))/pi*180;

        %check if line crosses top of pole and insert 90o point if so
        if sum(leftMapResults{j,3} == topInd)
            %if so find point on left of linkage
            tempInds = find(lineAzs < 0);
            [~,leftLink] = max(lineEls(tempInds));
            leftLink = tempInds(leftLink);
            %and then on right
            tempInds = find(lineAzs >= 0);
            [~,rightLink] = max(lineEls(tempInds));
            rightLink = tempInds(rightLink);

            if rightLink > leftLink
                lineEls = [lineEls(1:leftLink)' 90 lineEls(rightLink:end)'];    
                lineAzs = [lineAzs(1:leftLink)' 0 lineAzs(rightLink:end)'];
            else
                lineEls = [lineEls(1:rightLink)' 90 lineEls(leftLink:end)'];    
                lineAzs = [lineAzs(1:rightLink)' 0 lineAzs(leftLink:end)'];
            end
        end

        plot(lineAzs.*cos(lineEls/180*pi),lineEls,'-','color',beeCol','linewidth',2);
    end
end

%% plot profiles across entire mono fov

%integration of FULL FOV
azIntFullFOVF = figure; hold on
for i = [-180 -150 -120 -90 -30 -60 0 30 60 90 120 150 180]
   line([i i], [0 100], 'color', [0.75 0.75 0.75],'linewidth',1.5)  
end
ylabel('% Az. band in FOV');
xlabel('Azimuth band (^o)');
xlim([-180 180]); ylim([0 100]);
set(gca, 'YTick', [0 25 50 75 100]);
set(gca,'XTick',[-180 -120 -60 0 60 120 180],'TickDir','out');
set(gcf, 'Position', [50 50 1000 300]);

elIntFullFOVF = figure; hold on
for i =[-90 -60 -30 0 30 60 90]
   line([0 100], [i i], 'color', [0.75 0.75 0.75],'linewidth',1.5)  
end
set(gca, 'YTick', [-90 -60 -30 0 30 60 90]);
xlim([0 100]); ylim([-90 90]);
set(gca,'XTick',[0 20 50 75 100],'TickDir','out');
set(gcf, 'Position', [50 50 300 500]);

xlabel('% El. band in FOV');
ylabel('Elevation band (^o)');

for i = 1:8   
    if regGroup(i) == 1
        beeCol = 'k';
    else
        beeCol = volCols(round((eyeVols(i)-volRng(1))./(volRng(2)-volRng(1))*(size(volCols,1)-1))+1,:);
    end 
    
    pointEls = acos(sDat.refSphere.X(sDat.fullBeeData(i).inFOV,2))/pi*180-90;
    pointAzs = -atan2(sDat.refSphere.X(sDat.fullBeeData(i).inFOV,1),-sDat.refSphere.X(sDat.fullBeeData(i).inFOV,3))/pi*180;
    
    azSum = zeros(length(azBins)-1,1);
    for j = 1:length(azBins)-1
        azSum(j) = length(find(pointAzs > azBins(j) & pointAzs <= azBins(j+1)))./azSumT(j)*100;
    end
    elSum = zeros(length(elBins)-1,1);
    for j = 1:length(elBins)-1
        elSum(j) = length(find(pointEls > elBins(j) & pointEls <= elBins(j+1)))./elSumT(j)*100;
    end
    figure(azIntFullFOVF);
    plot(azBins(1:end-1)+angStep/2, azSum, '-','color',beeCol, 'linewidth',2);
    figure(elIntFullFOVF);
    plot(elSum,elBins(1:end-1)+angStep/2, '-','color',beeCol','linewidth',2);
end

%  integration across  bino fov
azIntBinoFOVF = figure; hold on
for i = [-180 -150 -120 -90 -30 -60 0 30 60 90 120 150 180]
   line([i i], [0 100], 'color', [0.75 0.75 0.75],'linewidth',1.5)  
end
ylabel('% Az. band in Binocular FOV');
xlabel('Azimuth band (^o)');
xlim([-180 180]); ylim([0 100]);
set(gca, 'YTick', [0 25 50 75 100]);
set(gca,'XTick',[-180 -120 -60 0 60 120 180],'TickDir','out');
set(gcf, 'Position', [50 50 1000 300]);

elIntBinoFOVF = figure; hold on
for i =[-90 -60 -30 0 30 60 90]
   line([0 100], [i i], 'color', [0.75 0.75 0.75],'linewidth',1.5)  
end
set(gca, 'YTick', [-90 -60 -30 0 30 60 90]);
xlim([0 100]); ylim([-90 90]);
set(gca,'XTick',[0 20 50 75 100],'TickDir','out');
set(gcf, 'Position', [50 50 300 500]);

xlabel('% El. band in Binocular FOV');
ylabel('Elevation band (^o)');

for i = 1:8   
    if regGroup(i) == 1
        beeCol = 'k';
    else
        beeCol = volCols(round((eyeVols(i)-volRng(1))./(volRng(2)-volRng(1))*(size(volCols,1)-1))+1,:);
    end 
    
    pointEls = acos(sDat.refSphere.X(sDat.fullBeeData(i).inFOVBino,2))/pi*180-90;
    pointAzs = -atan2(sDat.refSphere.X(sDat.fullBeeData(i).inFOVBino,1),-sDat.refSphere.X(sDat.fullBeeData(i).inFOVBino,3))/pi*180;
    
    azSum = zeros(length(azBins)-1,1);
    for j = 1:length(azBins)-1
        azSum(j) = length(find(pointAzs > azBins(j) & pointAzs <= azBins(j+1)))./azSumT(j)*100;
    end
    elSum = zeros(length(elBins)-1,1);
    for j = 1:length(elBins)-1
        elSum(j) = length(find(pointEls > elBins(j) & pointEls <= elBins(j+1)))./elSumT(j)*100;
    end
    figure(azIntBinoFOVF);
    plot(azBins(1:end-1)+angStep/2, azSum, '-','color',beeCol, 'linewidth',2);
    figure(elIntBinoFOVF);
    plot(elSum,elBins(1:end-1)+angStep/2, '-','color',beeCol','linewidth',2);
end

% plot avearge IO across field of view
azAvIOF = figure; hold on
for i = [-180 -150 -120 -90 -30 -60 0 30 60 90 120 150 180]
   line([i i], [0 4], 'color', [0.75 0.75 0.75],'linewidth',1.5)  
end
ylabel('Mean interommatidial angle across Az. band (^o)');
xlabel('Azimuth band (^o)');
xlim([-180 180]); ylim([0 4]);
set(gca, 'YTick', [0 1 2 3 4]);
set(gca,'XTick',[-180 -120 -60 0 60 120 180],'TickDir','out');
set(gcf, 'Position', [50 50 1000 300]);

elAvIOF = figure; hold on
for i =[-90 -60 -30 0 30 60 90]
   line([0 4], [i i], 'color', [0.75 0.75 0.75],'linewidth',1.5)  
end
set(gca, 'YTick', [-90 -60 -30 0 30 60 90]);
xlim([0 4]); ylim([-90 90]);
set(gca,'XTick',[0 1 2 3 4],'TickDir','out');
set(gcf, 'Position', [50 50 300 500]);

xlabel('Mean interommatidial angle across El. band (^o)');
ylabel('Elevation band (^o)');

for i = 1:8   
    if regGroup(i) == 1
        beeCol = 'k';
    else
        beeCol = volCols(round((eyeVols(i)-volRng(1))./(volRng(2)-volRng(1))*(size(volCols,1)-1))+1,:);
    end 
    
    pointEls = acos(sDat.refSphere.X(sDat.fullBeeData(i).inFOV,2))/pi*180-90;
    pointAzs = -atan2(sDat.refSphere.X(sDat.fullBeeData(i).inFOV,1),-sDat.refSphere.X(sDat.fullBeeData(i).inFOV,3))/pi*180;
    weights = projWeightInds{i};
    
    azSum = zeros(length(azBins)-1,1);
    for j = 1:length(azBins)-1
        inds = find(pointAzs > azBins(j) & pointAzs <= azBins(j+1));
        if ~isempty(inds)
            azSum(j) = wmean(sDat.fullBeeData(i).InterOOnSphere(inds),weights(inds))/pi*180;
        end
    end
    elSum = zeros(length(elBins)-1,1);
    for j = 1:length(elBins)-1
        inds = find(pointEls > elBins(j) & pointEls <= elBins(j+1));
        if ~isempty(inds)
            elSum(j) = wmean(sDat.fullBeeData(i).InterOOnSphere(inds),weights(inds))/pi*180;
        end
    end
    if remIrregularSpatialBands
        %set all after first drop out on az to zero
        %assumes that drop outs are just scattered points near poles...
        testD = azSum;
        testD(testD > 0) = 1;
        testD = diff(testD);
        testD = find(testD < 0);
        if ~isempty(testD)
           azSum(testD(1)+1:end) = NaN; 
        end
    end
    
    azSum(azSum == 0) = NaN;
    elSum(elSum == 0) = NaN;
    
    figure(azAvIOF);
    plot(azBins(1:end-1)+angStep/2, azSum, '-','color',beeCol, 'linewidth',2);
    figure(elAvIOF);
    plot(elSum,elBins(1:end-1)+angStep/2, '-','color',beeCol','linewidth',2);
end

% plot avearge facet diameter across field of view
azAvDF = figure; hold on
for i = [-180 -150 -120 -90 -30 -60 0 30 60 90 120 150 180]
   line([i i], [0 300], 'color', [0.75 0.75 0.75],'linewidth',1.5)  
end
ylabel('Mean facet diameter across Az. band (um)');
xlabel('Azimuth band (^o)');
xlim([-180 180]); ylim([0 30]);
set(gca, 'YTick', [0 10 20 30]);
set(gca,'XTick',[-180 -120 -60 0 60 120 180],'TickDir','out');
set(gcf, 'Position', [50 50 1000 300]);

elAvDF = figure; hold on
for i =[-90 -60 -30 0 30 60 90]
   line([0 300], [i i], 'color', [0.75 0.75 0.75],'linewidth',1.5)  
end
set(gca, 'YTick', [-90 -60 -30 0 30 60 90]);
xlim([0 30]); ylim([-90 90]);
set(gca,'XTick',[0 10 20 30],'TickDir','out');
set(gcf, 'Position', [50 50 300 500]);

xlabel('Mean facet diameter across El. band (um)');
ylabel('Elevation band (^o)');

for i = 1:8   
    if regGroup(i) == 1
        beeCol = 'k';
    else
        beeCol = volCols(round((eyeVols(i)-volRng(1))./(volRng(2)-volRng(1))*(size(volCols,1)-1))+1,:);
    end 
    
    pointEls = acos(sDat.refSphere.X(sDat.fullBeeData(i).inFOV,2))/pi*180-90;
    pointAzs = -atan2(sDat.refSphere.X(sDat.fullBeeData(i).inFOV,1),-sDat.refSphere.X(sDat.fullBeeData(i).inFOV,3))/pi*180;
    weights = projWeightInds{i};
    
    azSum = zeros(length(azBins)-1,1);
    for j = 1:length(azBins)-1
        inds = find(pointAzs > azBins(j) & pointAzs <= azBins(j+1));
        if ~isempty(inds)
            azSum(j) = wmean(sDat.fullBeeData(i).diameterOnSphere(inds),weights(inds));
        end
    end
    elSum = zeros(length(elBins)-1,1);
    for j = 1:length(elBins)-1
        inds = find(pointEls > elBins(j) & pointEls <= elBins(j+1));
        if ~isempty(inds)
            elSum(j) = wmean(sDat.fullBeeData(i).diameterOnSphere(inds),weights(inds));
        end
    end
    if remIrregularSpatialBands
        testD = azSum;
        testD(testD > 0) = 1;
        testD = diff(testD);
        testD = find(testD < 0);
        if ~isempty(testD)
           azSum(testD(1)+1:end) = NaN; 
        end
    end
    
    azSum(azSum == 0) = NaN;
    elSum(elSum == 0) = NaN;
    
    figure(azAvDF);
    plot(azBins(1:end-1)+angStep/2, azSum, '-','color',beeCol, 'linewidth',2);
    figure(elAvDF);
    plot(elSum,elBins(1:end-1)+angStep/2, '-','color',beeCol','linewidth',2);
end

% avearge ret thickness across field of view
azAvRetF = figure; hold on
for i = [-180 -150 -120 -90 -30 -60 0 30 60 90 120 150 180]
   line([i i], [0 800], 'color', [0.75 0.75 0.75],'linewidth',1.5)  
end
ylabel('Mean retina thickness across Az. band (um)');
xlabel('Azimuth band (^o)');
xlim([-180 180]); ylim([0 500]);
set(gca, 'YTick', [0 100 200 300 400 500]);
set(gca,'XTick',[-180 -120 -60 0 60 120 180],'TickDir','out');
set(gcf, 'Position', [50 50 1000 300]);

elAvRetF = figure; hold on
for i =[-90 -60 -30 0 30 60 90]
   line([0 800], [i i], 'color', [0.75 0.75 0.75],'linewidth',1.5)  
end
set(gca, 'YTick', [-90 -60 -30 0 30 60 90]);
xlim([0 500]); ylim([-90 90]);
set(gca,'XTick',[0 100 200 300 400 500],'TickDir','out');
set(gcf, 'Position', [50 50 300 500]);

xlabel('Mean retina thickness across El. band (um)');
ylabel('Elevation band (^o)');

for i = 1:8   
    if regGroup(i) == 1
        beeCol = 'k';
    else
        beeCol = volCols(round((eyeVols(i)-volRng(1))./(volRng(2)-volRng(1))*(size(volCols,1)-1))+1,:);
    end 
    
    pointEls = acos(sDat.refSphere.X(sDat.fullBeeData(i).inFOV,2))/pi*180-90;
    pointAzs = -atan2(sDat.refSphere.X(sDat.fullBeeData(i).inFOV,1),-sDat.refSphere.X(sDat.fullBeeData(i).inFOV,3))/pi*180;
    weights = projWeightInds{i};
    
    azSum = zeros(length(azBins)-1,1);
    for j = 1:length(azBins)-1
        inds = find(pointAzs > azBins(j) & pointAzs <= azBins(j+1));
        if ~isempty(inds)
            azSum(j) = wmean(sDat.fullBeeData(i).retTOnSphere(inds),weights(inds));
        end
    end
    elSum = zeros(length(elBins)-1,1);
    for j = 1:length(elBins)-1
        inds = find(pointEls > elBins(j) & pointEls <= elBins(j+1));
        if ~isempty(inds)
            elSum(j) = wmean(sDat.fullBeeData(i).retTOnSphere(inds),weights(inds));
        end
    end
    if remIrregularSpatialBands
        testD = azSum;
        testD(testD > 0) = 1;
        testD = diff(testD);
        testD = find(testD < 0);
        if ~isempty(testD)
           azSum(testD(1)+1:end) = NaN; 
        end
    end
    
    azSum(azSum == 0) = NaN;
    elSum(elSum == 0) = NaN;
    
    figure(azAvRetF);
    plot(azBins(1:end-1)+angStep/2, azSum, '-','color',beeCol, 'linewidth',2);
    figure(elAvRetF);
    plot(elSum,elBins(1:end-1)+angStep/2, '-','color',beeCol','linewidth',2);
end

% plot avearge CC thickness across field of view
azAvCCF = figure; hold on
for i = [-180 -150 -120 -90 -30 -60 0 30 60 90 120 150 180]
   line([i i], [0 800], 'color', [0.75 0.75 0.75],'linewidth',1.5)  
end
ylabel('Mean CC thickness across Az. band (um)');
xlabel('Azimuth band (^o)');
xlim([-180 180]); ylim([0 100]);
set(gca, 'YTick', [0 25 50 75 100]);
set(gca,'XTick',[-180 -120 -60 0 60 120 180],'TickDir','out');
set(gcf, 'Position', [50 50 1000 300]);

elAvCCF = figure; hold on
for i =[-90 -60 -30 0 30 60 90]
   line([0 800], [i i], 'color', [0.75 0.75 0.75],'linewidth',1.5)  
end
set(gca, 'YTick', [-90 -60 -30 0 30 60 90]);
xlim([0 100]); ylim([-90 90]);
set(gca,'XTick',[0 25 50 75 100],'TickDir','out');
set(gcf, 'Position', [50 50 300 500]);

xlabel('Mean retina thickness across El. band (um)');
ylabel('Elevation band (^o)');

for i = 1:8   
    if regGroup(i) == 1
        beeCol = 'k';
    else
        beeCol = volCols(round((eyeVols(i)-volRng(1))./(volRng(2)-volRng(1))*(size(volCols,1)-1))+1,:);
    end 
    
    pointEls = acos(sDat.refSphere.X(sDat.fullBeeData(i).inFOV,2))/pi*180-90;
    pointAzs = -atan2(sDat.refSphere.X(sDat.fullBeeData(i).inFOV,1),-sDat.refSphere.X(sDat.fullBeeData(i).inFOV,3))/pi*180;
    weights = projWeightInds{i};
    
    azSum = zeros(length(azBins)-1,1);
    for j = 1:length(azBins)-1
        inds = find(pointAzs > azBins(j) & pointAzs <= azBins(j+1));
        if ~isempty(inds)
            azSum(j) = wmean(sDat.fullBeeData(i).CCTOnSphere(inds),weights(inds));
        end
    end
    elSum = zeros(length(elBins)-1,1);
    for j = 1:length(elBins)-1
        inds = find(pointEls > elBins(j) & pointEls <= elBins(j+1));
        if ~isempty(inds)
            elSum(j) = wmean(sDat.fullBeeData(i).CCTOnSphere(inds),weights(inds));
        end
    end
    if remIrregularSpatialBands
        testD = azSum;
        testD(testD > 0) = 1;
        testD = diff(testD);
        testD = find(testD < 0);
        if ~isempty(testD)
           azSum(testD(1)+1:end) = NaN; 
        end
    end
    
    azSum(azSum == 0) = NaN;
    elSum(elSum == 0) = NaN;
    
    figure(azAvCCF);
    plot(azBins(1:end-1)+angStep/2, azSum, '-','color',beeCol, 'linewidth',2);
    figure(elAvCCF);
    plot(elSum,elBins(1:end-1)+angStep/2, '-','color',beeCol','linewidth',2);
end
% plot avearge lens thickness across field of view
azAvLensF = figure; hold on
for i = [-180 -150 -120 -90 -30 -60 0 30 60 90 120 150 180]
   line([i i], [0 800], 'color', [0.75 0.75 0.75],'linewidth',1.5)  
end
ylabel('Mean lens thickness across Az. band (um)');
xlabel('Azimuth band (^o)');
xlim([-180 180]); ylim([0 100]);
set(gca, 'YTick', [0 25 50 75 100]);
set(gca,'XTick',[-180 -120 -60 0 60 120 180],'TickDir','out');
set(gcf, 'Position', [50 50 1000 300]);

elAvLensF = figure; hold on
for i =[-90 -60 -30 0 30 60 90]
   line([0 800], [i i], 'color', [0.75 0.75 0.75],'linewidth',1.5)  
end
set(gca, 'YTick', [-90 -60 -30 0 30 60 90]);
xlim([0 100]); ylim([-90 90]);
set(gca,'XTick',[0 25 50 75 100],'TickDir','out');
set(gcf, 'Position', [50 50 300 500]);

xlabel('Mean lens thickness across El. band (um)');
ylabel('Elevation band (^o)');

for i = 1:8   
    if regGroup(i) == 1
        beeCol = 'k';
    else
        beeCol = volCols(round((eyeVols(i)-volRng(1))./(volRng(2)-volRng(1))*(size(volCols,1)-1))+1,:);
    end 
    
    pointEls = acos(sDat.refSphere.X(sDat.fullBeeData(i).inFOV,2))/pi*180-90;
    pointAzs = -atan2(sDat.refSphere.X(sDat.fullBeeData(i).inFOV,1),-sDat.refSphere.X(sDat.fullBeeData(i).inFOV,3))/pi*180;
    weights = projWeightInds{i};
    
    azSum = zeros(length(azBins)-1,1);
    for j = 1:length(azBins)-1
        inds = find(pointAzs > azBins(j) & pointAzs <= azBins(j+1));
        if ~isempty(inds)
            azSum(j) = wmean(sDat.fullBeeData(i).lensTOnSphere(inds),weights(inds));
        end
    end
    elSum = zeros(length(elBins)-1,1);
    for j = 1:length(elBins)-1
        inds = find(pointEls > elBins(j) & pointEls <= elBins(j+1));
        if ~isempty(inds)
            elSum(j) = wmean(sDat.fullBeeData(i).lensTOnSphere(inds),weights(inds));
        end
    end
    if remIrregularSpatialBands
        testD = azSum;
        testD(testD > 0) = 1;
        testD = diff(testD);
        testD = find(testD < 0);
        if ~isempty(testD)
           azSum(testD(1)+1:end) = NaN; 
        end
    end
    
    azSum(azSum == 0) = NaN;
    elSum(elSum == 0) = NaN;
    
    figure(azAvLensF);
    plot(azBins(1:end-1)+angStep/2, azSum, '-','color',beeCol, 'linewidth',2);
    figure(elAvLensF);
    plot(elSum,elBins(1:end-1)+angStep/2, '-','color',beeCol','linewidth',2);
end
% plot avearge curvature across field of view
azAvCurvF = figure; hold on
for i = [-180 -150 -120 -90 -30 -60 0 30 60 90 120 150 180]
   line([i i], [0 2000], 'color', [0.75 0.75 0.75],'linewidth',1.5)  
end
ylabel('Mean curvature across Az. band (um)');
xlabel('Azimuth band (^o)');
xlim([-180 180]); ylim([0 2000]);
set(gca, 'YTick', [0 500 1000 1500 2000]);
set(gca,'XTick',[-180 -120 -60 0 60 120 180],'TickDir','out');
set(gcf, 'Position', [50 50 1000 300]);

elAvCurvF = figure; hold on
for i =[-90 -60 -30 0 30 60 90]
   line([0 2000], [i i], 'color', [0.75 0.75 0.75],'linewidth',1.5)  
end
set(gca, 'YTick', [-90 -60 -30 0 30 60 90]);
xlim([0 2000]); ylim([-90 90]);
set(gca,'XTick',[0 500 1000 1500 2000],'TickDir','out');
set(gcf, 'Position', [50 50 300 500]);

xlabel('Mean curvature across El. band (um)');
ylabel('Elevation band (^o)');

for i = 1:8   
    if regGroup(i) == 1
        beeCol = 'k';
    else
        beeCol = volCols(round((eyeVols(i)-volRng(1))./(volRng(2)-volRng(1))*(size(volCols,1)-1))+1,:);
    end 
    
    pointEls = acos(sDat.refSphere.X(sDat.fullBeeData(i).inFOV,2))/pi*180-90;
    pointAzs = -atan2(sDat.refSphere.X(sDat.fullBeeData(i).inFOV,1),-sDat.refSphere.X(sDat.fullBeeData(i).inFOV,3))/pi*180;
    weights = projWeightInds{i};
    
    azSum = zeros(length(azBins)-1,1);
    for j = 1:length(azBins)-1
        inds = find(pointAzs > azBins(j) & pointAzs <= azBins(j+1));
        if ~isempty(inds)
            azSum(j) = wmean(sDat.fullBeeData(i).curvatureOnSphere(inds),weights(inds));
        end
    end
    elSum = zeros(length(elBins)-1,1);
    for j = 1:length(elBins)-1
        inds = find(pointEls > elBins(j) & pointEls <= elBins(j+1));
        if ~isempty(inds)
            elSum(j) = wmean(sDat.fullBeeData(i).curvatureOnSphere(inds),weights(inds));
        end
    end
    if remIrregularSpatialBands
        testD = azSum;
        testD(testD > 0) = 1;
        testD = diff(testD);
        testD = find(testD < 0);
        if ~isempty(testD)
           azSum(testD(1)+1:end) = NaN; 
        end
    end
    
    azSum(azSum == 0) = NaN;
    elSum(elSum == 0) = NaN;
    
    figure(azAvCurvF);
    plot(azBins(1:end-1)+angStep/2, azSum, '-','color',beeCol, 'linewidth',2);
    figure(elAvCurvF);
    plot(elSum,elBins(1:end-1)+angStep/2, '-','color',beeCol','linewidth',2);
end

%plot avearge eye parameter across field of view
azAvEyeParamF = figure; hold on
for i = [-180 -150 -120 -90 -30 -60 0 30 60 90 120 150 180]
   line([i i], [0 2], 'color', [0.75 0.75 0.75],'linewidth',1.5)  
end
ylabel('Mean eye parameter across Az. band (um.rad)');
xlabel('Azimuth band (^o)');
xlim([-180 180]); ylim([0 2]);
set(gca, 'YTick', [0 0.5 1 1.5 2]);
set(gca,'XTick',[-180 -120 -60 0 60 120 180],'TickDir','out');
set(gcf, 'Position', [50 50 1000 300]);

elAvEyeParamF = figure; hold on
for i =[-90 -60 -30 0 30 60 90]
   line([0 2], [i i], 'color', [0.75 0.75 0.75],'linewidth',1.5)  
end
set(gca, 'YTick', [-90 -60 -30 0 30 60 90]);
xlim([0 2]); ylim([-90 90]);
set(gca,'XTick',[0 0.5 1 1.5 2],'TickDir','out');
set(gcf, 'Position', [50 50 300 500]);

xlabel('Mean eye parameter across El. band (um.rad)');
ylabel('Elevation band (^o)');

for i = 1:8   
    if regGroup(i) == 1
        beeCol = 'k';
    else
        beeCol = volCols(round((eyeVols(i)-volRng(1))./(volRng(2)-volRng(1))*(size(volCols,1)-1))+1,:);
    end 
    
    pointEls = acos(sDat.refSphere.X(sDat.fullBeeData(i).inFOV,2))/pi*180-90;
    pointAzs = -atan2(sDat.refSphere.X(sDat.fullBeeData(i).inFOV,1),-sDat.refSphere.X(sDat.fullBeeData(i).inFOV,3))/pi*180;
    weights = projWeightInds{i};
    
    azSum = zeros(length(azBins)-1,1);
    for j = 1:length(azBins)-1
        inds = find(pointAzs > azBins(j) & pointAzs <= azBins(j+1));
        if ~isempty(inds)
            azSum(j) = wmean(sDat.fullBeeData(i).InterOOnSphere(inds).*sDat.fullBeeData(i).diameterOnSphere(inds),weights(inds));
        end
    end
    elSum = zeros(length(elBins)-1,1);
    for j = 1:length(elBins)-1
        inds = find(pointEls > elBins(j) & pointEls <= elBins(j+1));
        if ~isempty(inds)
            elSum(j) = wmean(sDat.fullBeeData(i).InterOOnSphere(inds).*sDat.fullBeeData(i).diameterOnSphere(inds),weights(inds));
        end
    end
    if remIrregularSpatialBands
        testD = azSum;
        testD(testD > 0) = 1;
        testD = diff(testD);
        testD = find(testD < 0);
        if ~isempty(testD)
           azSum(testD(1)+1:end) = NaN; 
        end
    end
    
    azSum(azSum == 0) = NaN;
    elSum(elSum == 0) = NaN;
    
    figure(azAvEyeParamF);
    plot(azBins(1:end-1)+angStep/2, azSum, '-','color',beeCol, 'linewidth',2);
    figure(elAvEyeParamF);
    plot(elSum,elBins(1:end-1)+angStep/2, '-','color',beeCol','linewidth',2);
end

%% make framework for plotting topology as image.

imAz = (-180:360/(360/imageConvRes):180)*-1;
imEl = (-90:180/(180/imageConvRes):90)*-1;

%create az/el mapping of image to be
imCoordsAz = zeros(length(imAz),length(imEl));
imCoordsEl = zeros(length(imAz),length(imEl));
for i = 1:length(imAz)
    imCoordsEl(i,:) = imEl;
end
for i = 1:length(imEl)
    imCoordsAz(:,i) = imAz;
end

%map image to sphere
imX = cos(imCoordsEl(:)/180*pi).*sin(imCoordsAz(:)/180*pi);
imZ = -cos(imCoordsEl(:)/180*pi).*cos(imCoordsAz(:)/180*pi);
imY = -sin(imCoordsEl(:)/180*pi);

%map original sphere to az and el
sphereEls = acos(sDat.refSphere.X(:,2))/pi*180-90;
sphereAzs = atan2(sDat.refSphere.X(:,1),-sDat.refSphere.X(:,3))/pi*180; %note that this is not negated to get correct az in image

gridAllocImage = zeros(length(imAz), length(imEl)); 
gridAllocImage(:) = assignSubsToVoroniDiagram(sDat.refSphere.X',imX, imY, imZ );

misingINDs = find(gridAllocImage(:) == 0);
goodINDs = find(gridAllocImage(:) > 0);
fovInterpolant = scatteredInterpolant(imCoordsAz(goodINDs), imCoordsEl(goodINDs), gridAllocImage(goodINDs), 'nearest', 'nearest');    
gridAllocImage(misingINDs) = round(fovInterpolant(imCoordsAz(misingINDs), imCoordsEl(misingINDs)));

uniqueFacetAllocs = unique(gridAllocImage(:));

%make new look up grid with cosine mapping
gridAllocImageSineProj = cell(length(uniqueFacetAllocs),1);
anyIndMap = gridAllocImage*0;
for i = 1:length(uniqueFacetAllocs)
    getInds = find(gridAllocImage == uniqueFacetAllocs(i));
    if ~isempty(getInds)
        newAz = imCoordsAz(getInds).*cos(imCoordsEl(getInds)/180*pi);
        %change to cell based approach
        gridAllocImageSineProj{i} = sub2ind(size(gridAllocImage),-round(newAz)+181,-imCoordsEl(getInds)+91);
        anyIndMap(gridAllocImageSineProj{i}) = 1;
    end
end
noInds = find(anyIndMap == 0);

% plot color map of IO 
ColMapToUse = IOColMap;
maxRange2Use = maxRangeIO;
minRange2Use = minRangeIO;

resCMapF = figure;
colormap(ColMapToUse);
hCB = colorbar;
colSteps = (([0 25 50 75 100]/100*(maxRange2Use-minRange2Use) + minRange2Use)/pi*180);
colLabs = {sprintf('%.1f', colSteps(1)), sprintf('%.1f', colSteps(2)), sprintf('%.1f', colSteps(3)),sprintf('%.1f', colSteps(4)),sprintf('%.1f', colSteps(5))};
set(hCB, 'Ticks', [0 25 50 75 100]/100, 'TickLabels',colLabs,'YAxisLocation','right','TickDirection','out');
title('Interommatidial angle (^o)')
IOIndFs = [];

for i = toPlotMaps
    if regGroup(i) == 1
        beeCol = 'k'; 
    else
        beeCol = volCols(round((eyeVols(i)-volRng(1))./(volRng(2)-volRng(1))*(size(volCols,1)-1))+1,:);
    end 
    
    IOIndFs(end+1) = figure;
    
    tempIm = gridAllocImage*0;
    tempNum = gridAllocImage*0;
     for j = 1:length(uniqueFacetAllocs)   
        if uniqueFacetAllocs(j) ~= 0
            linkInd = find(sDat.fullBeeData(i).inFOV == uniqueFacetAllocs(j));
            if ~isempty(linkInd)
                %average overlap
                inds = gridAllocImageSineProj{j};
                tempIm(inds) = tempIm(inds) + sDat.fullBeeData(i).InterOOnSphere(linkInd);
                tempNum(inds) = tempNum(inds) + 1;
            end
        end
     end
     tempIm = tempIm./tempNum; %correct for multiple additions
     tempIm = round((tempIm-minRange2Use)./(maxRange2Use-minRange2Use)*99+1);
     tempIm(tempIm < 1 & ~isnan(tempIm) & ~isinf(tempIm)) = 1; tempIm(tempIm > 100 & ~isnan(tempIm) & ~isinf(tempIm)) = 100;
     tempIm(noInds) = Inf; %point in sphere with not image will go to nan from div by zero
     
     colIm = zeros(size(tempIm,1), size(tempIm,2),3);
     for j = 1:size(tempIm,1)
         for k = 1:size(tempIm,2)
             if isnan(tempIm(j,k))
                colIm(j,k,:) = [0.25 0.25 0.25]; %[255 153 153]/255;
             elseif isinf(tempIm(j,k))
                colIm(j,k,:) = 1;
             else    
                 colIm(j,k,:) = ColMapToUse(tempIm(j,k),:);
             end
         end
     end
     
     %plot border of fov as well.
     imshow_w_by_h(colIm); hold on
     
     if plotFOVBorderOnMap
         %plot the fov border
         tempSubs = [circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,1)',smoothRngFOV), circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,2)',smoothRngFOV), circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,3)',smoothRngFOV)];
        tempSubs = vertcat(tempSubs,tempSubs(1,:));
        lineEls = acos(tempSubs(:,2))/pi*180-90;
        lineAzs = -atan2(tempSubs(:,1),-tempSubs(:,3))/pi*180;

        %check if line crosses top of pole and insert 90o point if so
        if sum(sDat.fullBeeData(i).inFOV == topInd)
            %if so find point on left of linkage
            tempInds = find(lineAzs < 0);
            [~,leftLink] = max(lineEls(tempInds));
            leftLink = tempInds(leftLink);
            %and then on right
            tempInds = find(lineAzs >= 0);
            [~,rightLink] = max(lineEls(tempInds));
            rightLink = tempInds(rightLink);

            if rightLink > leftLink
                lineEls = [lineEls(1:leftLink)' 90 lineEls(rightLink:end)'];    
                lineAzs = [lineAzs(1:leftLink)' 0 lineAzs(rightLink:end)'];
            else
                lineEls = [lineEls(1:rightLink)' 90 lineEls(leftLink:end)'];    
                lineAzs = [lineAzs(1:rightLink)' 0 lineAzs(leftLink:end)'];
            end
        end

        plot(lineAzs.*cos(lineEls/180*pi)+181,-lineEls+91,'-','color',beeCol','linewidth',2);
     end
     
     if plotBinoBorderOnMap
        leftMapResults = mapResultsForBees{i};
        numLeftReg = numRegionsForBees(i);
    
         for j = 1:numLeftReg
            tempSubs = leftMapResults{j,1};            
            [~,indsToClip] = intersect(tempSubs,sDat.fullBeeData(i).inFOVBorderSubs,'rows');
            tempSubs = [circularSmooth(tempSubs(1:end-1,1)',smoothRngFOV), circularSmooth(tempSubs(1:end-1,2)',smoothRngFOV), circularSmooth(tempSubs(1:end-1,3)',smoothRngFOV)];
            tempSubs = vertcat(tempSubs,tempSubs(1,:));
            lineEls = acos(tempSubs(:,2))/pi*180-90;
            lineAzs = -atan2(tempSubs(:,1),-tempSubs(:,3))/pi*180;
            lineAzs(indsToClip) = NaN;
            lineEls(indsToClip) = NaN;
            
            plot(lineAzs.*cos(lineEls/180*pi)+181,-lineEls+91,':','color','w','linewidth',2);
         end
     end
     %draw some lines - ignore bits that cross where visual field is
     for j = [-120 -60 0 60 120]
         lineAzs = (j*ones(1,181)).*cos((-90:90)/180*pi)+181;
         lineEls = -(-90:90)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end
     
     for j = [-90 -60 -30 0 30 60 90]
         lineAzs = (-180:180).*cos(j/180*pi)+181;
         lineEls = -j*ones(361,1)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end

     %do +-180 az bounds seperatley so we never impinge on fov
     plot(-180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
     plot(180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
end

% plot color map for facet diameter
ColMapToUse = DColMap;
maxRange2Use = maxRangeD;
minRange2Use = minRangeD;

facetDCMapF = figure;
colormap(ColMapToUse);
hCB = colorbar;
colSteps = (([0 25 50 75 100]/100*(maxRange2Use-minRange2Use) + minRange2Use));
colLabs = {sprintf('%.1f', colSteps(1)), sprintf('%.1f', colSteps(2)), sprintf('%.1f', colSteps(3)),sprintf('%.1f', colSteps(4)),sprintf('%.1f', colSteps(5))};
set(hCB, 'Ticks', [0 25 50 75 100]/100, 'TickLabels',colLabs,'YAxisLocation','right','TickDirection','out');
title('Lens diameter (um)')
facetDIndFs = [];
for i = toPlotMaps
    if regGroup(i) == 1
        beeCol = 'k';
    else
        beeCol = volCols(round((eyeVols(i)-volRng(1))./(volRng(2)-volRng(1))*(size(volCols,1)-1))+1,:);
    end 
    
    facetDIndFs(end+1) = figure;
    tempIm = gridAllocImage*0;
    tempNum = gridAllocImage*0;
     for j = 1:length(uniqueFacetAllocs)   
        if uniqueFacetAllocs(j) ~= 0
            linkInd = find(sDat.fullBeeData(i).inFOV == uniqueFacetAllocs(j));
            if ~isempty(linkInd)
                %average overlap
                inds = gridAllocImageSineProj{j};
                tempIm(inds) = tempIm(inds) + sDat.fullBeeData(i).diameterOnSphere(linkInd);
                tempNum(inds) = tempNum(inds) + 1;
            end
        end
     end
     tempIm = tempIm./tempNum; %correct for multiple additions
     tempIm = round((tempIm-minRange2Use)./(maxRange2Use-minRange2Use)*99+1);
     tempIm(tempIm < 1 & ~isnan(tempIm) & ~isinf(tempIm)) = 1; tempIm(tempIm > 100 & ~isnan(tempIm) & ~isinf(tempIm)) = 100;
     tempIm(noInds) = Inf; %point in sphere with not image will go to nan from div by zero
     
     colIm = zeros(size(tempIm,1), size(tempIm,2),3);
     for j = 1:size(tempIm,1)
         for k = 1:size(tempIm,2)
             if isnan(tempIm(j,k))
                colIm(j,k,:) = [0.25 0.25 0.25]; %[255 153 153]/255;
             elseif isinf(tempIm(j,k))
                colIm(j,k,:) = 1;
             else    
                 colIm(j,k,:) = ColMapToUse(tempIm(j,k),:);
             end
         end
     end
     
     imshow_w_by_h(colIm); hold on
     
     if plotFOVBorderOnMap
         %plot the fov border
         tempSubs = [circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,1)',smoothRngFOV), circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,2)',smoothRngFOV), circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,3)',smoothRngFOV)];
        tempSubs = vertcat(tempSubs,tempSubs(1,:));
        lineEls = acos(tempSubs(:,2))/pi*180-90;
        lineAzs = -atan2(tempSubs(:,1),-tempSubs(:,3))/pi*180;

        %check if line crosses top of pole and insert 90o point if so
        if sum(sDat.fullBeeData(i).inFOV == topInd)
            %if so find point on left of linkage
            tempInds = find(lineAzs < 0);
            [~,leftLink] = max(lineEls(tempInds));
            leftLink = tempInds(leftLink);
            %and then on right
            tempInds = find(lineAzs >= 0);
            [~,rightLink] = max(lineEls(tempInds));
            rightLink = tempInds(rightLink);

            if rightLink > leftLink
                lineEls = [lineEls(1:leftLink)' 90 lineEls(rightLink:end)'];    
                lineAzs = [lineAzs(1:leftLink)' 0 lineAzs(rightLink:end)'];
            else
                lineEls = [lineEls(1:rightLink)' 90 lineEls(leftLink:end)'];    
                lineAzs = [lineAzs(1:rightLink)' 0 lineAzs(leftLink:end)'];
            end
        end

        plot(lineAzs.*cos(lineEls/180*pi)+181,-lineEls+91,'-','color',beeCol','linewidth',2);
     end
     
     if plotBinoBorderOnMap
        leftMapResults = mapResultsForBees{i};
        numLeftReg = numRegionsForBees(i);
    
         for j = 1:numLeftReg
            tempSubs = leftMapResults{j,1};            
            [~,indsToClip] = intersect(tempSubs,sDat.fullBeeData(i).inFOVBorderSubs,'rows');
            tempSubs = [circularSmooth(tempSubs(1:end-1,1)',smoothRngFOV), circularSmooth(tempSubs(1:end-1,2)',smoothRngFOV), circularSmooth(tempSubs(1:end-1,3)',smoothRngFOV)];
            tempSubs = vertcat(tempSubs,tempSubs(1,:));
            lineEls = acos(tempSubs(:,2))/pi*180-90;
            lineAzs = -atan2(tempSubs(:,1),-tempSubs(:,3))/pi*180;
            lineAzs(indsToClip) = NaN;
            lineEls(indsToClip) = NaN;
            
            plot(lineAzs.*cos(lineEls/180*pi)+181,-lineEls+91,':','color','w','linewidth',2);
         end
     end
     %draw some lines - ignore bits that cross where visual field is
     %and text - will have to move text later
     for j = [-120 -60 0 60 120]
         lineAzs = (j*ones(1,181)).*cos((-90:90)/180*pi)+181;
         lineEls = -(-90:90)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end

     for j = [-90 -60 -30 0 30 60 90]
         lineAzs = (-180:180).*cos(j/180*pi)+181;
         lineEls = -j*ones(361,1)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end
     %do +-180 az bounds seperatley so we never impinge on fov
     plot(-180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
     plot(180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
end

% plot color map of rhabdom thickness
ColMapToUse = RLColMap;
maxRange2Use = maxRangeRL;
minRange2Use = minRangeRL;

rhabLengthCMap = figure;
colormap(ColMapToUse);
hCB = colorbar;
colSteps = (([0 25 50 75 100]/100*(maxRange2Use-minRange2Use) + minRange2Use));
colLabs = {sprintf('%.1f', colSteps(1)), sprintf('%.1f', colSteps(2)), sprintf('%.1f', colSteps(3)),sprintf('%.1f', colSteps(4)),sprintf('%.1f', colSteps(5))};
set(hCB, 'Ticks', [0 25 50 75 100]/100, 'TickLabels',colLabs,'YAxisLocation','right','TickDirection','out');
title('Rhabdom length (um)')
rhabLIndFs = [];

for i = toPlotMaps
    if regGroup(i) == 1
        beeCol = 'k';
    else
        beeCol = volCols(round((eyeVols(i)-volRng(1))./(volRng(2)-volRng(1))*(size(volCols,1)-1))+1,:);
    end 
    
    rhabLIndFs(end+1) = figure;
    tempIm = gridAllocImage*0;
    tempNum = gridAllocImage*0;
     for j = 1:length(uniqueFacetAllocs)   
        if uniqueFacetAllocs(j) ~= 0
            linkInd = find(sDat.fullBeeData(i).inFOV == uniqueFacetAllocs(j));
            if ~isempty(linkInd)
                %average overlap
                inds = gridAllocImageSineProj{j};
                tempIm(inds) = tempIm(inds) + sDat.fullBeeData(i).retTOnSphere(linkInd);
                tempNum(inds) = tempNum(inds) + 1;
            end
        end
     end
     tempIm = tempIm./tempNum; %correct for multiple additions
     tempIm = round((tempIm-minRange2Use)./(maxRange2Use-minRange2Use)*99+1);
     tempIm(tempIm < 1 & ~isnan(tempIm) & ~isinf(tempIm)) = 1; tempIm(tempIm > 100 & ~isnan(tempIm) & ~isinf(tempIm)) = 100;
     tempIm(noInds) = Inf; %point in sphere with not image will go to nan from div by zero
     
     colIm = zeros(size(tempIm,1), size(tempIm,2),3);
     for j = 1:size(tempIm,1)
         for k = 1:size(tempIm,2)
             if isnan(tempIm(j,k))
                colIm(j,k,:) = [0.25 0.25 0.25]; %[255 153 153]/255;
             elseif isinf(tempIm(j,k))
                colIm(j,k,:) = 1;
             else    
                 colIm(j,k,:) = ColMapToUse(tempIm(j,k),:);
             end
         end
     end
     
     %plot border of fov as well.
     imshow_w_by_h(colIm); hold on
     
      if plotFOVBorderOnMap
         %plot the fov border
         tempSubs = [circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,1)',smoothRngFOV), circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,2)',smoothRngFOV), circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,3)',smoothRngFOV)];
        tempSubs = vertcat(tempSubs,tempSubs(1,:));
        lineEls = acos(tempSubs(:,2))/pi*180-90;
        lineAzs = -atan2(tempSubs(:,1),-tempSubs(:,3))/pi*180;

        %check if line crosses top of pole and insert 90o point if so
        if sum(sDat.fullBeeData(i).inFOV == topInd)
            %if so find point on left of linkage
            tempInds = find(lineAzs < 0);
            [~,leftLink] = max(lineEls(tempInds));
            leftLink = tempInds(leftLink);
            %and then on right
            tempInds = find(lineAzs >= 0);
            [~,rightLink] = max(lineEls(tempInds));
            rightLink = tempInds(rightLink);

            if rightLink > leftLink
                lineEls = [lineEls(1:leftLink)' 90 lineEls(rightLink:end)'];    
                lineAzs = [lineAzs(1:leftLink)' 0 lineAzs(rightLink:end)'];
            else
                lineEls = [lineEls(1:rightLink)' 90 lineEls(leftLink:end)'];    
                lineAzs = [lineAzs(1:rightLink)' 0 lineAzs(leftLink:end)'];
            end
        end

        plot(lineAzs.*cos(lineEls/180*pi)+181,-lineEls+91,'-','color',beeCol','linewidth',2);
      end
     
     if plotBinoBorderOnMap
            leftMapResults = mapResultsForBees{i};
            numLeftReg = numRegionsForBees(i);
    
          for j = 1:numLeftReg
            tempSubs = leftMapResults{j,1};            
            [~,indsToClip] = intersect(tempSubs,sDat.fullBeeData(i).inFOVBorderSubs,'rows');
            tempSubs = [circularSmooth(tempSubs(1:end-1,1)',smoothRngFOV), circularSmooth(tempSubs(1:end-1,2)',smoothRngFOV), circularSmooth(tempSubs(1:end-1,3)',smoothRngFOV)];
            tempSubs = vertcat(tempSubs,tempSubs(1,:));
            lineEls = acos(tempSubs(:,2))/pi*180-90;
            lineAzs = -atan2(tempSubs(:,1),-tempSubs(:,3))/pi*180;
            lineAzs(indsToClip) = NaN;
            lineEls(indsToClip) = NaN;

            plot(lineAzs.*cos(lineEls/180*pi)+181,-lineEls+91,':','color','w','linewidth',2);
         end
     end
     %draw some lines - ignore bits that cross where visual field is
     for j = [-120 -60 0 60 120]
         lineAzs = (j*ones(1,181)).*cos((-90:90)/180*pi)+181;
         lineEls = -(-90:90)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end
   
     for j = [-90 -60 -30 0 30 60 90]
         lineAzs = (-180:180).*cos(j/180*pi)+181;
         lineEls = -j*ones(361,1)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end
     %do +-180 az bounds seperatley so we never impinge on fov
     plot(-180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
     plot(180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
end

% plot color map of eye parameter
ColMapToUse = EyeParaColMap;
maxRange2Use = maxRangeEyePara;
minRange2Use = minRangeEyePara;

eyeParaCMapF = figure;
colormap(ColMapToUse);
hCB = colorbar;
colSteps = (([0 25 50 75 100]/100*(maxRange2Use-minRange2Use) + minRange2Use));
colLabs = {sprintf('%.1f', colSteps(1)), sprintf('%.1f', colSteps(2)), sprintf('%.1f', colSteps(3)),sprintf('%.1f', colSteps(4)),sprintf('%.1f', colSteps(5))};
set(hCB, 'Ticks', [0 25 50 75 100]/100, 'TickLabels',colLabs,'YAxisLocation','right','TickDirection','out');
title('Eye parameter (um.rad)');
eyeParaIndFs = [];

for i = toPlotMaps
    if regGroup(i) == 1
        beeCol = 'k';
    else
        beeCol = volCols(round((eyeVols(i)-volRng(1))./(volRng(2)-volRng(1))*(size(volCols,1)-1))+1,:);
    end 
    
    eyeParaIndFs(end+1) = figure;
    
    tempIm = gridAllocImage*0;
    tempNum = gridAllocImage*0;
     for j = 1:length(uniqueFacetAllocs)   
        if uniqueFacetAllocs(j) ~= 0
            linkInd = find(sDat.fullBeeData(i).inFOV == uniqueFacetAllocs(j));
            if ~isempty(linkInd)
                %average overlap
                inds = gridAllocImageSineProj{j};
                tempIm(inds) = tempIm(inds) + sDat.fullBeeData(i).InterOOnSphere(linkInd).*sDat.fullBeeData(i).diameterOnSphere(linkInd);
                tempNum(inds) = tempNum(inds) + 1;
            end
        end
     end
     tempIm = tempIm./tempNum; %correct for multiple additions
     tempIm = round((tempIm-minRange2Use)./(maxRange2Use-minRange2Use)*99+1);
     tempIm(tempIm < 1 & ~isnan(tempIm) & ~isinf(tempIm)) = 1; tempIm(tempIm > 100 & ~isnan(tempIm) & ~isinf(tempIm)) = 100;
     tempIm(noInds) = Inf; %point in sphere with not image will go to nan from div by zero
     
     colIm = zeros(size(tempIm,1), size(tempIm,2),3);
     for j = 1:size(tempIm,1)
         for k = 1:size(tempIm,2)
             if isnan(tempIm(j,k))
                colIm(j,k,:) = [0.25 0.25 0.25]; %[255 153 153]/255;
             elseif isinf(tempIm(j,k))
                colIm(j,k,:) = 1;
             else    
                 colIm(j,k,:) = ColMapToUse(tempIm(j,k),:);
             end
         end
     end
     
     %plot border of fov as well.
     imshow_w_by_h(colIm); hold on
     
     if plotFOVBorderOnMap
         %plot the fov border
         tempSubs = [circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,1)',smoothRngFOV), circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,2)',smoothRngFOV), circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,3)',smoothRngFOV)];
        tempSubs = vertcat(tempSubs,tempSubs(1,:));
        lineEls = acos(tempSubs(:,2))/pi*180-90;
        lineAzs = -atan2(tempSubs(:,1),-tempSubs(:,3))/pi*180;

        %check if line crosses top of pole and insert 90o point if so
        if sum(sDat.fullBeeData(i).inFOV == topInd)
            %if so find point on left of linkage
            tempInds = find(lineAzs < 0);
            [~,leftLink] = max(lineEls(tempInds));
            leftLink = tempInds(leftLink);
            %and then on right
            tempInds = find(lineAzs >= 0);
            [~,rightLink] = max(lineEls(tempInds));
            rightLink = tempInds(rightLink);

            if rightLink > leftLink
                lineEls = [lineEls(1:leftLink)' 90 lineEls(rightLink:end)'];    
                lineAzs = [lineAzs(1:leftLink)' 0 lineAzs(rightLink:end)'];
            else
                lineEls = [lineEls(1:rightLink)' 90 lineEls(leftLink:end)'];    
                lineAzs = [lineAzs(1:rightLink)' 0 lineAzs(leftLink:end)'];
            end
        end

        plot(lineAzs.*cos(lineEls/180*pi)+181,-lineEls+91,'-','color',beeCol','linewidth',2);
     end
    if plotBinoBorderOnMap
            leftMapResults = mapResultsForBees{i};
            numLeftReg = numRegionsForBees(i);
     
         for j = 1:numLeftReg
                tempSubs = leftMapResults{j,1};            
                [~,indsToClip] = intersect(tempSubs,sDat.fullBeeData(i).inFOVBorderSubs,'rows');
                tempSubs = [circularSmooth(tempSubs(1:end-1,1)',smoothRngFOV), circularSmooth(tempSubs(1:end-1,2)',smoothRngFOV), circularSmooth(tempSubs(1:end-1,3)',smoothRngFOV)];
                tempSubs = vertcat(tempSubs,tempSubs(1,:));
                lineEls = acos(tempSubs(:,2))/pi*180-90;
                lineAzs = -atan2(tempSubs(:,1),-tempSubs(:,3))/pi*180;
                lineAzs(indsToClip) = NaN;
                lineEls(indsToClip) = NaN;

                plot(lineAzs.*cos(lineEls/180*pi)+181,-lineEls+91,':','color','w','linewidth',2);
         end
    end 
     %draw some lines - ignore bits that cross where visual field is
     for j = [-120 -60 0 60 120]
         lineAzs = (j*ones(1,181)).*cos((-90:90)/180*pi)+181;
         lineEls = -(-90:90)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end
     
     for j = [-90 -60 -30 0 30 60 90]
         lineAzs = (-180:180).*cos(j/180*pi)+181;
         lineEls = -j*ones(361,1)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end

     %do +-180 az bounds seperatley so we never impinge on fov
     plot(-180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
     plot(180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
end

% plot color map of lens thickness
ColMapToUse = LensTColMap;
maxRange2Use = maxRangeLensT;
minRange2Use = minRangeLensT;

lensTCMapF = figure;
colormap(ColMapToUse);
hCB = colorbar;
colSteps = (([0 25 50 75 100]/100*(maxRange2Use-minRange2Use) + minRange2Use));
colLabs = {sprintf('%.1f', colSteps(1)), sprintf('%.1f', colSteps(2)), sprintf('%.1f', colSteps(3)),sprintf('%.1f', colSteps(4)),sprintf('%.1f', colSteps(5))};
set(hCB, 'Ticks', [0 25 50 75 100]/100, 'TickLabels',colLabs,'YAxisLocation','right','TickDirection','out');
title('Lens Thickness (um)');
lensTIndFs = [];

for i = toPlotMaps
    if regGroup(i) == 1
        beeCol = 'k';
    else
        beeCol = volCols(round((eyeVols(i)-volRng(1))./(volRng(2)-volRng(1))*(size(volCols,1)-1))+1,:);
    end 
    
    lensTIndFs(end+1) = figure;
    
    tempIm = gridAllocImage*0;
    tempNum = gridAllocImage*0;
     for j = 1:length(uniqueFacetAllocs)   
        if uniqueFacetAllocs(j) ~= 0
            linkInd = find(sDat.fullBeeData(i).inFOV == uniqueFacetAllocs(j));
            if ~isempty(linkInd)
                %average overlap
                inds = gridAllocImageSineProj{j};
                tempIm(inds) = tempIm(inds) + sDat.fullBeeData(i).lensTOnSphere(linkInd);
                tempNum(inds) = tempNum(inds) + 1;
            end
        end
     end
     tempIm = tempIm./tempNum; %correct for multiple additions
     tempIm = round((tempIm-minRange2Use)./(maxRange2Use-minRange2Use)*99+1);
     tempIm(tempIm < 1 & ~isnan(tempIm) & ~isinf(tempIm)) = 1; tempIm(tempIm > 100 & ~isnan(tempIm) & ~isinf(tempIm)) = 100;
     tempIm(noInds) = Inf; %point in sphere with not image will go to nan from div by zero
     
     colIm = zeros(size(tempIm,1), size(tempIm,2),3);
     for j = 1:size(tempIm,1)
         for k = 1:size(tempIm,2)
             if isnan(tempIm(j,k))
                colIm(j,k,:) = [0.25 0.25 0.25]; %[255 153 153]/255;
             elseif isinf(tempIm(j,k))
                colIm(j,k,:) = 1;
             else    
                 colIm(j,k,:) = ColMapToUse(tempIm(j,k),:);
             end
         end
     end
     
     %plot border of fov as well.
     imshow_w_by_h(colIm); hold on
     
     if plotFOVBorderOnMap
         %plot the fov border
         tempSubs = [circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,1)',smoothRngFOV), circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,2)',smoothRngFOV), circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,3)',smoothRngFOV)];
        tempSubs = vertcat(tempSubs,tempSubs(1,:));
        lineEls = acos(tempSubs(:,2))/pi*180-90;
        lineAzs = -atan2(tempSubs(:,1),-tempSubs(:,3))/pi*180;

        %check if line crosses top of pole and insert 90o point if so
        if sum(sDat.fullBeeData(i).inFOV == topInd)
            %if so find point on left of linkage
            tempInds = find(lineAzs < 0);
            [~,leftLink] = max(lineEls(tempInds));
            leftLink = tempInds(leftLink);
            %and then on right
            tempInds = find(lineAzs >= 0);
            [~,rightLink] = max(lineEls(tempInds));
            rightLink = tempInds(rightLink);

            if rightLink > leftLink
                lineEls = [lineEls(1:leftLink)' 90 lineEls(rightLink:end)'];    
                lineAzs = [lineAzs(1:leftLink)' 0 lineAzs(rightLink:end)'];
            else
                lineEls = [lineEls(1:rightLink)' 90 lineEls(leftLink:end)'];    
                lineAzs = [lineAzs(1:rightLink)' 0 lineAzs(leftLink:end)'];
            end
        end

        plot(lineAzs.*cos(lineEls/180*pi)+181,-lineEls+91,'-','color',beeCol','linewidth',2);
     end
     if plotBinoBorderOnMap
            leftMapResults = mapResultsForBees{i};
            numLeftReg = numRegionsForBees(i);
    
         for j = 1:numLeftReg
                tempSubs = leftMapResults{j,1};            
                [~,indsToClip] = intersect(tempSubs,sDat.fullBeeData(i).inFOVBorderSubs,'rows');
                tempSubs = [circularSmooth(tempSubs(1:end-1,1)',smoothRngFOV), circularSmooth(tempSubs(1:end-1,2)',smoothRngFOV), circularSmooth(tempSubs(1:end-1,3)',smoothRngFOV)];
                tempSubs = vertcat(tempSubs,tempSubs(1,:));
                lineEls = acos(tempSubs(:,2))/pi*180-90;
                lineAzs = -atan2(tempSubs(:,1),-tempSubs(:,3))/pi*180;
                lineAzs(indsToClip) = NaN;
                lineEls(indsToClip) = NaN;

                plot(lineAzs.*cos(lineEls/180*pi)+181,-lineEls+91,':','color','w','linewidth',2);
         end
     end
     %draw some lines - ignore bits that cross where visual field is
     for j = [-120 -60 0 60 120]
         lineAzs = (j*ones(1,181)).*cos((-90:90)/180*pi)+181;
         lineEls = -(-90:90)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end
     
     for j = [-90 -60 -30 0 30 60 90]
         lineAzs = (-180:180).*cos(j/180*pi)+181;
         lineEls = -j*ones(361,1)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end

     %do +-180 az bounds seperatley so we never impinge on fov
     plot(-180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
     plot(180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
end

% plot CC thicknes on col maps
ColMapToUse = CCTColMap;
maxRange2Use = maxRangeCCT;
minRange2Use = minRangeCCT;

CCTCMapF = figure;
colormap(ColMapToUse);
hCB = colorbar;
colSteps = (([0 25 50 75 100]/100*(maxRange2Use-minRange2Use) + minRange2Use));
colLabs = {sprintf('%.1f', colSteps(1)), sprintf('%.1f', colSteps(2)), sprintf('%.1f', colSteps(3)),sprintf('%.1f', colSteps(4)),sprintf('%.1f', colSteps(5))};
set(hCB, 'Ticks', [0 25 50 75 100]/100, 'TickLabels',colLabs,'YAxisLocation','right','TickDirection','out');
title('CC Thickness (um)');
CCIndFs = [];

for i = toPlotMaps
    if regGroup(i) == 1
        beeCol = 'k';
    else
        beeCol = volCols(round((eyeVols(i)-volRng(1))./(volRng(2)-volRng(1))*(size(volCols,1)-1))+1,:);
    end 
    
    CCIndFs(end+1) = figure;
    
    tempIm = gridAllocImage*0;
    tempNum = gridAllocImage*0;
     for j = 1:length(uniqueFacetAllocs)   
        if uniqueFacetAllocs(j) ~= 0
            linkInd = find(sDat.fullBeeData(i).inFOV == uniqueFacetAllocs(j));
            if ~isempty(linkInd)
                %average overlap
                inds = gridAllocImageSineProj{j};
                tempIm(inds) = tempIm(inds) + sDat.fullBeeData(i).CCTOnSphere(linkInd);
                tempNum(inds) = tempNum(inds) + 1;
            end
        end
     end
     tempIm = tempIm./tempNum; %correct for multiple additions
     tempIm = round((tempIm-minRange2Use)./(maxRange2Use-minRange2Use)*99+1);
     tempIm(tempIm < 1 & ~isnan(tempIm) & ~isinf(tempIm)) = 1; tempIm(tempIm > 100 & ~isnan(tempIm) & ~isinf(tempIm)) = 100;
     tempIm(noInds) = Inf; %point in sphere with not image will go to nan from div by zero
     
     colIm = zeros(size(tempIm,1), size(tempIm,2),3);
     for j = 1:size(tempIm,1)
         for k = 1:size(tempIm,2)
             if isnan(tempIm(j,k))
                colIm(j,k,:) = [0.25 0.25 0.25]; %[255 153 153]/255;
             elseif isinf(tempIm(j,k))
                colIm(j,k,:) = 1;
             else    
                 colIm(j,k,:) = ColMapToUse(tempIm(j,k),:);
             end
         end
     end
     
     %plot border of fov as well.
     imshow_w_by_h(colIm); hold on
     
     if plotFOVBorderOnMap
         %plot the fov border
         tempSubs = [circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,1)',smoothRngFOV), circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,2)',smoothRngFOV), circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,3)',smoothRngFOV)];
        tempSubs = vertcat(tempSubs,tempSubs(1,:));
        lineEls = acos(tempSubs(:,2))/pi*180-90;
        lineAzs = -atan2(tempSubs(:,1),-tempSubs(:,3))/pi*180;

        %check if line crosses top of pole and insert 90o point if so
        if sum(sDat.fullBeeData(i).inFOV == topInd)
            %if so find point on left of linkage
            tempInds = find(lineAzs < 0);
            [~,leftLink] = max(lineEls(tempInds));
            leftLink = tempInds(leftLink);
            %and then on right
            tempInds = find(lineAzs >= 0);
            [~,rightLink] = max(lineEls(tempInds));
            rightLink = tempInds(rightLink);

            if rightLink > leftLink
                lineEls = [lineEls(1:leftLink)' 90 lineEls(rightLink:end)'];    
                lineAzs = [lineAzs(1:leftLink)' 0 lineAzs(rightLink:end)'];
            else
                lineEls = [lineEls(1:rightLink)' 90 lineEls(leftLink:end)'];    
                lineAzs = [lineAzs(1:rightLink)' 0 lineAzs(leftLink:end)'];
            end
        end

        plot(lineAzs.*cos(lineEls/180*pi)+181,-lineEls+91,'-','color',beeCol','linewidth',2);
     end
     if plotBinoBorderOnMap
            leftMapResults = mapResultsForBees{i};
            numLeftReg = numRegionsForBees(i);
    
         for j = 1:numLeftReg
                tempSubs = leftMapResults{j,1};            
                [~,indsToClip] = intersect(tempSubs,sDat.fullBeeData(i).inFOVBorderSubs,'rows');
                tempSubs = [circularSmooth(tempSubs(1:end-1,1)',smoothRngFOV), circularSmooth(tempSubs(1:end-1,2)',smoothRngFOV), circularSmooth(tempSubs(1:end-1,3)',smoothRngFOV)];
                tempSubs = vertcat(tempSubs,tempSubs(1,:));
                lineEls = acos(tempSubs(:,2))/pi*180-90;
                lineAzs = -atan2(tempSubs(:,1),-tempSubs(:,3))/pi*180;
                lineAzs(indsToClip) = NaN;
                lineEls(indsToClip) = NaN;

                plot(lineAzs.*cos(lineEls/180*pi)+181,-lineEls+91,':','color','w','linewidth',2);
         end
     end
     %draw some lines - ignore bits that cross where visual field is
     for j = [-120 -60 0 60 120]
         lineAzs = (j*ones(1,181)).*cos((-90:90)/180*pi)+181;
         lineEls = -(-90:90)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end
     
     for j = [-90 -60 -30 0 30 60 90]
         lineAzs = (-180:180).*cos(j/180*pi)+181;
         lineEls = -j*ones(361,1)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end

     %do +-180 az bounds seperatley so we never impinge on fov
     plot(-180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
     plot(180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
end

% plot curvature on col maps
ColMapToUse = CurvColMap;
maxRange2Use = maxRangeCurv;
minRange2Use = minRangeCurv;

CurvCMapF = figure;
colormap(ColMapToUse);
hCB = colorbar;
colSteps = (([0 25 50 75 100]/100*(maxRange2Use-minRange2Use) + minRange2Use));
colLabs = {sprintf('%.1f', colSteps(1)), sprintf('%.1f', colSteps(2)), sprintf('%.1f', colSteps(3)),sprintf('%.1f', colSteps(4)),sprintf('%.1f', colSteps(5))};
set(hCB, 'Ticks', [0 25 50 75 100]/100, 'TickLabels',colLabs,'YAxisLocation','right','TickDirection','out');
title('Curvature (um)');
CurvIndFs = [];

for i = toPlotMaps
    if regGroup(i) == 1
        beeCol = 'k';
    else
        beeCol = volCols(round((eyeVols(i)-volRng(1))./(volRng(2)-volRng(1))*(size(volCols,1)-1))+1,:);
    end 
    
    CurvIndFs(end+1) = figure;
    
    tempIm = gridAllocImage*0;
    tempNum = gridAllocImage*0;
     for j = 1:length(uniqueFacetAllocs)   
        if uniqueFacetAllocs(j) ~= 0
            linkInd = find(sDat.fullBeeData(i).inFOV == uniqueFacetAllocs(j));
            if ~isempty(linkInd)
                %average overlap
                inds = gridAllocImageSineProj{j};
                tempIm(inds) = tempIm(inds) + sDat.fullBeeData(i).curvatureOnSphere(linkInd);
                tempNum(inds) = tempNum(inds) + 1;
            end
        end
     end
     tempIm = tempIm./tempNum; %correct for multiple additions
     tempIm = round((tempIm-minRange2Use)./(maxRange2Use-minRange2Use)*99+1);
     tempIm(tempIm < 1 & ~isnan(tempIm) & ~isinf(tempIm)) = 1; tempIm(tempIm > 100 & ~isnan(tempIm) & ~isinf(tempIm)) = 100;
     tempIm(noInds) = Inf; %point in sphere with not image will go to nan from div by zero
     
     colIm = zeros(size(tempIm,1), size(tempIm,2),3);
     for j = 1:size(tempIm,1)
         for k = 1:size(tempIm,2)
             if isnan(tempIm(j,k))
                colIm(j,k,:) = [0.25 0.25 0.25]; %[255 153 153]/255;
             elseif isinf(tempIm(j,k))
                colIm(j,k,:) = 1;
             else    
                 colIm(j,k,:) = ColMapToUse(tempIm(j,k),:);
             end
         end
     end
     
     %plot border of fov as well.
     imshow_w_by_h(colIm); hold on
     
     if plotFOVBorderOnMap
         %plot the fov border
         tempSubs = [circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,1)',smoothRngFOV), circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,2)',smoothRngFOV), circularSmooth(sDat.fullBeeData(i).inFOVBorderSubs(1:end-1,3)',smoothRngFOV)];
        tempSubs = vertcat(tempSubs,tempSubs(1,:));
        lineEls = acos(tempSubs(:,2))/pi*180-90;
        lineAzs = -atan2(tempSubs(:,1),-tempSubs(:,3))/pi*180;

        %check if line crosses top of pole and insert 90o point if so
        if sum(sDat.fullBeeData(i).inFOV == topInd)
            %if so find point on left of linkage
            tempInds = find(lineAzs < 0);
            [~,leftLink] = max(lineEls(tempInds));
            leftLink = tempInds(leftLink);
            %and then on right
            tempInds = find(lineAzs >= 0);
            [~,rightLink] = max(lineEls(tempInds));
            rightLink = tempInds(rightLink);

            if rightLink > leftLink
                lineEls = [lineEls(1:leftLink)' 90 lineEls(rightLink:end)'];    
                lineAzs = [lineAzs(1:leftLink)' 0 lineAzs(rightLink:end)'];
            else
                lineEls = [lineEls(1:rightLink)' 90 lineEls(leftLink:end)'];    
                lineAzs = [lineAzs(1:rightLink)' 0 lineAzs(leftLink:end)'];
            end
        end

        plot(lineAzs.*cos(lineEls/180*pi)+181,-lineEls+91,'-','color',beeCol','linewidth',2);
     end
     
     if plotBinoBorderOnMap
        leftMapResults = mapResultsForBees{i};
        numLeftReg = numRegionsForBees(i);
    
         for j = 1:numLeftReg
            tempSubs = leftMapResults{j,1};
            tempSubs = [circularSmooth(tempSubs(1:end-1,1)',smoothRngFOV), circularSmooth(tempSubs(1:end-1,2)',smoothRngFOV), circularSmooth(tempSubs(1:end-1,3)',smoothRngFOV)];
            tempSubs = vertcat(tempSubs,tempSubs(1,:));
            lineEls = acos(tempSubs(:,2))/pi*180-90;
            lineAzs = -atan2(tempSubs(:,1),-tempSubs(:,3))/pi*180;

            %check if line crosses top of pole and insert 90o point if so
            if sum(leftMapResults{j,3} == topInd)
                %if so find point on left of linkage
                tempInds = find(lineAzs < 0);
                [~,leftLink] = max(lineEls(tempInds));
                leftLink = tempInds(leftLink);
                %and then on right
                tempInds = find(lineAzs >= 0);
                [~,rightLink] = max(lineEls(tempInds));
                rightLink = tempInds(rightLink);

                if rightLink > leftLink
                    lineEls = [lineEls(1:leftLink)' 90 lineEls(rightLink:end)'];    
                    lineAzs = [lineAzs(1:leftLink)' 0 lineAzs(rightLink:end)'];
                else
                    lineEls = [lineEls(1:rightLink)' 90 lineEls(leftLink:end)'];    
                    lineAzs = [lineAzs(1:rightLink)' 0 lineAzs(leftLink:end)'];
                end
            end
            plot(lineAzs.*cos(lineEls/180*pi)+181,-lineEls+91,':','color',beeCol','linewidth',2);
         end
     end
     %draw some lines - ignore bits that cross where visual field is
     for j = [-120 -60 0 60 120]
         lineAzs = (j*ones(1,181)).*cos((-90:90)/180*pi)+181;
         lineEls = -(-90:90)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end
     
     for j = [-90 -60 -30 0 30 60 90]
         lineAzs = (-180:180).*cos(j/180*pi)+181;
         lineEls = -j*ones(361,1)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end

     %do +-180 az bounds seperatley so we never impinge on fov
     plot(-180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
     plot(180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
end

% plot color leged for volumes
volCMapF = figure;
colormap(volCols);
hCB = colorbar;
colSteps = (([0 25 50 75 100]/100*(volRng(2)-volRng(1)) + volRng(1)));
colLabs = {sprintf('%.1f', colSteps(1)), sprintf('%.1f', colSteps(2)), sprintf('%.1f', colSteps(3)),sprintf('%.1f', colSteps(4)),sprintf('%.1f', colSteps(5))};
set(hCB, 'Ticks', [0 25 50 75 100]/100, 'TickLabels',colLabs,'YAxisLocation','right','TickDirection','out');
title('Volume^0^.^3 (um)')

%% do visual field scaling mapping
indsGr = find(regGroup == 2); %just for bb group

pvalCutOff = 1;
minBeeNum = 4;

%sort by size
tempVols = eyeVols(indsGr);
[~, indS] = sort(tempVols);
indsGr = indsGr(indS);

targetEyeVols = log10(([sDat.fullBeeData(indsGr).EyeVolume]/10^9).^(1/3));

%store resulting values 
%always, 1st scaling exponent, 2nd zero intercept, 3rd r^2 val, 4th pval
eyeParamScaling = zeros(length(uniqueFacetAllocs),4)*NaN;
facetDScaling = zeros(length(uniqueFacetAllocs),4)*NaN;
interOScaling = zeros(length(uniqueFacetAllocs),4)*NaN;
lensTScaling = zeros(length(uniqueFacetAllocs),4)*NaN;
retTScaling = zeros(length(uniqueFacetAllocs),4)*NaN;
CCTScaling = zeros(length(uniqueFacetAllocs),4)*NaN;
curvScaling = zeros(length(uniqueFacetAllocs),4)*NaN;

numberPoints = zeros(length(uniqueFacetAllocs),1);
numDone = 0;

%go across the unique sphere allocations and check for ones that have at least two vals
for i = 1:length(uniqueFacetAllocs) 
    %get points on bee eyes pointing this way
    indsToUse = zeros(length(indsGr),1);
    
    eyePara2Use = zeros(length(indsGr),1)*NaN;
    facetD2Use = zeros(length(indsGr),1)*NaN;
    interO2Use = zeros(length(indsGr),1)*NaN;
    lensT2Use = zeros(length(indsGr),1)*NaN;
    CCT2Use = zeros(length(indsGr),1)*NaN;
    retT2Use = zeros(length(indsGr),1)*NaN;
    curv2Use = zeros(length(indsGr),1)*NaN;
    
    for j = 1:length(indsGr)
        linkInd = find(sDat.fullBeeData(indsGr(j)).inFOV == uniqueFacetAllocs(i));
        if ~isempty(linkInd)
            indsToUse(j) = linkInd; %pretty sure this can only ever have one linkage...
            
            eyePara2Use(j) = log10((sDat.fullBeeData(indsGr(j)).InterOOnSphere(linkInd).*sDat.fullBeeData(indsGr(j)).diameterOnSphere(linkInd)));
            facetD2Use(j) = log10(sDat.fullBeeData(indsGr(j)).diameterOnSphere(linkInd));
            interO2Use(j) = log10(sDat.fullBeeData(indsGr(j)).InterOOnSphere(linkInd));
            lensT2Use(j) = log10(sDat.fullBeeData(indsGr(j)).lensTOnSphere(linkInd));
            retT2Use(j) = log10(sDat.fullBeeData(indsGr(j)).retTOnSphere(linkInd));
            CCT2Use(j) = log10(sDat.fullBeeData(indsGr(j)).CCTOnSphere(linkInd));
            curv2Use(j) = log10(sDat.fullBeeData(indsGr(j)).curvatureOnSphere(linkInd));
        end
    end
    beesToUse = find(indsToUse);
    numberPoints(i) = length(beesToUse);
    
    %need 3 or more points to calcaulte a p val
    %also fitting a line to two points is, well, just a line not a fit right?
    if length(beesToUse) >= minBeeNum       
        tempRes = polyfit(targetEyeVols(beesToUse),eyePara2Use(beesToUse)',1);
        tempRes(2) = 10.^tempRes(2);
        [R,P] = corrcoef(targetEyeVols(beesToUse),eyePara2Use(beesToUse));
        eyeParamScaling(i,:) = [tempRes, R(2)^2, P(2)];
        
        tempRes = polyfit(targetEyeVols(beesToUse),facetD2Use(beesToUse)',1);
        tempRes(2) = 10.^tempRes(2);
        [R,P] = corrcoef(targetEyeVols(beesToUse),facetD2Use(beesToUse));
        facetDScaling(i,:) = [tempRes, R(2)^2, P(2)];
        
        tempRes = polyfit(targetEyeVols(beesToUse),interO2Use(beesToUse)',1);
        tempRes(2) = 10.^tempRes(2);
        [R,P] = corrcoef(targetEyeVols(beesToUse),interO2Use(beesToUse));
        interOScaling(i,:) = [tempRes, R(2)^2, P(2)];
        
        tempRes = polyfit(targetEyeVols(beesToUse),lensT2Use(beesToUse)',1);
        tempRes(2) = 10.^tempRes(2);
        [R,P] = corrcoef(targetEyeVols(beesToUse),lensT2Use(beesToUse));
        lensTScaling(i,:) = [tempRes, R(2)^2, P(2)];
        
        tempRes = polyfit(targetEyeVols(beesToUse),retT2Use(beesToUse)',1);
        tempRes(2) = 10.^tempRes(2);
        [R,P] = corrcoef(targetEyeVols(beesToUse),retT2Use(beesToUse));
        retTScaling(i,:) = [tempRes, R(2)^2, P(2)];
        
        tempRes = polyfit(targetEyeVols(beesToUse),CCT2Use(beesToUse)',1);
        tempRes(2) = 10.^tempRes(2);
        [R,P] = corrcoef(targetEyeVols(beesToUse),CCT2Use(beesToUse));
        CCTScaling(i,:) = [tempRes, R(2)^2, P(2)];
        
        tempRes = polyfit(targetEyeVols(beesToUse),curv2Use(beesToUse)',1);
        tempRes(2) = 10.^tempRes(2);
        [R,P] = corrcoef(targetEyeVols(beesToUse),curv2Use(beesToUse));
        curvScaling(i,:) = [tempRes, R(2)^2, P(2)];
        
        numDone = numDone + 1;
    end
end

min(interOScaling(:,1))

%clean up to delete ones with low P values...
eyeParamScaling(eyeParamScaling(:,4) > pvalCutOff) = NaN;
facetDScaling(facetDScaling(:,4) > pvalCutOff) = NaN;
interOScaling(interOScaling(:,4) > pvalCutOff) = NaN;
lensTScaling(lensTScaling(:,4) > pvalCutOff) = NaN;
retTScaling(retTScaling(:,4) > pvalCutOff) = NaN;
CCTScaling(CCTScaling(:,4) > pvalCutOff) = NaN;
curvScaling(curvScaling(:,4) > pvalCutOff) = NaN;

ColMapToUse = ScaleColMap;
maxRange2Use = maxRangeScale;
minRange2Use = minRangeScale;

scaleMapF = figure;
colormap(ColMapToUse);
hCB = colorbar;
colSteps = (([0 75/4 75/2 75*3/4 75]/size(ScaleColMap,1)*(maxRange2Use-minRange2Use) + minRange2Use));
colLabs = {sprintf('%.1f', colSteps(1)), sprintf('%.1f', colSteps(2)), sprintf('%.1f', colSteps(3)),sprintf('%.1f', colSteps(4)),sprintf('%.1f', colSteps(5))};
set(hCB, 'Ticks', [0 75/4 75/2 75*3/4 75]/size(ScaleColMap,1), 'TickLabels',colLabs,'YAxisLocation','right','TickDirection','out');
title('Scaling exponent');
 
% plot the eye param scaling
eyePScalFig = figure;
tempIm = gridAllocImage*0;
tempNum = gridAllocImage*0;
for i = 1:length(uniqueFacetAllocs)   
    inds = gridAllocImageSineProj{i};
    if ~isnan(eyeParamScaling(i,1))
        %average overlap
        %if was missing, remove
        for j = 1:length(inds)
            if tempIm(inds(j)) == -999
                tempIm(inds(j)) = 0;
            end
        end
        %to average overlap
        tempIm(inds) = tempIm(inds) + eyeParamScaling(i,1);
        tempNum(inds) = tempNum(inds) + 1;
    elseif  numberPoints(i) > 0
        %label missing
        for j = 1:length(inds)
            if tempIm(inds(j)) == 0
                tempIm(inds(j)) = -999;
            end
        end
    end
end
 greyOutInds = find(tempIm == -999);
 tempIm(greyOutInds) = 0;
 tempIm = tempIm./tempNum; %correct for multiple additions
tempIm = round((tempIm-minRange2Use)./(maxRange2Use-minRange2Use)*(size(ColMapToUse,1)-1)+1);
 tempIm(tempIm < 1 & ~isnan(tempIm) & ~isinf(tempIm)) = 1; tempIm(tempIm > size(ColMapToUse,1) & ~isnan(tempIm) & ~isinf(tempIm)) = size(ColMapToUse,1);
 tempIm(noInds) = Inf; %point in sphere with not image will go to nan from div by zero
 tempIm(greyOutInds) = -1;
 
colIm = zeros(size(tempIm,1), size(tempIm,2),3);
 for j = 1:size(tempIm,1)
     for k = 1:size(tempIm,2)
         if isnan(tempIm(j,k))
            colIm(j,k,:) = [0.25 0.25 0.25]; %[255 153 153]/255;
         elseif isinf(tempIm(j,k))
            colIm(j,k,:) = 1;
         elseif tempIm(j,k) < 0
             colIm(j,k,:) = [0.75 0.75 0.75];
         else
             colIm(j,k,:) = ColMapToUse(tempIm(j,k),:);
         end
     end
 end
 
 imshow_w_by_h(colIm); hold on
 %draw some lines - ignore bits that cross where visual field is
     for j = [-120 -60 0 60 120]
         lineAzs = (j*ones(1,181)).*cos((-90:90)/180*pi)+181;
         lineEls = -(-90:90)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end
     
     for j = [-90 -60 -30 0 30 60 90]
         lineAzs = (-180:180).*cos(j/180*pi)+181;
         lineEls = -j*ones(361,1)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end

     %do +-180 az bounds seperatley so we never impinge on fov
     plot(-180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
     plot(180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
set(gcf, 'Position', [50 50 750 375]); axis off

 % plot the interO scaling
interOScalFig = figure;
tempIm = gridAllocImage*0;
tempNum = gridAllocImage*0;
for i = 1:length(uniqueFacetAllocs)   
    inds = gridAllocImageSineProj{i};
    if ~isnan(interOScaling(i,1))
        %if was missing, remove
        for j = 1:length(inds)
            if tempIm(inds(j)) == -999
                tempIm(inds(j)) = 0;
            end
        end
        %to average overlap
        tempIm(inds) = tempIm(inds) + interOScaling(i,1);
        tempNum(inds) = tempNum(inds) + 1;
    elseif  numberPoints(i) > 0
        %label missing
        for j = 1:length(inds)
            if tempIm(inds(j)) == 0
                tempIm(inds(j)) = -999;
            end
        end
    end
end
 greyOutInds = find(tempIm == -999);
 tempIm(greyOutInds) = 0;
 tempIm = tempIm./tempNum; %correct for multiple additions
tempIm = round((tempIm-minRange2Use)./(maxRange2Use-minRange2Use)*(size(ColMapToUse,1)-1)+1);
 tempIm(tempIm < 1 & ~isnan(tempIm) & ~isinf(tempIm)) = 1; tempIm(tempIm > size(ColMapToUse,1) & ~isnan(tempIm) & ~isinf(tempIm)) = size(ColMapToUse,1);
 tempIm(noInds) = Inf; %point in sphere with not image will go to nan from div by zero
 tempIm(greyOutInds) = -1;
 
colIm = zeros(size(tempIm,1), size(tempIm,2),3);
 for j = 1:size(tempIm,1)
     for k = 1:size(tempIm,2)
         if isnan(tempIm(j,k))
            colIm(j,k,:) = [0.25 0.25 0.25]; %[255 153 153]/255;
         elseif isinf(tempIm(j,k))
            colIm(j,k,:) = 1;
         elseif tempIm(j,k) < 0
             colIm(j,k,:) = [0.75 0.75 0.75];
         else
             colIm(j,k,:) = ColMapToUse(tempIm(j,k),:);
         end
     end
 end
 
 imshow_w_by_h(colIm); hold on
 %draw some lines - ignore bits that cross where visual field is
     for j = [-120 -60 0 60 120]
         lineAzs = (j*ones(1,181)).*cos((-90:90)/180*pi)+181;
         lineEls = -(-90:90)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end
     
     for j = [-90 -60 -30 0 30 60 90]
         lineAzs = (-180:180).*cos(j/180*pi)+181;
         lineEls = -j*ones(361,1)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end

     %do +-180 az bounds seperatley so we never impinge on fov
     plot(-180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
     plot(180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
set(gcf, 'Position', [50 50 750 375]); axis off

% plot the facet diameter scaling
facetDScalFig = figure;
tempIm = gridAllocImage*0;
tempNum = gridAllocImage*0;
for i = 1:length(uniqueFacetAllocs)   
    inds = gridAllocImageSineProj{i};
    if ~isnan(facetDScaling(i,1))
        %if was missing, remove
        for j = 1:length(inds)
            if tempIm(inds(j)) == -999
                tempIm(inds(j)) = 0;
            end
        end
        %to average overlap
        tempIm(inds) = tempIm(inds) + facetDScaling(i,1);
        tempNum(inds) = tempNum(inds) + 1;
    elseif  numberPoints(i) > 0
        %label missing
        for j = 1:length(inds)
            if tempIm(inds(j)) == 0
                tempIm(inds(j)) = -999;
            end
        end
    end
end
 greyOutInds = find(tempIm == -999);
 tempIm(greyOutInds) = 0;
 tempIm = tempIm./tempNum; %correct for multiple additions
tempIm = round((tempIm-minRange2Use)./(maxRange2Use-minRange2Use)*(size(ColMapToUse,1)-1)+1);
 tempIm(tempIm < 1 & ~isnan(tempIm) & ~isinf(tempIm)) = 1; tempIm(tempIm > size(ColMapToUse,1) & ~isnan(tempIm) & ~isinf(tempIm)) = size(ColMapToUse,1);
 tempIm(noInds) = Inf; %point in sphere with not image will go to nan from div by zero
 tempIm(greyOutInds) = -1;
 
colIm = zeros(size(tempIm,1), size(tempIm,2),3);
 for j = 1:size(tempIm,1)
     for k = 1:size(tempIm,2)
         if isnan(tempIm(j,k))
            colIm(j,k,:) = [0.25 0.25 0.25]; %[255 153 153]/255;
         elseif isinf(tempIm(j,k))
            colIm(j,k,:) = 1;
         elseif tempIm(j,k) < 0
             colIm(j,k,:) = [0.75 0.75 0.75];
         else
             colIm(j,k,:) = ColMapToUse(tempIm(j,k),:);
         end
     end
 end
 
 imshow_w_by_h(colIm); hold on
 %draw some lines - ignore bits that cross where visual field is
     for j = [-120 -60 0 60 120]
         lineAzs = (j*ones(1,181)).*cos((-90:90)/180*pi)+181;
         lineEls = -(-90:90)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end
     
     for j = [-90 -60 -30 0 30 60 90]
         lineAzs = (-180:180).*cos(j/180*pi)+181;
         lineEls = -j*ones(361,1)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end

     %do +-180 az bounds seperatley so we never impinge on fov
     plot(-180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
     plot(180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
set(gcf, 'Position', [50 50 750 375]); axis off

 % plot the lens Thickness scaling
lensTScalFig = figure;
tempIm = gridAllocImage*0;
tempNum = gridAllocImage*0;
for i = 1:length(uniqueFacetAllocs)   
    inds = gridAllocImageSineProj{i};
    if ~isnan(lensTScaling(i,1))
        %if was missing, remove
        for j = 1:length(inds)
            if tempIm(inds(j)) == -999
                tempIm(inds(j)) = 0;
            end
        end
        %to average overlap
        tempIm(inds) = tempIm(inds) + lensTScaling(i,1);
        tempNum(inds) = tempNum(inds) + 1;
    elseif  numberPoints(i) > 0
        %label missing
        for j = 1:length(inds)
            if tempIm(inds(j)) == 0
                tempIm(inds(j)) = -999;
            end
        end
    end
end
 greyOutInds = find(tempIm == -999);
 tempIm(greyOutInds) = 0;
 tempIm = tempIm./tempNum; %correct for multiple additions
tempIm = round((tempIm-minRange2Use)./(maxRange2Use-minRange2Use)*(size(ColMapToUse,1)-1)+1);
 tempIm(tempIm < 1 & ~isnan(tempIm) & ~isinf(tempIm)) = 1; tempIm(tempIm > size(ColMapToUse,1) & ~isnan(tempIm) & ~isinf(tempIm)) = size(ColMapToUse,1);
 tempIm(noInds) = Inf; %point in sphere with not image will go to nan from div by zero
 tempIm(greyOutInds) = -1;
 
colIm = zeros(size(tempIm,1), size(tempIm,2),3);
 for j = 1:size(tempIm,1)
     for k = 1:size(tempIm,2)
         if isnan(tempIm(j,k))
            colIm(j,k,:) = [0.25 0.25 0.25]; %[255 153 153]/255;
         elseif isinf(tempIm(j,k))
            colIm(j,k,:) = 1;
         elseif tempIm(j,k) < 0
             colIm(j,k,:) = [0.75 0.75 0.75];
         else
             colIm(j,k,:) = ColMapToUse(tempIm(j,k),:);
         end
     end
 end
 
 imshow_w_by_h(colIm); hold on
%draw some lines - ignore bits that cross where visual field is
     for j = [-120 -60 0 60 120]
         lineAzs = (j*ones(1,181)).*cos((-90:90)/180*pi)+181;
         lineEls = -(-90:90)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end
     
     for j = [-90 -60 -30 0 30 60 90]
         lineAzs = (-180:180).*cos(j/180*pi)+181;
         lineEls = -j*ones(361,1)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end

     %do +-180 az bounds seperatley so we never impinge on fov
     plot(-180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
     plot(180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
set(gcf, 'Position', [50 50 750 375]); axis off

% plot the retina thickness scaling
retTScalFig = figure;
tempIm = gridAllocImage*0;
tempNum = gridAllocImage*0;
for i = 1:length(uniqueFacetAllocs)   
    inds = gridAllocImageSineProj{i};
    if ~isnan(retTScaling(i,1))
        %if was missing, remove
        for j = 1:length(inds)
            if tempIm(inds(j)) == -999
                tempIm(inds(j)) = 0;
            end
        end
        %to average overlap
        tempIm(inds) = tempIm(inds) + retTScaling(i,1);
        tempNum(inds) = tempNum(inds) + 1;
    elseif  numberPoints(i) > 0
        %label missing
        for j = 1:length(inds)
            if tempIm(inds(j)) == 0
                tempIm(inds(j)) = -999;
            end
        end
    end
end
 greyOutInds = find(tempIm == -999);
 tempIm(greyOutInds) = 0;
 tempIm = tempIm./tempNum; %correct for multiple additions
tempIm = round((tempIm-minRange2Use)./(maxRange2Use-minRange2Use)*(size(ColMapToUse,1)-1)+1);
 tempIm(tempIm < 1 & ~isnan(tempIm) & ~isinf(tempIm)) = 1; tempIm(tempIm > size(ColMapToUse,1) & ~isnan(tempIm) & ~isinf(tempIm)) = size(ColMapToUse,1);
 tempIm(noInds) = Inf; %point in sphere with not image will go to nan from div by zero
 tempIm(greyOutInds) = -1;
 
colIm = zeros(size(tempIm,1), size(tempIm,2),3);
 for j = 1:size(tempIm,1)
     for k = 1:size(tempIm,2)
         if isnan(tempIm(j,k))
            colIm(j,k,:) = [0.25 0.25 0.25]; %[255 153 153]/255;
         elseif isinf(tempIm(j,k))
            colIm(j,k,:) = 1;
         elseif tempIm(j,k) < 0
             colIm(j,k,:) = [0.75 0.75 0.75];
         else
             colIm(j,k,:) = ColMapToUse(tempIm(j,k),:);
         end
     end
 end
 
 imshow_w_by_h(colIm); hold on
%draw some lines - ignore bits that cross where visual field is
     for j = [-120 -60 0 60 120]
         lineAzs = (j*ones(1,181)).*cos((-90:90)/180*pi)+181;
         lineEls = -(-90:90)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end
     
     for j = [-90 -60 -30 0 30 60 90]
         lineAzs = (-180:180).*cos(j/180*pi)+181;
         lineEls = -j*ones(361,1)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end

     %do +-180 az bounds seperatley so we never impinge on fov
     plot(-180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
     plot(180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
set(gcf, 'Position', [50 50 750 375]); axis off

% plot the CC thickness scaling
CCTScalFig = figure;
tempIm = gridAllocImage*0;
tempNum = gridAllocImage*0;
for i = 1:length(uniqueFacetAllocs)   
    inds = gridAllocImageSineProj{i};
    if ~isnan(CCTScaling(i,1))
        %if was missing, remove
        for j = 1:length(inds)
            if tempIm(inds(j)) == -999
                tempIm(inds(j)) = 0;
            end
        end
        %to average overlap
        tempIm(inds) = tempIm(inds) + CCTScaling(i,1);
        tempNum(inds) = tempNum(inds) + 1;
    elseif  numberPoints(i) > 0
        %label missing
        for j = 1:length(inds)
            if tempIm(inds(j)) == 0
                tempIm(inds(j)) = -999;
            end
        end
    end
end
 greyOutInds = find(tempIm == -999);
 tempIm(greyOutInds) = 0;
 tempIm = tempIm./tempNum; %correct for multiple additions
tempIm = round((tempIm-minRange2Use)./(maxRange2Use-minRange2Use)*(size(ColMapToUse,1)-1)+1);
 tempIm(tempIm < 1 & ~isnan(tempIm) & ~isinf(tempIm)) = 1; tempIm(tempIm > size(ColMapToUse,1) & ~isnan(tempIm) & ~isinf(tempIm)) = size(ColMapToUse,1);
 tempIm(noInds) = Inf; %point in sphere with not image will go to nan from div by zero
 tempIm(greyOutInds) = -1;
 
colIm = zeros(size(tempIm,1), size(tempIm,2),3);
 for j = 1:size(tempIm,1)
     for k = 1:size(tempIm,2)
         if isnan(tempIm(j,k))
            colIm(j,k,:) = [0.25 0.25 0.25]; %[255 153 153]/255;
         elseif isinf(tempIm(j,k))
            colIm(j,k,:) = 1;
         elseif tempIm(j,k) < 0
             colIm(j,k,:) = [0.75 0.75 0.75];
         else
             colIm(j,k,:) = ColMapToUse(tempIm(j,k),:);
         end
     end
 end
 
 imshow_w_by_h(colIm); hold on
%draw some lines - ignore bits that cross where visual field is
     for j = [-120 -60 0 60 120]
         lineAzs = (j*ones(1,181)).*cos((-90:90)/180*pi)+181;
         lineEls = -(-90:90)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end
     
     for j = [-90 -60 -30 0 30 60 90]
         lineAzs = (-180:180).*cos(j/180*pi)+181;
         lineEls = -j*ones(361,1)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end

     %do +-180 az bounds seperatley so we never impinge on fov
     plot(-180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
     plot(180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
set(gcf, 'Position', [50 50 750 375]); axis off

 % plot the curvature scaling
curvScalFig = figure;
tempIm = gridAllocImage*0;
tempNum = gridAllocImage*0;
for i = 1:length(uniqueFacetAllocs)   
     inds = gridAllocImageSineProj{i};
    if ~isnan(curvScaling(i,1))
        %if was missing, remove
        for j = 1:length(inds)
            if tempIm(inds(j)) == -999
                tempIm(inds(j)) = 0;
            end
        end
        %to average overlap
        tempIm(inds) = tempIm(inds) + curvScaling(i,1);
        tempNum(inds) = tempNum(inds) + 1;
    elseif  numberPoints(i) > 0
        %label missing
        for j = 1:length(inds)
            if tempIm(inds(j)) == 0
                tempIm(inds(j)) = -999;
            end
        end
    end
end
 greyOutInds = find(tempIm == -999);
 tempIm(greyOutInds) = 0;
 tempIm = tempIm./tempNum; %correct for multiple additions
tempIm = round((tempIm-minRange2Use)./(maxRange2Use-minRange2Use)*(size(ColMapToUse,1)-1)+1);
 tempIm(tempIm < 1 & ~isnan(tempIm) & ~isinf(tempIm)) = 1; tempIm(tempIm > size(ColMapToUse,1) & ~isnan(tempIm) & ~isinf(tempIm)) = size(ColMapToUse,1);
 tempIm(noInds) = Inf; %point in sphere with not image will go to nan from div by zero
 tempIm(greyOutInds) = -1;
 
colIm = zeros(size(tempIm,1), size(tempIm,2),3);
 for j = 1:size(tempIm,1)
     for k = 1:size(tempIm,2)
         if isnan(tempIm(j,k))
            colIm(j,k,:) = [0.25 0.25 0.25]; %[255 153 153]/255;
         elseif isinf(tempIm(j,k))
            colIm(j,k,:) = 1;
         elseif tempIm(j,k) < 0
             colIm(j,k,:) = [0.75 0.75 0.75];
         else
             colIm(j,k,:) = ColMapToUse(tempIm(j,k),:);
         end
     end
 end
 
 imshow_w_by_h(colIm); hold on
 %draw some lines - ignore bits that cross where visual field is
     for j = [-120 -60 0 60 120]
         lineAzs = (j*ones(1,181)).*cos((-90:90)/180*pi)+181;
         lineEls = -(-90:90)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end
     
     for j = [-90 -60 -30 0 30 60 90]
         lineAzs = (-180:180).*cos(j/180*pi)+181;
         lineEls = -j*ones(361,1)+91;
         for k = 1:length(lineAzs)
             if ~isnan(tempIm(round(lineAzs(k)),round(lineEls(k)))) & ~isinf(tempIm(round(lineAzs(k)),round(lineEls(k))))
                 lineAzs(k) = NaN; lineEls(k) = NaN;
             end
         end
         plot(lineAzs, lineEls, '-','color', [0.75 0.75 0.75],'linewidth',1.5)  
     end

     %do +-180 az bounds seperatley so we never impinge on fov
     plot(-180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
     plot(180.*cos((-90:90)/180*pi)+181, (-90:90)+91, '-','color', [0.75 0.75 0.75],'linewidth',1.5)
set(gcf, 'Position', [50 50 750 375]); axis off

%% calculate correlations and plot
superCoreF = figure; hold on
for i = 1:length(regressionGroups)
    tempX = ([sDat.fullBeeData(regGroup == i).EyeVolume]/10^9).^(1/3);
    inds = find(regGroup == i);
    
    lensDvsIO = zeros(length(inds),1)';
    lensDvsRet = zeros(length(inds),1)';
    lensDvsLens = zeros(length(inds),1)';
    lensDvsCC = zeros(length(inds),1)';
    lensDvsCurv = zeros(length(inds),1)';
    
    IOvsRet = zeros(length(inds),1)';
    IOvsLens = zeros(length(inds),1)';
    IOvsCC = zeros(length(inds),1)';
    IOvsCurv = zeros(length(inds),1)';
    
    RetvsLens = zeros(length(inds),1)';
    RetvsCC = zeros(length(inds),1)';
    RetvsCurv = zeros(length(inds),1)';
    
    LensvsCC = zeros(length(inds),1)';
    LensvsCurv = zeros(length(inds),1)';
    
    CCvsCurv = zeros(length(inds),1)';
    
    for j = 1:length(inds)
        %note d and IO lines swapped!
        temp = weightedcorrs([sDat.fullBeeData(inds(j)).calInterO,sDat.fullBeeData(inds(j)).avgDiameter], facetWeightInds{inds(j)});
        lensDvsIO(j) = temp(2);  
        temp = weightedcorrs([sDat.fullBeeData(inds(j)).calInterO,sDat.fullBeeData(inds(j)).RetThickness], facetWeightInds{inds(j)});
        lensDvsRet(j) = temp(2);  
        temp = weightedcorrs([sDat.fullBeeData(inds(j)).calInterO,sDat.fullBeeData(inds(j)).lensThickness], facetWeightInds{inds(j)});
        lensDvsLens(j) = temp(2);  
        temp = weightedcorrs([sDat.fullBeeData(inds(j)).calInterO,sDat.fullBeeData(inds(j)).CCThickness], facetWeightInds{inds(j)});
        lensDvsCC(j) = temp(2);  
        temp = weightedcorrs([sDat.fullBeeData(inds(j)).calInterO,sDat.fullBeeData(inds(j)).avgRadius], facetWeightInds{inds(j)});
        lensDvsCurv(j) = temp(2);     
    
        temp = weightedcorrs([sDat.fullBeeData(inds(j)).avgDiameter,sDat.fullBeeData(inds(j)).RetThickness], facetWeightInds{inds(j)});
        IOvsRet(j) = temp(2);  
        temp = weightedcorrs([sDat.fullBeeData(inds(j)).avgDiameter,sDat.fullBeeData(inds(j)).lensThickness], facetWeightInds{inds(j)});
        IOvsLens(j) = temp(2);  
        temp = weightedcorrs([sDat.fullBeeData(inds(j)).avgDiameter,sDat.fullBeeData(inds(j)).CCThickness], facetWeightInds{inds(j)});
        IOvsCC(j) = temp(2);  
        temp = weightedcorrs([sDat.fullBeeData(inds(j)).avgDiameter,sDat.fullBeeData(inds(j)).avgRadius], facetWeightInds{inds(j)});
        IOvsCurv(j) = temp(2);
        
        temp = weightedcorrs([sDat.fullBeeData(inds(j)).RetThickness,sDat.fullBeeData(inds(j)).lensThickness], facetWeightInds{inds(j)});
        RetvsLens(j) = temp(2);  
        temp = weightedcorrs([sDat.fullBeeData(inds(j)).RetThickness,sDat.fullBeeData(inds(j)).CCThickness], facetWeightInds{inds(j)});
        RetvsCC(j) = temp(2);  
        temp = weightedcorrs([sDat.fullBeeData(inds(j)).RetThickness,sDat.fullBeeData(inds(j)).avgRadius], facetWeightInds{inds(j)});
        RetvsCurv(j) = temp(2);
        
        temp = weightedcorrs([sDat.fullBeeData(inds(j)).lensThickness,sDat.fullBeeData(inds(j)).CCThickness], facetWeightInds{inds(j)});
        LensvsCC(j) = temp(2);  
        temp = weightedcorrs([sDat.fullBeeData(inds(j)).lensThickness,sDat.fullBeeData(inds(j)).avgRadius], facetWeightInds{inds(j)});
        LensvsCurv(j) = temp(2);
        
        temp = weightedcorrs([sDat.fullBeeData(inds(j)).CCThickness,sDat.fullBeeData(inds(j)).avgRadius], facetWeightInds{inds(j)});
        CCvsCurv(j) = temp(2);
    end
    subplot(5,5,1); hold on
    plot(tempX,lensDvsIO, groupMarkerStyles{i}, 'color', facetCol,'markerfacecolor','w','markersize',10);
    if i==2; text(0.7,-0.6,sprintf('%.2f', mean(lensDvsIO)));end
    subplot(5,5,2); hold on
    plot(tempX,lensDvsRet, groupMarkerStyles{i}, 'color', facetCol,'markerfacecolor','w','markersize',10);
    if i==2; text(0.7,-0.6,sprintf('%.2f', mean(lensDvsRet)));end
    subplot(5,5,3); hold on
    plot(tempX,lensDvsLens, groupMarkerStyles{i}, 'color', facetCol,'markerfacecolor','w','markersize',10);
    if i==2; text(0.7,-0.6,sprintf('%.2f', mean(lensDvsLens)));end
    subplot(5,5,4); hold on
    plot(tempX,lensDvsCC, groupMarkerStyles{i}, 'color', facetCol,'markerfacecolor','w','markersize',10);
    if i==2; text(0.7,-0.6,sprintf('%.2f', mean(lensDvsCC)));end
    subplot(5,5,5); hold on
    plot(tempX,lensDvsCurv, groupMarkerStyles{i}, 'color', facetCol,'markerfacecolor','w','markersize',10);
    if i==2; text(0.7,0.6,sprintf('%.2f', mean(lensDvsCurv)));end
    
    subplot(5,5,7); hold on
    plot(tempX,IOvsRet, groupMarkerStyles{i}, 'color', facetCol,'markerfacecolor','w','markersize',10);
    if i==2; text(0.7,-0.6,sprintf('%.2f', mean(IOvsRet)));end
    subplot(5,5,8); hold on
    if i==2; text(0.7,-0.6,sprintf('%.2f', mean(IOvsLens)));end
    plot(tempX,IOvsLens, groupMarkerStyles{i}, 'color', facetCol,'markerfacecolor','w','markersize',10);
    subplot(5,5,9); hold on
    if i==2; text(0.7,-0.6,sprintf('%.2f', mean(IOvsCC)));end
    plot(tempX,IOvsCC, groupMarkerStyles{i}, 'color', facetCol,'markerfacecolor','w','markersize',10);
    subplot(5,5,10); hold on
    if i==2; text(0.7,-0.6,sprintf('%.2f', mean(IOvsCurv)));end
    plot(tempX,IOvsCurv, groupMarkerStyles{i}, 'color', facetCol,'markerfacecolor','w','markersize',10);
    
    subplot(5,5,13); hold on
    plot(tempX,RetvsLens, groupMarkerStyles{i}, 'color', facetCol,'markerfacecolor','w','markersize',10);
    if i==2; text(0.7,-0.6,sprintf('%.2f', mean(RetvsLens)));end
    subplot(5,5,14); hold on
    if i==2; text(0.7,-0.6,sprintf('%.2f', mean(RetvsCC)));end
    plot(tempX,RetvsCC, groupMarkerStyles{i}, 'color', facetCol,'markerfacecolor','w','markersize',10);
    subplot(5,5,15); hold on
    if i==2; text(0.7,-0.6,sprintf('%.2f', mean(RetvsCurv)));end
    plot(tempX,RetvsCurv, groupMarkerStyles{i}, 'color', facetCol,'markerfacecolor','w','markersize',10);
    
    subplot(5,5,19); hold on
    if i==2; text(0.7,-0.6,sprintf('%.2f', mean(LensvsCC)));end
    plot(tempX,LensvsCC, groupMarkerStyles{i}, 'color', facetCol,'markerfacecolor','w','markersize',10);
    subplot(5,5,20); hold on
    if i==2; text(0.7,-0.6,sprintf('%.2f', mean(LensvsCurv)));end
    plot(tempX,LensvsCurv, groupMarkerStyles{i}, 'color', facetCol,'markerfacecolor','w','markersize',10);
    
    subplot(5,5,25); hold on
    plot(tempX,CCvsCurv, groupMarkerStyles{i}, 'color', facetCol,'markerfacecolor','w','markersize',10);
    if i==2; text(0.7,-0.6,sprintf('%.2f', mean(CCvsCurv)));end
end

for i = 1:5
    for j = 1:5
       if j >= i
            subplot(5,5,(i-1)*5+j);
           ylim([-1 1]); xlim([0.5 0.9]);
           set(gca, 'YTick', [-1 0 1]);
           set(gca,'XTick',[0.5 0.7 0.9],'TickDir','out');
       end
    end
end
set(gcf, 'Position', [50 50 1600 1000]);

%% save all the stuff
if saveAllThePlots
    cd(saveFolder);
    
    saveas(volAllomF, 'volAllomFig.svg');
    saveas(facetAllomF, 'facetAllomF.svg');
    saveas(fovRatioAndEyeP, 'fovRatioAndEyeP.svg');
    saveas(fovAllomF, 'fovAllomF.pdf');
    saveas(sensAllomF, 'sensAllomF.pdf');
    
    saveas(IOHistsF, 'IOHistsF.svg');
    saveas(DiamHistsF, 'DiamHistsF.svg');
    saveas(RhabdomsHistsF, 'RhabdomsHistsF.svg');
    saveas(corneaHistF, 'corneaHistF.svg');
    saveas(CCHistF, 'CCHistF.svg');
    saveas(curvHistF, 'curvHistF.svg');
    saveas(eyeParamHistF, 'eyeParamHistF.svg');
    
    %has transparent patch - matlab can't save a nice vector graphic, fml
    saveas(worldSphereF, 'worldSphereF.tif');
    saveas(projFullFov, 'projFullFov.svg');
    saveas(projBinoFOV, 'projBinoFOV.svg');
    
    saveas(azIntFullFOVF, 'azIntFullFOVF.svg');
    saveas(elIntFullFOVF, 'elIntFullFOVF.svg');
    saveas(azIntBinoFOVF, 'azIntBinoFOVF.svg');
    saveas(elIntBinoFOVF, 'elIntBinoFOVF.svg');
    
    saveas(azAvIOF, 'azAvIOF.svg');
    saveas(elAvIOF, 'elAvIOF.svg');
    saveas(azAvDF, 'azAvDF.svg');
    saveas(elAvDF, 'elAvDF.svg');
    saveas(azAvRetF, 'azAvRetF.svg');
    saveas(elAvRetF, 'elAvRetF.svg');
    saveas(azAvCCF, 'azAvCCF.svg');
    saveas(elAvCCF, 'elAvCCF.svg');
    saveas(azAvLensF, 'azAvLensF.svg');
    saveas(elAvLensF, 'elAvLensF.svg');
    saveas(azAvCurvF, 'azAvCurvF.svg');
    saveas(elAvCurvF, 'elAvCurvF.svg');
    saveas(azAvEyeParamF, 'azAvEyeParamF.svg');
    saveas(elAvEyeParamF, 'elAvEyeParamF.svg');

    
    %only half of c bar opens in inkscape when saved as svg
    %... but pdf leads to unformated text
    saveas(resCMapF, 'resCMapF.pdf');
    for i = 1:length(toPlotMaps)
        saveas(IOIndFs(i), sprintf('IOIndF_%s.pdf', sDat.fullBeeData(toPlotMaps(i)).BeeID));
    end
    saveas(facetDCMapF, 'facetDCMapF.pdf');
    for i = 1:length(toPlotMaps)
        saveas(facetDIndFs(i), sprintf('facetDIndFs_%s.pdf', sDat.fullBeeData(toPlotMaps(i)).BeeID));
    end
    saveas(rhabLengthCMap, 'rhabLengthCMap.pdf');
    for i = 1:length(toPlotMaps)
        saveas(rhabLIndFs(i), sprintf('rhabLIndFs_%s.pdf', sDat.fullBeeData(toPlotMaps(i)).BeeID));
    end
    saveas(eyeParaCMapF, 'eyeParaCMap.pdf');
    for i = 1:length(toPlotMaps)
        saveas(eyeParaIndFs(i), sprintf('eyeParaIndFs_%s.pdf', sDat.fullBeeData(toPlotMaps(i)).BeeID));
    end
    saveas(lensTCMapF, 'lensTCMapF.pdf');
    for i = 1:length(toPlotMaps)
        saveas(lensTIndFs(i), sprintf('lensTIndFs_%s.pdf', sDat.fullBeeData(toPlotMaps(i)).BeeID));
    end
    saveas(CCTCMapF, 'CCTCMapF.pdf');
    for i = 1:length(toPlotMaps)
        saveas(CCIndFs(i), sprintf('CCTIndFs_%s.pdf', sDat.fullBeeData(toPlotMaps(i)).BeeID));
    end
    saveas(CurvCMapF, 'CurvCMapF.pdf');
    for i = 1:length(toPlotMaps)
        saveas(CurvIndFs(i), sprintf('CurvIndFs_%s.pdf', sDat.fullBeeData(toPlotMaps(i)).BeeID));
    end
    saveas(volCMapF, 'volCMapF.pdf');
    
    saveas(scaleMapF, 'scaleMapF.pdf');
    saveas(eyePScalFig, 'eyePScalFig.pdf');
    saveas(interOScalFig, 'interOScalFig.pdf');
    saveas(facetDScalFig, 'facetDScalFig.pdf');
    saveas(lensTScalFig, 'lensTScalFig.pdf');
    saveas(retTScalFig, 'retTScalFig.pdf');
    saveas(CCTScalFig, 'CCTScalFig.pdf');
    saveas(curvScalFig, 'curvScalFig.pdf');
    
    saveas(superCoreF, 'superCoreF.svg');
    
    display('all the plots are saved');
end
