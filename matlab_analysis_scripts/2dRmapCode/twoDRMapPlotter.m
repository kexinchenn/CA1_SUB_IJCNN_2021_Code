function map = twoDRMapPlotter(neuronStruct, recStruct, iNeu, data2PlotName, plotData, mapToUse, filterName, filterValues, sampleRate, caxisMax)
%twoDRMapPlotter One line description of what the function or script performs (H1 line)
%		% SubCA1 Plotter - Cleans out the bad spots, plots the data or returns the map
%
% 	map = twoDRMapPlotter(neuronStruct, recStruct, iNeu, data2PlotName, plotData, filterName, filterValues)
%
% INPUTS:
%		neuronStruct - Description
%		recStruct - Description
%		iNeu - index of neuron in neuronStruct to use
%		data2PlotName (optional - make empty '' or [] to bypass) - field to map and plot
%		plotData (optional - make empty '' or [] to bypass) - logical true/false if you want the
%           plot created or suppressed
%		mapToUse (optional - make empty '' or [] to bypass) - colormap to use
%		filterName (optional - make empty '' or [] to bypass) - field to use as a filter
%		filterValues (optional - make empty '' or [] to bypass) - value(s) to use as a filter
%       sampleRate (optional - make empy '' or [] to bypass) - value to transform map from per-sample to per-second 
%
% OUTPUTS:
%		map - 2D array of values plotted.
%
%
% EXAMPLES:
%    twoDRMapPlotter(neuronStruct, recStruct, iNeu, 'sampleNSpikes', true, mapToUse);
%
%	 twoDRMapPlotter(neuronStruct, recStruct, iNeu, 'sampleNSpikes', true, mapToUse,...
%            'samplePaths',thisPath);
%
%
% REMARKS
%
% OTHER M-FILES REQUIRED: 
%   mapVar2D
%   filter2DMatrices
%   sc
% SUBFUNCTIONS:
% MAT-FILES REQUIRED:
%
% SEE ALSO:
%
%
% CREATED BY: Jacob M Olson
% EMAIL: jolson1129@gmail.com
% WEBSITE: http://www.jmolson.com
% CREATED ON: 07-May-2020
% LAST MODIFIED BY:
% LAST MODIFIED ON:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize Constants, handle inputs
XMIN = 50;
XMAX = 600;
YMIN = 50;
YMAX = 800;

% max bin values can change w/out being considered a jump and removed
JUMP_MAX = 35;

if ~exist('plotData','var')
    plotData = true;
elseif  isempty(plotData) || ~islogical(plotData)
    plotData = true;
end

if ~exist('mapToUse','var')
    mapToUse = parula;
elseif  isempty(mapToUse)
    plotData = parula;
end

if ~exist('filterName','var')
    filterName = '';
elseif ~exist('filterValues','var')
    filterName = '';
    disp('No filter values - ignoring filter named');
end

if ~exist('sampleRate','var')
    sampleRate = 60; %Default to sampling rate of 60Hz
elseif  isempty(sampleRate)
    sampleRate = 60; %Default to sampling rate of 60Hz
end

if ~exist('caxisMax','var')
    caxisMax = '';
elseif  isempty(caxisMax)
    caxisMax = '';
end

%%
recIndex = neuronStruct.recStructIndex(iNeu);

startEpochSample = recStruct.sessionTimeStamps{recIndex}(3);
endEpochSample = recStruct.sessionTimeStamps{recIndex}(4);
tracking = recStruct.Behavior.posXY{recIndex}(startEpochSample:endEpochSample,:);
data2Plot = neuronStruct.(data2PlotName){iNeu}(startEpochSample:endEpochSample);

% Clean up tracking for plotting
if ~isempty(filterName)
    filterData = (recStruct.Behavior.(filterName){recIndex}(startEpochSample:endEpochSample));
    xs = tracking(filterData == filterValues,1);
    ys = tracking(filterData == filterValues,2);
else
    xs = tracking(:,1);
    ys = tracking(:,2);
end
%take out rows with nans
badSpotsNans = isnan(xs) | isnan(ys);
% remove off track from plotting
badSpotsOffTrack = xs<XMIN | xs>XMAX | ys<YMIN | ys>YMAX;
% remove jumps from glare from plotting.
badSpotsJumps = [0; abs(diff(xs)) > JUMP_MAX | abs(diff(ys)) > JUMP_MAX ];

anyBadSpots = badSpotsNans | badSpotsOffTrack | badSpotsJumps;
xsClean = xs(~anyBadSpots);
ysClean = ys(~anyBadSpots);

map = mapVar2D(data2Plot(~anyBadSpots), xsClean, ysClean, [XMAX,YMAX]);
map = map.*sampleRate; %make map convert from val/sample to val/sec
filtNormMap = filter2DMatrices(map, false);
% max of colorbar which = mean(velocity)+2*std(velocity)

%Define color axis
colormax = colorAxisCalculator(filtNormMap, caxisMax);

if plotData
    fig = figure('Name',['Neuron ',num2str(iNeu),'. ','Color Axis Max = ',num2str(colormax)]);
    sc(filtNormMap,[0,colormax],mapToUse,'w',isnan(filtNormMap));   

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of Code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


