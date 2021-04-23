function [twoDRMaps] = twoDPathsRateMapper( DataStruct, NeuronStruct, ...
    binSizePixels, xPixelSize, yPixelSize, neuronList)
%twoDPathsRateMapper 
%   Takes a neuron list and plots 2D rate maps of paths using the bin and
%   pixelSize parameters.
%
%   Identical to twoDRateMapper but for paths.
%
%   Written by Jake Olson, December 2015 updated april 2017

% 6.666 pixels is ~2cm square bins for platform - RECOMMENDED
% probably 6 pixels would be about the same distance for the track, since
% it is about 1 brick lower. Try 6 pixels - RECOMMENDED
% For the platform, 300 by 300 pixelSize variables gives a 45 bin by 45 bin square.

%For regular datasets of TTT maze in Nitz lab @UCSD 
maxX = 650;
maxY = 500;

% 1 as an input value means to use all of the data.
if xPixelSize == 1
	xPixelSize = maxX;  
end
if yPixelSize == 1
    yPixelSize = maxY;
end

nPaths = size(DataStruct.Behavior.pathList,1);

% initialize maps
xSize = ceil(xPixelSize/binSizePixels);
ySize = ceil(yPixelSize/binSizePixels);
twoDRMaps = nan(xSize,ySize,nPaths,length(neuronList));

for iNeuron = 1:length(neuronList)
    % data for the neuron of interest
    thisNeuron = neuronList(iNeuron);
    iRec = NeuronStruct.recStructIndex(thisNeuron);
    sampleRate = DataStruct.Behavior.trackingSampleRate_Hz(iRec);
    pathList = DataStruct.Behavior.pathList{iRec};
    pathRunsLineMarkers = DataStruct.Behavior.pathRunsLineMarkers{iRec};
    
    % location tracking samples, and spikes at each tracking sample
    iNeuNSpikes = NeuronStruct.sampleNSpikes{thisNeuron};
    %remake pixelDVT due to changes in data processing 
    newDVT(:,1) = DataStruct.Behavior.sampleIndices{iRec};
    newDVT(:,2) = DataStruct.Behavior.posXY{iRec}(:,1);
    newDVT(:,3) = DataStruct.Behavior.posXY{iRec}(:,2);
    pos = newDVT;
    
    
    for iPath = 1:numel(pathList)
        thisPath = pathList(iPath);
        
        % Find all samples and corresponding spikes that were labeled as on a given path
        posThisPath = zeros(0,3);
        posNSpikesThisPath = zeros(0,1);
        for iRun = 1:numel(pathRunsLineMarkers{iPath})
            posThisPath = [posThisPath;...
                pos(pathRunsLineMarkers{iPath,1}(iRun): pathRunsLineMarkers{iPath,2}(iRun),:)];
            
            posNSpikesThisPath = [posNSpikesThisPath;iNeuNSpikes(pathRunsLineMarkers{iPath,1}(iRun): pathRunsLineMarkers{iPath,2}(iRun),:)];
        end
    
        % Binning the twoD bins - perhaps downsampling.
        if xPixelSize == maxX
            xOffset = 0;
        else
            plotMinX = nanmean(posThisPath(posThisPath(:,2)>1,2))-xPixelSize/2;
            xOffset = ceil(plotMinX/binSizePixels)-1;
        end
        if yPixelSize == maxY
            yOffset = 0;
        else
            plotMinY = nanmean(posThisPath(posThisPath(:,2)>1,3))-yPixelSize/2;
            yOffset = floor(plotMinY/binSizePixels)-1;
        end
        binX = ceil(posThisPath(:,2)/binSizePixels)-xOffset;
        binY = ceil(posThisPath(:,3)/binSizePixels)-yOffset;
    
        % Throwing out points where it would be outside the grid.
        sampleOutOfGrid = binX <= 0 | binY <= 0 | isnan(binX) | isnan(binY) |...
            binX > xPixelSize/binSizePixels | binY > yPixelSize/binSizePixels;
        binXGood = binX(~sampleOutOfGrid);
        binYGood = binY(~sampleOutOfGrid);
        posNSpikesGood = posNSpikesThisPath(~sampleOutOfGrid);
    
        % Create 2D map from path specific (downsampled and cleaned) data
        spikes = zeros(xSize,ySize);
        occs = spikes;
        isOcc = false(size(occs));
        rates = nan(size(occs));
        for iSample = 1:length(posNSpikesGood)
            occs(binXGood(iSample),binYGood(iSample)) = ...
                occs(binXGood(iSample),binYGood(iSample))+1;
            
            spikes(binXGood(iSample),binYGood(iSample)) = ...
                spikes(binXGood(iSample),binYGood(iSample)) +...
                posNSpikesGood(iSample);
        end
        isOcc = occs > 0;
        rates(isOcc) = (spikes(isOcc)./occs(isOcc)).*sampleRate;
        twoDRMaps(:,:,thisPath,iNeuron) = rates;
        
        % I don't report spikes and occs unless troubleshooting - dataset
        % becomes very large (3x)
        %             twoDRMaps(:,:,iNeuron,1) = rates;
        %             twoDRMaps(:,:,iNeuron,2) = spikes;
        %             twoDRMaps(:,:,iNeuron,3) = occs;
    end
end

