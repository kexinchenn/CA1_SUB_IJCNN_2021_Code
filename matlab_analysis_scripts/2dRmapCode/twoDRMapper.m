function [spikemap,occmap,normmap,filtnormmap,spiketimes] = twoDRMapper(cellnum, pixelDvt, allTfiles, timeStamps)

% Written by Alex Johnson

FRAMERATE = 60;

occmap  = zeros(700,500);
spikemap= zeros(700,500);

switch nargin
    case 3
            timeStamps(1,3) = find(abs(pixelDvt(:,11))<11, 1 );
            timeStamps(1,4) = find(abs(pixelDvt(:,11))<11, 1, 'last' );
    case 4
    otherwise
        disp('Wrong number of inputs - need 2 or 3');
        return;
end

lightXCol = 9;
lightYCol = 10;
timeStampCol = 2;


tracking(:,1) = pixelDvt(timeStamps(3):timeStamps(4),lightXCol);
tracking(:,2) = pixelDvt(timeStamps(3):timeStamps(4),lightYCol);
tracking(:,3) = pixelDvt(timeStamps(3):timeStamps(4),timeStampCol);
%velocity(:,1) = vel(timeStamps(3):timeStamps(4),1,1);

spiketimes = allTfiles{cellnum};

% get rid of samples outside tracking area
isGoodSample = tracking(:,2) > 0 & tracking(:,2) < 500 & ...
         tracking(:,1) > 0 & tracking(:,1) < 700;
for iSample = 1:length(tracking)
    if isGoodSample(iSample)
        occmap(tracking(iSample,1),tracking(iSample,2)) = ...
            occmap(tracking(iSample,1),tracking(iSample,2))+1;
    end
end

isGoodSpike = spiketimes<max(tracking(:,3)) & spiketimes>min(tracking(:,3)); 
for iSpike = 1:length(spiketimes)
    if isGoodSpike(iSpike)
        spiketime = spiketimes(iSpike);
        [~,nearestTrackingSampleInd] = min(abs(spiketime-tracking(:,3)));
    if isGoodSample(nearestTrackingSampleInd)
        spikemap(tracking(nearestTrackingSampleInd,1),tracking(nearestTrackingSampleInd,2)) = ...
            spikemap(tracking(nearestTrackingSampleInd,1),tracking(nearestTrackingSampleInd,2))+1;
    end
    end
end
%occmap
%spikemap
occmap(occmap == 0) = nan;
spikemap(1,1) = nan;
occmap(1,1) = nan;

normmap=spikemap./occmap;
normmap=normmap.*FRAMERATE;

filtnormmap = filter2DMatrices(normmap, 0);

end
