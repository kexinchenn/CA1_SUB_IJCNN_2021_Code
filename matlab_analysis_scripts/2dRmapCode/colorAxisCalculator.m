function [colormax] = colorAxisCalculator(twoDMap, style)

% Defines color axis max 
%Input for 'style' is a string defining how to calculate the color map max:
%   'Max'      - Colormap max is max value in data
%   'Mean'     - Colormap max is Mean value 
%   'MeanXSD'  - *Default X2*  7 character string Colormap max is Mean value + X Standard Deviations character 5 is an integar
%   'PcntXX'   - 6 character string first 4 'Pcnt' identify, second 2 are integer value to define percent of max to use 

%Set default colormap calculator 
if ~exist('style','var')
    style = 'Mean2SD';
elseif  isempty(style)
    style = 'Mean2SD';
end

%Pull Max, Mean, Standard Dev. Initialize colormax
colormax = 0;
valMax = max(max(twoDMap));
valMean = nanmean(nanmean(twoDMap));
valSTD = nanstd(twoDMap(:));

%different if-statements to define different inputs 
%ADD MORE TO END IF WANTED 

if strcmp(style,'Max')
    colormax = valMax;
end

if strcmp(style,'Mean')
    colormax = valMean;
end

if length(style)==7
    if strcmp(style(1:4),'Mean') && strcmp(style(6:7),'SD') %double check
        cscale = str2double(style(5)); %turn digit 5 into a number of SD to scale by
        colormax = valMean+(cscale*valSTD);
    end
end

if length(style)==6;
    if strcmp(style(1:4),'Pcnt') %double check
        cscale = ((str2double(style(5:6))/100)); %turn last 2 digits into a number to scale by
        colormax = valMax*cscale;
    end
end

%old hardcoded percent 
% if strcmp(style,'Pcnt90');
%     colormax = valMax*0.9;
% end
% 
% if strcmp(style,'Pcnt66');
%     colormax = valMax*(2/3);
% end
% 
% if strcmp(style,'Pcnt50');
%     colormax = valMax/2;
% end

end