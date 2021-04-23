function [finalMap, occMap] = mapVar2D(varToPlot, xs, ys, mapSize, isCircular)
    %SSSSSS One line description of what the function or script performs (H1 line)
    %		Optional file header info (to give more details about the function than in the H1 line)  
    %
    % 	[output1, output2] = ssssss(input1, input2)
    %
    % INPUTS: 
    %		varToPlot - Description
    %		tracking - Description
    %		input3 - Description
    %		                      
    %
    % OUTPUTS: 
    %		output1 - Description
    %		                       
    %
    % EXAMPLES: 
    %		Line 1 of example
    %		Line 2 of example
    %		Line 3 of example
    %		                   
    %
    % REMARKS:
    %
    %       Does not scale - [makes values per sample not values per second]
    %
    % OTHER M-FILES REQUIRED: 
    % SUBFUNCTIONS: 
    % MAT-FILES REQUIRED: 
    %
    % SEE ALSO:%
    % CREATED WITH MATLAB VERSION: 9.6.0.1150989 (R2019a) Update 4 on Microsoft Windows 10 Enterprise Version 10.0 (Build 18363)
    %
    % CREATED BY: Jacob M Olson
    % EMAIL: jolson1129@gmail.com
    % WEBSITE: http://www.jmolson.com
    % CREATED ON: 28-Apr-2020
    % LAST MODIFIED BY: 
    % LAST MODIFIED ON: 
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%rea
    if exist('mapSize','var')
        if ~isempty(mapSize)
            occMap = zeros(mapSize(1),mapSize(2));
        else
            occMap = zeros(700,500);
        end
    else
        occMap = zeros(700,500);
    end
    valueMap = occMap;
    finalMap = occMap;
    circleMap = cell(size(occMap));

    for i = 1:length(xs)
        valueMap(xs(i),ys(i)) = valueMap(xs(i),ys(i))+varToPlot(i);
        circleMap{xs(i),ys(i)} = [circleMap{xs(i),ys(i)}, varToPlot(i)];
        occMap(xs(i),ys(i)) = occMap(xs(i),ys(i))+1;
    end

    if exist('isCircular','var')
        for i = 1:size(occMap,1)
            for j = 1:size(occMap,2)
                if ~isempty(circleMap{i,j})
                    finalMap(i,j) = circ_mean(circleMap{i,j},[],2);
                else
                    finalMap(i,j) = NaN;
                end
            end
        end
    else
        finalMap = valueMap./occMap;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of Code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


