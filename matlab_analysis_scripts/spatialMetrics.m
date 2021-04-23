function SummaryStats = spatialMetrics(ClassicSpatialMeasuresHPC, ClassicSpatialMeasuresSUB)
% Classic Spatial Stats Wrapper
% Main fn in classicSpatialMeasurs.m
%
% Written by Jacob Olson October 2020
%

% [ ClassicSpatialMeasuresHPC ] = convertTo2dMaps(HpcRec, HpcNeuron, hpcNeurIds);
% [ ClassicSpatialMeasuresSUB ] = convertTo2dMaps(SubRec, SubNeuron, subNeurIds);

AllFieldNames = fieldnames(ClassicSpatialMeasuresHPC);
meanFieldValsHPC = structfun(@mean, ClassicSpatialMeasuresHPC);
meanFieldValsSUB = structfun(@mean, ClassicSpatialMeasuresSUB);
stdFieldValsHPC = structfun(@nanstd, ClassicSpatialMeasuresHPC);
stdFieldValsSUB = structfun(@nanstd, ClassicSpatialMeasuresSUB);

SummaryStats = array2table([meanFieldValsHPC,stdFieldValsHPC,meanFieldValsSUB,stdFieldValsSUB],...
    'VariableNames',{'meanHPC', 'stdHPC','meanSUB','stdSUB'});
SummaryStats.FieldNames = AllFieldNames;
SummaryStats = movevars(SummaryStats,'FieldNames','Before',1);

% tail right means first var bigger is the alt hypothesis
% creating table - some measures we hypoth and would test hpc bigger, others sub. Just grabbing the
% relevant value out of the table.
for iField= 1:length(SummaryStats.FieldNames)
    iFieldName = SummaryStats.FieldNames{iField};
    areDifferent(iField) = ranksum(ClassicSpatialMeasuresHPC.(iFieldName),ClassicSpatialMeasuresSUB.(iFieldName));
    hpcBigger(iField) = ranksum(ClassicSpatialMeasuresHPC.(iFieldName),ClassicSpatialMeasuresSUB.(iFieldName),...
        'tail','right');
    subBigger(iField) = ranksum(ClassicSpatialMeasuresHPC.(iFieldName),ClassicSpatialMeasuresSUB.(iFieldName),...
        'tail','left');
end
SummaryStats.areDifferentRankSumP = areDifferent';
SummaryStats.HPCBiggerRankSumP = hpcBigger';
SummaryStats.SUBBiggerRankSumP = subBigger';

end