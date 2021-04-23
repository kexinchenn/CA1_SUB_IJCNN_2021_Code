function figure2_SUBAN(SubNeuron, HpcNeuron, SubRec, HpcRec,ClusterQuality_forUsed, caxisMax)
%% Figure 2 script - 2D Ratemap examples

%% Load in Isolation Distances from Previously Compiled Files 
% 2B - Isolation Distances in log,log 
figure(1)
    histogram(log10(ClusterQuality_forUsed{2,1}),'binwidth',0.1)
    title("logxlog plot of CA1 neurons Isolation Distances")
    set(gca,'YScale','log')
    xlim([0 5])
    ylim([0 100])
    figure(2)
    histogram(log10(ClusterQuality_forUsed{2,2}),'binwidth',0.1)
    title("logxlog plot of SUB neurons Isolation Distances")
    set(gca,'YScale','log')
    xlim([0 5])
    ylim([0 100])


%% 2 C - D 2d ratemaps 

% HPC Neurons
% BL02 rec 23 8 - 62 in array as of 8/28/2020
% NS23 rec 18 5 - 208
% NS23 rec 19 8 - 225
% JL01 rec 31 2 - 273

%OLD_NEURON_LIST_HPC = [56, 180, 190, 200];
NEURON_LIST_HPC = [208, 225, 273, 62]; %right cells for 8/28/20

% SUB Neurons
% NS23 rec 21 7  - 9 in array as of 4/29/2020
% NS23 rec 22 2 - 13 
% NS15 rec 14 14 - 143
% NS15 rec 16 3  -149
% NS15 rec 19 4  - 162
% NS15 rec 23 6  - 210
% NS15 rec 29 11 - 314
% NS16 rec 06 1  - 352

NEURON_LIST_SUB = [9, 13, 143, 149, 162, 210, 314, 352];

%% Plot all examples
mapToUse = parula;

%Sub neurons
neuronStruct = SubNeuron;
recStruct = SubRec;
for iNeu = NEURON_LIST_SUB
    twoDRMapPlotter(neuronStruct, recStruct, iNeu, 'sampleNSpikes', true, mapToUse,[],[],[],caxisMax);
end

%Hippocampus neurons
neuronStruct = HpcNeuron;
recStruct = HpcRec;
for iNeu = NEURON_LIST_HPC
    twoDRMapPlotter(neuronStruct, recStruct, iNeu, 'sampleNSpikes', true, mapToUse,[],[],[],caxisMax);
end

end



