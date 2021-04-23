set(0,'defaultAxesFontSize',20)
%%
c1 = [243 169 114]./255;    % Orange theme
c2 = [128 193 219]./255;    % Blue theme
colors = [c1;c2];
%% FIGURE 3: Network Fitness
addpath stdshade.m
resRootDir = '~/Downloads/IJCNN_2021_Code/matlab_analysis_scripts/data/';
% resRootDir = 'data/';
names = {'hpc', 'sub'};
fitFile = 'UnlesionedFitnessScores.csv';
numI = 15;
numGen = 50;
fitnessMaxAll = cell(2,1);

for r = 1:length(names)
    name = names{r};
    root = [resRootDir, name];
    result_dir = dir(root);  % get a list of all files and folders
    folders = {result_dir.name};  % put the names into cell array
    folders = folders([result_dir.isdir] & startsWith(folders, name)); 

    numRuns = length(folders);
    fitnessMaxAll{r} = zeros(numRuns, numGen);
    for n=1:numRuns
        d = folders{n};
        fitness = csvread(fullfile(root, d, fitFile));
        fitness = fitness(2:2:end); % first number: max mean FR; second number: fitness
        fitnessGrouped = reshape(fitness, numI, numGen);
        fitnessMax = max(fitnessGrouped, [], 1);
        % find best so far fitness
        for g = 2:numGen
            if fitnessMax(g-1) > fitnessMax(g)
                fitnessMax(g) = fitnessMax(g-1);
            end
        end 
        fitnessMaxAll{r}(n,:) = fitnessMax;
    end
end

figure
[aline(1), aFill(1)] = stdshade(fitnessMaxAll{1},0.2,c1);
hold on
[aline(2), aFill(2)] = stdshade(fitnessMaxAll{2},0.2,c2);
hold off
legend(aline, {'CA1', 'SUB'})
set(gca, 'YGrid', 'on', 'XGrid', 'off')
set(gca, 'FontSize', 18)
yticks(90:10:220)
ylim([90 220])

xlabel('Generation','fontweight','bold')
ylabel('Fitness', 'fontweight','bold')
title('Best-So-Far Fitness','fontweight','bold') 
%% FIGURE 4: Rate Maps
addpath 2dRmapCode/
addpath ojwoodford-sc-21df4f5/

% CA1,1_Driver/Unlesioned
% neurons = [34,112,113,217,239,306,379,380,391,451,592]; 
% SUB,3_Driver/Unlesioned
neurInds = [3,1,2;4,3,1;5,4,1;6,2,5;8,5,4;8,5,5;8,4,4;12,1,1;12,4,5];
neurons = (neurInds(:,1)-1)*25+(neurInds(:,2)-1)*5+neurInds(:,3)-1; 
nNeur = length(neurons);

nTrials = 100;
x = simStruct.Behavior.PosX; %  numTrials x numBins
y = simStruct.Behavior.PosY; %  numTrials x numBins

for n = 1:nNeur
    neur = neurons(n);
    ratemapSim = zeros(700,500);
    maxRSim = 0;

    for t=1:nTrials
        pos = [floor(x(t,:));floor(y(t,:))];
        rSim = simStruct.LinearRates{neur}(t,:)';
        rSim(rSim<1) = 0.01;
        I = sub2ind(size(ratemapSim),pos(1,:),pos(2,:));
        ratemapSim(I)=ratemapSim(I)'+rSim;
        maxRSim = max(maxRSim, max(rSim));
    end

    fig = figure('Name',['Neuron ',num2str(neur),'. ','maxFR ', num2str(ceil(maxRSim))]);
    ratemapSim(ratemapSim == 0) = nan;
    ratemapSim = filter2DMatrices(ratemapSim,true);
    
    colormax = colorAxisCalculator(ratemapSim, '');
    ratemapSim = ratemapSim';

    sc(ratemapSim,[0,colormax],parula,[0.5 0.5 0.5],isnan(ratemapSim));   
    shading flat
%     saveas(gcf, fullfile(saveDir,['Neuron ',num2str(neur),'. ','Color Axis Max = ',num2str(colormax),'.png']))
end

%% FIGURE 5: Connection Weights
addpath ~/CARLsim5/tools/offline_analysis_toolbox/
nAngVel = 12;
nHD = 8;
nPos = 25*18;
nLinVel = 12;

angVelInd = 1:nAngVel;
HDInd = angVelInd(end)+[1:nHD];
posInd = HDInd(end)+[1:nPos];
linVelInd = posInd(end)+[1:nLinVel];

regions = {'CA1', 'SUB'};
conn = 'inp_exc';
figure
for r = 1:2
    simDir = [resRootDir, regions{r}, '_'];
    if r == 1
        indi = 2;
    elseif r == 2
        indi = 4;
    end
    netName = [num2str(indi-1) '_Driver/'];    
    weightFile = ['conn_' conn '.dat'];
    CR = ConnectionReader([simDir, netName,weightFile]);
    [allTimestamps, allWeights] = CR.readWeights();
    allWeights = reshape(allWeights(1,:), CR.getNumNeuronsPost(), CR.getNumNeuronsPre());

    wtGrouped = cell(4,1);
    wtGrouped{1} = allWeights(:, angVelInd);
    wtGrouped{3} = allWeights(:, HDInd);
    wtGrouped{4} = allWeights(:, posInd);
    wtGrouped{2} = allWeights(:, linVelInd);

    inNames = ["AV", "LV", "HD", "Pos"];
    for i=1:4
        subplot(2,4,i+(r-1)*4)
        h1 = bar(hist(wtGrouped{i}(:)) ./ sum(hist(wtGrouped{i}(:))),...
        'FaceColor', colors(r,:),'EdgeColor', [0.5 0.5 0.5],'FaceAlpha', 1.0);
        ax = gca;
        ax.LineWidth = 1.0;

        xlim([0.5 10.5])
        ylim([0 1])
        if r == 1
            title(inNames(i))
            xticks([])
        elseif r == 2
            xticks([0.5 10])
            xticklabels({'0', 'Max'})
        end

        if i == 1
            yticks([0 1])
        else
            yticks([])
        end

        if i == 1
            legend(regions{r}, 'Location','northwest')
        end
    end
end
%% FIGURE 6: STDP
paramFile = [resRootDir, 'evolvedParametersCA1SUB.mat'];
load(paramFile)

W_m = cell(2,1);
W_p = cell(2,1);
modes = ["EE" "EI" "IE"];

nTrials = [5,5];
figure()
for m = 1:length(modes)
    mode = modes(m);
    subplot(1,length(modes),m)
    t_neg = 0:1:100;
    t_pos = 0:1:100;

    switch mode
        case 'EE'
            cs = [1 4 7 8];
        case 'EI'
            cs = [2 5 9 10];
        case 'IE'
            cs = [6 3 11 12];
    end
    
    for r = 1:2
        W_p{r} = zeros(nTrials(r),length(t_pos)); 
        W_m{r} = zeros(nTrials(r),length(t_neg));
        for i = 1:nTrials(r)
            p = params{r}(i,cs);
            A_p = p(1);
            A_m = p(2);
            tao_p = p(3);
            tao_m = p(4);

            w_p = A_p*exp(-t_pos/tao_p);
            w_m = A_m*exp((-t_neg)/tao_m);

            W_m{r}(i,:) = w_m; 
            W_p{r}(i,:) = w_p;
        end
    end

    CA1_indi = 1:5;
    SUB_indi = 1:5;
    [aline(1), aFill(1)] = stdshade([flip(W_m{1}(CA1_indi,:),2),W_p{1}(CA1_indi,:)],...
        0.2,c1,[flip(-t_neg),t_pos]);
    hold on; 
    [aline(2), aFill(2)] = stdshade([flip(W_m{2}(SUB_indi,:),2),W_p{2}(SUB_indi,:)],...
        0.2,c2,[flip(-t_neg),t_pos]);
    hold off

    xlimit = get(gca, 'xlim');
    xlabel('t_{post-pre}')

    ymax = 4*10^(-3);
    ylimit = [-ymax ymax];
    ylim(ylimit)
    ylabel('\Delta w')
    refline(0,0)
    line([0 0], ylimit)
    
    if m == length(modes)
        f=get(gca,'Children');
        legend([f(5),f(3)],'CA1','SUB')
    end
    
    if mode == "EE"
        title('E-STDP on Exc (EE)')
    elseif mode == "EI"
        title('E-STDP on Inh (EI)')
    elseif mode == "IE"
        title('I-STDP on Exc (IE)')
    end 
end
%% FIGURE 7: Correlation Matrices
numBins = 954;
evenfrFile = 'Even.csv';
oddfrFile = 'Odd.csv';
figure()
for r = 1:2
    simDir = [resRootDir, regions{r}, '_'];
    if r == 1
        indi = 2;
    elseif r == 2
        indi = 4;
    end
    netName = [num2str(indi-1) '_Driver/'];    
    evenfr = csvread([simDir, netName, evenfrFile]); 
    evenfr = evenfr(:,1:numBins);
    oddfr = csvread([simDir, netName, oddfrFile]); 
    oddfr = oddfr(:,1:numBins);

    subplot(1,2,r)
    mat = corr(evenfr,oddfr);
    imagesc(mat)
    title("Odd/Even Correlation Matrix")
    ylabel("position bin (even runs)")
    xlabel("position bin (odd runs)")
    caxis([0 1]) 
end
%% FIGURE 8: Ablation Studies
ablationDataFile = [resRootDir, 'AblationStudies.mat'];
load(ablationDataFile)

lesionsCA1 = {'Unlesioned', 'AV', 'HD', 'Pos', 'LV', ...
'HD_Pos', 'AV_HD_LV'};
lesionsSUB = {'Unlesioned', 'AV', 'HD', 'Pos', 'LV', ...
'HD_Pos', 'AV_HD_LV'};

meanStdLesionScores = cell(2,1);
for r = 1:2
    meanStdLesionScores{r} = cellfun(@(x) [mean(x), nanstd(x)],lesionScores{r}(:,1),'UniformOutput',false);
    meanStdLesionScores{r} = cell2mat(meanStdLesionScores{r});
end
regions = {'CA1','SUB'};
for r = 1:2
    if r == 1
        lesions = lesionsCA1;
    elseif r == 2
        lesions = lesionsSUB;
    end
    subplot(2,1,r)
    [sortedScores, sortInd] = sort(meanStdLesionScores{r}(:,1), 'descend');
    h = bar(sortedScores);
    h.FaceColor = colors(r,:);
    h.EdgeColor = [0.5 0.5 0.5];
    
    hold on
    stdwidth = meanStdLesionScores{r}(:,2);
    stdwidth = stdwidth(sortInd);
    er = errorbar(1:length(lesions), sortedScores,...
        stdwidth,stdwidth);    
    er.Color = [0.5 0.5 0.5];                            
    er.LineStyle = 'none';  
    er.LineWidth = 1;

    hold off
    
    xticklabels(lesions(sortInd))
    ylim([0 0.8])
    if (r==2)
        xlabel('Lesion Model', 'fontweight','bold')
    end
    ylabel('Similarity Score', 'fontweight','bold')
    set(gca,'TickLabelInterpreter','none')
    legend(regions{r})
end
%% TABLE 1: Spatial Analyses
spatialDataFile = [resRootDir, 'hpcSubSpatialStats.mat'];
load(spatialDataFile)
SummaryStats = spatialMetrics(ClassicSpatialMeasuresHPC, ClassicSpatialMeasuresSUB);
disp(SummaryStats)