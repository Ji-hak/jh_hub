function plot_small_model_app_results(results,plotHorizon,TeXfiguresDir,useSlideFormat)
% plot_small_model_app_results This function plots charts for sim results
% for the small model application section.
%   Inputs are a results structure, plot horizon scalar and folder name for
%   saving output. Third (optional) input controls whether charts are saved
%   in a format better suited for slideshow presentations.

%% SOME OPTIONS SET WITHIN THIS FUNCTION
plotVarMnems = {'pie';'y';'rAnn'};
plotVarNames = {'Quarterly inflation (%)';'Output gap (%)';'Policy rate (annualised, %)'};
% Common y-axis limits for the imperfect credibility comparison
yLims = {[-0.1 7];[-1 7];[0 7]};
% Option to control whether ex post paths in LFL sims are plotted
plotExPostPaths = false;

%% WORK OUT OTHER OPTIONS HERE
if nargin < 3
    saveTeXfigures = false;
else
    saveTeXfigures = true;
end
figStub = '';
if nargin < 4
    useSlideFormat = false;
elseif useSlideFormat
    figStub = '_ForSlideShow';
end
% Handle sensitivity case
if isfield(results.LFLk,'HighCred')
    figStubs = {'ImpCred';'HighCred';'LowCred'};
    experimentFields = {'ImpCred';'HighCred';'LowCred'};
else
    figStubs = {'ImpCred'};
    experimentFields = {'ImpCred'};
end
nExperiments = numel(experimentFields);

%% ENSURE FIGURES DIRECTORY EXISTS
if saveTeXfigures
    if ~exist(TeXfiguresDir,'dir')
        mkdir(TeXfiguresDir);
    end
end

%% PLOT EXPOSITION CHART
nLFLvarsToPlot = numel(plotVarMnems);
baseSimFigHandle = figure;
for iVar = 1:nLFLvarsToPlot
    iMnem = plotVarMnems{iVar};
    [iVarType,iVarIndex] = lookup_LSS_model_variable_type_and_index_number(results.Model,iMnem);
    iBasePath = results.Baseline.Paths.(iVarType)(iVarIndex,1:plotHorizon,end);
    iFullCred = results.LFLk.FullCred.Paths.(iVarType)(iVarIndex,1:plotHorizon,end);
    iFullCredK1 = results.LFLkPlus1.FullCred.Paths.(iVarType)(iVarIndex,1:plotHorizon,end);
    subplot(1,3,iVar); hold on;
    plot(1:plotHorizon,iBasePath,'color',0.5*ones(3,1),'linewidth',2.5);
    plot(1:plotHorizon,iFullCred,'-bs','linewidth',1);
    plot(1:plotHorizon,iFullCredK1,'-ro','linewidth',1);
    plot(1:plotHorizon,zeros(1,plotHorizon),'k','linewidth',0.5);
    if iVar==nLFLvarsToPlot
        legend({'Baseline', ...
            ['Liftoff in period ' num2str(results.nLFLperiods+1)], ...
            ['Liftoff in period ' num2str(results.nLFLperiods+2)]},'Location','NorthWest');
        legend boxoff;
    end
    xlim([1 plotHorizon]); ylim(yLims{iVar});
    title(plotVarNames{iVar});
end

%% PLOT LOSS COMPARISON CHART
fullCredLossCompHandle = figure;
subplot(1,2,1); hold on;
plot(1:results.nLFLperiods+1,[results.Baseline.Losses(1:results.nLFLperiods) nan],'color',0.5*ones(3,1),'linewidth',2.5);
plot(1:results.nLFLperiods+1,[results.LFLk.FullCred.Losses(1,1:results.nLFLperiods,end) nan],'-bs','linewidth',1);
plot(1:results.nLFLperiods+1,results.LFLkPlus1.FullCred.Losses(1,1:results.nLFLperiods+1,end),'-ro','linewidth',1);
xlim([1 results.nLFLperiods+1]);
legend({'Baseline', ...
    ['Liftoff in period ' num2str(results.nLFLperiods+1)], ...
    ['Liftoff in period ' num2str(results.nLFLperiods+2)]},'Location','NorthEast'); legend boxoff;
title('Output-equivalent losses');
subplot(1,2,2); hold on;
net1 = -[results.Baseline.Losses(1:results.nLFLperiods) nan] + [results.LFLk.FullCred.Losses(1,1:results.nLFLperiods,end) nan];
net2 = -[results.Baseline.Losses(1:results.nLFLperiods) nan] + results.LFLkPlus1.FullCred.Losses(1,1:results.nLFLperiods+1,end);
plot(1:results.nLFLperiods+1,net1,'-bs','linewidth',1);
plot(1:results.nLFLperiods+1,net2,'-ro','linewidth',1);
plot(1:results.nLFLperiods+1,zeros(1,results.nLFLperiods+1),'k','linewidth',0.5);
xlim([1 results.nLFLperiods+1]); legend({ ...
    ['K minus baseline'],['K+1 minus baseline']},'Location','SouthEast'); legend boxoff;
title('Net losses');

%% PLOT IMPERFECT CREDIBILITY SIMS (3x3 PANELS)
figHandles = cell(nExperiments,1);
for iExpt = 1:nExperiments
    iField = experimentFields{iExpt};
    figHandles{iExpt} = figure;
    for iVar = 1:nLFLvarsToPlot
        [iVarType,iVarIndex] = lookup_LSS_model_variable_type_and_index_number(results.Model,plotVarMnems{iVar});
        base = results.Baseline.Paths.(iVarType)(iVarIndex,1:plotHorizon,end);
        fullK = results.LFLk.FullCred.Paths.(iVarType)(iVarIndex,1:plotHorizon,end);
        impK = results.LFLk.(iField).ExpectedPaths.(iVarType)(iVarIndex,1:plotHorizon,1);
        fullK1 = results.LFLkPlus1.FullCred.Paths.(iVarType)(iVarIndex,1:plotHorizon,end);
        impK1 = results.LFLkPlus1.(iField).ExpectedPaths.(iVarType)(iVarIndex,1:plotHorizon,1);
        subplot(3,3,iVar); hold on;
        plot(1:plotHorizon,base,'color',0.5*ones(3,1),'linewidth',2.5);
        plot(1:plotHorizon,fullK,'--','color',0.5*ones(3,1),'linewidth',2.5);
        plot(1:plotHorizon,impK,'-bs','linewidth',1);
        plot(1:plotHorizon,zeros(1,plotHorizon),'k','linewidth',0.5);
        xlim([1 plotHorizon]); ylim(yLims{iVar}); title(plotVarNames{iVar});
        subplot(3,3,iVar+3); hold on;
        plot(1:plotHorizon,base,'color',0.5*ones(3,1),'linewidth',2.5);
        plot(1:plotHorizon,fullK1,'--','color',0.5*ones(3,1),'linewidth',2.5);
        plot(1:plotHorizon,impK1,'-ro','linewidth',1);
        plot(1:plotHorizon,zeros(1,plotHorizon),'k','linewidth',0.5);
        xlim([1 plotHorizon]); ylim(yLims{iVar}); title(plotVarNames{iVar});
    end
    subplot(3,3,7); hold on;
    plot(results.LFLk.(iField).NetLosses,'-bs','linewidth',1);
    plot(results.LFLkPlus1.(iField).NetLosses,'-ro','linewidth',1);
    title('Output equivalent net benefits from reneging'); xlim([1 results.nLFLperiods+1]);
    subplot(3,3,8); hold on;
    plot(1:results.nLFLperiods+1,[results.LFLk.(iField).ExAnteExitProbs(1:results.nLFLperiods) 0],'-bs','linewidth',1);
    plot(1:results.nLFLperiods+1,results.LFLkPlus1.(iField).ExAnteExitProbs,'-ro','linewidth',1);
    xlim([1 results.nLFLperiods+1]); title('Ex ante renege probabilities');
    subplot(3,3,9); hold on;
    plot(1:results.nLFLperiods+2,[results.LFLk.(iField).SumLiftOffDateProbs nan],'-bs','linewidth',1);
    plot(1:results.nLFLperiods+2,results.LFLkPlus1.(iField).SumLiftOffDateProbs,'-ro','linewidth',1);
    xlim([1 results.nLFLperiods+2]); ylim([0 1]); title('Liftoff date probabilities');
    legend({'Liftoff K','Liftoff K+1'},'Location','NorthWest'); legend boxoff;
end

%% SAVE CHARTS
if saveTeXfigures
    % Baseline
    figName = 'baseline';
    figFullPathName = fullfile(TeXfiguresDir,[figName,figStub,'.eps']);
    set(baseSimFigHandle,'PaperPositionMode','auto');
    set(baseSimFigHandle,'PaperUnits',get(baseSimFigHandle,'Units'));
    pp = get(baseSimFigHandle,'PaperPosition'); set(baseSimFigHandle,'PaperSize',pp(3:4));
    print(baseSimFigHandle,'-depsc',figFullPathName);
    % Loss comparison
    figName = 'fullCredLossComparison';
    figFullPathName = fullfile(TeXfiguresDir,[figName,figStub,'.eps']);
    set(fullCredLossCompHandle,'PaperPositionMode','auto');
    set(fullCredLossCompHandle,'PaperUnits',get(fullCredLossCompHandle,'Units'));
    pp = get(fullCredLossCompHandle,'PaperPosition'); set(fullCredLossCompHandle,'PaperSize',pp(3:4));
    print(fullCredLossCompHandle,'-depsc',figFullPathName);
    % Imperfect credibility
    for iCase=1:nExperiments
        figName=['LFLCompareHorizon',figStubs{iCase}];
        figFullPathName=fullfile(TeXfiguresDir,[figName,figStub,'.eps']);
        h=figHandles{iCase};
        set(h,'PaperPositionMode','auto');
        set(h,'PaperUnits',get(h,'Units')); pp=get(h,'PaperPosition'); set(h,'PaperSize',pp(3:4));
        if useSlideFormat, set(h,'PaperPosition',[0 0 16 8]); else, set(h,'PaperPosition',[0 0 10 8]); end
        print(h,'-depsc',figFullPathName);
    end
end
end
