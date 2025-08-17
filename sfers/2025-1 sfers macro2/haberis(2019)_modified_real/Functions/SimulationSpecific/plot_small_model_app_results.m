function  plot_small_model_app_results(results,plotHorizon,...
    TeXfiguresDir,useSlideFormat)
% plot_small_model_app_results This function plots charts for sim results
% for the small model application section.
%   Inputs are a results structure, plot horizon scalar and folder name for
%   saving output. empty, no saving occurs. Third (optional) input controls
%   whether charts are saved in a format better suited for slideshow 
%   presentations.

%% SOME OPTIONS SET WITHIN THIS FUNCTION
plotVarMnems = {'pie';'y';'rAnn'};
plotVarNames = {'Quarterly inflation (%)';'Output gap (%)';...
    'Policy rate (annualised, %)'};
% Common y-axis limits for the imperfect credibility comparison
yLims = {[-0.5 0.15];[-7 1.25];[0 3]};
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
% Handle sensitivity case ... i.e. if sensitivity is in results set or not
if isfield(results.LFLk,'HighCred')
    figStubs = {'ImpCred';'HighCred';'LowCred'};
    experimentFields = {'ImpCred';'HighCred';'LowCred'};
else
    figStubs = {'ImpCred'};
    experimentFields = {'ImpCred'};
end
nExperiments = size(experimentFields,1);

%% EXPOSITION CHART (1x3 PANELS) ILLUSTRATING TEMPTATION TO RENEGE
nLFLvarsToPlot = size(plotVarMnems,1);
baseSimFigHandle = figure;
baseYlims = {[-0.5 0.15];[-7 1.25];[0 3]};
for iVar = 1:nLFLvarsToPlot
    iMnem = plotVarMnems{iVar};
    [iVarType,iVarIndex] = ...
        lookup_LSS_model_variable_type_and_index_number(results.Model,iMnem);
    % Extract the paths
    iBasePath = ...
        results.Baseline.Paths.(iVarType)(iVarIndex,1:plotHorizon,end);
    % Extract perfect credibility paths too
    iFullCredPath = ...
        results.LFLk.FullCred.Paths.(iVarType)(iVarIndex,...
        1:plotHorizon,end);
    % ... K+1 period case, full cred ...
    iFullCredPathKplus1 = ...
        results.LFLkPlus1.FullCred.Paths.(iVarType)(iVarIndex,...
        1:plotHorizon,end);
    % Plot the paths
    subplot(1,3,iVar);
    hold on;
    plot(1:plotHorizon,iBasePath,'color',0.5*ones(3,1),'linewidth',2.5);
    plot(1:plotHorizon,iFullCredPath,'-bs','linewidth',1);%2.5);
    plot(1:plotHorizon,iFullCredPathKplus1,'-ro','linewidth',1);%,2.5);
    plot(1:plotHorizon,zeros(1,plotHorizon),'k','linewidth',0.5);
    % Add legend (TBC)
    if iVar == nLFLvarsToPlot
        legendStrings = cell(3,1);
        legendStrings{1,1} = 'Baseline';
        legendStrings{2,1} = ['Liftoff in period ' ...
            num2str(results.nLFLperiods+1)];
        legendStrings{3,1} = ['Liftoff in period ' ...
            num2str(results.nLFLperiods+2)];
        legend(legendStrings,'Location','NorthWest');
        legend boxoff;
    end
    xlim([1 plotHorizon]);
    ylim(baseYlims{iVar});
    title(plotVarNames{iVar});
end

%% CHART FOR LOSS COMPARISON OF FULLY CREDIBLE LFL POLICIES
%  Updated to implement two panels (including loss differences)
fullCredLossCompHandle = figure;
subplot(1,2,1);
hold on;
plot(1:results.nLFLperiods+1,...
    [results.Baseline.Losses(1:results.nLFLperiods) nan],...
    'color',0.5*ones(3,1),'linewidth',2.5);
plot(1:results.nLFLperiods+1,...
    [results.LFLk.FullCred.Losses(1,1:results.nLFLperiods,end) nan],'-bs',...
    'linewidth',1);
plot(1:results.nLFLperiods+1,...
    results.LFLkPlus1.FullCred.Losses(1,1:results.nLFLperiods+1,end),'-ro',...
    'linewidth',1);
xlim([1 results.nLFLperiods+1]);
legendStrings = cell(3,1);
legendStrings{1,1} = 'Baseline';
legendStrings{2,1} = ['Liftoff in period ' ...
    num2str(results.nLFLperiods+1)];
legendStrings{3,1} = ['Liftoff in period ' ...
    num2str(results.nLFLperiods+2)];
legend(legendStrings,'Location','NorthEast');
legend boxoff;
title('Output-equivalent losses');
%
subplot(1,2,2);
hold on;
plot(1:results.nLFLperiods+1,...
    -[results.Baseline.Losses(1:results.nLFLperiods) nan]+ ...
    [results.LFLk.FullCred.Losses(1,1:results.nLFLperiods,end) nan],'-bs',...
    'linewidth',1);%2.5);
plot(1:results.nLFLperiods+1,...
    -[results.Baseline.Losses(1:results.nLFLperiods) nan]+ ...
    results.LFLkPlus1.FullCred.Losses(1,1:results.nLFLperiods+1,end),'-ro',...
    'linewidth',1);%2.5);
plot(1:results.nLFLperiods+1,zeros(1,results.nLFLperiods+1),'k',...
    'linewidth',0.5);
xlim([1 results.nLFLperiods+1]);
legendStrings = cell(2,1);
legendStrings{1,1} = ['Liftoff in period ' ...
    num2str(results.nLFLperiods+1) ' minus baseline'];
legendStrings{2,1} = ['Liftoff in period ' ...
    num2str(results.nLFLperiods+2) ' minus baseline'];
legend(legendStrings,'Location','SouthEast');
legend boxoff;
title('Net losses');

%% PLOT SIMULATIONS OF LFL POLICY UNDER ENDOGENOUS IMPERFECT CREDIBILITY
figHandles = cell(nExperiments,1);
for iExpt = 1:nExperiments
    % TBC collect figure handles in cell array to use later when saving
    iField = experimentFields{iExpt};
    % K+1 vs K period LFL experiment
    kPlusOneFigHandle = figure;
    figHandles{iExpt,1} = kPlusOneFigHandle;
    for iVar = 1:nLFLvarsToPlot
        iMnem = plotVarMnems{iVar};
        [iVarType,iVarIndex] = ...
            lookup_LSS_model_variable_type_and_index_number(...
            results.Model,iMnem);
        iYlim = yLims{iVar,1};
        % Extract the paths ...
        % ... baseline simulation ...
        iBasePath = ...
            results.Baseline.Paths.(iVarType)(iVarIndex,...
            1:plotHorizon,end);
        % .... K period case, full and imp cred ...
        iFullCredPath = ...
            results.LFLk.FullCred.Paths.(iVarType)(iVarIndex,...
            1:plotHorizon,end);
        iImpCredExpected = ...
            results.LFLk.(iField).ExpectedPaths.(iVarType)(iVarIndex,...
            1:plotHorizon,1);
        iImpCredExPostPath = ...
            results.LFLk.(iField).Paths.(iVarType)(iVarIndex,...
            1:plotHorizon,end);                                         %#ok<NASGU>
        % .... K+1 period case, full and imp cred ...
        iFullCredPathKplus1 = ...
            results.LFLkPlus1.FullCred.Paths.(iVarType)(iVarIndex,...
            1:plotHorizon,end);
        iImpCredExpectedKplus1 = ...
            results.LFLkPlus1.(iField).ExpectedPaths.(iVarType)(iVarIndex,...
            1:plotHorizon,1);
        iImpCredExPostPathKplus1 = ...
            results.LFLkPlus1.(iField).Paths.(iVarType)(iVarIndex,...
            1:plotHorizon,end);                                         %#ok<NASGU>
        % Plot the K period case
        subplot(3,3,iVar);
        hold on;
        plot(1:plotHorizon,iBasePath,...
            'color',0.5*ones(3,1),'linewidth',2.5);
        plot(1:plotHorizon,iFullCredPath,'--',...
            'color',0.5*ones(3,1),'linewidth',2.5);
        plot(1:plotHorizon,iImpCredExpected,'-bs','linewidth',1);
        if plotExPostPaths
            plot(1:plotHorizon,iImpCredExPostPath,'gx','linewidth',1, ...
                'MarkerSize',9);
        end
        plot(1:plotHorizon,zeros(1,plotHorizon),'k','linewidth',0.5);
        xlim([1 plotHorizon]);
        ylim(iYlim);
        title(plotVarNames{iVar});
        % Plot the K+1 period case
        subplot(3,3,iVar+3);
        hold on;
        plot(1:plotHorizon,iBasePath,...
            'color',0.5*ones(3,1),'linewidth',2.5);
        plot(1:plotHorizon,iFullCredPathKplus1,'--',...
            'color',0.5*ones(3,1),'linewidth',2.5);
        plot(1:plotHorizon,iImpCredExpectedKplus1,'-ro',...
            'linewidth',1);
        if plotExPostPaths
            plot(1:plotHorizon,iImpCredExPostPathKplus1,'gx',...
                'linewidth',1,'MarkerSize',9);
        end
        plot(1:plotHorizon,zeros(1,plotHorizon),'k','linewidth',0.5);
        xlim([1 plotHorizon]);
        ylim(iYlim);
        title(plotVarNames{iVar});
        % Add legend (TBC)
    end
    subplot(3,3,7);
    hold on;
    plot(results.LFLk.(iField).NetLosses,'-bs','linewidth',1);
    plot(results.LFLkPlus1.(iField).NetLosses,'-ro','linewidth',1);
    title('Output equivalent net benefits from reneging');
    xlim([1 results.nLFLperiods+1]);
    % Plot renege probabilities
    subplot(3,3,8);
    hold on;
    %plot(1:nLFLperiods,LFLexitProbs,'b--','linewidth',2);
    plot(1:results.nLFLperiods+1,...
        [results.LFLk.(iField).ExAnteExitProbs(1:results.nLFLperiods) 0],...
        '-bs','linewidth',1);
    plot(1:results.nLFLperiods+1,...
        results.LFLkPlus1.(iField).ExAnteExitProbs(1:results.nLFLperiods+1),...
        '-ro','linewidth',1);
    xlim([1 results.nLFLperiods+1]);
    title('Ex ante renege probabilities');
    % Plot liftoff distribution
    subplot(3,3,9);
    hold on;
    plot(1:results.nLFLperiods+2,...
        [results.LFLk.(iField).SumLiftOffDateProbs nan],'-bs',...
        'linewidth',1);
    plot(1:results.nLFLperiods+2,...
        results.LFLkPlus1.(iField).SumLiftOffDateProbs,'-ro',...
        'linewidth',1);
    xlim([1 results.nLFLperiods+2]);
    ylim([0 1]);
    title('Liftoff date probabilities');
    legendStrings = cell(2,1);
    legendStrings{1,1} = ['Liftoff in period ' num2str(results.nLFLperiods+1)];
    legendStrings{2,1} = ['Liftoff in period ' num2str(results.nLFLperiods+2)];
    legend(legendStrings,'Location','NorthWest');
    legend boxoff;
end

%% SAVE CHARTS
% Re-size and save to designated folder
if saveTeXfigures
    % Baseline simulation case
    figName = 'baseline';
    figFullPathName = [TeXfiguresDir,'\',figName,figStub,'.eps'];
    set(baseSimFigHandle,'PaperUnits','inches',...
        'PaperPosition',[0 0 10 2]);
    print(baseSimFigHandle,'-depsc',figFullPathName);
    % Exposition just showing loss comparison
    figName = 'fullCredLossComparison';
    figFullPathName = [TeXfiguresDir,'\',figName,figStub,'.eps'];
    set(fullCredLossCompHandle,'PaperUnits','inches',...
        'PaperPosition',[0 0 10 3]);
    print(fullCredLossCompHandle,'-depsc',figFullPathName);
    for iCase = 1:nExperiments
        iStub = figStubs{iCase};
        % K versus K+1 comparison chart
        figName = ['LFLCompareHorizon' iStub];
        figFullPathName = [TeXfiguresDir,'\',figName,figStub,'.eps'];
        iFigHandle = figHandles{iCase,1};
        if useSlideFormat
            set(iFigHandle,'PaperUnits','inches',...
                'PaperPosition',[0 0 16 8]);
        else
            set(iFigHandle,'PaperUnits','inches',...
                'PaperPosition',[0 0 10 8]);
        end
        print(iFigHandle,'-depsc',figFullPathName);
    end
end

end