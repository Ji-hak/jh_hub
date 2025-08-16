%% THIS SCRIPT PRODUCES FIGURES FOR SECTION 4 OF UNCERTAIN POLICY PROMISES
% This script loads results for the FRBUS TBFG application section 4
% of the paper, plots charts and saves figures for TeX.

%% HOUSEKEEPING
close all;
%restoredefaultpath;
clear variables;
delete('*.asv');
clc; 

%% ADD DIRECTORIES
fullPathNameForThisFile = mfilename('fullpath');
fullPathFileNameSplitByFinalBackSlash = regexp(...
    fullPathNameForThisFile,'[^\\]+$','split');
thisDir = fullPathFileNameSplitByFinalBackSlash{1};
addpath(genpath(thisDir));

%% SET OPTIONS FOR SCRIPT
saveTexFigures = true;
formatFigsForSlides = false;

%% DEFINE CHARTS TO PRODUCE
% 1st column is name of figure
% Cell in 2nd column is info about first experiment 
%  -- 1st column is file name of 1st set of results to plot
%  -- 2nd column is legend name of 1st set of results
%  -- 3rd column is line plot style
% Cell in 3rd column is info about second experiment with same format
chartsToMake = {...
    'TBFGusingFRBUSuThreshold5point75triggerVersus6opendEnded',...
        {'TBFGusingFRBUSuThreshold5point75trigger',...
        '$\bar{u} = 5.75\%$ (``Trigger'''')','r--'},...
        {'TBFGusingFRBUSuThreshold6openEnded',...
        '$\bar{u} = 6.0\%$ (``Open-ended'''')','g-.'};
    'TBFGusingFRBUSuThreshold5point75triggerVersus5point75opendEnded',...
        {'TBFGusingFRBUSuThreshold5point75trigger',...
        '$\bar{u} = 5.75\%$ (``Trigger'''')','r--'},...
        {'TBFGusingFRBUSuThreshold5point75openEnded',...
        '$\bar{u} = 5.75\%$ (``Open-ended'''')','g-.'};   
    'TBFGusingFRBUSuThreshold6triggerVersus6opendEnded',...
        {'TBFGusingFRBUSuThreshold6trigger',...
        '$\bar{u} = 6.0\%$ (``Trigger'''')','r--'},...
        {'TBFGusingFRBUSuThreshold6openEnded',...
        '$\bar{u} = 6.0\%$ (``Open-ended'''')','g-.'}; 
        };    

%% CHARTING OPTIONS
plotHorizon = 20;
varsToPlot = ...
    {'annualPCEinf';'core';'annGDPgrowth';'fedFunds';'urate'};
varNames = {'Annual PCE inflation, %';...
    'Annual core PCE inflation, %';...
    'Annual GDP growth, %';...
    'Federal funds rate, %';...
    'Unemployment rate, %'};
addVerticalLine = true;
addZeroLines = [false;false;false;true;false;true];
nRows = 2;
nBackDataPeriods = 4;
tickMarkGap = 4;
plotToAddLegendTo = 'annualPCEinf';
legendFontSize = 9;
legendLocation = 'NorthEast';
dataStartYear = 1985;
dataStartQuarter = 1;

%% DEFINE DIRECTORIES
TeXfiguresDir = [thisDir,'Figures\'];
resultsDir = [thisDir,'Results\'];

%% MAKE CHARTS
nChartsToMake = size(chartsToMake,1);
for iChart = 1:nChartsToMake
    % Extract info about the two policies to plot and load results
    figName = chartsToMake{iChart,1};
    firstPolicyResultsFileName = [...
        resultsDir,chartsToMake{iChart,2}{1},'.mat'];
    firstPolicyLegendName = chartsToMake{iChart,2}{2};
    firstPolicyLineStyle = chartsToMake{iChart,2}{3};
    secondPolicyResultsFileName = [...
        resultsDir,chartsToMake{iChart,3}{1},'.mat'];
    secondPolicyLegendName = chartsToMake{iChart,3}{2};
    secondPolicyLineStyle = chartsToMake{iChart,3}{3};
    [BaselineForecast,Model,FirstPolicyExpectedPaths,...
        firstPolicyPathProbabilities,firstPolicyUnempThreshold] = ...
        load_variables_from_mat_file(firstPolicyResultsFileName,...
        'baselineForecast','Model','TBFGexpectedPaths','TBFGpathProbs',...
        'uThresholdPct');
    [SecondPolicyExpectedPaths,secondPolicyPathProbabilities,...
        secondPolicyUnempThreshold] = ...
        load_variables_from_mat_file(secondPolicyResultsFileName,...
        'TBFGexpectedPaths','TBFGpathProbs','uThresholdPct');
    % Compute dimesnions of data and plot
    nPastPeriods = size(BaselineForecast.Past.rawObservables,2);
    nForecastPeriods = ...
        size(BaselineForecast.Forecast.rawObservables,2);
    nVarsToPlot = size(varsToPlot,1);
    nCols = ceil(nVarsToPlot/nRows);
    % Axis data handling
    dateStringVector = ...
        create_quarterly_date_strings(dataStartYear,dataStartQuarter,...
        nPastPeriods+nForecastPeriods);
    dateStringVector(:,1:2)=[]; % Strip out '20's
    plotDateStrings = ...
        dateStringVector(end-nBackDataPeriods-nForecastPeriods+1:end,:);
    nPlotPeriods = nBackDataPeriods + plotHorizon;
    plotDateIndices = 1:nPlotPeriods;
    datesTicksToInclude = (1:tickMarkGap:nPlotPeriods);
    % Lookupm info about variables to plot in model
    varInfo = lookup_LSS_mnemonic_type_and_position(Model,varsToPlot);
    FigHandle = figure;
    for iVar = 1:nVarsToPlot
        subplot(nRows,nCols,iVar);
        % Plot baseline projection
        iBaseBackData = ...
            BaselineForecast.Past.(varInfo{iVar,2})(varInfo{iVar,3},:);
        iBaseForecastData = ...
            BaselineForecast.Forecast.(varInfo{iVar,2})(varInfo{iVar,3},:);
        iBaseBackData = iBaseBackData(1,end-nBackDataPeriods+1:end);
        nanBackData = [nan(1,length(iBaseBackData)-1) iBaseBackData(end)];
        iDataToPlot = [iBaseBackData iBaseForecastData(1:plotHorizon)];
        plot(plotDateIndices,iDataToPlot,'color',0.5*ones(1,3),...
            'linewidth',2.5);
        hold on;
        % Plot expected paths from 1st policy
        iFirstPolicyExpectation = FirstPolicyExpectedPaths...
            .(varInfo{iVar,2})(varInfo{iVar,3},:,1);
        iSecondPolicyExpectationToPlot = ...
            [nanBackData iFirstPolicyExpectation(1:plotHorizon)];
        plot(plotDateIndices,iSecondPolicyExpectationToPlot,...
            firstPolicyLineStyle,'linewidth',1.5);        
        % Plot empirical TBFG
        iSecondPolicyExpectation = SecondPolicyExpectedPaths...
            .(varInfo{iVar,2})(varInfo{iVar,3},:,1);
        iSecondPolicyExpectationToPlot = ...
            [nanBackData iSecondPolicyExpectation(1:plotHorizon)];
        plot(plotDateIndices,iSecondPolicyExpectationToPlot,...
            secondPolicyLineStyle,'linewidth',1.5);       
        % Add lines for threshold levels (if plot is of unemployment rate)
        if strcmpi(varsToPlot{iVar},'urate')
            plot(plotDateIndices,...
                0*plotDateIndices+firstPolicyUnempThreshold,...
                firstPolicyLineStyle,'linewidth',1.5);
            plot(plotDateIndices,...
                0*plotDateIndices+secondPolicyUnempThreshold,...
                secondPolicyLineStyle,'linewidth',1.5);
        end
        % Add zero line
        if addZeroLines(iVar)
            plot(plotDateIndices,0*plotDateIndices,'k','linewidth',1);
        end
        % Sort out dates
        AxisHandle = get(FigHandle,'CurrentAxes');
        set(AxisHandle,'xLim',[1 nPlotPeriods]);
        set(AxisHandle,'XTick',datesTicksToInclude);
        set(AxisHandle,'xTickLabel',plotDateStrings(datesTicksToInclude,:));
        % Add vertical line
        iAxis = axis;
        plot([nBackDataPeriods nBackDataPeriods],[iAxis(3) iAxis(4)],...
            'k','linewidth',0.5);
        set(AxisHandle,'yLim',[iAxis(3) iAxis(4)]);
        % Add title
        title(varNames{iVar});
        % Add legend to plot
        if strcmpi(varsToPlot{iVar},plotToAddLegendTo)
            TBFGlegendStrings = {...
                'Baseline forecast',...
                firstPolicyLegendName,...
                secondPolicyLegendName};           
            legend(TBFGlegendStrings,'Location',legendLocation,...
                'Fontsize',legendFontSize,'interpreter','latex');
            legend boxoff;
        end
    end
    % Plot probabilities
    nProbPers = size(firstPolicyPathProbabilities,3);
    firstPolicyPathProbabilitiesToPlot = ...
        reshape(firstPolicyPathProbabilities(1,1,:),1,nProbPers);
    secondPolicyPathProbabilitiesToPlot = ...
        reshape(secondPolicyPathProbabilities(1,1,:),1,nProbPers);
    subplot(nRows,nCols,nVarsToPlot+1);
    hold on;
    plot(1:nProbPers,firstPolicyPathProbabilitiesToPlot,...
        [firstPolicyLineStyle,'d'],'linewidth',1.5);
    plot(1:nProbPers,secondPolicyPathProbabilitiesToPlot,...
        [secondPolicyLineStyle,'s'],'linewidth',1.5);
    endTick = nProbPers;
    probTicks = 1:tickMarkGap:endTick;
    probDateStrings = plotDateStrings(nBackDataPeriods+1:endTick + ...
        nBackDataPeriods,:);
    set(gca,'XTick',probTicks,'XTickLabel',probDateStrings(probTicks,:),...
        'XLim',[1 endTick]);
    title('Liftoff probabilities');  
    % Save figure
    if saveTexFigures
        if formatFigsForSlides
            figDimensions = [0 0 10 5];                                     %#ok<UNRCH>
        else
            figDimensions = [0 0 10 6];
        end
        set(FigHandle,'PaperUnits','inches','PaperPosition',figDimensions);
        print(FigHandle,'-depsc',[TeXfiguresDir,figName,'.eps']);
    end
end