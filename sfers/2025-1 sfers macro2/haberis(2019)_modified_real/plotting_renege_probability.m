%% plotting_renege_probability.m
% 
% Creates a plot comparing the probability of reneging for simulations
% with and without financial frictions.
% 
% Based on plotFigureB2.m from the original authors.
% 

%% 1. PRELIMINARIES

% Load results from simulations
load withoutfinancial.mat;
SimOutput_noFin = SimOutput;
clear SimOutput;

load withfinancial.mat;
SimOutput_Fin = SimOutput;
clear SimOutput;

% Define plot settings
recession_shade_color = [0.85,0.85,0.85];
start_plot_date = 2009.0;
end_plot_date = 2019.75;
fig_size = [2.5,2.5,6,4.5];
font_size = 10;
line_width = 1.5;
legend_location = 'NorthEast';

% Create date strings for x-axis
[qtr_dates,qtr_labels] = ...
  create_quarterly_date_strings(start_plot_date,end_plot_date,4);


%% 2. CREATE PLOT

figure;

% Plot renege probabilities
plot(qtr_dates,SimOutput_noFin.EXANTE.renegeProbs(1,:),"b-",...
     qtr_dates,SimOutput_Fin.EXANTE.renegeProbs(1,:),"r--",...
     'LineWidth',line_width);

% Add labels and title
xlabel('Year','FontSize',font_size);
ylabel('Probability','FontSize',font_size);
title('Probability of Reneging','FontSize',font_size);

% Set axis limits and ticks
set(gca,'XLim',[start_plot_date-0.25,end_plot_date],
    'XTick',qtr_dates(1:4:end),
    'XTickLabel',qtr_labels(1:4:end),
    'FontSize',font_size);
ylim([0 1]);

% Add legend
legend({'Without Financial Frictions','With Financial Frictions'},
       'Location',legend_location);
legend boxoff;

% Save the figure
fig_name = 'compare_renege_prob.pdf';
set(gcf,'PaperUnits','inches','PaperPosition',fig_size);
saveas(gcf,['Figures/',fig_name]);

fprintf('Figure saved as %s\n', fig_name);