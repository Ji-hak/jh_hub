%% COMPARE Endogenous Credibility (No-FS vs +FS) — with legend fix & peg shading
% This script compares endogenous-credibility expected paths from two
% result files:
%   - withoutfinancial (no FS terms in rule/loss)
%   - withfinancial    (+ FS terms)
% and plots three panels (inflation, output gap, policy rate).
% MODS vs previous canvas version:
%   (1) Legend now binds to LINE HANDLES explicitly (fixes color mismatch)
%   (2) Peg-shading patch excluded from legend via 'HandleVisibility','off'
%   (3) Shading uses peg horizon from the loaded results

%% HOUSEKEEPING
close all; clear variables; delete('*.asv'); clc;

%% PATHS (edit if needed)
file_noFS = [thisDir,'Results\withoutfinancial.mat'];
file_FS   = [thisDir,'Results\withfinancial.mat'];

%% ADD DIRECTORIES
fullPathNameForThisFile = mfilename('fullpath');
fullPathFileNameSplitByFinalBackSlash = regexp(fullPathNameForThisFile,'[^\]+$','split');
thisDir = fullPathFileNameSplitByFinalBackSlash{1};
addpath(genpath(thisDir));
addpath(genpath('Models')); addpath(genpath('Functions'));
addpath(genpath('Figures')); addpath(genpath('Data'));
addpath(genpath('MAPSlite')); addpath(genpath('Results'));

%% OPTIONS
plotHorizon   = 14;        % first H quarters
saveTeXfigures = true;
useSlideFormat = false; %#ok<NASGU>

%% FIGURES OUT DIR
thisDir       = fileparts(mfilename('fullpath'));
TeXfiguresDir = fullfile(thisDir, 'Figures');
if ~exist(TeXfiguresDir, 'dir'); mkdir(TeXfiguresDir); end

%% LOAD RESULTS
S1 = load(file_noFS, 'results'); results_noFS = S1.results;
S2 = load(file_FS,   'results'); results_FS   = S2.results;

%% INDICES (assume same mnemonics order across models)
rAnnIndex = unpack_model_metadata_and_lookup_index_numbers(results_noFS.Model,'Ymnems','rAnn');
piAnnIndex = unpack_model_metadata_and_lookup_index_numbers(results_noFS.Model,'Ymnems','piAnn');
yIndex    = unpack_model_metadata_and_lookup_index_numbers(results_noFS.Model,'xMnems','y');

%% Extract endogenous-credibility expected paths from each results file
% Prefer LFLk (liftoff in period k); if absent, fall back to LFLkPlus1
if isfield(results_noFS,'LFLk') && isfield(results_noFS.LFLk,'ImpCred')
    Exp_noFS = results_noFS.LFLk.ImpCred.ExpectedPaths;  % K
    pegEnd_noFS = results_noFS.nLFLperiods;
elseif isfield(results_noFS,'LFLkPlus1') && isfield(results_noFS.LFLkPlus1,'ImpCred')
    Exp_noFS = results_noFS.LFLkPlus1.ImpCred.ExpectedPaths; % K+1
    pegEnd_noFS = results_noFS.nLFLperiods + 1;
else
    error('No endogenous-credibility paths found in withoutfinancial file.');
end

if isfield(results_FS,'LFLk') && isfield(results_FS.LFLk,'ImpCred')
    Exp_FS = results_FS.LFLk.ImpCred.ExpectedPaths;   % K
    pegEnd_FS = results_FS.nLFLperiods;
elseif isfield(results_FS,'LFLkPlus1') && isfield(results_FS.LFLkPlus1,'ImpCred')
    Exp_FS = results_FS.LFLkPlus1.ImpCred.ExpectedPaths; % K+1
    pegEnd_FS = results_FS.nLFLperiods + 1;
else
    error('No endogenous-credibility paths found in withfinancial file.');
end

% Choose a common shading horizon (use no-FS as reference)
pegEnd = pegEnd_noFS;

%% Series (Horizon crop)
H = plotHorizon; t = 1:H;

pi_noFS = squeeze(Exp_noFS.modelObservables(piAnnIndex,t,1));
pi_FS   = squeeze(Exp_FS.modelObservables(piAnnIndex,t,1));

y_noFS  = squeeze(Exp_noFS.modelVariables(yIndex,t,1));
y_FS    = squeeze(Exp_FS.modelVariables(yIndex,t,1));

r_noFS  = squeeze(Exp_noFS.modelObservables(rAnnIndex,t,1));
r_FS    = squeeze(Exp_FS.modelObservables(rAnnIndex,t,1));

%% Colors & styles
coBlue = [0.00 0.45 0.74];      % no-FS
coRed  = [0.85 0.33 0.10];      % +FS
lw     = 1.9; ms = 5;

fig = figure('Color','w','Position',[100 100 980 360]);

% ================= (1) INFLATION =================
subplot(1,3,1); hold on; box on;
% Lines first (so ylim is data-driven)
hNoFS_pi = plot(t, pi_noFS, 's-', 'Color', coBlue, 'LineWidth', lw, 'MarkerSize', ms);
hFS_pi   = plot(t, pi_FS,   'o-', 'Color', coRed,  'LineWidth', lw, 'MarkerSize', ms);
yline(0,'k-','LineWidth',0.5);
% Shading for peg [1..pegEnd], excluded from legend
yl = ylim; ph = patch([1 pegEnd pegEnd 1], [yl(1) yl(1) yl(2) yl(2)], [0.5 0.5 0.5], ...
    'FaceAlpha',0.08, 'EdgeColor','none', 'HandleVisibility','off'); %#ok<NASGU>
% Put patch behind lines
uistack(hNoFS_pi,'top'); uistack(hFS_pi,'top');
% Legend with explicit handles (legend color fix)
legend([hNoFS_pi hFS_pi], {'Endog cred (no FS)','Endog cred (+FS)'}, 'Location','southeast');

title('Inflation (annualised, %)');

% ================= (2) OUTPUT GAP =================
subplot(1,3,2); hold on; box on;
hNoFS_y = plot(t, y_noFS, 's-', 'Color', coBlue, 'LineWidth', lw, 'MarkerSize', ms);
hFS_y   = plot(t, y_FS,   'o-', 'Color', coRed,  'LineWidth', lw, 'MarkerSize', ms);
yline(0,'k-','LineWidth',0.5);
yl = ylim; ph = patch([1 pegEnd pegEnd 1], [yl(1) yl(1) yl(2) yl(2)], [0.5 0.5 0.5], ...
    'FaceAlpha',0.08, 'EdgeColor','none', 'HandleVisibility','off'); %#ok<NASGU>
uistack(hNoFS_y,'top'); uistack(hFS_y,'top');
legend([hNoFS_y hFS_y], {'Endog cred (no FS)','Endog cred (+FS)'}, 'Location','southeast');

title('Output gap (%)');

% ================= (3) POLICY RATE =================
subplot(1,3,3); hold on; box on;
hNoFS_r = plot(t, r_noFS, 's-', 'Color', coBlue, 'LineWidth', lw, 'MarkerSize', ms);
hFS_r   = plot(t, r_FS,   'o-', 'Color', coRed,  'LineWidth', lw, 'MarkerSize', ms);
yline(0,'k-','LineWidth',0.5);
yl = ylim; ph = patch([1 pegEnd pegEnd 1], [yl(1) yl(1) yl(2) yl(2)], [0.5 0.5 0.5], ...
    'FaceAlpha',0.08, 'EdgeColor','none', 'HandleVisibility','off'); %#ok<NASGU>
uistack(hNoFS_r,'top'); uistack(hFS_r,'top');
legend([hNoFS_r hFS_r], {'Endog cred (no FS)','Endog cred (+FS)'}, 'Location','southeast');

title('Policy rate (annualised, %)');

% ===== SAVE =====
if saveTeXfigures
    outName = fullfile(TeXfiguresDir, 'compare_endogcred_noFS_vs_FS.pdf');
    exportgraphics(fig, outName, 'ContentType','vector');
end
