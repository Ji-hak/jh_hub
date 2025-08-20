%% plotting_liftoff_probability.m
% Creates a plot comparing the liftoff date probabilities for simulations
% with and without financial frictions.

%% 1. PRELIMINARIES
clear; clc;

% Load results from simulations
load withoutfinancial.mat;    % must contain SimOutput.EXANTE.dateNums & liftoffProbs
SimOutput_noFin = SimOutput;  clear SimOutput;

load withfinancial.mat;
SimOutput_Fin = SimOutput;    clear SimOutput;

% Define plot settings (years as decimals; e.g., 2013.00, 2013.25, ...)
start_plot_date = 2013.00;
end_plot_date   = 2017.75;

fig_size   = [2.5, 2.5, 6, 4.5]; % [left bottom width height] in inches
font_size  = 10;
line_width = 1.5;
legend_location = 'northeast'; % case-insensitive

% Ensure identical date grids across the two simulations (safer than ==)
dates_noFin = SimOutput_noFin.EXANTE.dateNums(:);
dates_Fin   = SimOutput_Fin.EXANTE.dateNums(:);

% 보통 두 결과의 날짜축이 동일해야 합니다. 혹시 다르면 교집합을 취합니다.
common_dates = intersect(dates_noFin, dates_Fin, 'stable');

% 관심 구간 필터 (부동소수 허용오차 포함)
epsTol = 1e-10;
in_window = (common_dates >= start_plot_date - epsTol) & (common_dates <= end_plot_date + epsTol);
plot_dates = common_dates(in_window);

if isempty(plot_dates)
    error('선택한 기간에 해당하는 날짜가 없습니다. start/end_plot_date 값을 점검하세요.');
end

% 각 시뮬레이션의 liftoffProbs를 공통 날짜에 맞춰 정렬/슬라이스
% (liftoffProbs 가 T×1 또는 1×T라고 가정; 차원 안전화)
% 먼저 noFin 쪽 인덱스 매핑
[~, idx_noFin] = ismember(plot_dates, dates_noFin);
[~, idx_Fin]   = ismember(plot_dates, dates_Fin);

if any(idx_noFin==0) || any(idx_Fin==0)
    error('공통 날짜 매핑에 실패했습니다. dateNums가 서로 완전히 일치하지 않습니다.');
end

y_noFin_full = SimOutput_noFin.EXANTE.liftoffProbs;
y_Fin_full   = SimOutput_Fin.EXANTE.liftoffProbs;

% liftoffProbs가 행벡터/열벡터 여부와 1행 선택 가정에 대비
y_noFin = y_noFin_full(:);
y_Fin   = y_Fin_full(:);

% 혹시 liftoffProbs가 2차원(예: Nseries×T)이면 1행만 사용하려 했다면 아래처럼 선택:
% y_noFin = squeeze(y_noFin_full(1, :)).';
% y_Fin   = squeeze(y_Fin_full(1, :)).';

y_noFin = y_noFin(idx_noFin);
y_Fin   = y_Fin(idx_Fin);

% x축 라벨(연/분기) 생성
[qtr_dates, qtr_labels] = create_quarterly_date_strings(plot_dates(1), plot_dates(end), 4);

% qtr_dates는 분기 간격의 숫자형 벡터여야 하며 plot_dates와 길이가 같도록 맞춥니다.
% 간혹 create_quarterly_date_strings가 양끝 포함 규칙이 달라 길이가 어긋나면 보간/재생성 필요
if numel(qtr_dates) ~= numel(plot_dates)
    % qtr_dates를 plot_dates로 대체 (라벨만 재활용)
    qtr_dates = plot_dates;
    % qtr_labels 길이도 맞춥니다 (간단히 연도.분기 형식으로 다시 만듦)
    qtr_labels = arrayfun(@(x) sprintf('%dQ%d', floor(x), round(mod(x,1)*4)+1), qtr_dates, 'UniformOutput', false);
end

%% 2. CREATE PLOT
f = figure('Units','inches','Position',fig_size);

p1 = plot(qtr_dates, y_noFin, 'b-', 'LineWidth', line_width, 'DisplayName','Without Financial Frictions'); hold on;
p2 = plot(qtr_dates, y_Fin,   'r--', 'LineWidth', line_width, 'DisplayName','With Financial Frictions');

xlabel('Year',       'FontSize', font_size);
ylabel('Probability','FontSize', font_size);
title('Liftoff Date Probabilities','FontSize', font_size);

% X축 눈금: 1년 간격(=4분기 간격)으로 라벨링
set(gca, 'XLim', [plot_dates(1)-0.25, plot_dates(end)], ...
         'XTick', qtr_dates(1:4:end), ...
         'XTickLabel', qtr_labels(1:4:end), ...
         'FontSize', font_size);

legend('Location', legend_location); legend boxoff;

% 저장 경로 보장
outdir = 'Figures';
if ~exist(outdir,'dir'); mkdir(outdir); end

% PDF로 저장 (exportgraphics 권장)
fig_name = 'compare_liftoff_prob.pdf';
exportgraphics(gca, fullfile(outdir, fig_name), 'ContentType','vector', 'BackgroundColor','white');

fprintf('Figure saved as %s\n', fullfile(outdir, fig_name));
