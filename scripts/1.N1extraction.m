%% =====================================================================
%  GROUP‑LEVEL ERP ANALYSIS  &  N1 EXTRACTION
%  ---------------------------------------------------------------------
%  Author : Suong Welp
%  Date   : 2025-07-24
%  Loads all *_ERP.mat files created by preprocessing2ERP_final.m, 
%  builds grand‑average ERPs, extracts the N1 peak (50–150 ms) at Cz for 
%  every participant & condition, writes out tidy .csv/.mat tables, 
%  and produces Cz waveforms & topoplots.
%  =====================================================================

% ----------------------------- SET‑UP ----------------------------------
clc; close all; clear all;
addpath(genpath(fullfile(pwd,'scripts')));        % for chanlocs.mat, etc.

cfg          = struct();
cfg.erpDir   = fullfile(pwd,'data','ERP');
cfg.outCsv   = fullfile(cfg.erpDir,'Data_N1.csv');
cfg.outMat   = fullfile(cfg.erpDir,'Data_N1.mat');
cfg.chanCz   = 18;                 % Cz index 
cfg.winN1    = [50 150];           % ms
cfg.sRate    = 500;                % Hz – final sampling rate
cfg.tVec     = -250:1000/cfg.sRate:748;  % time‑vector (matches Section 1)

condLabels   = {'veridical_speak','stranger_match','stranger_listen', ...
                'stranger_mis','veridical_listen'};

% Pre‑allocate containers (subjects × conditions × time × channels)
maxSubs      = 60;                % safe upper‑bound for mem alloc
chanCount    = 59;                 % after eog and emg channel removal
ptCount      = numel(cfg.tVec);
C            = nan(chanCount, ptCount, maxSubs, numel(condLabels));
subIDs       = nan(1,maxSubs);

%% ---------------------------- LOAD FILES ------------------------------
files = dir(fullfile(cfg.erpDir,'P*_ERP.mat'));
subIdx = 0;
for f = 1:numel(files)
    sid = str2double(files(f).name(2:4));  % parse "P###"
    subIdx = subIdx + 1;  subIDs(subIdx) = sid;

    S = load(fullfile(files(f).folder, files(f).name), 'allERP');
    if ~isfield(S,'allERP'), warning('File %s has no allERP – skipping',files(f).name); continue; end
    for c = 1:numel(condLabels)
        C(:,:,subIdx,c) = S.allERP(:,:,c);
    end
end

% Trim unused pre‑allocation
C      = C(:,:,1:subIdx,:);
subIDs = subIDs(1:subIdx);
subCt  = numel(subIDs);

%% ----------------------- N1 PEAK EXTRACTION ---------------------------
N1   = nan(subCt, numel(condLabels));
idxN1= find(cfg.tVec>=cfg.winN1(1) & cfg.tVec<=cfg.winN1(2));

for c = 1:numel(condLabels)
    for s = 1:subCt
        trace      = squeeze(C(cfg.chanCz, idxN1, s, c));
        [N1(s,c),~] = min(trace);   % negative peak
    end
end

% ---------------------- WRITE RESULTS TABLE ---------------------------
T = array2table(N1,'VariableNames',condLabels,'RowNames',cellstr(string(subIDs.')));

writetable(T,cfg.outCsv,'WriteRowNames',true);
save(cfg.outMat,'N1','T','cfg');

fprintf('\nSaved N1 amplitudes for %d subjects → %s\n',subCt,cfg.outCsv);

%% ------------------------- GRAND AVERAGE ------------------------------
G = squeeze(mean(C,3,'omitnan'));  % channels × time × conditions
figure('Name','Grand‑Average Cz ERPs','Color','w','Position',[200 200 1000 500]);
colors = [ ...
    68   1  84 ;   % veridical speak   – purple
    33 145 140 ;   % stranger match    – teal
    50 120  70 ;   % stranger listen   – green
   133  94  66 ;   % stranger mismatch – brown
    59  82 139 ];  % veridical listen  – indigo
colors = colors/255;
lineStyles = {'-','-.','--',':','-'};

hold on;
text(-0.06, 1.07, 'a) N1 Event Related Potentials & Topography in different condition', 'Units', 'normalized', 'FontSize', 17);
for c = [1,5,2,4,3] %order condition to match labels
    plot(cfg.tVec, G(cfg.chanCz,:,c), ...
        'LineWidth',1.5, ...
        'Color', colors(c,:), ...
        'LineStyle', lineStyles{c});
end
legend({'veridical_speak', 'veridical_listen','stranger_match','stranger_mismatch','stranger_listen'}, 'Location', 'best','Interpreter','none');
hold off;lg=legend(); lg.AutoUpdate='off'; 
xlim([-250 750]); ylim([-12 10]);box on;
xlabel('Time (ms)'); ylabel('N1 Amplitude (\muV)');
line([0 0],ylim,'Color',[0 0 0 0.4],'LineStyle',':');  % stim onset
line([100 100],ylim,'Color',[0 0 0 0.4],'LineStyle','--'); % typical N1

erpC1 = G(:,:,1);erpC5 = G(:,:,5);

% Find peak latencies (indices) in N1 window
[~, lat1c] = min(erpC1(cfg.chanCz, idxN1));
[~, lat5c] = min(erpC5(cfg.chanCz, idxN1));
latIdx1c = idxN1(lat1c);
latIdx5c = idxN1(lat5c);
colorRange = [-10, 10];
% First topoplot
axes('Position', [0.5, 0.15, 0.15, 0.15]); % Adjust as needed
topoplot(erpC1(:, latIdx1c), chanlocs, 'electrodes', 'off');
caxis(colorRange);
title('N1 veridical speak');

% Second topoplot
axes('Position', [0.7, 0.15, 0.15, 0.15]); % Adjust as needed
topoplot(erpC5(:, latIdx5c), chanlocs, 'electrodes', 'off');
caxis(colorRange);
title('N1 veridical listen');

% Add a colorbar and set its title
cb = colorbar('Position', [0.92, 0.11, 0.02, 0.3]); % Adjust position as needed
cb.Title.String = 'N1 Amplitude [µV]';
cb.Title.FontSize = 7; % Adjust font size as needed

% Manually setting the colorbar ticks and tick labels
cb.Ticks = [-10, -5, 0, 5, 10];
cb.TickLabels = {'-10', '-5', '0', '5', '10'};

print(gcf, 'GrandERP', '-dsvg');

%% ------------------------- TOPOPLOTS ----------------------------------
load(fullfile(pwd,'scripts','chanlocs.mat'));  % supplies 'chanlocs'
clim = [-10 10];
idxN1pk = nan(numel(condLabels),1);
for c = 1:numel(condLabels)
    [~,lat]  = min(G(cfg.chanCz,idxN1,c));
    idxN1pk(c) = idxN1(lat);
end

figure('Name','N1 Topographies','Color','w');
for c = 1:numel(condLabels)
    subplot(2,3,c);
    topoplot(G(:,idxN1pk(c),c),chanlocs,'electrodes','off');
    caxis(clim); title(condLabels{c},'Interpreter','none');
end
colormap('parula');
c = colorbar('Position',[0.92 0.15 0.02 0.7]); c.Label.String = '\muV';



%% -------------------- SUPPLEMENT: SPQ GROUP ERPs Veridical ------------------
% The plots for the PDI cutoff can be produced the same way, with
% pdi_threshold = 5
% Load SPQ data and group indices
demo = readtable(fullfile(pwd, 'data', 'Dat0_Demo_SPQ_PDI.csv'));
spq_threshold = 22;
control_idx = demo.SPQ < spq_threshold;
high_idx    = demo.SPQ >= spq_threshold;

% Time vector and N1 time window
time = cfg.tVec;
idxN1 = find(time >= cfg.winN1(1) & time <= cfg.winN1(2));

% ERP data extraction function
erpCz = @(condIdx, groupIdx) mean(C(:,:,groupIdx,condIdx), 3);
erpC1con  = erpCz(1, control_idx); % veridical speak (control)
erpC5con  = erpCz(5, control_idx); % veridical listen (control)
erpC1high = erpCz(1, high_idx);    % veridical speak (high SPQ)
erpC5high = erpCz(5, high_idx);    % veridical listen (high SPQ)

% Plot setup
figure('Color','w','Position',[200 200 1000 400]);

% --- Subplot 1: Control group
subplot(1,2,1);hold on;text(-0.15, 1.05, 'a)', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 12);
h1=plot(time, erpC1con(cfg.chanCz,:), '--', 'LineWidth', 1.5, 'Color', [0.5 0 0.5]); 
h2=plot(time, erpC5con(cfg.chanCz,:), '-', 'LineWidth', 1.5, 'Color', [0 0.6 0.6]);  
legend([h1 h2], {'veridical_speak', 'veridical_listen'}, 'Location', 'best','Interpreter','none');
hold off;lg=legend(); lg.AutoUpdate='off'; 
xlim([-250 748]); ylim([-12 12]);
title('Control Group: N1 in Veridical Conditions');
xlabel('Time [ms]'); ylabel('Amplitude [\muV]');
line([0 0], ylim, 'Color', [0 0 0 0.3], 'LineStyle', ':');
line([100 100], ylim, 'Color', [0 0 0 0.3], 'LineStyle', '--');

% --- Subplot 2: High SPQ group
subplot(1,2,2);hold on;text(-0.15, 1.05, 'b)', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 12);
h3=plot(time, erpC1high(cfg.chanCz,:), '--', 'LineWidth', 1.5, 'Color', [0.5 0 0.5]); 
h4=plot(time, erpC5high(cfg.chanCz,:), '-', 'LineWidth', 1.5, 'Color', [0 0.6 0.6]);  
hold off;lg=legend(); lg.AutoUpdate='off'; 
legend([h3 h4], {'veridical_speak', 'veridical_listen'}, 'Location', 'best','Interpreter','none');
xlim([-250 748]); ylim([-12 12]);
title('High SPQ Group: N1 in Veridical Conditions');
xlabel('Time [ms]'); ylabel('Amplitude [\muV]');
legend([h3 h4], {'veridical speak', 'veridical listen'}, 'Location', 'best');
line([0 0], ylim, 'Color', [0 0 0 0.3], 'LineStyle', ':');
line([100 100], ylim, 'Color', [0 0 0 0.3], 'LineStyle', '--');

% ----------------------- TOPOPLOTS: Control group -----------------------
% Get subplot position 
pos = get(subplot(1,2,1), 'Position');
topoWidth = 0.08; topoHeight = 0.10;
topoBottom = pos(2) + 0.03;
topoLefts = [pos(1) + pos(3) - 2*topoWidth - 0.015, pos(1) + pos(3) - topoWidth];

% Find peak latencies (indices) in N1 window
[~, lat1c] = min(erpC1con(cfg.chanCz, idxN1));
[~, lat5c] = min(erpC5con(cfg.chanCz, idxN1));
latIdx1c = idxN1(lat1c);
latIdx5c = idxN1(lat5c);

% Topoplots for control group
axes('Position', [topoLefts(1), topoBottom, topoWidth, topoHeight]);
topoplot(erpC1con(:, latIdx1c), chanlocs, 'electrodes', 'off', 'headrad', 0.5);
caxis([-10 10]); title('N1 veridical speak', 'FontSize', 8);

axes('Position', [topoLefts(2), topoBottom, topoWidth, topoHeight]);
topoplot(erpC5con(:, latIdx5c), chanlocs, 'electrodes', 'off');
caxis([-10 10]); title('N1 veridical listen', 'FontSize', 8);

% ---------------------- TOPOPLOTS: High SPQ group -----------------------
pos2 = get(subplot(1,2,2), 'Position');
topoLefts2 = [pos2(1) + pos2(3) - 2*topoWidth - 0.015, pos2(1) + pos2(3) - topoWidth];

% Find peak latencies
[~, lat1h] = min(erpC1high(cfg.chanCz, idxN1));
[~, lat5h] = min(erpC5high(cfg.chanCz, idxN1));
latIdx1h = idxN1(lat1h);
latIdx5h = idxN1(lat5h);

% Topoplots for high SPQ group
axes('Position', [topoLefts2(1), topoBottom, topoWidth, topoHeight]);
topoplot(erpC1high(:, latIdx1h), chanlocs, 'electrodes', 'off', 'headrad', 0.5);
caxis([-10 10]); title('N1 veridical speak', 'FontSize', 8);

axes('Position', [topoLefts2(2), topoBottom, topoWidth, topoHeight]);
topoplot(erpC5high(:, latIdx5h), chanlocs, 'electrodes', 'off');
caxis([-10 10]); title('N1 veridical listen', 'FontSize', 8);

% ------------------------------ Colorbar -------------------------------
cb = colorbar('Position', [topoLefts2(2) + topoWidth + 0.01, topoBottom, 0.01, topoHeight]);
cb.Title.String = '\muV'; cb.Title.FontSize = 7;
cb.Ticks = -10:5:10; cb.TickLabels = arrayfun(@num2str, cb.Ticks, 'UniformOutput', false);

print(gcf, 'SPQ_VeridicalConditions', '-dsvg');

%% ----------- SUPPLEMENT: SPQ GROUP ERPs – Stranger Conditions ----------
% Define condition indices
cond_match  = 2;  % stranger match
cond_listen = 3;  % stranger listen
cond_mis    = 4;  % stranger mismatch

% Average ERP per group and condition
erpCz = @(condIdx, groupIdx) mean(C(:,:,groupIdx,condIdx), 3);
erpMatchCon = erpCz(cond_match, control_idx);
erpMisCon   = erpCz(cond_mis, control_idx);
erpListCon  = erpCz(cond_listen, control_idx);

erpMatchHigh = erpCz(cond_match, high_idx);
erpMisHigh   = erpCz(cond_mis, high_idx);
erpListHigh  = erpCz(cond_listen, high_idx);

% Plot setup
figure('Color','w','Position',[100 200 1100 420]);
time = cfg.tVec;

% Define new color/style mapping
colors = [0.5 0 0.5;    % purple – match
          0 0.6 0.6;    % teal – mismatch
          0.5 0.3 0.1]; % brown – listen
styles = {'--', ':', '-'};

% --- Subplot 1: Control group
subplot(1,2,1); hold on;text(-0.15, 1.05, 'c)', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 12);
h1 = plot(time, erpMatchCon(cfg.chanCz,:), styles{1}, 'Color', colors(1,:), 'LineWidth', 1.5);
h2 = plot(time, erpMisCon(cfg.chanCz,:), styles{2}, 'Color', colors(2,:), 'LineWidth', 1.5);
h3 = plot(time, erpListCon(cfg.chanCz,:), styles{3}, 'Color', colors(3,:), 'LineWidth', 1.5);
legend([h1 h2 h3], {'stranger_match', 'stranger_mismatch', 'stranger_listen'}, ...
       'Location', 'best', 'Interpreter', 'none');
lg = legend(); lg.AutoUpdate = 'off';
hold off;
xlim([-250 748]); ylim([-12 10]);
title('Control Group: N1 in Stranger Conditions');
xlabel('Time [ms]'); ylabel('Amplitude [\muV]');
line([0 0], ylim, 'Color', [0 0 0 0.3], 'LineStyle', ':');
line([100 100], ylim, 'Color', [0 0 0 0.3], 'LineStyle', '--');

% --- Subplot 2: High SPQ group
subplot(1,2,2); hold on;text(-0.15, 1.05, 'd)', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 12);
h4 = plot(time, erpMatchHigh(cfg.chanCz,:), styles{1}, 'Color', colors(1,:), 'LineWidth', 1.5);
h5 = plot(time, erpMisHigh(cfg.chanCz,:), styles{2}, 'Color', colors(2,:), 'LineWidth', 1.5);
h6 = plot(time, erpListHigh(cfg.chanCz,:), styles{3}, 'Color', colors(3,:), 'LineWidth', 1.5);
legend([h4 h5 h6], {'stranger_match', 'stranger_mismatch', 'stranger_listen'}, ...
       'Location', 'best', 'Interpreter', 'none');
lg = legend(); lg.AutoUpdate = 'off';
hold off;
xlim([-250 748]); ylim([-12 12]);
title('High SPQ Group: N1 in Stranger Conditions');
xlabel('Time [ms]'); ylabel('Amplitude [\muV]');
line([0 0], ylim, 'Color', [0 0 0 0.3], 'LineStyle', ':');
line([100 100], ylim, 'Color', [0 0 0 0.3], 'LineStyle', '--');

% ------------------------ TOPOPLOTS: Control group -----------------------
pos1 = get(subplot(1,2,1), 'Position');
topoWidth = 0.07; topoHeight = 0.09;
topoBottom = pos1(2) + 0.02;
topoLefts = pos1(1) + pos1(3) - (3:-1:1)*(topoWidth + 0.003);

% N1 latencies
[~, latMatchC] = min(erpMatchCon(cfg.chanCz, idxN1)); latIdxMatchC = idxN1(latMatchC);
[~, latMisC]   = min(erpMisCon(cfg.chanCz, idxN1));   latIdxMisC   = idxN1(latMisC);
[~, latListC]  = min(erpListCon(cfg.chanCz, idxN1));  latIdxListC  = idxN1(latListC);

axes('Position', [topoLefts(1), topoBottom, topoWidth, topoHeight]);
topoplot(erpMatchCon(:, latIdxMatchC), chanlocs, 'electrodes', 'off'); caxis([-10 10]); title('Match','FontSize',8);
axes('Position', [topoLefts(2), topoBottom, topoWidth, topoHeight]);
topoplot(erpMisCon(:, latIdxMisC), chanlocs, 'electrodes', 'off'); caxis([-10 10]); title('Mismatch','FontSize',8);
axes('Position', [topoLefts(3), topoBottom, topoWidth, topoHeight]);
topoplot(erpListCon(:, latIdxListC), chanlocs, 'electrodes', 'off'); caxis([-10 10]); title('Listen','FontSize',8);

% ------------------------ TOPOPLOTS: High SPQ group ----------------------
pos2 = get(subplot(1,2,2), 'Position');
topoLefts2 = pos2(1) + pos2(3) - (3:-1:1)*(topoWidth + 0.003);

[~, latMatchH] = min(erpMatchHigh(cfg.chanCz, idxN1)); latIdxMatchH = idxN1(latMatchH);
[~, latMisH]   = min(erpMisHigh(cfg.chanCz, idxN1));   latIdxMisH   = idxN1(latMisH);
[~, latListH]  = min(erpListHigh(cfg.chanCz, idxN1));  latIdxListH  = idxN1(latListH);

axes('Position', [topoLefts2(1), topoBottom, topoWidth, topoHeight]);
topoplot(erpMatchHigh(:, latIdxMatchH), chanlocs, 'electrodes', 'off'); caxis([-10 10]); title('Match','FontSize',8);
axes('Position', [topoLefts2(2), topoBottom, topoWidth, topoHeight]);
topoplot(erpMisHigh(:, latIdxMisH), chanlocs, 'electrodes', 'off'); caxis([-10 10]); title('Mismatch','FontSize',8);
axes('Position', [topoLefts2(3), topoBottom, topoWidth, topoHeight]);
topoplot(erpListHigh(:, latIdxListH), chanlocs, 'electrodes', 'off'); caxis([-10 10]); title('Listen','FontSize',8);

% ----------------------------- Colorbar ---------------------------------
cb = colorbar('Position', [topoLefts2(3) + topoWidth + 0.01, topoBottom, 0.01, topoHeight]);
cb.Title.String = '\muV'; cb.Title.FontSize = 7;
cb.Ticks = -10:5:10; cb.TickLabels = arrayfun(@num2str, cb.Ticks, 'UniformOutput', false);

print(gcf, 'SPQ_StrangerConditions', '-dsvg');