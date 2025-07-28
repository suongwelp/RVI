%% =====================================================================
%  EEG PRE-PROCESSING AND ERP EXTRACTION PIPELINE (EEGLAB) 
%  ---------------------------------------------------------------------
%  Author : Suong Welp
%  Date   : 2025-07-24
%  Purpose: pipeline for cleaning raw EEG,removing ocular / speaking artifacts, 
%           extracting condition‑specific ERPs for the Rubber Voice Illusion experiment.
%  ---------------------------------------------------------------------
%  Dependencies
%  * EEGLAB 2023+  (https://eeglab.org)
%  * CleanRawdata plugin (for ASR)
%  * ICLabel plugin   (for automated IC annotation)
%  =====================================================================
clc; clear all; close all;
%% --------------------------- USER INPUT --------------------------------
subID        = 130;                 % <<< change or iterate over a list
rawDir       = fullfile(pwd,'data','raw');
outDir       = fullfile(pwd,'data','ERP');
mkdir(outDir);

%% ---------------------- ANALYSIS PARAMETERS ----------------------------
% All parameters gathered in a single structure for readability
P = struct();

% Filters (linear‑phase FIR, Hamming window)
P.hpf_ica     = 1;      % [Hz] high‑pass cut‑off *before ICA*
P.lpf_ica     = 40;     % [Hz] low‑pass  cut‑off *before ICA*
P.hpf_final   = 1;      % [Hz] high‑pass for final data
P.lpf_final   = 30;     % [Hz] low‑pass  for final data
P.fir_order   = struct('HP',414,'LP',166); % pre‑computed orders (fs=250)

% Resampling
P.fs_ica      = 250;    % [Hz] – lower rate speeds ICA / ASR
P.fs_final    = 500;    % [Hz] – uniform rate for ERP alignment

% Epoching & Baseline
P.largeEpoch  = [-2 1];         % [s] big epoch around speech onset (S1,S2,S4)
P.finalEpoch  = [-0.25 0.75];   % [s] ERP window of interest
P.baseWin     = [-0.25 0];      % [s] 250‑ms pre‑stimulus baseline

% Artifact rejection
P.rejectAmp   = 100;     % [µV] +/- amplitude threshold 
P.icLabelThr  = 0.8;    % probability threshold for EOG/Muscle ICs
P.asr         = struct('flat',5,'chan',0.8,'line',4,'burst',40,...
                       'window',0.25,'distance','Euclidian');

% Condition codes (trigger strings in the .set file)
%  col1 = auditory stimulus trigger, col2 = fixation cross trigger (if any), col3 = condition label
conds = { ...
    'S 1','S 12','veridical_talk';
    'S 2','S 22','stranger_match';
    'S 3','',    'stranger_listen';    %no fixation cross prompting speech
    'S 4','S 42','stranger_mis';
    'S 5','',    'veridical_listen'};  %no fixation cross prompting speech

%% -------------------------- INITIALISE EEGLAB --------------------------
% --------------------------- LOAD DATA ---------------------------------
[ALLEEG,EEG,CURRENTSET] = eeglab('nogui');

setName  = sprintf('P%d.set',subID);
rawPath  = fullfile(rawDir,setName);
if ~isfile(rawPath)
    error('Cannot find file: %s',rawPath);
end
EEG = pop_loadset('filename', setName,'filepath', rawDir);
[ALLEEG,EEG] = eeg_store(ALLEEG, EEG, 1);

%% ----------------- BASIC RE‑REFERENCING OF EOG CHANNELS ---------------
% Convert monopolar VEOG/HEOG to bipolar derivations 
veog         = EEG.data(2,:)  - EEG.data(20,:); % IO –> IO‑Fp2
heog         = EEG.data(29,:) - EEG.data(30,:); % HEOG‑L –> HEOG‑L – HEOG‑R
EEG.data(20,:) = veog;   % overwrite channel 20 with bipolar VEOG
EEG.data(29,:) = heog;   % overwrite channel 29 with bipolar HEOG
rawEEG         = EEG;    % keep a copy of the continuous raw data

%% ------------------------- ICA PREPARATION -----------------------------
EEG = pop_firws(EEG,'fcutoff',P.lpf_ica,'ftype','lowpass','forder',P.fir_order.LP,'wtype','hamming');
EEG = pop_resample(EEG,P.fs_ica);
EEG = pop_firws(EEG,'fcutoff',P.hpf_ica,'ftype','highpass','forder',P.fir_order.HP,'wtype','hamming');

% Remove continuous flat segments (pauses)
EEG = eeg_regepochs(EEG);

% Automatic high‑variance segment removal (joint prob.)
EEG = pop_jointprob(EEG,1,1:EEG.nbchan,3,3,0,1);

% Run ICA (RunICA algorithm)
EEG = pop_runica(EEG,'icatype','runica','extended',1,'interrupt','on','PCA',EEG.nbchan);

% Save ICA matrices for re‑application on raw data
icaWeights = EEG.icaweights; icaSphere = EEG.icasphere; icaWinv = EEG.icawinv;

%% ---------------- RE‑APPLY ICA TO ORIGINAL RAW CONTINUOUS --------------
EEG = rawEEG;
EEG.icaweights = icaWeights; EEG.icasphere = icaSphere; EEG.icawinv = icaWinv;
EEG = eeg_checkset(EEG);

% ICLabel + marking of bad ICs
EEG  = pop_iclabel(EEG,'default');
EEG  = pop_icflag(EEG,[NaN NaN; NaN NaN; P.icLabelThr 1; ... % Brain
                         P.icLabelThr 1; NaN NaN; 0.9 1; NaN NaN]);
EEG.badcomps.eog    = find(EEG.reject.gcompreject==1)'; % Auto EOG components

% Manual visual inspection (interactive) -------------------------------
pop_viewprops(EEG,0,1:size(EEG.icaweights,1),{'freqrange' [2 50]});
uiwait(msgbox('Inspect ICs and note indices of speech / muscle / bad‑chan components.','ICA Inspection','modal'));

% --- Manually enter additional ICs to remove
speechICs  = input('Enter speech artifact ICs (e.g., [4 5]): ');
muscleICs  = input('Enter muscle  artifact ICs (e.g., [12]): ');
badChanICs = input('Enter bad‑channel ICs (e.g., []): ');
EEG.badcomps.speech  = speechICs;
EEG.badcomps.muscle  = muscleICs;
EEG.badcomps.channel = badChanICs;

% Remove marked ICs -----------------------------------------------------
allRemove = unique([EEG.badcomps.eog, speechICs, muscleICs, badChanICs]);
EEG       = pop_subcomp(EEG, allRemove, 0);

%% -------------------- CLEANRAW (ASR) AND INTERPOLATION ---------------
EEG = pop_select(EEG,'rmchannel',{'IO','heogl','heogr','lip'});  % remove externals
EEG.chaninfo.removedchans = [];
EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',P.asr.flat, ...
                                'ChannelCriterion',P.asr.chan, ...
                                'LineNoiseCriterion',P.asr.line, ...
                                'Highpass',[0.25 0.75], ...
                                'BurstCriterion',P.asr.burst, ...
                                'WindowCriterion',P.asr.window, ...
                                'BurstRejection','on', ...
                                'Distance',P.asr.distance, ...
                                'WindowCriterionTolerances',[-Inf 10]);
EEG = eeg_interp(EEG,EEG.chaninfo.removedchans,'spherical');

%% ---------------- FINAL FILTERING & RE‑REFERENCE ----------------------
EEG = pop_eegfiltnew(EEG,'locutoff',P.hpf_final,'hicutoff',P.lpf_final);
EEG = pop_reref(EEG,[]);          % common average reference
EEG = pop_resample(EEG,P.fs_final);

%% --------------------- ERP EXTRACTION PER CONDITION -------------------
chanCz = find(strcmpi({EEG.chanlocs.labels},'Cz'));
numCon = size(conds,1);
allERP = zeros(EEG.nbchan, diff(P.finalEpoch)*P.fs_final, numCon);

for c = 1:numCon
    condTrig  = conds{c,1};
    baseTrig  = conds{c,2};

    % --- (A) Large epoch around fixation cross --------------
    if ~isempty(baseTrig)
        EEGc = pop_epoch(EEG,{condTrig},P.largeEpoch,'epochinfo','yes');

        % baseline subtraction using 250‑ms pre‑baseline‑event
        for ep = 1:EEGc.trials
            evTypes   = EEGc.epoch(ep).eventtype;
            evLat     = EEGc.epoch(ep).eventlatency;
            baseIdx   = find(strcmp(evTypes, baseTrig),1);
            if isempty(baseIdx)
                warning('Baseline trigger (%s) not found in epoch %d. Using condition trigger.',baseTrig,ep);
                baseIdx = find(strcmp(evTypes,condTrig),1);
            end
            t0        = evLat{baseIdx};
            baseMask  = (EEGc.times >= (t0+P.baseWin(1))) & (EEGc.times <= (t0+P.baseWin(2)));
            if any(baseMask)
                baseVal = mean(EEGc.data(:,baseMask,ep),2);
                EEGc.data(:,:,ep) = EEGc.data(:,:,ep) - baseVal;
            end
        end
    else
        EEGc = EEG; % no special baseline trigger
    end

    % --- (B) Re‑epoch to final ERP window & reject artefacts ------------
    EEGc = pop_epoch(EEGc,{condTrig},P.finalEpoch,'epochinfo','yes');
    EEGc = pop_eegthresh(EEGc,1,[], -P.rejectAmp, P.rejectAmp, P.finalEpoch(1), P.finalEpoch(2),2,0);
    EEGc = pop_jointprob(EEGc,1,1:EEGc.nbchan,3,3,0,1);

    % --- Store grand‑average ERP (µV) ----------------------------------
    allERP(:,:,c) = mean(EEGc.data,3);

    % --- Quick QC plot of Cz ERP & topography --------------------------
    figTitle = sprintf('P%d – %s',subID,conds{c,3});
    figure('Name',figTitle);
    subplot(1,2,1);
    plot(EEGc.times,mean(EEGc.data(chanCz,:,:),3),'LineWidth',1.2);
    title(['Cz ERP – ' conds{c,3}],'Interp','none'); xlabel('Time (ms)'); ylabel('\muV'); 
    grid on;xlim([-250 750]);


    % N1 peak detection
    n1Win = find(EEGc.times>=50 & EEGc.times<=150);
    [amp,lat] = min(mean(EEGc.data(chanCz,n1Win,:),3));
    subplot(1,2,2);
    topoplot(mean(EEGc.data(:,n1Win(1)+lat-1,:),3), EEGc.chanlocs);
    colorbar; title(sprintf('N1 @ %d ms', EEGc.times(n1Win(1)+lat-1)));
end

%% ------------------------- SAVE OUTPUT ---------------------------------
matName = sprintf('P%d_ERP.mat',subID);
save(fullfile(outDir,matName),'allERP','P','conds','-v7.3');

fprintf('\n=== Finished subject %d – data saved to %s ===\n',subID,matName);
