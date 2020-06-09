%% Ketamine Pilot Study
% KJS pilot experiment: January 17-29, 2019
% 
% Subject IDs: A201, A202 = male SD rats
%   Electrode implants in vHPC, dHPC, mPFC-PL,and mPFC-IL
%
% Aim:
%   Determine how acute low-dose ketamine (10 mg/kg, i.p.) affects HPC-mPFC
%   signaling in a familiar arena at +10min, +3hr, and +24hr after
%   injection. Determine effects relative to saline injection (same time
%   points).
% 
% Calls upon scripts:
%   - KETrexlist
%   - ImportCSC2
%   - PowSpec
%   - HilbKJS
%   - ImportVTBL
%   - vline
%   - ...
% 
% 
% 
%% SETUP (hard-coded, user must modify)
subjs = {'A201' 'A202'}; %animal IDs
tx = {'SAL' 'KET10'}; %treatments
timepts = {'BL' '1MIN' '10MIN' '3HR' '24HR'}; %timepoints recorded
    % SAL: timepts{1,3,4,5} (no 1MIN)    KET: timepts{1,2,3,4,5}

% Root input directory
rootdrIn = [uigetdir(pwd,'Select root NLX input directory holding male data: KETPilotExperiment') filesep];  %data input directory
	% 'K:\Datasets Storage\Neuralynx\[user]\M\KETPilotExperiment\';

% Root output directory
rootdrOut = [uigetdir(pwd, 'Select root output directory for RawEEG: KETPilot') filesep];
	% 'K:\Personal Folders\[user]\Neuralynx\DATA\REVAMPED\dat\RawEEG\KETPilot\';

% Add file paths containing necessary functions
addpath(uigetdir(pwd, 'Select folder containing functions called'))
	% 'K:\Personal Folders\[user]\MATLAB\gitRepo\m\code'
	% 'K:\Personal Folders\[user]\MATLAB\gitRepo\toolboxes\removePLI'
	% 'K:\Personal Folders\[user]\MATLAB\gitRepo\m\sandbox'

cd(rootdrIn)

%% Import and preprocess data
for si = 1:length(subjs) %loop thru subjects
    cd(subjs{si})
    subjID = subjs{si}; %subject ID (ex: 'A201')
    for ri = 1:length(tx) %loop thru treatments
        cd(tx{ri})
        
        % Get session ID listing
            files=dir;
            if strcmp(tx{ri},'SAL') %SAL recordings
                [BL,TenMin,ThreeHr,TwoFourHr] = KETrexlist(files); 
                sessID = {BL; TenMin; ThreeHr; TwoFourHr};
            else %KET recordings
                [BL,TenMin,ThreeHr,TwoFourHr,OneMin] = KETrexlist(files); 
                sessID = {BL; OneMin; TenMin; ThreeHr; TwoFourHr};
            end
            clear files
        
        %% Import data
            if strcmp(tx{ri},'SAL')
                for rexi = 1:length(sessID) %loop thru recording sessions
                    % Identify directory holding 'VT1_PrcFields_OrdCor_SmoothedNoBins.nvt'  	*** MUST BE HARD-CODED ***
                    switch sessID{rexi}(end-3:end)
                        case 'L_BL'
                            % drInvid = ['P:\Personal Folders\[user]\Neuralynx\DATA\REVAMPED\dat\VT\KETPilot\NLX\' subjID filesep 'SAL\SAL_BL\'];
                        case '0MIN'
                            % drInvid = ['P:\Personal Folders\[user]\Neuralynx\DATA\REVAMPED\dat\VT\KETPilot\NLX\' subjID filesep 'SAL\SAL2_10MIN\'];
                        case '_3HR'
                            % drInvid = ['P:\Personal Folders\[user]\Neuralynx\DATA\REVAMPED\dat\VT\KETPilot\NLX\' subjID filesep 'SAL\SAL3_3HR\'];
                        case '24HR'
                            % drInvid = ['P:\Personal Folders\[user]\Neuralynx\DATA\REVAMPED\dat\VT\KETPilot\NLX\' subjID filesep 'SAL\SAL4_24HR\'];
                    end 
                    
                    cd(sessID{rexi})
                    drIn = [pwd filesep];
                    %% 1. Import CSC (EEG/LFP) files
                    [EEG,thetadata,AllDat,Fs]= ImportCSC2(pwd);                
                    %% 2. LFP Power Spectra - all velocities considered
                    [fourierCoefs,frex] = PowSpec(AllDat,Fs);
                    %% 3. Hilbert transform on thetadata
                    [thetadata,itpc,prefAngle] = HilbKJS(thetadata);
                    %% 4. Import, trim, and adjust VT data (video tracker)
                    eegtsec = EEG.t*10^-4; %convert LFP ts to second units
                    [pos,ExpKeys] = ImportVTBL(drInvid,eegtsec);
                    %% SAVE entire dataset
                    cd(rootdrOut)
                    cd(subjID)
                    cd SAL
                    savefile = [subjID '_' sessID{rexi} '_AllDat.mat'];
                    disp('Saving large dataset..')
                    save(savefile,'AllDat','d*','E*','eegtsec','f*','Fs','itpc','pos','prefAngle','savefile','subjID','thetadata','rexi','-v7.3')
                    disp('Data saved.')
                    
                    %reset workspace for next recording
                    clear EEG thetadata AllDat fourierCoefs ExpKeys pos itpc prefAngle eegtsec save* frex
                    cd(drIn)
                    cd ..
                end   
                cd ..
                clear BL drI* rexi TenMin ThreeHr TwoFourHr
                
            else %KET
                for rexi = 1:length(sessID)
                     % Identify directory holding 'VT1_PrcFields_OrdCor_SmoothedNoBins.nvt'		*** MUST BE HARD-CODED ***
                    switch sessID{rexi}(end-3:end)
                        case 'T_BL'
                            % drInvid = ['P:\Personal Folders\[user]\Neuralynx\DATA\REVAMPED\dat\VT\KETPilot\NLX\' subjID filesep 'KET\KET_BL\'];
                        case '1MIN'
                            % drInvid = ['P:\Personal Folders\[user]\Neuralynx\DATA\REVAMPED\dat\VT\KETPilot\NLX\' subjID filesep 'KET\KET1_1MIN\'];
                        case '0MIN'
                            % drInvid = ['P:\Personal Folders\[user]\Neuralynx\DATA\REVAMPED\dat\VT\KETPilot\NLX\' subjID filesep 'KET\KET2_10MIN\'];
                        case '_3HR'
                            % drInvid = ['P:\Personal Folders\[user]\Neuralynx\DATA\REVAMPED\dat\VT\KETPilot\NLX\' subjID filesep 'KET\KET3_3HR\'];
                        case '24HR'
                            % drInvid = ['P:\Personal Folders\[user]\Neuralynx\DATA\REVAMPED\dat\VT\KETPilot\NLX\' subjID filesep 'KET\KET4_24HR\'];
                    end 
                    
                    cd(sessID{rexi})
                    drIn = [pwd filesep];
                    %% 1. Import CSC (EEG/LFP) files
                    [EEG,thetadata,AllDat,Fs]= ImportCSC2(drIn);                
                    %% 2. LFP Power Spectra - all velocities considered
                    [fourierCoefs,frex] = PowSpec(AllDat,Fs);
                    %% 3. Hilbert transform on thetadata
                    [thetadata,itpc,prefAngle] = HilbKJS(thetadata);
                    %% 4. Import, trim, and adjust VT data (video tracker)
                    eegtsec = EEG.t*10^-4; %convert LFP ts to second units
                    [pos,ExpKeys] = ImportVTBL(drInvid,eegtsec);
                    %% SAVE entire dataset
                    cd(rootdrOut)
                    cd(subjID)
                    cd KET
                    savefile = [subjID '_' sessID{rexi} '_AllDat.mat'];
                    disp('Saving large dataset..')
                    save(savefile,'AllDat','d*','E*','eegtsec','f*','Fs','itpc','pos','prefAngle','savefile','subjID','thetadata','rexi','-v7.3')
                    disp('Data saved.')
                    
                    %reset workspace for next recording
                    clear EEG thetadata AllDat fourierCoefs ExpKeys pos itpc prefAngle eegtsec save* frex
                    cd(drIn)
                    cd ..
                end  %loop thru sessions 
                cd ..
                clear BL drI* rexi TenMin ThreeHr TwoFourHr
            end %if SAL/KET
    end %loop thru treatments
    cd(rootdrIn)    
end %subjects

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Reducing channel count
% Set input directory (1 per animal)
subjID = input('Input subject ID: ','s');
rootdrIn = [uigetdir(pwd) filesep];
	% ['P:\Personal Folders\[user]\Neuralynx\DATA\REVAMPED\dat\RawEEG\KETPilot\' subjID];

% Set Output directory
drOut = [uigetdir(pwd) filesep];
	% ['P:\Personal Folders\[user]\Neuralynx\DATA\REVAMPED\dat\ReducedEEG\KETPilot\' subjID '\'];

cd(rootdrIn)

%% Load SAL dataset
cd SAL
fname = [subjID '_SAL_AllDat.mat'];
load(fname,'SALAllDat','BL','TenMin','ThreeHr','TwoFourHr')

% Load pos structs
    fn = [subjID '_' BL '_AllDat.mat'];
    Pos{1} = load(fn,'pos'); %BL
    fn = [subjID '_' TenMin '_AllDat.mat'];
    Pos{2} = load(fn,'pos'); %10min
    fn = [subjID '_' ThreeHr '_AllDat.mat'];
    Pos{3} = load(fn,'pos'); %3hr
    fn = [subjID '_' TwoFourHr '_AllDat.mat'];
    Pos{4} = load(fn,'pos'); %24hr
    clear fn
cd ..
%% Load KET dataset
cd KET
fname = [subjID '_KET_AllDat.mat'];
load(fname,'KETAllDat','Fs','BL','TenMin','ThreeHr','TwoFourHr','OneMin')

% Load pos structs
    fn = [subjID '_' BL '_AllDat.mat'];
    PosK{1} = load(fn,'pos'); %BL
    fn = [subjID '_' OneMin '_AllDat.mat'];
    PosK{2} = load(fn,'pos'); %1min
    fn = [subjID '_' TenMin '_AllDat.mat'];
    PosK{3} = load(fn,'pos'); %10min
    fn = [subjID '_' ThreeHr '_AllDat.mat'];
    PosK{4} = load(fn,'pos'); %3hr
    fn = [subjID '_' TwoFourHr '_AllDat.mat'];
    PosK{5} = load(fn,'pos'); %24hr
    clear fn
cd ..

%% Velocity data: total distance
for i = 1:length(Pos)  %SAL
    [VelSAL(i)] = deal(Pos{1,i}.pos.totaldist);
end

for i = 1:length(PosK) %KET
    [VelKET(i)] = deal(PosK{1,i}.pos.totaldist);
end
clear i

% Also get all BL velocities for this animal
    % This file is created by 'VTPosAdjust.m' (lines 155-197)
lfpdrIn = [uigetdir(pwd), filesep];
	% ['P:\Personal Folders\[user]\Neuralynx\DATA\REVAMPED\dat\VT\BL\MAT\' subjID filesep];
fn = [lfpdrIn subjID '_BLTotalDistance.mat'];
load(fn,'d')
BLVel = d; clear d

AllVel = [BLVel; VelSAL'; VelKET'];

% Plot velocity over sessions
dh = figure('position',[491 529 1075 420]);
set(dh,'DefaultAxesFontSize',14)
    scatter(1:length(BLVel),BLVel,'ko','markerfacecolor','k')
    hold on
    scatter(length(BLVel)+1:length(BLVel)+4,VelSAL,'bo','markerfacecolor','b') %SAL
    scatter(length(BLVel)+5:length(AllVel),VelKET,'ro','markerfacecolor','r') %SAL
    xlabel('sessionID')
    ylabel('meters (10min)')
    title([subjID ' Familiar arena locomotion'])
    ax=gca;
    set(ax,'ticklabelinterpreter','none','xticklabelmode','manual','xtickmode','manual')
    box off
    xticks(1:length(AllVel))
    xticklabels({1:length(BLVel),'BL', 'SAL10m', 'SAL3h', 'SAL24h','BL','KET1m','KET10m','KET3h','KET24h'})
    xtickangle(45)
    vline(length(BLVel)+0.5,'b')
    vline(length(BLVel)+length(VelSAL)+0.5,'r')
    set(ax,'fontsize',12,'TitleFontSizeMultiplier',1.5)
    ax.YLabel.FontSize = 16;
    ax.XLabel.FontSize = 16;
    legend FamiliarArena SAL KET
    if strcmp(subjID,'A202')
        legend('location','northwest')
    end
% save figure
    fd = [uigetdir(pwd) filesep] ; %figure output directory
		% 'K:\Personal Folders\[user]\Neuralynx\DATA\REVAMPED\figs\VT\KETPilot\';
    fn = [fd subjID '_TotalDistTraveled.tif'];
    saveas(dh,fn)
    fn = [fd subjID '_TotalDistTraveled.fig'];
    saveas(dh,fn)
    close(dh)
    clear dh ax fd fn


%% Save compiled data: AllDat and Pos
clear BL OneMin TenMin ThreeHr TwoFourHr
fd = [uigetdir(pwd) filesep] ; %figure output directory
	% 'K:\Personal Folders\[user]\Neuralynx\DATA\REVAMPED\dat\VT\KETPilot\MAT\';
fn = [fd subjID '_CompiledVelocity-AllDat16ch.mat'];
save(fn)
disp('Data saved.')
clear fn


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Combined velocity/time figure for A201 and A202
clear AllVel BLVel f* KET* lfpd* Pos* SAL* Vel* subj*
cd(fd) %(from line 265)
	% 'K:\Personal Folders\[user]\Neuralynx\DATA\REVAMPED\dat\VT\KETPilot\MAT\')
load('A201_CompiledVelocity-AllDat16ch.mat','AllVel')
A201vel = AllVel;
clear AllVel
load('A202_CompiledVelocity-AllDat16ch.mat','AllVel')
A202vel = AllVel;
clear AllVel
AllVel = [A201vel A202vel]; %concatenate

dh = figure('position',[491 529 1075 420]);
set(dh,'DefaultAxesFontSize',14)
    scatter(1:27,BLVel,'ko','markerfacecolor','k')
    hold on
    scatter(length(BLVel)+1:length(BLVel)+4,VelSAL,'bo','markerfacecolor','b') %SAL
    scatter(length(BLVel)+5:length(AllVel),VelKET,'ro','markerfacecolor','r') %SAL
    xlabel('sessionID')
    ylabel('meters (10min)')
    title([subjID ' Familiar arena locomotion'])
    ax=gca;
    set(ax,'ticklabelinterpreter','none','xticklabelmode','manual','xtickmode','manual')
    box off
    xticks(1:length(AllVel))
    xticklabels({1:length(BLVel),'BL', 'SAL10m', 'SAL3h', 'SAL24h','BL','KET1m','KET10m','KET3h','KET24h'})
    xtickangle(45)
    vline(length(BLVel)+0.5,'b')
    vline(length(BLVel)+length(VelSAL)+0.5,'r')
    set(ax,'fontsize',12,'TitleFontSizeMultiplier',1.5)
    ax.YLabel.FontSize = 16;
    ax.XLabel.FontSize = 16;
    legend FamiliarArena SAL KET

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure: representative traces (A201)
h = figure('units','normalized','position',[0 0 1 1]);
subplot 421
    plot(KET_IL(10000:16000,1),'k','linew',1)
    ax=gca;
    set(ax,'fontsize',12,'TitleFontSizeMultiplier',1)
    axis tight
    set(gca,'xticklabel',{[]},'xlabel',[],'ylim',[-15000 15000])
    ylabel('\muV')
    box off
    grid off
    title([subjID ' mPFC-IL: BL'])
subplot 422
    plot(KET_IL(10000:16000,3),'k','linew',1)
    ax=gca;
    set(ax,'fontsize',12,'TitleFontSizeMultiplier',1)
    axis tight
    set(gca,'xticklabel',{[]},'xlabel',[],'ylim',[-15000 15000])
    ylabel('\muV')
    box off
    grid off
    title([subjID ' mPFC-IL: KET +10min'])
subplot 423
    plot(KET_PL(10000:16000,1),'k','linew',1)
    ax=gca;
    set(ax,'fontsize',12,'TitleFontSizeMultiplier',1)
    axis tight
    set(gca,'xticklabel',{[]},'xlabel',[],'ylim',[-15000 15000])
    ylabel('\muV')
    box off
    grid off
    title([subjID ' mPFC-PL: BL'])
subplot 424
    plot(KET_PL(10000:16000,3),'k','linew',1)
    ax=gca;
    set(ax,'fontsize',12,'TitleFontSizeMultiplier',1)
    axis tight
    set(gca,'xticklabel',{[]},'xlabel',[],'ylim',[-15000 15000])
    ylabel('\muV')
    box off
    grid off
    title([subjID ' mPFC-PL: KET +10min'])
subplot 425
    plot(KET_DHIP(10000:16000,1),'k','linew',1)
    ax=gca;
    set(ax,'fontsize',12,'TitleFontSizeMultiplier',1)
    axis tight
    set(gca,'xticklabel',{[]},'xlabel',[],'ylim',[-15000 15000])
    ylabel('\muV')
    box off
    grid off
    title([subjID ' dHPC: BL'])
subplot 426
    plot(KET_DHIP(10000:16000,3),'k','linew',1)
    ax=gca;
    set(ax,'fontsize',12,'TitleFontSizeMultiplier',1)
    axis tight
    set(gca,'xticklabel',{[]},'xlabel',[],'ylim',[-15000 15000])
    ylabel('\muV')
    box off
    grid off
    title([subjID ' dHPC: KET +10min'])
subplot 427
    plot(KET_VHIP(10000:16000,1),'k','linew',1)
    ax=gca;
    set(ax,'fontsize',12,'TitleFontSizeMultiplier',1)
    axis tight
    set(gca,'xticklabel',{[1:6]},'ylim',[-15000 15000])
    xlabel('time (s)')
     ylabel('\muV')
    box off
    grid off
    title([subjID ' vHPC: BL'])
subplot 428
    plot(KET_VHIP(10000:16000,3),'k','linew',1)
    ax=gca;
    set(ax,'fontsize',12,'TitleFontSizeMultiplier',1)
    axis tight
    set(ax,'xticklabel',{[1:6]},'ylim',[-15000 15000])
    xlabel('time (s)')
    ylabel('\muV')
    box off
    grid off
    title([subjID ' vHPC: KET +10min'])
    
% Save figure
 fd = [uigetdir(pwd) filesep]; %figure output directory
  % 'K:\Personal Folders\Kristin Schoepfer\Lab Meetings\20190227\ket\';
    fn = [fd subjID '_representativeKETtraces']; %figure name
    saveas(h,[fn '.png'])
	saveas(h,[fn '.fig'])
    close(h)
    clear dh ax fd fn


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% %%%%%%%%%%% Power spectra figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cd('K:\Personal Folders\[user]\Neuralynx\DATA\REVAMPED\dat\ReducedEEG\KETPilot\A201\')
load('A201_SALvsKET_AllDat-ReducedChan-Cohere.mat')
set(0,'DefaultAxesFontSize',18)
% cd('K:\Personal Folders\[user\Lab Meetings\20190227\ket\powspec\')
% cd A202

 %%%%% SALINE FIGURES %%%%%
% SAL IL
fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
plot(10*log10(welchdat.SAL_IL(:,1)),'k','linew',3)
hold on
plot(10*log10(welchdat.SAL_IL(:,2)),'b','linew',3)
plot(10*log10(welchdat.SAL_IL(:,3)),'c','linew',3)
plot(10*log10(welchdat.SAL_IL(:,4)),'m','linew',3)
plot(10*log10(welchdat.KET_IL(:,1)),'g','linew',3)
xlim([1 100])
ylim([40 70])
set(gca,'fontsize',30,'TitleFontSizeMultiplier',1.5)
axis square
box off
legend Baseline 10MIN 3HR 24HR KET-BL
xlabel('Frequency (Hz)')
set(gca,'xtick',(10:10:100))
ylabel('Power (dB/Hz)')
title([subjID ' mPFC-IL: Saline'])
    % Save figure
    saveas(fig1,'SAL_IL.tif')
    saveas(fig1,'SAL_IL.fig')
    close(fig1)

% SAL PL
fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
plot(10*log10(welchdat.SAL_PL(:,1)),'k','linew',3)
hold on
plot(10*log10(welchdat.SAL_PL(:,2)),'b','linew',3)
plot(10*log10(welchdat.SAL_PL(:,3)),'c','linew',3)
plot(10*log10(welchdat.SAL_PL(:,4)),'m','linew',3)
plot(10*log10(welchdat.KET_PL(:,1)),'g','linew',3)
xlim([1 100])
ylim([40 70])
set(gca,'fontsize',30,'TitleFontSizeMultiplier',1.5)
axis square
box off
legend Baseline 10MIN 3HR 24HR KET-BL
xlabel('Frequency (Hz)')
set(gca,'xtick',(10:10:100))
ylabel('Power (dB/Hz)')
title([subjID ' mPFC-PL: Saline'])
    saveas(fig1,'SAL_PL.tif')
    saveas(fig1,'SAL_PL.fig')
    close(fig1) 

% SAL DHIP
fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
plot(10*log10(welchdat.SAL_DHIP(:,1)),'k','linew',3)
hold on
plot(10*log10(welchdat.SAL_DHIP(:,2)),'b','linew',3)
plot(10*log10(welchdat.SAL_DHIP(:,3)),'c','linew',3)
plot(10*log10(welchdat.SAL_DHIP(:,4)),'m','linew',3)
plot(10*log10(welchdat.KET_DHIP(:,1)),'g','linew',3)
xlim([1 100])
ylim([40 70])
set(gca,'fontsize',30,'TitleFontSizeMultiplier',1.5)
axis square
box off
legend Baseline 10MIN 3HR 24HR KET-BL
xlabel('Frequency (Hz)')
set(gca,'xtick',(10:10:100))
ylabel('Power (dB/Hz)')
title([subjID ' dHPC: Saline'])
    saveas(fig1,'SAL_dhpc.tif')
    saveas(fig1,'SAL_dhpc.fig')
    close(fig1)
  
% SAL VHIP
fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
plot(10*log10(welchdat.SAL_VHIP(:,1)),'k','linew',3)
hold on
plot(10*log10(welchdat.SAL_VHIP(:,2)),'b','linew',3)
plot(10*log10(welchdat.SAL_VHIP(:,3)),'c','linew',3)
plot(10*log10(welchdat.SAL_VHIP(:,4)),'m','linew',3)
plot(10*log10(welchdat.KET_VHIP(:,1)),'g','linew',3)
xlim([1 100])
ylim([40 70])
set(gca,'fontsize',30,'TitleFontSizeMultiplier',1.5)
axis square
box off
legend Baseline 10MIN 3HR 24HR KET-BL
xlabel('Frequency (Hz)')
set(gca,'xtick',(10:10:100))
ylabel('Power (dB/Hz)')
title([subjID ' vHPC: Saline'])
    saveas(fig1,'SAL_vhpc.tif')
    saveas(fig1,'SAL_vhpc.fig')
    close(fig1)
    
 %%%%% KET FIGURES %%%%%
%KET IL
fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
plot(10*log10(welchdat.KET_IL(:,1)),'k','linew',3)
hold on
plot(10*log10(welchdat.KET_IL(:,2)),'b','linew',3)
% plot(10*log10(welchdat.KET_IL(:,3)),'c','linew',3)
plot(10*log10(welchdat.KET_IL(:,3)),'g','linew',3)
plot(10*log10(welchdat.KET_IL(:,4)),'m','linew',3)
xlim([1 100])
ylim([40 70])
set(gca,'fontsize',30,'TitleFontSizeMultiplier',1.5)
axis square
box off
% legend Baseline 10MIN 3HR 24HR 
legend Baseline +10MIN +3HR +24HR 
xlabel('Frequency (Hz)')
set(gca,'xtick',(10:10:100))
ylabel('Power (dB/Hz)')
title([subjID ' mPFC-IL: KET'])
    saveas(fig1,'KET_IL.tif')
    saveas(fig1,'KET_IL.fig')
    close(fig1)

% KET PL
fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
plot(10*log10(welchdat.KET_PL(:,1)),'k','linew',3)
hold on
plot(10*log10(welchdat.KET_PL(:,2)),'b','linew',3)
% plot(10*log10(welchdat.KET_PL(:,3)),'c','linew',3)
plot(10*log10(welchdat.KET_PL(:,3)),'g','linew',3)
plot(10*log10(welchdat.KET_PL(:,4)),'m','linew',3)
xlim([1 100])
ylim([40 70])
set(gca,'fontsize',30,'TitleFontSizeMultiplier',1.5)
axis square
box off
% legend Baseline 10MIN 3HR 24HR 
legend Baseline +10MIN +3HR +24HR 
xlabel('Frequency (Hz)')
set(gca,'xtick',(10:10:100))
ylabel('Power (dB/Hz)')
title([subjID ' mPFC-PL: KET'])
    saveas(fig1,'KET_PL.tif')
    saveas(fig1,'KET_PL.fig')
    close(fig1)
    
% KET DHIP
fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
plot(10*log10(welchdat.KET_DHIP(:,1)),'k','linew',3)
hold on
plot(10*log10(welchdat.KET_DHIP(:,2)),'b','linew',3)
% plot(10*log10(welchdat.KET_DHIP(:,3)),'c','linew',3)
plot(10*log10(welchdat.KET_DHIP(:,3)),'g','linew',3)
plot(10*log10(welchdat.KET_DHIP(:,4)),'m','linew',3)
xlim([1 100])
ylim([40 70])
set(gca,'fontsize',30,'TitleFontSizeMultiplier',1.5)
axis square
box off
% legend Baseline 10MIN 3HR 24HR 
legend Baseline +10MIN +3HR +24HR 
xlabel('Frequency (Hz)')
set(gca,'xtick',(10:10:100))
ylabel('Power (dB/Hz)')
title([subjID ' dHPC: KET'])
    saveas(fig1,'KET_dhip.tif')
    saveas(fig1,'KET_dhip.fig')
    close(fig1)

% KET VHIP
fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
plot(10*log10(welchdat.KET_VHIP(:,1)),'k','linew',3)
hold on
plot(10*log10(welchdat.KET_VHIP(:,2)),'b','linew',3)
% plot(10*log10(welchdat.KET_VHIP(:,3)),'c','linew',3)
plot(10*log10(welchdat.KET_VHIP(:,3)),'g','linew',3)
plot(10*log10(welchdat.KET_VHIP(:,4)),'m','linew',3)
xlim([1 100])
ylim([40 70])
set(gca,'fontsize',30,'TitleFontSizeMultiplier',1.5)
axis square
box off
% legend Baseline 10MIN 3HR 24HR 
legend Baseline +10MIN +3HR +24HR 
xlabel('Frequency (Hz)')
set(gca,'xtick',(10:10:100))
ylabel('Power (dB/Hz)')
title([subjID ' vHPC: KET'])
    saveas(fig1,'KET_vhip.tif')
    saveas(fig1,'KET_vhip.fig')
    close(fig1)
% end PowSpec figures



%% %%%%%%%%%%%%%%%%%% COHERENCE FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kjs edit 2019/03/12

% cd('K:\Personal Folders\[user]\Research talks - Seminars - Annual reports\Annual Progress Report 2018 - MK R01 KET\KETpilot\')
    % You may need to adjust drOut for later use!
    
% IL-vHPC
fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
    plot(KET_coh.ilvh(:,1),'k','linew',3)
    hold on
    plot(KET_coh.ilvh(:,2),'b','linew',3)
    plot(KET_coh.ilvh(:,3),'g','linew',3)
    plot(KET_coh.ilvh(:,4),'m','linew',3)
    xlim([1 100])
    ylim([0 1])
    set(gca,'fontsize',30,'TitleFontSizeMultiplier',1.5)
    axis square
    box off
%     legend Baseline +10MIN +3HR +24HR
    title('mPFCil-vHPC')
    ylabel('Coherence')
    xlabel('Frequency (Hz)')
    set(gca,'xtick',(10:10:100))
    set(gca,'ytick',(0:0.2:1))
    saveas(fig1,'KETcohere_ilvhip.tif')
    saveas(fig1,'KETcohere_ilvhip.fig')
    close(fig1)
	
% IL-dHPC
fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
    plot(KET_coh.ildh(:,1),'k','linew',3)
    hold on
    plot(KET_coh.ildh(:,2),'b','linew',3)
    plot(KET_coh.ildh(:,3),'g','linew',3)
    plot(KET_coh.ildh(:,4),'m','linew',3)
    xlim([1 100])
    ylim([0 1])
    set(gca,'fontsize',30,'TitleFontSizeMultiplier',1.5)
    axis square
    box off
    legend Baseline +10MIN +3HR +24HR
    title('mPFCil-dHPC')
    ylabel('Coherence')
    xlabel('Frequency (Hz)')
    set(gca,'xtick',(10:10:100))
    set(gca,'ytick',(0:0.2:1))
    saveas(fig1,'KETcohere_ildhip.tif')
    saveas(fig1,'KETcohere_ildhip.fig')
    close(fig1)
	
% PL-vHPC
fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
    plot(KET_coh.vhpl(:,1),'k','linew',3)
    hold on
    plot(KET_coh.vhpl(:,2),'b','linew',3)
    plot(KET_coh.vhpl(:,3),'g','linew',3)
    plot(KET_coh.vhpl(:,4),'m','linew',3)
    xlim([1 100])
    ylim([0 1])
    set(gca,'fontsize',30,'TitleFontSizeMultiplier',1.5)
    axis square
    box off
%     legend Baseline +10MIN +3HR +24HR
    title('mPFCpl-vHPC')
    ylabel('Coherence')
    xlabel('Frequency (Hz)')
    set(gca,'xtick',(10:10:100))
    set(gca,'ytick',(0:0.2:1))
    saveas(fig1,'KETcohere_plvhip.tif')
    saveas(fig1,'KETcohere_plvhip.fig')
    close(fig1)
	
% PL-dHPC
fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
    plot(KET_coh.dhpl(:,1),'k','linew',3)
    hold on
    plot(KET_coh.dhpl(:,2),'b','linew',3)
    plot(KET_coh.dhpl(:,3),'g','linew',3)
    plot(KET_coh.dhpl(:,4),'m','linew',3)
    xlim([1 100])
    ylim([0 1])
    set(gca,'fontsize',30,'TitleFontSizeMultiplier',1.5)
    axis square
    box off
%     legend Baseline +10MIN +3HR +24HR
    title('mPFCpl-dHPC')
    ylabel('Coherence')
    xlabel('Frequency (Hz)')
    set(gca,'xtick',(10:10:100))
    set(gca,'ytick',(0:0.2:1))
    saveas(fig1,'KETcohere_pldhip.tif')
    saveas(fig1,'KETcohere_pldhip.fig')
    close(fig1)