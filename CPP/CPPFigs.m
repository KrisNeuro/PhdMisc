%% CPP movement figures
% For Devin: KET vs NLX experiment

%% SETUP
% EXPT struct: Experimental conditions
Expt.Pretx = {'sal' 'nlx'}; %expt. pretreatment drug options
Expt.Tx = {'sal' 'ket'}; % expt. treatment drug options
% Expt.FeatNames = {'Control ','Treatment'};	%names of relevant features to be analyzed

% Load Trial and TimeCourse from file made with CPPMvmtCounts.m
[fname,fpath] = uigetfile;
cd(fpath)
load(fname,'Trial','TimeCourse');

%% Get means and std for each group per conditioning session
SS = [];
SK = [];
NS = [];
NK = [];
for prei = 1:length(Expt.Pretx) %loop thru pretx conds
	pretxname = Expt.Pretx{prei};
	for txi = 1:length(Expt.Tx) %loop thru tx conds
		txname = Expt.Tx{txi};
		switch [pretxname(1) txname(1)]
			case 'ss'
				for sesi = 1:4 %loop thru sessions per thing
					[SS.Conmeans(sesi)] = mean([Trial.(pretxname).(txname).(['Con' num2str(sesi)])].');
					[SS.Constds(sesi)] = std([Trial.(pretxname).(txname).(['Con' num2str(sesi)])].');
					[SS.Txmeans(sesi)] = mean([Trial.(pretxname).(txname).(['Tx' num2str(sesi)])].');
					[SS.Txstds(sesi)] = std([Trial.(pretxname).(txname).(['Tx' num2str(sesi)])].');
				end
			case 'sk'
				for sesi = 1:4 %loop thru sessions per thing
					[SK.Conmeans(sesi)] = mean([Trial.(pretxname).(txname).(['Con' num2str(sesi)])].');
					[SK.Constds(sesi)] = std([Trial.(pretxname).(txname).(['Con' num2str(sesi)])].');
					[SK.Txmeans(sesi)] = mean([Trial.(pretxname).(txname).(['Tx' num2str(sesi)])].');	%
					[SK.Txstds(sesi)] = std([Trial.(pretxname).(txname).(['Tx' num2str(sesi)])].');	%
				end
			case 'ns'
				for sesi = 1:4 %loop thru sessions per thing
					[NS.Conmeans(sesi)] = mean([Trial.(pretxname).(txname).(['Con' num2str(sesi)])].');
					[NS.Constds(sesi)] = std([Trial.(pretxname).(txname).(['Con' num2str(sesi)])].');
					[NS.Txmeans(sesi)] = mean([Trial.(pretxname).(txname).(['Tx' num2str(sesi)])].');
					[NS.Txstds(sesi)] = std([Trial.(pretxname).(txname).(['Tx' num2str(sesi)])].');
				end
			case 'nk'
				for sesi = 1:4 %loop thru sessions per thing
					[NK.Conmeans(sesi)] = mean([Trial.(pretxname).(txname).(['Con' num2str(sesi)])].');
					[NK.Constds(sesi)] = std([Trial.(pretxname).(txname).(['Con' num2str(sesi)])].');
					[NK.Txmeans(sesi)] = mean([Trial.(pretxname).(txname).(['Tx' num2str(sesi)])].');
					[NK.Txstds(sesi)] = std([Trial.(pretxname).(txname).(['Tx' num2str(sesi)])].');
				end
		end %switch
	end %txi for
end %prei for
clear txname pretxname prei txi sesi

%% Standard error calcs
    %SAL/SAL
        [SS.ConErr] = SS.Constds/(sqrt(length(Trial.sal.sal)));
        [SS.TxErr] = SS.Txstds/(sqrt(length(Trial.sal.sal)));
    %SAL/KET
        [SK.ConErr] = SK.Constds/(sqrt(length(Trial.sal.ket)));
        [SK.TxErr] = SK.Txstds/(sqrt(length(Trial.sal.ket)));
    %NLX/SAL
        [NS.ConErr] = NS.Constds/(sqrt(length(Trial.nlx.sal)));
        [NS.TxErr] = NS.Txstds/(sqrt(length(Trial.nlx.sal)));
    %NLX/KET
        [NK.ConErr] = NK.Constds/(sqrt(length(Trial.nlx.ket)));
        [NK.TxErr] = NK.Txstds/(sqrt(length(Trial.nlx.ket)));

%% Total movement count figs, 1 fig per group
% change defaults for plotting
	set(0,'DefaultAxesFontSize',16)
%% SAL/SAL
fig1 = figure;
subplot(1,2,1)
	bar([1 2 3 4], [SS.Conmeans]);
	hold on;
	errorbar([1 2 3 4], [SS.Conmeans], [SS.ConErr],'.');
	xlabel('Conditioning session')
	ylabel('Mvmt counts (45 mins)')
	title('SAL/SAL')
	legend '+SAL (Con)'
    ylim([0 3500])
subplot(1,2,2)
	bar([1 2 3 4], [SS.Txmeans]);
	hold on;
	errorbar([1 2 3 4], [SS.Txmeans], [SS.TxErr],'.');
	xlabel('Conditioning session')
	ylabel('Mvmt counts (45 mins)')
	legend '+SAL (Tx)'
    ylim([0 3500])

%% SAL/KET
fig2 = figure;
subplot(1,2,1)
	bar([1 2 3 4], [SK.Conmeans]);
	hold on;
	errorbar([1 2 3 4], [SK.Conmeans], [SK.ConErr],'.');
	xlabel('Conditioning session')
	ylabel('Mvmt counts (45 mins)')
	title('SAL/KET')
	legend '+SAL/KET (Con)'
	ylim([0 3500])
subplot(1,2,2)
	bar([1 2 3 4], [SK.Txmeans]);
	hold on;
	errorbar([1 2 3 4], [SK.Txmeans], [SK.TxErr],'.');
	xlabel('Conditioning session')
	ylabel('Mvmt counts (45 mins)')
	legend '+SAL/KET (Tx)'
    ylim([0 3500])

%% NLX/SAL
fig3 = figure;
subplot(1,2,1)
	bar([1 2 3 4], [NS.Conmeans]);
	hold on;
	errorbar([1 2 3 4], [NS.Conmeans], [NS.ConErr],'.');
	xlabel('Conditioning session')
	ylabel('Mvmt counts (45 mins)')
	title('NLX/SAL')
	% ylim([0 3000])
	ylim([0 3500])
	legend '+NLX/SAL (Con)'
subplot(1,2,2)
	bar([1 2 3 4], [NS.Txmeans]);
	hold on;
	errorbar([1 2 3 4], [NS.Txmeans], [NS.TxErr],'.');
	xlabel('Conditioning session')
	ylabel('Mvmt counts (45 mins)')
	legend '+NLX/SAL (Tx)'
	ylim([0 3500])
	
%% NLX/KET
fig4 = figure;
subplot(1,2,1)
	bar([1 2 3 4], [NK.Conmeans]);
	hold on;
	errorbar([1 2 3 4], [NK.Conmeans], [NK.ConErr],'.');
	xlabel('Conditioning session')
	ylabel('Mvmt counts (45 mins)')
	title('NLX/KET')
	ylim([0 3500])
	legend '+NLX/KET (Con)'
subplot(1,2,2)
	bar([1 2 3 4], [NK.Txmeans]);
	hold on;
	errorbar([1 2 3 4], [NK.Txmeans], [NK.TxErr],'.');
	xlabel('Conditioning session')
	ylabel('Mvmt counts (45 mins)')
	legend '+NLX/KET (Tx)'
	ylim([0 3500])

%% TIME COURSE FIGS
t = [1:45]; %time bins (mins)

%% SAL/SAL
%legend fig
    figure; plot(t,zeros(length(t)),'Linewidth',7);
    legend A39 A46 A49 A56 A59 A60 A65 A70
%timepoint fig
sstime = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(4,2,1)
    title('SAL/SAL: Con 1')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    yticks([0:100:500])
    xlim([0 45])
    hold on;
    for i = 1:length(TimeCourse.sal.sal) 
        plot(t,TimeCourse.sal.sal(i).Con1,'Linewidth',2)
    end
subplot(4,2,2)
    title('SAL/SAL: Tx 1')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    yticks([0:100:500])
    xlim([0 45])
    hold on;
    for i = 1:length(TimeCourse.sal.sal) 
        plot(t,TimeCourse.sal.sal(i).Tx1,'Linewidth',2)
    end
subplot(4,2,3)
    title('SAL/SAL: Con 2')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    yticks([0:100:500])
    xlim([0 45])
    hold on;
    for i = 1:length(TimeCourse.sal.sal) 
        plot(t,TimeCourse.sal.sal(i).Con2,'Linewidth',2)
    end
subplot(4,2,4)
    title('SAL/SAL: Tx 2')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    yticks([0:100:500])
    xlim([0 45])
    hold on;
    for i = 1:length(TimeCourse.sal.sal) 
        plot(t,TimeCourse.sal.sal(i).Tx2,'Linewidth',2)
    end
subplot(4,2,5)
    title('SAL/SAL: Con 3')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    yticks([0:100:500])
    xlim([0 45])
    hold on;
    for i = 1:length(TimeCourse.sal.sal) 
        plot(t,TimeCourse.sal.sal(i).Con3,'Linewidth',2)
    end
subplot(4,2,6)
    title('SAL/SAL: Tx 3')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    yticks([0:100:500])
    xlim([0 45])
    hold on;
    for i = 1:length(TimeCourse.sal.sal) 
        plot(t,TimeCourse.sal.sal(i).Tx3,'Linewidth',2)
    end
subplot(4,2,7)
    title('SAL/SAL: Con 4')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    yticks([0:100:500])
    xlim([0 45])
    hold on;
    for i = 1:length(TimeCourse.sal.sal) 
        plot(t,TimeCourse.sal.sal(i).Con4,'Linewidth',2)
    end
subplot(4,2,8)
    title('SAL/SAL: Tx 4')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    yticks([0:100:500])
    xlim([0 45])
    hold on;
    for i = 1:length(TimeCourse.sal.sal) 
        plot(t,TimeCourse.sal.sal(i).Tx4,'Linewidth',2)
    end
    
%% NLX/SAL
    %legend fig
figure; plot(t,zeros(length(t)),'Linewidth',7);
    legend A42 A43 A48 A50 A51 A52 A64 A66

nstime = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(4,2,1)
    title('NLX/SAL: Con 1')
    xlabel('time (mins)')
    ylabel('Movement counts')
    %     ylim([0 200])
    ylim([0 500])
    yticks([0:100:500])
    xlim([0 45])
    hold on;
    for i = 1:length(TimeCourse.nlx.sal)
        plot(t,TimeCourse.nlx.sal(i).Con1,'Linewidth',2)
    end
subplot(4,2,2)
    title('NLX/SAL: Tx 1')
    xlabel('time (mins)')
    ylabel('Movement counts')
    %     ylim([0 200])
    ylim([0 500])
    yticks([0:100:500])
    xlim([0 45])
    hold on;
for i = [1 2 4 5 6 7 8] %skip animal 3
        plot(t,TimeCourse.nlx.sal(i).Tx1,'Linewidth',2)
    end
subplot(4,2,3)
    title('NLX/SAL: Con 2')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    yticks([0:100:500])
    xlim([0 45])
    hold on;
    for i = 1:length(TimeCourse.nlx.sal) 
        plot(t,TimeCourse.nlx.sal(i).Con2,'Linewidth',2)
    end
subplot(4,2,4)
    title('NLX/SAL: Tx 2')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    yticks([0:100:500])
    xlim([0 45])
    hold on;
    for i = 1:length(TimeCourse.nlx.sal) 
        plot(t,TimeCourse.nlx.sal(i).Tx2,'Linewidth',2)
    end
subplot(4,2,5)
    title('NLX/SAL: Con 3')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    yticks([0:100:500])
    xlim([0 45])
    hold on;
    for i = 1:length(TimeCourse.nlx.sal) 
        plot(t,TimeCourse.nlx.sal(i).Con3,'Linewidth',2)
    end
subplot(4,2,6)
    title('NLX/SAL: Tx 3')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    yticks([0:100:500])
    xlim([0 45])
    hold on;
    for i = 1:length(TimeCourse.nlx.sal) 
        plot(t,TimeCourse.nlx.sal(i).Tx3,'Linewidth',2)
    end
subplot(4,2,7)
    title('NLX/SAL: Con 4')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    yticks([0:100:500])
    xlim([0 45])
    hold on;
    for i = 1:length(TimeCourse.nlx.sal)
        plot(t,TimeCourse.nlx.sal(i).Con4,'Linewidth',2)
    end
subplot(4,2,8)
    title('NLX/SAL: Tx 4')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    yticks([0:100:500])
    xlim([0 45])
    hold on;
    for i = 1:length(TimeCourse.nlx.sal)
        plot(t,TimeCourse.nlx.sal(i).Tx4,'Linewidth',2)
    end
    
    
%% NLX/KET
%legend 
    figure;
    plot(t,zeros(length(t)),'Linewidth',7);
    legend A41 A44 A45 A47 A53 A58 A61 A67
%timecourse fig
nktime = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(4,2,1)
    title('NLX/KET: Con 1')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    xlim([0 45])
    yticks([0:100:500])
    hold on;
    for i = 1:length(TimeCourse.nlx.ket) 
        plot(t,TimeCourse.nlx.ket(i).Con1,'Linewidth',2)
    end
subplot(4,2,2)
    title('NLX/KET: Tx 1')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    xlim([0 45])
    yticks([0:100:500])
    hold on;
    for i = 1:length(TimeCourse.nlx.ket) 
        plot(t,TimeCourse.nlx.ket(i).Tx1,'Linewidth',2)
    end
subplot(4,2,3)
    title('NLX/KET: Con 2')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    xlim([0 45])
    yticks([0:100:500])
    hold on;
    for i = 1:length(TimeCourse.nlx.ket)
        plot(t,TimeCourse.nlx.ket(i).Con2,'Linewidth',2)
    end
subplot(4,2,4)
    title('NLX/KET: Tx 2')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    xlim([0 45])
    yticks([0:100:500])
    hold on;
    for i = 1:length(TimeCourse.nlx.ket)
        plot(t,TimeCourse.nlx.ket(i).Tx2,'Linewidth',2)
    end
subplot(4,2,5)
    title('NLX/KET: Con 3')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    xlim([0 45])
    yticks([0:100:500])
    hold on;
    for i = 1:length(TimeCourse.nlx.ket)
        plot(t,TimeCourse.nlx.ket(i).Con3,'Linewidth',2)
    end
subplot(4,2,6)
    title('NLX/KET: Tx 3')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    xlim([0 45])
    yticks([0:100:500])
    hold on;
    for i = 1:length(TimeCourse.nlx.ket)
        plot(t,TimeCourse.nlx.ket(i).Tx3,'Linewidth',2)
    end
subplot(4,2,7)
    title('NLX/KET: Con 4')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    xlim([0 45])
    yticks([0:100:500])
    hold on;
    for i = 1:length(TimeCourse.nlx.ket) 
        plot(t,TimeCourse.nlx.ket(i).Con4,'Linewidth',2)
    end
subplot(4,2,8)
    title('NLX/KET: Tx 4')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    xlim([0 45])
    yticks([0:100:500])
    hold on;
    for i = 1:length(TimeCourse.nlx.ket)
        plot(t,TimeCourse.nlx.ket(i).Tx4,'Linewidth',2)
    end
    
    
%% SAL/KET
%legend
    figure; plot(t,zeros(length(t)),'Linewidth',7);
    legend A40 A54 A55 A57 A62 A63 A68 A69
%timecourse
sktime = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(4,2,1)
    title('SAL/KET: Con 1')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    xlim([0 45])
    yticks([0:100:500])
    hold on;
    for i = 1:length(TimeCourse.sal.ket) 
        plot(t,TimeCourse.sal.ket(i).Con1,'Linewidth',2)
    end
subplot(4,2,2)
    title('SAL/KET: Tx 1')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    xlim([0 45])
    yticks([0:100:500])
    hold on;
    for i = 1:length(TimeCourse.sal.ket) 
        plot(t,TimeCourse.sal.ket(i).Tx1,'Linewidth',2)
    end
subplot(4,2,3)
    title('SAL/KET: Con 2')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    xlim([0 45])
    yticks([0:100:500])
    hold on;
    for i = 1:length(TimeCourse.sal.ket) 
        plot(t,TimeCourse.sal.ket(i).Con2,'Linewidth',2)
    end
subplot(4,2,4)
    title('SAL/KET: Tx 2')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    xlim([0 45])
    yticks([0:100:500])
    hold on;
    for i = 1:length(TimeCourse.sal.ket) 
        plot(t,TimeCourse.sal.ket(i).Tx2,'Linewidth',2)
    end
subplot(4,2,5)
    title('SAL/KET: Con 3')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    xlim([0 45])
    yticks([0:100:500])
    hold on;
    for i = 1:length(TimeCourse.sal.ket) 
        plot(t,TimeCourse.sal.ket(i).Con3,'Linewidth',2)
    end
subplot(4,2,6)
    title('SAL/KET: Tx 3')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    xlim([0 45])
    yticks([0:100:500])
    hold on;
    for i = 1:length(TimeCourse.sal.ket) 
        plot(t,TimeCourse.sal.ket(i).Tx3,'Linewidth',2)
    end
subplot(4,2,7)
    title('SAL/KET: Con 4')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    xlim([0 45])
    yticks([0:100:500])
    hold on;
    for i = 1:length(TimeCourse.sal.ket) 
        plot(t,TimeCourse.sal.ket(i).Con4,'Linewidth',2)
    end
subplot(4,2,8)
    title('SAL/KET: Tx 4')
    xlabel('time (mins)')
    ylabel('Movement counts')
    ylim([0 500])
    xlim([0 45])
    yticks([0:100:500])
    hold on;
    for i = 2:8 
    %     for i = 1:length(TimeCourse.sal.ket) 
        plot(t,TimeCourse.sal.ket(i).Tx4,'Linewidth',2)
    end
