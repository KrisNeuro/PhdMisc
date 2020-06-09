    function [traceresults, traceout, outfig, outfig2, outfig3, outfig4] = traceanalysisscript_yesfigs(filename)

%% Threshold Parameters
 spikedetectthreshold = -30; % Threshold for spike detection (mV)
 minspikedistance=10; % The minimum spike distance in samples, used as find peaks filter.
 minspikeheight=5;
 maxpeakwidth=3000;
 htcstart=1; %Default 1 
 htctime=100;
 memcaptime=50; %Default 50  Membrane capacitance(?) 
 
%% READ DATA FILE
    [traceout, si] = abf2load(filename); %load abf file recorded from pclamp 10. si is in microseconds
 
    [totalrows,~,numberoftraces] = size(traceout);
    [~,name,~] =fileparts(filename);
    
    timefactor=1000/si; % # of samples per ms
    
    % time vector
    time(:,1) = (1:totalrows) / timefactor; %in fractions of ms
    
    % voltage traces matrix
    voltagetraces(:,:) = traceout(:,1,:);  %2D matrix of unfiltered voltage: rows=time, colums=sweeps
     
    % Voltage Trace Filter
    sk = 3; % filter parameter
    sf = 81; % filter parameter
    voltagetraces = sgolayfilt(voltagetraces,sk,sf); % voltage traces filtered
        % rows = timestamps; colums = sweeps

    % current traces matrix 
    currenttraces(:,:) = traceout(:,2,:);
    [CTrows,~] = size(currenttraces);
     
%% Calculate Current Injection on and off times in samples and ms 
 % filtered current trace
    sk = 10; % filter parameter 
    sf=81;
    sf = 10*timefactor+1; % filter parameter
    maxIcol= find(max(abs(currenttraces))==max(max(abs(currenttraces))));
    filteredcurrenttiming(:,1) = sgolayfilt(currenttraces(:,maxIcol(1)),sk,sf); % trace filter
    
 % Calculate Derivative of current trace 
 currenttraces_1 = zeros (CTrows,1); % preallocate
 for i_h=2:CTrows
    currenttraces_1(i_h,1) = abs((filteredcurrenttiming(i_h,1) - filteredcurrenttiming(i_h-1,1)) / (time(i_h,1) - time(i_h-1,1)));
 end
 
currentonthreshold = (max(currenttraces_1)*0.8); % Assumes current on derivitive is within 80% of current off derivitve
[~, currentpkslocs]= findpeaks(currenttraces_1,'MinPeakHeight',currentonthreshold, 'MinPeakDistance',10*timefactor); % assumes current injection is greater than 10ms
currentons =(currentpkslocs(1)); % current on time in samples
currentoffs=(currentpkslocs(2)); % current off time in samples
currentpkslocsms= currentpkslocs /timefactor;        %converting from samples to ms
currentinjectiontime = (currentoffs - currentons) / timefactor;  % Duration of current injection in samples
currentonms = (currentpkslocsms(1)) ;    % current on time in ms
currentoffms = (currentpkslocsms(2));   % current off time in ms
traceend=length(time(:,1)); %end of sweep, in timestamps

% Normalize Current Injection Values
    currenttracesnormalize = mean(mean(currenttraces((1):(currentons-(round(currentons/3))),:)));
    currenttraces = currenttraces -currenttracesnormalize;
    % Calculate current values in picoamps
    currentvalues= zeros (1, numberoftraces); % preallocate
    for i=1:numberoftraces %loop thru sweeps
        currentvalues(i)=mean(currenttraces(currentons:currentoffs,i));       
    end

% To be used for setting axes limits in figures
minvoltage=min(min(voltagetraces));
maxvoltage=max(max(voltagetraces));
mincurrent=min(min(currenttraces));
maxcurrent=max(max(currenttraces));
currentonoffcheck = (zeros(CTrows,1)); % Values to plot current on and off times
currentonoffcheck(currentons)=(currentvalues(maxIcol));
currentonoffcheck(currentoffs)=(currentvalues(maxIcol));

%% resting potential: average from first 100ms of all sweeps
restingpotential = mean(mean(voltagetraces(1:1000,:)));
 
%% FIND SPIKES/BURSTS
%preallocate cells
spikepks=cell(1);
spikelocs=cell(1);
halfpromwidth=cell(1);
proms=cell(1);

spikepkscit=cell(1);
spikelocscit=cell(1);
halfpromwidthcit=cell(1);
promscit=cell(1);

spikepksmin=cell(1);
spikelocsmin=cell(1);
halfpromwidthmin=cell(1);
promsmin=cell(1);

spikepksfilt=cell(1);
spikelocsfilt=cell(1);
halfpromwidthfilt=cell(1);
promsfilt=cell(1);


for i_h=1:numberoftraces
 % Total Number of spikes across whole trace
[spikepks{:,i_h}, spikelocs{:,i_h}] = find(voltagetraces(:,i_h)>spikedetectthreshold);
    if ~isempty(spikepks{:,i_h}(:)) %if it contains data
        [spikepks{:,i_h}, spikelocs{:,i_h}, halfpromwidth{:,i_h}, proms{:,i_h}]=findpeaks(voltagetraces(:,i_h),'MinPeakHeight',spikedetectthreshold, 'MinPeakDistance', minspikedistance, 'MinPeakProminence',minspikeheight, 'Annotate','extents');
    
         % Total Number of spikes during curent injection only
         [spikepkscit{:,i_h}, spikelocscit{:,i_h}, halfpromwidthcit{:,i_h}, promscit{:,i_h}]=findpeaks(voltagetraces(currentons:currentoffs,i_h),'MinPeakHeight',spikedetectthreshold, 'MinPeakDistance', minspikedistance, 'MinPeakProminence',minspikeheight, 'Annotate','extents');
            %cit = during current injection time
         [spikepksmin{:,i_h}, spikelocsmin{:,i_h}, halfpromwidthmin{:,i_h}, promsmin{:,i_h}]=findpeaks(voltagetraces(currentons:currentoffs,i_h),'MinPeakHeight',spikedetectthreshold, 'MinPeakDistance', minspikedistance, 'MinPeakProminence',minspikeheight, 'Annotate','extents');

         % Spike widths. Stronger filtering of spikes
         [spikepksfilt{:,i_h}, spikelocsfilt{:,i_h}, halfpromwidthfilt{:,i_h}, promsfilt{:,i_h}]=findpeaks(voltagetraces(:,i_h),'MinPeakHeight',spikedetectthreshold, 'MinPeakDistance', minspikedistance, 'MinPeakProminence',minspikeheight, 'MaxPeakWidth', maxpeakwidth, 'Annotate','extents');
         
    else % if there are no spikes for that sweep
        [spikepks{i_h}]=0;
        [spikelocs{:,i_h}] =0;
        [halfpromwidth{:,i_h} ]=0;
        [proms{:,i_h}] = 0;
        [spikepkscit{i_h}]=0;
        [spikelocscit{i_h}]=0;
        [halfpromwidthcit{:,i_h}]=0;
        [promscit{:,i_h}]=0;
        [spikepksmin{i_h}]=0; 
        [spikelocsmin{i_h}]=0;
        [halfpromwidthmin{:,i_h}]=0;
        [promsmin{:,i_h}]=0;
        [spikepksfilt{i_h}]=0;
        [spikelocsfilt{i_h}]=0;
        [halfpromwidthfilt{:,i_h}]=0;
        [promsfilt{:,i_h}]=0;
    end
end


for i=1:numberoftraces
    if ~isempty(spikepks(i)) %if there are spikes
        totalspikes(i)=numel(cell2mat(spikepks(i)));  %total spikes for whole sweep
        totalspikescit(i)=numel(cell2mat(spikepkscit(i))); %total spikes during current injection
    else %if no spikes
        totalspikes(i) = 0;
        totalspikescit(i) = 0;
    end
end

%preallocate space
spikepksmatrix=NaN(max(totalspikes),numberoftraces);
spikelocsmatrix=NaN(max(totalspikes),numberoftraces);
spikelocsmatrixcit=NaN(max(totalspikescit),numberoftraces);
minmatrix=NaN(max(totalspikescit),numberoftraces);
spikewidthmatrix=NaN(max(totalspikes),numberoftraces);
promsmatrix=NaN(max(totalspikes),numberoftraces);

% Repackage all spikes into concatenated matrices; one for spike peaks, one
% for spike peak times
for i=1:numberoftraces %loop thru sweeps
    for i_a=1:(numel(cell2mat(spikepks(i))))    %loop thru spikes per sweep
        spikepksmatrix(i_a,i)=spikepks{i}(i_a);
    end
    for i_a=1:(numel(cell2mat(spikelocs(i))))
        spikelocsmatrix(i_a,i)=spikelocs{i}(i_a);
    end
% For spikes during current injection
    for i_a=1:(numel(cell2mat(spikelocscit(i))))
        spikelocsmatrixcit(i_a,i)=spikelocscit{i}(i_a)+currentons;
    end
%for times, spike widths, spike prominence
    for i_a=1:(numel(cell2mat(spikelocsmin(i))))
        minmatrix(i_a,i)= spikelocsmin{i}(i_a)+currentons;
    end

    for i_a=1:(numel(cell2mat(halfpromwidthfilt(i))))
        spikewidthmatrix(i_a,i)=(halfpromwidthfilt{i}(i_a))/timefactor;
    end

    for i_a=1:(numel(cell2mat(proms(i))))       
        promsmatrix(i_a,i)=proms{i}(i_a);
    end
end

interspikeint=diff(spikelocsmatrixcit)./timefactor; % Only measured during current injection time

spikelocsplot= cellfun(@(v) v./timefactor, spikelocs, 'UniformOutput', false);
 
% Preallocation
    adaptrate_spikewidth = 0;  % Array to save the adaptation rate of the spike durations
    adaptrate_interspikeint = NaN(1,numberoftraces);  % Array to save the adaptation rate of the interspikes intervals
    adaptrate_spikeamp = 0;  % Array to save the adaptation rate of the spikes amplitudes 
    spikefrequency=zeros(1,numberoftraces);
    ave_spikear = NaN(1,numberoftraces);
    ave_isi = NaN(1,numberoftraces);  
    std_isi = NaN(1,numberoftraces);     
    ave_sw = NaN(1,numberoftraces);  
    std_sw =  NaN(1,numberoftraces);     
    ave_amps = NaN(1,numberoftraces);
    first_amp = NaN(1,numberoftraces);
    first_sw = NaN(1,numberoftraces); 
    std_amps =  NaN(1,numberoftraces);   
    sag_vmin= NaN(1,numberoftraces);
    sag_mintime1= NaN(1,numberoftraces);
    sag_mintime2= NaN(1,numberoftraces);
    ave_sag_amp= NaN(1,numberoftraces);
    sag_amp= NaN(1,numberoftraces);
    sag_ratio= NaN(1,numberoftraces);
    vdrop= NaN(1,numberoftraces);
    sag_vend= NaN(1,numberoftraces);
    reboundlatency=NaN(1,numberoftraces);
    reboundspike=NaN(1,numberoftraces);
    spikelatency=NaN(1,numberoftraces);
    spikelatencytemp=NaN(1,numberoftraces);
    relatencytemp=NaN(1,numberoftraces);
    rebounddepol= NaN(1,numberoftraces);
    fitcoeffs= [0,0,0,0];
    inputresistance=[];
    timeconstant=NaN;
    membranecapacitance=[];
    
%Voltage during current pulse (time x sweeps)
    voltagecurrenttime(:,:) = voltagetraces(currentons:currentoffs,:);
   
    sagavetimeon=0; 
    sagavetimeoff=0;  
    sagavetimeend=0;
    
    vmintimeon=0; 
    vmintimeoff=0; 

    
   if (traceend-currentoffs) < (300*timefactor)
        afterdepolchecktime= traceend;
   else
       afterdepolchecktime =currentoffs+(300*timefactor);
   end
       
   
for i_d=1:numberoftraces
    % Look for the sag if the applied current is negative    
    if currentvalues(i_d) < -5;
           [sag_vmin(i_d), sag_mintime1(i_d)] = min(voltagecurrenttime(:,i_d));    % Minimum Voltage in the interval when the current is applied
           sag_mintime2(i_d)= sag_mintime1(i_d)+currentons ;
           ave_sag_amp(i_d) = mean(voltagetraces((sag_mintime2(i_d)-(sagavetimeon*timefactor)):(sag_mintime2(i_d)+(sagavetimeoff*timefactor)),i_d));
           vmin(i_d) = mean(voltagetraces((sag_mintime2(i_d)-(vmintimeon*timefactor)):(sag_mintime2(i_d)+(vmintimeoff*timefactor)),i_d));
           vdrop(i_d) = voltagetraces(currentons-10, i_d)-vmin(i_d);
           sag_vend(i_d) = mean(voltagetraces(((currentoffs-(sagavetimeend*timefactor)):currentoffs),i_d)); % Volatage at end of current injection
           sag_amp(i_d) = abs(ave_sag_amp(i_d)-voltagetraces(currentoffs,i_d)); % Amplitude of sag
           sag_ratio(i_d) = sag_amp(i_d)/abs(ave_sag_amp(i_d)); % sag ratio
           rebounddepol(i_d) = (abs(max(voltagetraces(currentoffs:afterdepolchecktime,i_d))-(voltagetraces(currentons-10,i_d))));
    end

    % CALCULATE ADAPTATION RATES
    if totalspikescit(i_d)>1
        isicol = interspikeint(:,i_d);
        isicol2 = isicol(~isnan(isicol));
        adaptrate_interspikeint(i_d) = (isicol2(end)) / (isicol2(1));
        for i_h=1: numel(interspikeint(:,i_d))
            adaptrate_interspikeint2(i_h,i_d) = (interspikeint(i_h,i_d)) / (interspikeint(1,i_d));
        end

        for i_h=1: numel(totalspikescit(:,i_d))
             adaptrate_spikewidth(i_h,i_d) = spikewidthmatrix(i_h,i_d)/( spikewidthmatrix(1,i_d)); 
             adaptrate_spikeamp(i_h,i_d) = promsmatrix(i_h,i_d) / promsmatrix(1,i_d); 
        end    

        % CALCULATE AVERAGE AND STANDARD DEVIATIONS of interspike interval
        ave_spikear(i_d) = nanmean(adaptrate_interspikeint(:,i_d));
        ave_isi(i_d) = nanmean(interspikeint(:,i_d));  % Calculate average interspikes intervals 
        std_isi(i_d) = nanstd(interspikeint(:,i_d));   % Calculate standard deviation interspikes intervals 
    end

    if totalspikescit(i_d)>0
        ave_sw(i_d) = nanmean(spikewidthmatrix(:,i_d));   % Calculate average width of spikes
        std_sw(i_d) =  nanstd(spikewidthmatrix(:,i_d));   % Calculate standard deviation of spike width
        first_sw(i_d) = spikewidthmatrix(1,i_d);

        ave_amps(i_d) = nanmean(promsmatrix(:,i_d));    % Calculate average spikes amplitudes
        std_amps(i_d) =  nanstd(promsmatrix(:,i_d));    % Calculate standard deviation of spikes amplitudes 
        if minmatrix(1,i_d) > 0
            first_amp(i_d) = voltagetraces(spikelocsmatrixcit(1,i_d),i_d) - voltagetraces(minmatrix(1,i_d),i_d);
        end
    end

    % Calculate spike latency after current injection
    if currentvalues(i_d) > 5 && totalspikes(i_d) > 0;
        for i_h=1: numel(spikelocsmatrix(i_d))
                if spikelocsmatrix(i_h,i_d) > currentons
                    spikelatencytemp(i_h,i_d) = spikelocsmatrix(i_h,i_d); 
                end
        end
        spikelatency(i_d) =  ((min(spikelatencytemp(:,i_d))) - currentons) /timefactor ; % Calculate spike latency  
    end

    % Calculate the rebound latency
    if currentvalues(i_d) < -5 && totalspikes(i_d) > 0;
        for i_h=1: numel(spikelocsmatrix(i_d))
            if spikelocsmatrix(i_h,i_d)> currentoffs
                relatencytemp(i_h,i_d) = spikelocsmatrix(i_h,i_d);
            end
        end
        reboundlatency(i_d) = ((min(relatencytemp(:,i_d))) - currentoffs) /timefactor; % Calculate rebound latency
    end

    if currentvalues(i_d) < -5 && totalspikes(i_d) > 0;
        reboundspike(i_d) = 1;
    elseif currentvalues(i_d) < -5 && totalspikes(i_d) == 0;
        reboundspike(i_d) = 0;
    end

    % calculate spike frequency
    if totalspikescit(1,i_d) > 0
        spikefrequency(i_d)= (totalspikescit(i_d) * 1000)/ currentinjectiontime; 

    else
        spikefrequency(1,i_d) = NaN;    
    end
end 

%% Passive Membrane properties
% May need to be adjusted for passive step protocols
for i=1:numberoftraces
    if currentvalues(i) < -7  && currentvalues(i) > -13   % figure out the right voltage trace to check input resistance. first negative current injection 
        pmcol10(i) = i;
    else
        pmcol10(i)=0;
    end 
    
    if currentvalues(i) < -23  && currentvalues(i) > -27   % figure out the right voltage trace to check input resistance. first negative current injection 
        pmcol25(i) = i;
    else
        pmcol25(i)=0;
    end 
end

pmcol10= max(pmcol10);
pmcol25= max(pmcol25);

if pmcol10 > 0
    pmcol=pmcol10;
else
    pmcol=pmcol10;
end

if pmcol > 0
    fittime=time(currentons:currentons+(memcaptime*timefactor),1)-time(currentons);
    fitvoltage= (voltagetraces(currentons:(currentons+(memcaptime*timefactor)),pmcol)) - (min(voltagetraces(currentons:(currentons+(memcaptime*timefactor)),pmcol)));
    [fit1]= fit(fittime,fitvoltage,'exp2'); % Fit of two-term exponential to first 10 ms of voltage trace following current injection. f2(x) = a*exp(b*x) + c*exp(d*x) 
    fitcoeffs= coeffvalues(fit1); % Values of fit. [a, b] b is decay constant
    if abs(fitcoeffs(2)) < 0.15
        timeconstant= abs(((1 / (fitcoeffs(4))) *10^-3)); % Time constant tau (s)
        inputresistance = ((min(voltagetraces(:,pmcol))-restingpotential)*10^-3) / ((currentvalues(pmcol))*10^-12);  % Membrane resistance in Ohms R(ohms) = E(volts)/I(amps)  
        membranecapacitance= (timeconstant/inputresistance) *10^12;  % Membrane capacitance in picoFarads Cm = tau / Rm
        if membranecapacitance > 500
            membranecapacitance = [];
        end
        inputresistance= inputresistance *10^-6;  % ohms to megaohms
        timeconstant= timeconstant *10^3; % seconds to miliseconds
    end  
end
   
%% Spike and Rebound Threshold
rethresholdall=ones(1,numberoftraces);
spikethresholdall=ones(1,numberoftraces);
firstspikeindex=0;
for i=1:numberoftraces
    if reboundlatency(i)>1
        rethresholdall(i)=currentvalues(i);
    else
        rethresholdall(i)=NaN;
    end
    reboundthreshold=max(rethresholdall);  
   
    if spikelatency(i)>1
        spikethresholdall(i)=currentvalues(i);
        firstspikeindex=i;
    else
        spikethresholdall(i)=NaN;
    end
    spikethreshold=min(spikethresholdall);
end   

% "Final" values
sagratio200 = [];
sagratio300 = [];
vdrop200 = [];
htctimeconstant= [];
reboundlatencyfinal = [];
spikelatencyfinal= [];
arfinal= [];
spikefrequencyfinal= [];
ave_isifinal= [];
std_isifinal= [];
afterdpfinal= [];
hypcol=0;
depcol=0;

maxrise=1;
maxfall=1;
apreliability=1;

for i=1:numberoftraces
    if currentvalues(i) < -190 && currentvalues(i) > -210
        hypcol=i;
        sagratio200= sag_ratio(i);
        vdrop200= vdrop(i);
        reboundlatencyfinal= reboundlatency(i);
        htcfittime=time(currentons+(htcstart*timefactor):currentons+(htctime*timefactor),1)-time(currentons);
        htcfitvoltage= (voltagetraces(currentons+(htcstart*timefactor):(currentons+(htctime*timefactor)),i)) - (min(voltagetraces(currentons+(htcstart*timefactor):(currentons+(htctime*timefactor)),i)));
        [htcfit]= fit(htcfittime,htcfitvoltage,'exp2'); % Fit of single exponential to first 10 ms of voltage trace following current injection. y=ae^bx 
        htcfitcoeffs= coeffvalues(htcfit); % Values of fit. [a, b] b is decay constant 
        htctimeconstant= abs(((1 / (htcfitcoeffs(4))) *10^-3)); % Time constant tau (s)
        htctimeconstant= htctimeconstant *10^3; % seconds to miliseconds
        afterdpfinal= rebounddepol(i);
    end

    if currentvalues(i) < -295 && currentvalues(i) > -305
        sagratio300= sag_ratio(i);
    end

    if ~isempty(strfind(filename,'x')) 
        if currentvalues(i) > 90 && currentvalues(i) < 110
            spikelatencyfinal= spikelatency(i);
        end
    end

    if ~isempty(strfind(filename,'int')) 
        if currentvalues(i) > 90 && currentvalues(i) < 110
            spikelatencyfinal= spikelatency(i);
        end
    end

    if ~isempty(strfind(filename,'ra'))
        if currentvalues(i) > 190 && currentvalues(i) < 210
            spikelatencyfinal= spikelatency(i);
        end
    end

    if currentvalues(i) > 195 && currentvalues(i) < 210
        depcol=i;
        spikefrequencyfinal=spikefrequency(i);
        if isnan(spikefrequencyfinal)
            spikefrequencyfinal= 0.1;
        end
        arfinal= adaptrate_interspikeint(i);
        ave_isifinal=ave_isi(i); 
        std_isifinal= std_isi(i);

    end
end

% Clear out NaNs
if isnan(timeconstant)
    timeconstant= [];
end

if isnan(spikelatencyfinal)
    spikelatencyfinal= [];
end

if isnan(reboundlatencyfinal)
    reboundlatencyfinal= [];
end

if isnan(reboundthreshold)
    reboundthreshold= [];
end  

if isnan(spikethreshold)
    spikethreshold= [];
end   

if isnan(spikefrequencyfinal)
    spikefrequencyfinal= [];
end  

if isnan(arfinal)
    arfinal= [];
end  

if isnan(ave_isifinal)
    ave_isifinal= [];
end 

if isnan(std_isifinal)
    std_isifinal= [];
end

finalavespikewidth= nanmean(first_sw); 
if finalavespikewidth== 0
    finalavespikewidth= [];
end 

if isnan(finalavespikewidth)
    finalavespikewidth= [];
end  

finalavespikeamplitude= nanmean(first_amp);
if finalavespikeamplitude== 0
    finalavespikeamplitude= [];
end

if isnan(finalavespikeamplitude)
    finalavespikeamplitude= [];
end  

%-------------------------------------------------------------------------
% Reorder data based on sort(currentvalues)
% used in case current injections are not in sequence

[currentvalues_ordered, cv_order] = sort(currentvalues);
spikefrequency_ordered= spikefrequency(cv_order);
reboundlatency_ordered= reboundlatency(cv_order);
sag_ratio_ordered= sag_ratio(cv_order);
ave_spikear_ordered= ave_spikear(cv_order);
ave_isi_ordered= ave_isi(cv_order);
sag_vend_ordered= sag_vend(cv_order);
sag_vmin_ordered= sag_vmin(cv_order);
vmin_amp_ordered= ave_sag_amp(cv_order);
std_amps_ordered= std_amps(cv_order);
ave_amps_ordered= ave_amps(cv_order);
std_sw_ordered= std_sw(cv_order);
ave_sw_ordered= ave_sw(cv_order);
spikelatency_ordered= spikelatency(cv_order);
std_isi_ordered= std_isi(cv_order);
totalspikes_ordered= totalspikes(cv_order);
vdrop_ordered= vdrop(cv_order);
reboundspike_ordered= reboundspike(cv_order);
rebounddepol_ordered= rebounddepol(cv_order);

%% Output: traceresults
% Setup and Headers
traceresults = cell(15,numberoftraces);
% Col 1 Names
    % Membrane Properties
    traceresults{1,1}=   'Filename'; % filename 
    traceresults{2,1}=   'Resting Potential'; % resting potential
    traceresults{3,1}=   'Time Constant'; % tau
    traceresults{4,1}=   'Input Resistance'; %input resistance
    traceresults{5,1}=   'Membrane Capacitance'; % membrane capacitace
    traceresults{6,1}=   'Fit Coeffs'; % membrane capacitace

    % Spike Timing
    traceresults{7,1}=   'Spike Frequency'; % 
    traceresults{8,1}=   'Adaptation Rate'; % 
    traceresults{9,1}=   'Average Interspike Interval'; % 
    traceresults{10,1}=  'Standard Deviation Interspike Interval'; % 
    
% Col 3 Names
    % Action Potential Properties
    traceresults{2,3}=   'Spike Threshold'; % 
    traceresults{3,3}=   'Spike Latency'; % 
    traceresults{4,3}=   'Max Rise'; % 
    traceresults{5,3}=   'Max Fall'; % 
    traceresults{6,3}=   'Spike Width'; % 
    traceresults{7,3}=   'Spike Reilability'; % 
    traceresults{8,3}=   'Spike Amplitude'; % 
    
% Col 5 Names
    % Hyperpolarizing Response
    traceresults{2,5}=   'Sag Ratio -200'; % Sag ratio
    traceresults{3,5}=   'Hyper Time Constant'; % Sag ratio
    traceresults{4,5}=   'rethresh'; % Rebound spike threshold current
    traceresults{5,5}=   'Rebound Spike Timing'; % sag ratio at -200
    traceresults{6,5}=   'After Depolarization'; % After Depolarization 
    traceresults{7,5}=   'Voltage Drop';
   
% Full Results   
    traceresults{13,1}=  'currentvalues';
    traceresults{14,1}=  'Total Spikes';
    traceresults{15,1}=  'Spike Frequency';
    traceresults{16,1}=  'Average Frequency Adaptation Rate';
    traceresults{17,1}=  'Average Interspike Interval';
    traceresults{18,1}=  'Standard Deviation Interspike Interval';
    traceresults{19,1}=  'Spike Latency';
    
    traceresults{22,1}=  'Average Spike Width'; 
    traceresults{23,1}=  'Standard Deviation Spike Width';
   
    traceresults{25,1}=  'Average Spike Amplitude';
    traceresults{26,1}=  'Standard Deviation Spike Amplitude';
    
% Sag
    traceresults{27,1}=  'Sag Ratio';
    traceresults{28,1}=  'Sag Amplitude';
    traceresults{29,1}=  'Rebound Latency';
    traceresults{30,1}=  'Rebound Spike?';
    traceresults{31,1}=  'Voltage Drop';
    traceresults{32,1}=  'Rebound Depolarization';

%% Assign Main Results to traceresults
% Col 2 main results
    % Membrane Properties
    traceresults{1,2}=   filename; % filename 
    traceresults{2,2}=   restingpotential; % resting potential
    traceresults{3,2}=   timeconstant; % tau
    traceresults{4,2}=   inputresistance; %input resistance
    traceresults{5,2}=   membranecapacitance; % cell capacitace
    traceresults{6,2}=   num2str(fitcoeffs);
    % Spike Timing
    traceresults{7,2}=   spikefrequencyfinal; % Spike Frequency
    traceresults{8,2}=   arfinal; % Adaptation Rate
    traceresults{9,2}=   ave_isifinal; % Average Interspike Interval
    traceresults{10,2}=  std_isifinal; % Standard Deviation Interspike Interval
    
% Col 4
    % Action Potential Properties
    traceresults{2,4}=   spikethreshold; % 
    traceresults{3,4}=   spikelatencyfinal; % 
    traceresults{4,4}=   maxrise; % 
    traceresults{5,4}=   maxfall; % 
    traceresults{6,4}=   finalavespikewidth; % 
    traceresults{7,4}=   apreliability; % 
    traceresults{8,4}=   finalavespikeamplitude; % 
    
% Col 6
    % Hyperpolarizing Response
    traceresults{2,6}=   sagratio200; %
    traceresults{3,6}=   htctimeconstant; %
    traceresults{4,6}=   reboundthreshold; % 
    traceresults{5,6}=   reboundlatencyfinal; %
    traceresults{6,6}=   afterdpfinal; %
    traceresults{7,6}=   vdrop200;
    
    % Full Results
for i_g=1:numberoftraces 
    traceresults{13,i_g+1}=  currentvalues_ordered(1,i_g);
    traceresults{14,i_g+1}=  totalspikes_ordered(1,i_g);
    traceresults{15,i_g+1}=  spikefrequency_ordered(1,i_g);
    traceresults{16,i_g+1}=  ave_spikear_ordered(1,i_g);
    traceresults{17,i_g+1}=  ave_isi_ordered(1,i_g);
    traceresults{18,i_g+1}=  std_isi_ordered(1,i_g);
    traceresults{19,i_g+1}=  spikelatency_ordered(1,i_g);
    
    traceresults{22,i_g+1}=  ave_sw_ordered(1,i_g); 
    traceresults{23,i_g+1}=  std_sw_ordered(1,i_g);
   
    traceresults{25,i_g+1}=  ave_amps_ordered(1,i_g);
    traceresults{26,i_g+1}=  std_amps_ordered(1,i_g);
    
    % Hyperpolarization
    traceresults{27,i_g+1}=  sag_ratio_ordered(1,i_g);
    traceresults{28,i_g+1}=  vmin_amp_ordered(1,i_g);
    traceresults{29,i_g+1}=  reboundlatency_ordered(1,i_g);
    traceresults{30,i_g+1}=  reboundspike_ordered(1,i_g);
    traceresults{31,i_g+1}=  vdrop_ordered(1,i_g);
    traceresults{32,i_g+1}=  rebounddepol_ordered(1,i_g); 
end


%--------------------------------------------------------------------------
%% Plot Figure Main Data
titlefontsize=20;
labelfontsize= 17;
axisfontsize= 11;
figfont='Arial';
figfontweight= 'bold';
    
%% Figure 1
outfig=figure('name',name,'position',[1 1 1474 744]);

subplot(2,2,1); % Voltage Traces
for i_h=1:numberoftraces
   plot(time(:,1),voltagetraces(:,i_h), 'color', [0.8 0.8 0.8], 'LineWidth',1)
   hold on
   
end
for i_h=1:numberoftraces
   if currentvalues(i_h) > 190 && currentvalues(i_h) < 210
       plot(time(:,1), voltagetraces(:,i_h), 'color', [56,83,163] ./255, 'LineWidth',2);
       %plot(cell2mat(spikelocsplot(:,i_h)),cell2mat(spikepks(:,i_h)),'o')
   end
   hold on
   if currentvalues(i_h) > -210 && currentvalues(i_h) < -190
       plot(time(:,1), voltagetraces(:,i_h), 'k', 'LineWidth',2)
       %plot(cell2mat(spikelocsplot(:,i_h)),cell2mat(spikepks(:,i_h)),'o')
   end
end
%plot(-50,(currentons-(round(currentons/2)):(currentons-(round(currentons/5)))),'color', 'r', 'LineWidth',5);
plot(time(currentons:currentoffs,1),voltagetraces(currentons:currentoffs,1),'color', 'r', 'LineWidth',5);

set(gca, 'FontSize', axisfontsize, 'FontName', figfont, 'FontWeight', figfontweight)
title('Voltage Traces','FontSize', titlefontsize, 'FontName', figfont, 'FontWeight', figfontweight)
xlabel('Time (ms)','FontSize', labelfontsize, 'FontName', figfont, 'FontWeight', figfontweight)
ylabel('Voltage (mV)','FontSize', labelfontsize, 'FontName', figfont, 'FontWeight', figfontweight)
%axis([currentonms-200,currentoffms+250, -150, 50])


subplot(2,2,2); % Spike Frequency
   [hAx,~,~]= plotyy(currentvalues_ordered(1,:),spikefrequency_ordered(1,:),currentvalues_ordered(1,:),ave_spikear_ordered(1,:));
   
set(hAx(1),'YLim',[0,425]);
set(hAx(1),'YTick',[0:50:425]);
set(hAx(2),'YLim',[0,4]);
set(hAx(2),'YTick',[0:.5:4]);
set(hAx(1),'XLim',[0,425]);
set(hAx(2),'XLim',[0,425]);
set(hAx(1),'box','off');

set(gca, 'FontSize', axisfontsize, 'FontName', figfont, 'FontWeight', figfontweight)
title('Average Spike Frequency and Adaptation Rate','FontSize', titlefontsize, 'FontName', figfont, 'FontWeight', figfontweight) ;
xlabel('Current Injection (pA)','FontSize', labelfontsize, 'FontName', figfont, 'FontWeight', figfontweight) ; 
ylabel(hAx(1),'Spike Frequency (Hz)','FontSize', labelfontsize, 'FontName', figfont, 'FontWeight', figfontweight) ; %left y-axis
ylabel(hAx(2),'Adaptation Rate','FontSize', labelfontsize, 'FontName', figfont, 'FontWeight', figfontweight) ; % right y-axis
 
subplot(2,2,3); % Sag
   hold on
   sagIend=plot(currentvalues_ordered(1,:),sag_vend_ordered(1,:),'r', 'LineWidth',2); 
   sagnadir=plot(currentvalues_ordered(1,:),sag_vmin_ordered(1,:),'color', [153 0 153] ./255, 'LineWidth',2);
set(gca, 'FontSize', axisfontsize, 'FontName', figfont, 'FontWeight', figfontweight)
title('Sag','FontSize', titlefontsize, 'FontName', figfont, 'FontWeight', figfontweight)
xlabel('Current Injection (pA)','FontSize', labelfontsize, 'FontName', figfont, 'FontWeight', figfontweight)
ylabel('Voltage (mV)','FontSize', labelfontsize, 'FontName', figfont, 'FontWeight', figfontweight)
axis([-260, 0, -130, -70])
legend([sagIend,sagnadir],{'Vend','Vmin'} ,'Location','northeast','FontSize', 15);
     
subplot(2,2,4);
[hAx,~,~]= plotyy(currentvalues_ordered(1,:),reboundlatency_ordered(1,:),currentvalues_ordered(1,:),sag_ratio_ordered(1,:));
set(hAx(1),'YLim',[0,250]);
set(hAx(1),'YTick',[0:25:250]);
set(hAx(2),'YLim',[0,.3]);
set(hAx(2),'YTick',[0:.05:.3]);
set(hAx(1),'XLim',[-400,0]);
set(hAx(2),'XLim',[-400,0]);
set(hAx(1),'box','off')

set(gca, 'FontSize', axisfontsize, 'FontName', figfont, 'FontWeight', figfontweight)
title('Rebound Timing and Sag Ratio vs. I applied', 'FontSize', titlefontsize, 'FontName', figfont, 'FontWeight', figfontweight);
ylabel(hAx(1),'Rebound Timing (ms)','FontSize', labelfontsize, 'FontName', figfont, 'FontWeight', figfontweight); % left y-axis 'color','k'
ylabel(hAx(2),'Sag Ratio','FontSize', labelfontsize, 'FontName', figfont, 'FontWeight', figfontweight); % right y-axis
xlabel('Current Injection (pA)','FontSize', labelfontsize, 'FontName', figfont, 'FontWeight', figfontweight)

%% Figure 2
outfig2=figure('name',name,'position',[1 1 1474 744]);

subplot(2,2,1); % Passive Membrane Fit 
if pmcol > 0
    plot(fit1,fittime, fitvoltage)
    hold on
end
if hypcol>0
    plot(htcfit,htcfittime,htcfitvoltage)
end
title('Passive Membrane Fit', 'FontSize', labelfontsize);    

subplot(2,2,2); % Current Injections Values
for i_h=1:numberoftraces
   plot(time(:,1),currenttraces(:,i_h), 'color', [0.8 0.8 0.8], 'LineWidth',1)
   hold on
   plot(time(:,1),currentonoffcheck, 'r', 'LineWidth',1)  
   plot(currentpkslocsms,(currentvalues(maxIcol)),'o')
end
axis([currentonms-50,currentoffms+50, mincurrent-50, maxcurrent+50]) 
title('Current Injection Values', 'FontSize', labelfontsize);    

subplot(2,2,3); % 
if pmcol > 0
plot(time(:,1),voltagetraces(:,pmcol), 'color', [0.8 0.8 0.8], 'LineWidth',1)
   hold on
   plot(time(currentons:currentons+(memcaptime*timefactor),1),voltagetraces(currentons:currentons+(memcaptime*timefactor),pmcol), 'r', 'LineWidth',2)
axis([currentonms-5,currentoffms+50, (min(voltagetraces(pmcol)))-15, restingpotential+2])
end 

if hypcol > 0
plot(time(:,1),voltagetraces(:,hypcol), 'color', [0.8 0.8 0.8], 'LineWidth',1)
   hold on
   plot(time(currentons+(htcstart*timefactor):currentons+(htctime*timefactor),1),voltagetraces(currentons+(htcstart*timefactor):currentons+(htctime*timefactor),hypcol), 'r', 'LineWidth',2)
axis([currentonms-5,currentonms+200, sag_vmin(hypcol)-5, restingpotential+2])
xlabel('time (ms)')
ylabel('mV')
title('Unknown title?','FontSize', labelfontsize)   %**********************UNKNOWN TITLE****************************
end 


subplot(2,2,4);
if firstspikeindex > 0
plot(time(:,1),voltagetraces(:,firstspikeindex), 'color', [0.8 0.8 0.8], 'LineWidth',1)
axis([(spikelocsmatrix(1,firstspikeindex))/timefactor-4,(spikelocsmatrix(1,firstspikeindex))/timefactor+4, (min(voltagetraces(:,firstspikeindex)))-2, (max(voltagetraces(:,firstspikeindex)))+2])
xlabel('time')
ylabel('mV')
title('Unknown title?','FontSize', labelfontsize)   %**********************UNKNOWN TITLE****************************
end

%% Figure 3: Voltage vs current injected, per sweep
outfig3=figure('name',name,'position',[1 1 1474 744]);

for i_h=1:numberoftraces
   subplot(3,round(numberoftraces/3)+1,i_h); % Voltage Traces
   findpeaks(voltagetraces(:,i_h),'MinPeakHeight',spikedetectthreshold, 'MinPeakDistance', minspikedistance, 'MinPeakProminence', minspikeheight,'MaxPeakWidth', maxpeakwidth,   'Annotate','extents')
   currentlegend = gca; 
   legend(currentlegend,'off');
   hold on
   title(currentvalues(i_h), 'FontSize', labelfontsize);
   axis([currentons-(50*timefactor), currentoffs+(100*timefactor), minvoltage-50, maxvoltage+50])
   
end

%% Figure 4: Voltage Traces
volscale=.25;
outfig4=figure('name',name,'position',[1 1 1474 744]);

for i=1:numberoftraces
   if currentvalues(i) > -210 && currentvalues(i) < -190
      lineneg200= plot(time(:,1), voltagetraces(:,i), 'k', 'LineWidth',4*volscale);
   end  
   hold on
   if currentvalues(i) > 190 && currentvalues(i) < 210
      linepos200= plot(time(:,1), voltagetraces(:,i), 'color', [56,83,163] ./255, 'LineWidth',4*volscale);
   end  
   
   if currentvalues(i) > -10 && currentvalues(i) < 10
       plot(time(:,1), voltagetraces(:,i), 'color', [0.8 0.8 0.8], 'LineWidth',2*volscale);
   end
   
   if currentvalues(i) > 40 && currentvalues(i) < 60
       plot(time(:,1), voltagetraces(:,i), 'color', [0.8 0.8 0.8], 'LineWidth',2*volscale);
   end
   if currentvalues(i) > 90 && currentvalues(i) < 110
       plot(time(:,1), voltagetraces(:,i), 'color', [0.8 0.8 0.8], 'LineWidth',2*volscale);
   end
   if currentvalues(i) > 140 && currentvalues(i) < 160
       plot(time(:,1), voltagetraces(:,i), 'color', [0.8 0.8 0.8], 'LineWidth',2*volscale);
   end
   if currentvalues(i) > 240 && currentvalues(i) < 260
       plot(time(:,1), voltagetraces(:,i), 'color', [0.8 0.8 0.8], 'LineWidth',2*volscale);
   end
    if currentvalues(i) > 295 && currentvalues(i) < 305
        plot(time(:,1), voltagetraces(:,i), 'color', [0.8 0.8 0.8], 'LineWidth',2*volscale);
    end
   
  % Negative Current Values
   
   if currentvalues(i) > -60 && currentvalues(i) < -40
       plot(time(:,1), voltagetraces(:,i), 'color', [0.8 0.8 0.8], 'LineWidth',2*volscale);
   end  
   if currentvalues(i) > -110 && currentvalues(i) < -90
       plot(time(:,1), voltagetraces(:,i), 'color', [0.8 0.8 0.8], 'LineWidth',2*volscale);
   end 
   if currentvalues(i) > -160 && currentvalues(i) < -140
       plot(time(:,1), voltagetraces(:,i), 'color', [0.8 0.8 0.8], 'LineWidth',1*volscale);
   end    
   if currentvalues(i) > -260 && currentvalues(i) < -240
       plot(time(:,1), voltagetraces(:,i), 'color', [0.8 0.8 0.8], 'LineWidth',2*volscale);
   end
   
end
if hypcol>0
uistack(lineneg200, 'top')
end
if depcol>0
uistack(linepos200, 'top')
end
set(gca, 'FontSize', axisfontsize, 'FontName', figfont, 'FontWeight', figfontweight)
title('Voltage Traces','FontSize', titlefontsize, 'FontName', figfont, 'FontWeight', figfontweight)
xlabel('Time (ms)','FontSize', labelfontsize, 'FontName', figfont, 'FontWeight', figfontweight)
ylabel('Voltage (mV)','FontSize', labelfontsize, 'FontName', figfont, 'FontWeight', figfontweight)
axis([currentonms-25,currentoffms+275, -200, 50])

    
end  %function