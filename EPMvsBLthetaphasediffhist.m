edges=-pi:0.2:pi; 
   
%Get animal data: EPM and same-day BL
cd('H:\Neuralynx\DATA-Backup\PreAnalysis\E108')
load('E108_EPMdat-BLdat.mat')
   
% setup for theta butter filter
Fs = 2000;
Fn = Fs/2;             % Nyquist Frequency
Wp = [4  12]/Fn;         % Theta Passband
Ws = [2.5 13.5]/(Fs/2);  %Theta StopBand
Rp = 3;                 %ripple pass
Rs = 40;                %ripple stop
[n,Wn] = buttord(Wp,Ws,Rp,Rs);
[z, p, k] = butter(n,Wn,'bandpass');
[sos,g] = zp2sos(z,p,k);
filt = dfilt.df2sos(sos,g);

%% Filter for theta
BLdat = filter(filt,BLdat);
EPMdat = filter(filt,EPMdat);
clear filt Fn g i k n p Rp Rs sos Wn Wp Ws z Fs

%% Get theta angles
for i = 1:size(BLdat,3) %loop thru recordings
[BL_IL_ang(:,i)] = angle(hilbert(BLdat(:,1,i)));
[BL_dHPC_ang(:,i)] = angle(hilbert(BLdat(:,2,i)));
[BL_vHPC_ang(:,i)] = angle(hilbert(BLdat(:,3,i)));
[BL_PL_ang(:,i)] = angle(hilbert(BLdat(:,4,i)));
[EPM_IL_ang(:,i)] = angle(hilbert(EPMdat(:,1,i)));
[EPM_dHPC_ang(:,i)] = angle(hilbert(EPMdat(:,2,i)));
[EPM_vHPC_ang(:,i)] = angle(hilbert(EPMdat(:,3,i)));
[EPM_PL_ang(:,i)] = angle(hilbert(EPMdat(:,4,i)));
end
clear BLdat EPMdat

%% Get theta lags for each pair
BLPhaseDiff.ILdHPC = BL_IL_ang - BL_dHPC_ang;
BLPhaseDiff.ILvHPC = BL_IL_ang - BL_vHPC_ang;
BLPhaseDiff.ILPL = BL_IL_ang - BL_PL_ang;
BLPhaseDiff.dHPCvHPC = BL_dHPC_ang - BL_vHPC_ang;
BLPhaseDiff.PLdHPC = BL_PL_ang - BL_dHPC_ang;
BLPhaseDiff.PLvHPC = BL_PL_ang - BL_vHPC_ang;

EPMPhaseDiff.ILdHPC = EPM_IL_ang - EPM_dHPC_ang;
EPMPhaseDiff.ILvHPC = EPM_IL_ang - EPM_vHPC_ang;
EPMPhaseDiff.ILPL = EPM_IL_ang - EPM_PL_ang;
EPMPhaseDiff.dHPCvHPC = BL_dHPC_ang - EPM_vHPC_ang;
EPMPhaseDiff.PLdHPC = EPM_PL_ang - EPM_dHPC_ang;
EPMPhaseDiff.PLvHPC = EPM_PL_ang - EPM_vHPC_ang;
EPMPhaseDiff.dHPCvHPC = EPM_dHPC_ang - EPM_vHPC_ang;

%% Limit lags from -pi:pi
% BL 
  for i = 1:size(BLPhaseDiff.ILvHPC,2) %loop through recordings
      disp(i)
   % IL-dHPC
    f=find(BLPhaseDiff.ILdHPC(:,i)>pi);%identifies phase differences larger than pi and subtracts pi from them.
    BLPhaseDiff.ILdHPC(f,i)=BLPhaseDiff.ILdHPC(f,i)-pi;
    g=find(BLPhaseDiff.ILdHPC(:,i)<-pi);%does the same thing, but for phase differences smaller than -pi
    BLPhaseDiff.ILdHPC(g,i)=BLPhaseDiff.ILdHPC(g,i)+pi;
    [BLPhaseDiffhist.ILdHPC(:,i)]=hist(BLPhaseDiff.ILdHPC(:,i),edges);%now il_dhip_phase_diff has only values in between -pi and +pi
    clear f g
  % IL-vHPC
    f=find(BLPhaseDiff.ILvHPC(:,i)>pi);
    BLPhaseDiff.ILvHPC(f,i)=BLPhaseDiff.ILvHPC(f,i)-pi;
    g=find(BLPhaseDiff.ILvHPC(:,i)<-pi);
    BLPhaseDiff.ILvHPC(g,i)=BLPhaseDiff.ILvHPC(g,i)+pi;
    [BLPhaseDiffhist.ILvHPC(:,i)]=hist(BLPhaseDiff.ILvHPC(:,i),edges);
    clear f g
   % IL-PL
    f=find(BLPhaseDiff.ILPL(:,i)>pi);
    BLPhaseDiff.ILPL(f,i)=BLPhaseDiff.ILPL(f,i)-pi;
    g=find(BLPhaseDiff.ILPL(:,i)<-pi);
    BLPhaseDiff.ILPL(g,i)=BLPhaseDiff.ILPL(g,i)+pi;
    [BLPhaseDiffhist.ILPL(:,i)]=hist(BLPhaseDiff.ILPL(:,i),edges);
    clear f g
   % dHPC-vHPC
    f=find(BLPhaseDiff.dHPCvHPC(:,i)>pi);
    BLPhaseDiff.dHPCvHPC(f,i)=BLPhaseDiff.dHPCvHPC(f,i)-pi;
    g=find(BLPhaseDiff.dHPCvHPC(:,i)<-pi);
    BLPhaseDiff.dHPCvHPC(g,i)=BLPhaseDiff.dHPCvHPC(g,i)+pi;
    [BLPhaseDiffhist.dHPCvHPC(:,i)]=hist(BLPhaseDiff.dHPCvHPC(:,i),edges);
    clear f g
   % PL-dHPC
    f=find(BLPhaseDiff.PLdHPC(:,i)>pi);
    BLPhaseDiff.PLdHPC(f,i)=BLPhaseDiff.PLdHPC(f,i)-pi;
    g=find(BLPhaseDiff.PLdHPC(:,i)<-pi);
    BLPhaseDiff.PLdHPC(g,i)=BLPhaseDiff.PLdHPC(g,i)+pi;
    [BLPhaseDiffhist.PLdHPC(:,i)]=hist(BLPhaseDiff.PLdHPC(:,i),edges);
    clear f g
   %PL-vHPC
    f=find(BLPhaseDiff.PLvHPC(:,i)>pi);
    BLPhaseDiff.PLvHPC(f,i)=BLPhaseDiff.PLvHPC(f,i)-pi;
    g=find(BLPhaseDiff.PLvHPC(:,i)<-pi);
    BLPhaseDiff.PLvHPC(g,i)=BLPhaseDiff.PLvHPC(g,i)+pi;
    [BLPhaseDiffhist.PLvHPC(:,i)]=hist(BLPhaseDiff.PLvHPC(:,i),edges);
    clear f g
  end
  
% EPM 
  for i = 1:size(EPMPhaseDiff.ILvHPC,2) %loop through recordings
      disp(i)
   % IL-dHPC
    f=find(EPMPhaseDiff.ILdHPC(:,i)>pi);%identifies phase differences larger than pi and subtracts pi from them.
    EPMPhaseDiff.ILdHPC(f,i)=EPMPhaseDiff.ILdHPC(f,i)-pi;
    g=find(EPMPhaseDiff.ILdHPC(:,i)<-pi);%does the same thing, but for phase differences smaller than -pi
    EPMPhaseDiff.ILdHPC(g,i)=EPMPhaseDiff.ILdHPC(g,i)+pi;
    [EPMPhaseDiffhist.ILdHPC(:,i)]=hist(EPMPhaseDiff.ILdHPC(:,i),edges);%now il_dhip_phase_diff has only values in between -pi and +pi
    clear f g
  % IL-vHPC
    f=find(EPMPhaseDiff.ILvHPC(:,i)>pi);
    EPMPhaseDiff.ILvHPC(f,i)=EPMPhaseDiff.ILvHPC(f,i)-pi;
    g=find(EPMPhaseDiff.ILvHPC(:,i)<-pi);
    EPMPhaseDiff.ILvHPC(g,i)=EPMPhaseDiff.ILvHPC(g,i)+pi;
    [EPMPhaseDiffhist.ILvHPC(:,i)]=hist(EPMPhaseDiff.ILvHPC(:,i),edges);
    clear f g
   % IL-PL
    f=find(EPMPhaseDiff.ILPL(:,i)>pi);
    EPMPhaseDiff.ILPL(f,i)=EPMPhaseDiff.ILPL(f,i)-pi;
    g=find(EPMPhaseDiff.ILPL(:,i)<-pi);
    EPMPhaseDiff.ILPL(g,i)=EPMPhaseDiff.ILPL(g,i)+pi;
    [EPMPhaseDiffhist.ILPL(:,i)]=hist(EPMPhaseDiff.ILPL(:,i),edges);
    clear f g
   % dHPC-vHPC
    f=find(EPMPhaseDiff.dHPCvHPC(:,i)>pi);
    EPMPhaseDiff.dHPCvHPC(f,i)=EPMPhaseDiff.dHPCvHPC(f,i)-pi;
    g=find(EPMPhaseDiff.dHPCvHPC(:,i)<-pi);
    EPMPhaseDiff.dHPCvHPC(g,i)=EPMPhaseDiff.dHPCvHPC(g,i)+pi;
    [EPMPhaseDiffhist.dHPCvHPC(:,i)]=hist(EPMPhaseDiff.dHPCvHPC(:,i),edges);
    clear f g
   % PL-dHPC
    f=find(EPMPhaseDiff.PLdHPC(:,i)>pi);
    EPMPhaseDiff.PLdHPC(f,i)=EPMPhaseDiff.PLdHPC(f,i)-pi;
    g=find(EPMPhaseDiff.PLdHPC(:,i)<-pi);
    EPMPhaseDiff.PLdHPC(g,i)=EPMPhaseDiff.PLdHPC(g,i)+pi;
    [EPMPhaseDiffhist.PLdHPC(:,i)]=hist(EPMPhaseDiff.PLdHPC(:,i),edges);
    clear f g
   %PL-vHPC
    f=find(EPMPhaseDiff.PLvHPC(:,i)>pi);
    EPMPhaseDiff.PLvHPC(f,i)=EPMPhaseDiff.PLvHPC(f,i)-pi;
    g=find(EPMPhaseDiff.PLvHPC(:,i)<-pi);
    EPMPhaseDiff.PLvHPC(g,i)=EPMPhaseDiff.PLvHPC(g,i)+pi;
    [EPMPhaseDiffhist.PLvHPC(:,i)]=hist(EPMPhaseDiff.PLvHPC(:,i),edges);
    clear f g
  end
  
  
  %% Half of maximum peak parameters
  % BL
for i = 1:size(BLPhaseDiffhist.ILdHPC,2) %loop through recordings
      disp(i)
   %IL-dHPC
    BL_halph_max_ILdHPC(i)=(max(BLPhaseDiffhist.ILdHPC(:,i)))/2;%identifies half of the peak height of the histogram
    BL_width_at_halph_max_ILdHPC(i)=length(find(BLPhaseDiffhist.ILdHPC(:,i)>BL_halph_max_ILdHPC(i)));%finds indices of the bins which have height higher than the peak_height/2, counts the number of bins
    BL_width_at_halph_max_ILdHPC(i)=BL_width_at_halph_max_ILdHPC(i)*mean(diff(edges));
   %IL-vHPC
    BL_halph_max_ILvHPC(i)=(max(BLPhaseDiffhist.ILvHPC(:,i)))/2;
    BL_width_at_halph_max_ILvHPC(i)=length(find(BLPhaseDiffhist.ILvHPC(:,i)>BL_halph_max_ILvHPC(i)));
    BL_width_at_halph_max_ILvHPC(i)=BL_width_at_halph_max_ILvHPC(i)*mean(diff(edges));
   % IL-PL
    BL_halph_max_ILPL(i)=(max(BLPhaseDiffhist.ILPL(:,i)))/2;
    BL_width_at_halph_max_ILPL(i)=length(find(BLPhaseDiffhist.ILPL(:,i)>BL_halph_max_ILPL(i)));
    BL_width_at_halph_max_ILPL(i)=BL_width_at_halph_max_ILPL(i)*mean(diff(edges));
   % dHPC-vHPC
    BL_halph_max_dHPCvHPC(i)=(max(BLPhaseDiffhist.dHPCvHPC(:,i)))/2;
    BL_width_at_halph_max_dHPCvHPC(i)=length(find(BLPhaseDiffhist.dHPCvHPC(:,i)>BL_halph_max_dHPCvHPC(i)));
    BL_width_at_halph_max_dHPCvHPC(i)=BL_width_at_halph_max_dHPCvHPC(i)*mean(diff(edges));
   % PL-dHPC
    BL_halph_max_PLdHPC(i)=(max(BLPhaseDiffhist.PLdHPC(:,i)))/2;
    BL_width_at_halph_max_PLdHPC(i)=length(find(BLPhaseDiffhist.PLdHPC(:,i)>BL_halph_max_PLdHPC(i)));
    BL_width_at_halph_max_PLdHPC(i)=BL_width_at_halph_max_PLdHPC(i)*mean(diff(edges));
   % PL-vHPC
    BL_halph_max_PLvHPC(i)=(max(BLPhaseDiffhist.PLvHPC(:,i)))/2;
    BL_width_at_halph_max_PLvHPC(i)=length(find(BLPhaseDiffhist.PLvHPC(:,i)>BL_halph_max_PLvHPC(i)));
    BL_width_at_halph_max_PLvHPC(i)=BL_width_at_halph_max_PLvHPC(i)*mean(diff(edges));
end
   
 % EPM
for i = 1:size(EPMPhaseDiffhist.ILdHPC,2) %loop through recordings
      disp(i)
   %IL-dHPC
    EPM_halph_max_ILdHPC(i)=(max(EPMPhaseDiffhist.ILdHPC(:,i)))/2;%identifies half of the peak height of the histogram
    EPM_width_at_halph_max_ILdHPC(i)=length(find(EPMPhaseDiffhist.ILdHPC(:,i)>EPM_halph_max_ILdHPC(i)));%finds indices of the bins which have height higher than the peak_height/2, counts the number of bins
    EPM_width_at_halph_max_ILdHPC(i)=EPM_width_at_halph_max_ILdHPC(i)*mean(diff(edges));
   %IL-vHPC
    EPM_halph_max_ILvHPC(i)=(max(EPMPhaseDiffhist.ILvHPC(:,i)))/2;
    EPM_width_at_halph_max_ILvHPC(i)=length(find(EPMPhaseDiffhist.ILvHPC(:,i)>EPM_halph_max_ILvHPC(i)));
    EPM_width_at_halph_max_ILvHPC(i)=EPM_width_at_halph_max_ILvHPC(i)*mean(diff(edges));
   % IL-PL
    EPM_halph_max_ILPL(i)=(max(EPMPhaseDiffhist.ILPL(:,i)))/2;
    EPM_width_at_halph_max_ILPL(i)=length(find(EPMPhaseDiffhist.ILPL(:,i)>EPM_halph_max_ILPL(i)));
    EPM_width_at_halph_max_ILPL(i)=EPM_width_at_halph_max_ILPL(i)*mean(diff(edges));
   % dHPC-vHPC
    EPM_halph_max_dHPCvHPC(i)=(max(EPMPhaseDiffhist.dHPCvHPC(:,i)))/2;
    EPM_width_at_halph_max_dHPCvHPC(i)=length(find(EPMPhaseDiffhist.dHPCvHPC(:,i)>EPM_halph_max_dHPCvHPC(i)));
    EPM_width_at_halph_max_dHPCvHPC(i)=EPM_width_at_halph_max_dHPCvHPC(i)*mean(diff(edges));
   % PL-dHPC
    EPM_halph_max_PLdHPC(i)=(max(EPMPhaseDiffhist.PLdHPC(:,i)))/2;
    EPM_width_at_halph_max_PLdHPC(i)=length(find(EPMPhaseDiffhist.PLdHPC(:,i)>EPM_halph_max_PLdHPC(i)));
    EPM_width_at_halph_max_PLdHPC(i)=EPM_width_at_halph_max_PLdHPC(i)*mean(diff(edges));
   % PL-vHPC
    EPM_halph_max_PLvHPC(i)=(max(EPMPhaseDiffhist.PLvHPC(:,i)))/2;
    EPM_width_at_halph_max_PLvHPC(i)=length(find(EPMPhaseDiffhist.PLvHPC(:,i)>EPM_halph_max_PLvHPC(i)));
    EPM_width_at_halph_max_PLvHPC(i)=EPM_width_at_halph_max_PLvHPC(i)*mean(diff(edges));
end

clear i ans

%% Save data
save('A202_BLvsEPM_thetaphasediffhist.mat')
disp('saved')


% A201 done: 2/27/2019
% A202 done: 2/27/2019