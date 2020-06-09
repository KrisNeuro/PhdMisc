%% SA_ActiveResponsesOverTime.m
%
% Standalone script for pulling active nose-poke responses from MEDPC
% output files
%
% Data file format:
%   - C: Infusions             K: Time series (10-min intervals, 12 values
%   - D: Active responses      L: Time series (10-min intervals, 12 values)
%   - E: Inactive responses    M: Time series (10-min intervals, 12 values
% 
% KJS init: 2019-08-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FR5_d5,E_d1,E_d10,R] = SA_ActiveResponsesOverTime()

%% SETUP
%% Subject ID listing
    % 'D*' = 0.5 mg/kg KET
    % 'E*' = 0.25 mg/kg KET
    % 'F*' = 0.125 mg/kg KET
    % 'S*' = 0.0 mg/kg KET (saline)
subjs = {
    'D1-1' 'D1-2' 'D2-1' 'D2-2' 'D3-2' 'D4-1' 'D4-2' 'D5-1' 'D5-2' 'D6-1' 'D6-2' 'D8-1' 'D8-2' 'D9-1' 'D10-1' 'D10-2' 'D11-1' 'D11-2' ...
    'E1-1' 'E1-2' 'E3-1' 'E3-2' 'E4-1' 'E4-2' 'E5-1' 'E6-1' 'E7-1' 'E7-2' 'E8-1' 'E8-2' 'E9-1' 'E9-2' 'E12-1' 'E12-2' 'E10-1' 'E10-2' 'E11-1' 'E11-2' ...
    'F1-1' 'F1-2' 'F2-1' 'F2-2' 'F3-1' 'F3-2' 'F4-1' 'F4-2' 'F5-1' 'F5-2' 'F6-1' 'F6-2' 'F7-1' 'F7-2' ...
    'S1-1' 'S1-2' 'S2-1' 'S2-2' 'S3-1' 'S3-2' 'S4-1' 'S4-2' 'S5-1' 'S5-2' 'S6-1' 'S6-2' 'S7-1' 'S7-2'}; %subject ID listing
    % n(D)=19    n(E)=20   n(F)=14   n(S)=14

%% Treatment group listing   [KET](mg/kg/inf)
Tx = [0.5 0.25 0.125 0];

%% Listing of testing date files -- USER MUST CHANGE IF DIFFERENT FROM DPH!!
FR5_dates = {'!2019-08-06' '!2019-08-16'}; % FR5 (day 5)
E1_dates = {}; % Extinction (day 1)
E10_dates = {};% Extinction (day 10)
R_dates = {}; % Reinstatement (day 1)

%% Create data tables (3) FR5_d5, E_d1, E_d10, R  (FR=Fixed ratio; E=Extinction; R=Reinstatement; d=day)
    % subjID, Sex, Dose, Session Date, Infusions, ActiveRes, InactiveRes
    % string  cat  cat   string         cell       cell       cell
FR5_d5 = table('Size',[length(subjs) 7],...
    'VariableTypes',{'string' 'categorical' 'categorical' 'string' 'cell' 'cell' 'cell'},...
    'VariableNames',{'subjID', 'Sex', 'Dose', 'Date', 'Infusions', 'AR', 'IR'});

E_d1 = table('Size',[length(subjs) 7],...
    'VariableTypes',{'string' 'categorical' 'categorical' 'string' 'cell' 'cell' 'cell'},...
    'VariableNames',{'subjID', 'Sex', 'Dose', 'Date', 'Infusions', 'AR', 'IR'});

E_d10 = table('Size',[length(subjs) 7],...
    'VariableTypes',{'string' 'categorical' 'categorical' 'string' 'cell' 'cell' 'cell'},...
    'VariableNames',{'subjID', 'Sex', 'Dose', 'Date', 'Infusions', 'AR', 'IR'});

R = table('Size',[length(subjs) 7],...
    'VariableTypes',{'string' 'categorical' 'categorical' 'string' 'cell' 'cell' 'cell'},...
    'VariableNames',{'subjID', 'Sex', 'Dose', 'Date', 'Infusions', 'AR', 'IR'});
    
% Deal subject IDs into tables
FR5_d5.subjID(:) = subjs;
E_d1.subjID(:) = subjs;
E_d10.subjID(:) = subjs;
R.subjID(:) = subjs;

%% Assign dose based on subject ID
    % D = 0.5 mg/kg KET
    % E = 0.25 mg/kg KET
    % F = 0.125 mg/kg KET
    % S = SAL
for i = 1:length(subjs)
    switch subjs{i}(1)
        case 'D'
            FR5_d5.dose(i) = Tx(1);
            E_d1.dose(i) = Tx(1);
            E_d10.dose(i) = Tx(1);
            R.dose(i) = Tx(1);
        case 'E'
            FR5_d5.dose(i) = Tx(2);
            E_d1.dose(i) = Tx(2);
            E_d10.dose(i) = Tx(2);
            R.dose(i) = Tx(2);
        case 'F'
            FR5_d5.dose(i) = Tx(3);
            E_d1.dose(i) = Tx(3);
            E_d10.dose(i) = Tx(3);
            R.dose(i) = Tx(3);
        case 'S'
            FR5_d5.dose(i) = Tx(4);
            E_d1.dose(i) = Tx(4);
            E_d10.dose(i) = Tx(4);
            R.dose(i) = Tx(4);
    end
end
clear i

%% User select data file input directory
drIn = uigetdir(pwd,'Select data input directory'); 
if drIn(end) ~= filesep %if file separator is not at the end of drIn
    drIn = [drIn filesep];  %append file separator
end
cd(drIn)

%%
%% OK TO HERE %%
%%
%% Data-scan parameters


%from fxn
%     delimiter = ' ';
%     startRow = [4,166,174,182];
%     endRow = [4,168,176,184];
%     formatSpec = '%*s%s%s%s%s%s%[^\n\r]';

% Text lines for raw data textscan
    subjRow = [5 193 381 569 757 945 1133 1321 1509 1697 1885 2073 2261 2449]; %subject ID
    subjRow = subjRow - 1;
    [formatSpec.subj] = '%*s%s%*s%*s%*s%*s%[^\n\r]'; %chars 10-13 will have subjID in these rows
% % Infusions
    [startRow.In] = [166,354,542,730,918,1106,1294,1482,1670,1858,2046,2234,2422,2610];
    [endRow.In] = startRow.In + 2;
    [formatSpec.In] = '%*s%f%f%f%f%f%[^\n\r]';
% % Active reseponses
%     [startRow.AR] = [175 363 551 739 927 1115 1303 1491 1679 1867 2055 2243 2431 2619];
%     [endRow.AR] = startRow.AR+2;
%     [formatSpec.AR] = '%*39s%f%[^\n\r]';
% % Inactive responses
%     [startRow.IR] = [183 371 559 747 935 1123 1311 1499 1687 1875 2063 2251 2439 2627];
%     [endRow.IR] = startRow.IR+2;
%     [formatSpec.IR] = '';

%% Get list of MedAssociates data files to extract
cd(drIn)
[files] = dir('!*'); %list data files in drIn
files = {files(:).name};

%% Extract data from txt files
%% FR5_d5
fr5fileidx = strcmp(FR5_dates,files); %find indices of listed FR5 testing dates within full data file list
fr5files = files(fr5fileidx); %cell array of FR5 recording date files

for fi = 1:length(fr5files) %loop thru data files
    % Open txt file
    session = fr5files{fi}; %txt file name
%     fid = fopen(session,'r'); %open txt file
    
%     Search all subjRow rows in text to see which subject IDs are here
    subi = {'D1-2' 'E1-1' 'D5-1'}; %temporary replacement for sandbox!
%     textscan(fid, '%[^\n\r]', subjRow(1)-1, 'WhiteSpace', '', 'ReturnOnError', false);
%     dataArray = textscan(fid, formatSpec.subj, subjRow(1)-subjRow(1)+1, 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
%     for block=2:length(subjRow)
%     frewind(fid);
%     textscan(fid, '%[^\n\r]', subjRow(block)-1, 'WhiteSpace', '', 'ReturnOnError', false);
%     dataArrayBlock = textscan(fid, formatSpec.subj, subjRow(block)-subjRow(block)+1, 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
%     dataArray{block} = [dataArray{1};dataArrayBlock{1}];
%     end
    %     for block=2:length(subjRow)
%         frewind(fid);
%         textscan(fid, '%[^\n\r]', subjRow(block)-1, 'WhiteSpace', '', 'ReturnOnError', false);
%         dataArrayBlock = textscan(fid, formatSpec.subj, subjRow(block)-subjRow(block)+1, 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
%         dataArray{block} = [dataArray{1};dataArrayBlock{1}];
%     end
%     dataArray = textscan(fid, formatSpec.subj, subjRow, 'Delimiter',' ', 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
%     s = textscan(fid);
%     [subi] = textscan(fid, formatSpec.subj, subjRow(1),'Delimiter', '', 'HeaderLines', subjRow(1)-1,'ReturnOnError', false, 'EndOfLine', '\r','TextType','char');
%     [subi] = textscan(fid,formatSpec.subj,subjRow,'Delimiter', '', 'WhiteSpace','','HeaderLines', subjRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');


    % {subi} = {1 x 14} cell array of char arrays
    for si = 1:length(subi) %loop thru subjects on this recording
        fid = fopen(session,'r'); %open file
        % Find row index in output tables for this subject
        sidx = strcmp(subi{si},FR5_d5.subjID(:));
        sidx = find(sidx==1); %index in table for this subject's data

        % INFUSIONS TIME SERIES DATA
        textscan(fid, '%[^\n\r]', startRow.In(si)-1, 'WhiteSpace', '', 'ReturnOnError', false);
        dataArray = textscan(fid, formatSpec.In, endRow.In(si)-startRow.In(si)+1, 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
        X = [dataArray{1:end-1}]';
        fclose(fid); %close file
        X = reshape(X,1,15);
        X = X(1:12); %Anything after the 12th spot is all 0s (test ends at +120, 10-min bins)
        FR5_d5.Infusions(sidx) = deal(mat2cell(X,1,12)); 
        clear fid X dataA*
    end
    

    
%     data = textscan(fid,formatSpec,endRow-startRow+1,'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n','TextType','char');
%     [Trial.(pretxname).(txname)(subi).(sessname)] = data{1,1}(1,1); %total counts for session
%     [TimeCourse.(pretxname).(txname)(subi).(sessname)] = data{1,1}(2:end); %minute bin counts
%     end %loop thru subjects
    
    fclose(fid); %close txt file
    clear fid session
end %loop thru data files
    
    
    
    
    
    
% 
% for prei = 1:length(Expt.Pretx) %loop thru pretx conds
% 	pretxname = Expt.Pretx{prei};
% 	prehere = cd(pretxname);
% 	for txi = 1:length(Expt.Tx) %loop thru tx conds
% 		txname = Expt.Tx{txi};
% 		txhere = cd(txname);
% 		idlist = dir; %list all files
% 		idlist = idlist(3:end); %exclude . and .. files
% 		for subi = 1:length(idlist)	%loop thru each animal data
% 			subj = idlist(subi).name;
% 			[Trial.(pretxname).(txname)(subi).AnimalID] = deal(subj); %each row = animal ID
% 			[TimeCourse.(pretxname).(txname)(subi).AnimalID] = deal(subj);
% 			cd(subj) %open animal's folder
% 			for typei = 1:length(Expt.FeatNames) %loop thru Control/Treatment files
% 				if typei==1
% 					sesslabel = [5:8];
% 					else
% 					sesslabel = [5:7];
% 				end
% 				typehere = cd(Expt.FeatNames{typei}); %open appropriate folder
% 				files = dir('*.txt'); %list txt files
% 				if length(files)==0 %check that files are in this dir
% 					fprintf('%s\n','Error! No data in folder:',typehere)
% 				end
% 				for sessi = 1:length(files) %loop thru sessions in this type
% 					session = files(sessi).name; %txt file name
% 					sessname = session(sesslabel);
% 					fid = fopen(session,'r'); %open txt file
% 					data = textscan(fid,formatSpec,endRow-startRow+1,'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
% 					[Trial.(pretxname).(txname)(subi).(sessname)] = data{1,1}(1,1); %total counts for session
% 					[TimeCourse.(pretxname).(txname)(subi).(sessname)] = data{1,1}(2:end); %minute bin counts
% 					fclose(fid); %close txt file
% 				end
% 				cd ..
% 			end
% 			cd ..
% 		end
% 		cd ..
% 	end
% 	cd ..
% end			
cd(drIn);
cd outputs

end

