%% SA_ActiveResponsesOverTime_V2
%
% Standalone script for pulling active nose-poke responses from MEDPC
% output files
%
% *** Bottom-up approach ***
% Read in all data files (loop) and input the values into cells, then make a table
% out of the data using cell2table
% Match subject IDs and dates post-hoc (?)
%
%
% Data file format:
%   - C: Infusions             K: Time series (10-min intervals, 12 values
%   - D: Active responses      L: Time series (10-min intervals, 12 values)
%   - E: Inactive responses    M: Time series (10-min intervals, 12 values
% 
% KJS init: 2019-09-11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [FR5_d5,E_d1,E_d10,R] = SA_ActiveResponsesOverTime_V2()
% Matches SA_ActiveResponsesOverTime from here until noted.
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

%% Treatment group listing   
    % [KET] mg/kg/inf
Tx = [0.5 0.25 0.125 0];

%% Listing of testing date files
FR5_dates = {'!2019-08-06' '!2019-08-07' '!2019-08-16'}; % FR5 (day 5)
E1_dates = {'!2019-08-07'}; % Extinction (day 1)
E10_dates = {};% Extinction (day 10)
R_dates = {}; % Reinstatement (day 1)

%% User select data file input directory
drIn = uigetdir(pwd,'Select data input directory'); 
if drIn(end) ~= filesep %if file separator is not at the end of drIn
    drIn = [drIn filesep];  %append file separator
end

%% Get list of MedAssociates data files to extract
cd(drIn)
[files] = dir('!*'); %list data files in drIn
files = {files(:).name};

    % Identify files for each test type
    FR5_files = files(find(contains(files,FR5_dates)));
    E1_files = files(find(contains(files,E1_dates)));
    E10_dates = files(find(contains(files,E10_dates)));
    R_files = files(find(contains(files,R_dates)));


%% Create data tables (3) FR5_d5, E_d1, E_d10, R  
    % FR=Fixed ratio
    % E=Extinction
    % R=Reinstatement   d=day
    
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


%% Extract data from FR5_d5 data files
for fi = 1:length(FR5_files) %loop thru data files
    fid = fopen(FR5_files{fi},'r'); %open txt file
    
    % Find row index in output tables for this subject
%         sidx = strcmp(subi{si},FR5_d5.subjID(:));
%         sidx = find(sidx==1); %index in table for this subject's data
% 
%         % INFUSIONS TIME SERIES DATA
%         textscan(fid, '%[^\n\r]', startRow.In(si)-1, 'WhiteSpace', '', 'ReturnOnError', false);
%         dataArray = textscan(fid, formatSpec.In, endRow.In(si)-startRow.In(si)+1, 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
%         X = [dataArray{1:end-1}]';

    fclose(fid); %close txt file
    
%         X = reshape(X,1,15);
%         X = X(1:12); %Anything after the 12th spot is all 0s (test ends at +120, 10-min bins)
%         FR5_d5.Infusions(sidx) = deal(mat2cell(X,1,12)); 
%         clear fid X dataA*
    end



end

