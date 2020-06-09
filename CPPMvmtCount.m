%% CPPMvmtCount.m
% CPP conditioning movement count extraction: MedAssociates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CPPMvmtCount.m extracts movement counts (beam breaks) from the
% MedAssociates CPP chambers for all conditioning sessions in an experiment.
%
% Use:
%   [Trial,TimeCourse] = CPPMvmtCount();
%
% INPUT:
% - .txt data files trimmed to contain only usable info, appropriately
% organized into subfolders within the root directory
% 
% OUTPUTS:
% - Trial       Total movement counts per animal, per CPP conditioning
% session (labeled as Con/control and Tx/drug-treated)
%       Trial outputs are automatically saved into labeled .xlsx
%       spreadsheets in a folder called 'outputs' within the root directory
% - TimeCourse  Time series data of movement counts per animal, per
% conditioning session, labeled similarly as in 'Trial'
%       Due to their multidimentional nature, TimeCourse data are not
%       auto-saved into .xlsx spreadsheets, but are relatively
%       user-friendly to navigate otherwise.
%
%
% Written by KJS, March 5, 2018
% please list any future edits and their dates below:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Trial,TimeCourse] = CPPMvmtCount();
%% SETUP
drIn = uigetdir;
%user selects the root directory containing appropriately formatted
%subfolders with trimmed CPP data
	if drIn(end) ~= filesep %if file separator is not at the end of drIn
        drIn = [drIn filesep];  %append file separator
    end
cd(drIn)
addpath(genpath(drIn)) %add all subfolders containing data

% EXPT struct: Experimental conditions
Expt.Pretx = {'sal' 'nlx'}; %expt. pretreatment drug options
Expt.Tx = {'sal' 'ket'}; % expt. treatment drug options
Expt.FeatNames = {'Control ','Treatment'};	%names of relevant features to be analyzed

%initialize holding places for final data
[Trial] = [];
[TimeCourse] = [];

% Variables for .txt textscan
startRow = 8;
endRow = 53;
formatSpec = '%*39s%f%[^\n\r]';

%% Extract data from txt files
for prei = 1:length(Expt.Pretx) %loop thru pretx conds
	pretxname = Expt.Pretx{prei};
	prehere = cd(pretxname);
	for txi = 1:length(Expt.Tx) %loop thru tx conds
		txname = Expt.Tx{txi};
		txhere = cd(txname);
		idlist = dir; %list all files
		idlist = idlist(3:end); %exclude . and .. files
		for subi = 1:length(idlist)	%loop thru each animal data
			subj = idlist(subi).name;
			[Trial.(pretxname).(txname)(subi).AnimalID] = deal(subj); %each row = animal ID
			[TimeCourse.(pretxname).(txname)(subi).AnimalID] = deal(subj);
			cd(subj) %open animal's folder
			for typei = 1:length(Expt.FeatNames) %loop thru Control/Treatment files
				if typei==1
					sesslabel = [5:8];
					else
					sesslabel = [5:7];
				end
				typehere = cd(Expt.FeatNames{typei}); %open appropriate folder
				files = dir('*.txt'); %list txt files
				if length(files)==0 %check that files are in this dir
					fprintf('%s\n','Error! No data in folder:',typehere)
				end
				for sessi = 1:length(files) %loop thru sessions in this type
					session = files(sessi).name; %txt file name
					sessname = session(sesslabel);
					fid = fopen(session,'r'); %open txt file
					data = textscan(fid,formatSpec,endRow-startRow+1,'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
					[Trial.(pretxname).(txname)(subi).(sessname)] = data{1,1}(1,1); %total counts for session
					[TimeCourse.(pretxname).(txname)(subi).(sessname)] = data{1,1}(2:end); %minute bin counts
					fclose(fid); %close txt file
				end
				cd ..
			end
			cd ..
		end
		cd ..
	end
	cd ..
end			
cd(drIn);
cd outputs

%% TOTAL COUNT OUTPUTS
% Reformat structs to cell arrays
salsal = struct2cell(Trial.sal.sal);
salket = struct2cell(Trial.sal.ket);
nlxsal = struct2cell(Trial.nlx.sal);
nlxket = struct2cell(Trial.nlx.ket);

% Save cell arrays as .xls documents
%SAL/SAL
xlfnamess = 'DevDat_SALSAL.xlsx'; %name of Excel sheet to be saved
xlswrite(xlfnamess,salsal(:,:));
fprintf('%s\n','SAL/SAL spreadsheet saved as:',xlfnamess)

%SAL/KET
xlfnamesk = 'DevDat_SALKET.xlsx'; %name of Excel sheet to be saved
xlswrite(xlfnamesk,salket(:,:));
fprintf('%s\n','SAL/KET spreadsheet saved as:',xlfnamesk)

%NLX/SAL
xlfnamens = 'DevDat_NLXSAL.xlsx'; %name of Excel sheet to be saved
xlswrite(xlfnamens,nlxsal(:,:));
fprintf('%s\n','NLX/SAL spreadsheet saved as:',xlfnamens);

%NLX/KET
xlfnamenk = 'DevDat_NLXKET.xlsx'; %name of Excel sheet to be saved
xlswrite(xlfnamenk,nlxket(:,:));
fprintf('%s\n','NLX/KET spreadsheet saved as:',xlfnamenk)

%% Save .mat data to \outputs
savefile = input('Name of mat file to be saved?: ','s');
savefile = [savefile '.mat']; %appdend .mat

save(savefile,'Trial','TimeCourse')
fprintf('%s\n','MATLAB data saved as: ',savefile)

%% WORKING SCRIPT ENDS HERE
%% 
%%%%%% SANDBOX START %%%%%%    
% Reformat structs to cell arrays
% salsaltime = struct2cell(TimeCourse.sal.sal);
% for pre = 1:length(
% for i = 1:length(TimeCourse.sal.sal) %loop thru subjects
% 	for ii = 1:size(TimeCourse.sal.sal
% salkettime = struct2cell(TimeCourse.sal.ket);
% nlxsaltime = struct2cell(TimeCourse.nlx.sal);
% nlxkettime = struct2cell(TimeCourse.nlx.ket);
% 	end
% end
% 
% % Save cell arrays as .xls documents
% %SAL/SAL
% xlfnamess = 'DevDat_SALSALtime.xlsx'; %name of Excel sheet to be saved
% xlswrite(xlfnamess,salsaltime{:,:});
% fprintf('%s\n','SAL/SAL timecourse spreadsheet saved as:',xlfnamess)
% 
% %SAL/KET
% xlfnamesk = 'DevDat_SALKETtime.xlsx'; %name of Excel sheet to be saved
% xlswrite(xlfnamesk,salket(:,:));
% fprintf('%s\n','SAL/KET timecourse spreadsheet saved as:',xlfnamesk)
% 
% %NLX/SAL
% xlfnamens = 'DevDat_NLXSALtime.xlsx'; %name of Excel sheet to be saved
% xlswrite(xlfnamens,nlxsal(:,:));
% fprintf('%s\n','NLX/SAL timecourse spreadsheet saved as:',xlfnamens);
% 
% %NLX/KET
% xlfnamenk = 'DevDat_NLXKETtime.xlsx'; %name of Excel sheet to be saved
% xlswrite(xlfnamenk,nlxket(:,:));
% fprintf('%s\n','NLX/KET timecourse spreadsheet saved as:',xlfnamenk)
% 
% %%%%%% SANDBOX END %%%%%%

end