%%GetLeonCases
function [Leondates,Leoncases,NewLeoncases,t,H1,H2] = GetLeonCases()
% [Leondates,Leoncases,NewLeoncases,t,H1,H2] = GetLeonCases();
% 
% Data source: https://github.com/nytimes/covid-19-data
% 
% KJS init: 2020-06-14 (as function)
% 

%% Set local data directory
if license == "731138"
    drIn = 'C:\Users\kristin.schoepfer\Desktop\rona\rona\covid-19-data\'; %local data input directory
else
    drIn = [uigetdir(pwd,'Select covid-19-data directory') filesep];
end

%% Import NY Times data
disp('Importing NY Times covid-19 data by county...')
filename = [drIn 'us-counties.csv']; %input file
delimiter = ',';
formatSpec = '%q%q%q%*q%q%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

% Converts text in the input cell array to numbers. Replaced non-numeric
% text with NaN.
rawData = dataArray{4};
for row=1:size(rawData, 1)
    % Create a regular expression to detect and remove non-numeric prefixes and suffixes.
    regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
    try
        result = regexp(rawData(row), regexstr, 'names');
        numbers = result.numbers;
        
        % Detected commas in non-thousand locations.
        invalidThousandsSeparator = false;
        if numbers.contains(',')
            thousandsRegExp = '^[-/+]*\d+?(\,\d{3})*\.{0,1}\d*$';
            if isempty(regexp(numbers, thousandsRegExp, 'once'))
                numbers = NaN;
                invalidThousandsSeparator = true;
            end
        end
        % Convert numeric text to numbers.
        if ~invalidThousandsSeparator
            numbers = textscan(char(strrep(numbers, ',', '')), '%f');
            numericData(row, 4) = numbers{1};
            raw{row, 4} = numbers{1};
        end
    catch
        raw{row, 4} = rawData{row};
    end
end

% Convert the contents of columns with dates to MATLAB datetimes using the specified date format.
try
    dates{1} = datetime(dataArray{1}, 'Format', 'yyyy-MM-dd', 'InputFormat', 'yyyy-MM-dd');
catch
    try
        % Handle dates surrounded by quotes
        dataArray{1} = cellfun(@(x) x(2:end-1), dataArray{1}, 'UniformOutput', false);
        dates{1} = datetime(dataArray{1}, 'Format', 'yyyy-MM-dd', 'InputFormat', 'yyyy-MM-dd');
    catch
        dates{1} = repmat(datetime([NaN NaN NaN]), size(dataArray{1}));
    end
end
dates = dates(:,1);
rawNumericColumns = raw(:, 4);
rawStringColumns = string(raw(:, [2,3]));

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Make sure any text containing <undefined> is properly converted to an <undefined> categorical
for catIdx = [1,2]
    idx = (rawStringColumns(:, catIdx) == "<undefined>");
    rawStringColumns(idx, catIdx) = "";
end

%% Create output variable
uscounties = table;
uscounties.date = dates{:, 1};
uscounties.county = categorical(rawStringColumns(:, 1));
uscounties.state = categorical(rawStringColumns(:, 2));
uscounties.cases = cell2mat(rawNumericColumns(:, 1));
uscounties(1,:) = [];

%% Index Leon County, FL data
FLidx = strcmp(string(uscounties.state),'Florida');
Leonidx = strcmp(string(uscounties.county),'Leon') ;
LeonFLidx = Leonidx & FLidx;

Leoncases = uscounties.cases(LeonFLidx); %cumulative cases
NewLeoncases = diff([0; Leoncases]); %new cases per day
Leondates = uscounties.date(LeonFLidx); %dates
t = [1:length(Leondates)]'; %#ok<NBRAK> %time (day) vector

%% Calculate moving averages
M7d = movmean(NewLeoncases,[6 0]); %7-day average
M14d = movmean(NewLeoncases, [13 0]); %14-day average

%% Plot Leon FL cases
H1 = figure('position',[1036 426 724 570]);
subplot(211); %New cases/day
    plot(Leondates,NewLeoncases,'k','linew',2)
    box off
    ylabel('# new positive tests/day')
    title('Leon county, FL: COVID-19 new cases')
    axis tight
subplot(212); %Cumulative cases
    plot(Leondates,Leoncases,'b','linew',2)
    box off
    ylabel('cumulative positive tests')
    title('Leon county, FL: Cumulative COVID-19 + cases')
    axis tight

%% Plot moving averages for new Leon cases/day
H2 = figure('position',[1036 41 724 310]);
    plot(Leondates, NewLeoncases,'k','linew',2)
    hold on
    plot(Leondates, M7d,'r','linew',1.5)
    plot(Leondates, M14d,'b','linew',1.5)
    box off
    axis tight
    legend RawVals 7dayAvg 14dayAvg
    legend('location','northwest') 
    title('New Leon cases/day')

end
