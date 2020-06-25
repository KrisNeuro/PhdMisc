%%GetCountyCases
function [uscounties,Dates,CountyCases,NewCountyCases,t,H1,H2] = GetCountyCases()
% [Dates,CountyCases,NewCountyCases,t,H1,H2] = GetCountyCases()
% 
% Imports NY Times COVID-19 data, filters daily cases by user-input state
% and county.
% Plots new cases/day, cumulative + cases, and new cases/day with 7-day and
%   14-day running averages overlaid.
% 
% Data source: https://github.com/nytimes/covid-19-data
% 
% KJS init: 2020-06-14 (as function)
% 

%% Set local data directory
drIn = [uigetdir(pwd,'Select covid-19-data directory (NYTimes master)') filesep];

redoflag = 1; %initialize flag for later

%% Import NY Times data (.csv)
disp('Importing NY Times COVID-19 data by county...')
filename = [drIn 'us-counties.csv']; %input file
delimiter = ',';
formatSpec = '%q%q%q%*q%q%[^\n\r]';
fileID = fopen(filename,'r');
    disp('[.     ]')
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
    disp('[..    ]')
fclose(fileID);
clear filename

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

disp('[...   ]')

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
disp('[....  ]')

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
disp('[..... ]')

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Make sure any text containing <undefined> is properly converted to an <undefined> categorical
for catIdx = [1,2]
    idx = (rawStringColumns(:, catIdx) == "<undefined>");
    rawStringColumns(idx, catIdx) = "";
end
disp('[......]')

%% Create data table
uscounties = table;
uscounties.date = dates{:, 1};
uscounties.county = categorical(rawStringColumns(:, 1));
uscounties.state = categorical(rawStringColumns(:, 2));
uscounties.cases = cell2mat(rawNumericColumns(:, 1));
uscounties(1,:) = [];
disp('Data loaded!')

% Clean up workspace
clearvars ans catIdx col dat* delimiter fileID formatSpec idx invalid* num* R raw* regexstr result row

%% User selects state & county; Index data within big table
states = unique(uscounties.state); %unique state names
% counties = unique(uscounties.county); %unique county names

while ~isempty(redoflag)
    
    % State
    [stindx] = listdlg('PromptString',{'Select a state.'},'SelectionMode','single','ListString',states);
    stname = char(states(stindx)); %state name
    clear stindx
    Stateidx = strcmp(string(uscounties.state),stname); %state index in big data table

    % County
    okcounties = unique(uscounties.county(Stateidx)); %counties within the selected state
    [ctyindx] = listdlg('PromptString',{'Select a county.'},'SelectionMode','single','ListString',okcounties);
    cname = char(okcounties(ctyindx)); %county name
    clear ctyindx
    Countyidx = strcmp(string(uscounties.county),cname); %county index in big table - **some county names are repeated b/t states!**

    CountyStateidx = Countyidx & Stateidx; %specific county within state
    fprintf('Getting COVID-19 new case data for %s County, %s...\n',cname, stname)

    CountyCases = uscounties.cases(CountyStateidx); %cumulative cases
    CountyCases =[0; CountyCases];
    NewCountyCases = diff([0; CountyCases]); %new cases per day
    Dates = uscounties.date(CountyStateidx); %dates
    Dates = [Dates(1)-1; Dates];
    t = [1:length(Dates)]'; %#ok<NBRAK> %time (day) vector

    %% Calculate moving averages
    M7d = movmean(NewCountyCases,[6 0]); %7-day average
    M14d = movmean(NewCountyCases, [13 0]); %14-day average

    %% Plot cases
    H1 = figure('position',[1036 426 724 570]);
    subplot(211); %New cases/day
        plot(Dates,NewCountyCases,'k','linew',2)
        box off
        ylabel('# new positive tests/day')
        title(sprintf('%s County, %s: COVID-19 new cases',cname,stname));
        axis tight
    subplot(212); %Cumulative cases
        plot(Dates,CountyCases,'b','linew',2)
        box off
        ylabel('cumulative positive tests')
        title(sprintf('%s County, %s: Cumulative COVID-19+ cases',cname,stname));
        axis tight

    %% Plot moving averages for new cases/day
    H2 = figure('position',[1036 41 724 310]);
        plot(Dates, NewCountyCases,'k','linew',2)
        hold on
        plot(Dates, M7d,'r','linew',1.5)
        plot(Dates, M14d,'b','linew',1.5)
        box off
        axis tight
        legend RawVals 7dayAvg 14dayAvg
        legend('location','northwest') 
        title('New Leon cases/day')
        title(sprintf('New %s County cases/day',cname));
    
        
    %% Option to re-do and view a new county
    goagain = [];
    while isempty(goagain)
        goagain = string(input('View another county? \n 0=No 1=Yes\n','s'));
        if strcmp(goagain,'0') %end
            disp('Plotting complete.')
            redoflag = []; %#ok<NASGU>
            return
        elseif ~strcmp(goagain,'1') && ~strcmp(goagain,'0') %typo
            disp('Typo. Try again.')
            goagain = [];
        else %repeat
            close(H1); close(H2);
            clear H1 H2 Dates cname County* M* NewCounty* ok* Stateidx stname t
        end
    end

end %while redoflag

end %function
