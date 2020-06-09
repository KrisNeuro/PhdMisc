function [BL,TenMin,ThreeHr,TwoFourHr,OneMin] = KETrexlist(files)
    %1. Baseline (BL)
    seek = 'BL';
    BLlist1 = {};
    for i = 3:length(files) %skip . and ..
        ii = i-2;
        id = files(i).name;
        check = strfind(id,seek);
            if isempty(check) %if an EPM file
                BLlist1{ii,1} = [];
            else %if a BL file
                BLlist1{ii,1} = deal(id);
            end
    end
    for i = 1:length(BLlist1)
        BLlist(i) = ~isempty(BLlist1{i});
    end
    BL = BLlist1(BLlist);
    BL = BL{:};
    
    %2. 10 mins (TenMin)
    seek = '10MIN';
    BLlist1 = {};
    for i = 3:length(files) %skip . and ..
        ii = i-2;
        id = files(i).name;
        check = strfind(id,seek);
            if isempty(check) %if an EPM file
                BLlist1{ii,1} = [];
            else %if a BL file
                BLlist1{ii,1} = deal(id);
            end
    end
    for i = 1:length(BLlist1)
        BLlist(i) = ~isempty(BLlist1{i});
    end
    TenMin = BLlist1(BLlist);
    TenMin = TenMin{:};
    
    %3. 3 hours (ThreeHr)
    seek = '3HR';
    BLlist1 = {};
    for i = 3:length(files) %skip . and ..
        ii = i-2;
        id = files(i).name;
        check = strfind(id,seek);
            if isempty(check) %if an EPM file
                BLlist1{ii,1} = [];
            else %if a BL file
                BLlist1{ii,1} = deal(id);
            end
    end
    for i = 1:length(BLlist1)
        BLlist(i) = ~isempty(BLlist1{i});
    end
    ThreeHr = BLlist1(BLlist);
    ThreeHr = ThreeHr{:};
    
    %4. 24 hours (TwoFourHr)
    seek = '24HR';
    BLlist1 = {};
    for i = 3:length(files) %skip . and ..
        ii = i-2;
        id = files(i).name;
        check = strfind(id,seek);
            if isempty(check) %if an EPM file
                BLlist1{ii,1} = [];
            else %if a BL file
                BLlist1{ii,1} = deal(id);
            end
    end
    for i = 1:length(BLlist1)
        BLlist(i) = ~isempty(BLlist1{i});
    end
    TwoFourHr = BLlist1(BLlist);
    TwoFourHr = TwoFourHr{:};
    
    %5. (ONLY FOR KET REX) 1 min (OneMin)
    if length(files)>6
        seek = '1MIN';
        BLlist1 = {};
        for i = 3:length(files) %skip . and ..
            ii = i-2;
            id = files(i).name;
            check = strfind(id,seek);
                if isempty(check) %if an EPM file
                    BLlist1{ii,1} = [];
                else %if a BL file
                    BLlist1{ii,1} = deal(id);
                end
        end
            for i = 1:length(BLlist1)
                BLlist(i) = ~isempty(BLlist1{i});
            end
        OneMin = BLlist1(BLlist);
        OneMin = OneMin{:};
    end
end