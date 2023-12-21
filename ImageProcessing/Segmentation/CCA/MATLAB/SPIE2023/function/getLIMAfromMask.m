function [LI_AUTO,MA_AUTO]=getLIMAfromMask(out)
    bwper = bwperim(out);
    [rows,cols,~] = find(bwper==1);
    
    %% there are two points for each column, except for first and last
    %% we need to avoid them and make two distinct profiles

    % take first element
    [cols_first,uni_1] = unique(cols);
    rows_first = rows(uni_1);
    
    % avoid first element
    avoid = true(length(cols),1); avoid(uni_1)=false;
    cols_not = cols(avoid);
    rows_not = rows(avoid);
    
    % take last element
    [cols_second,uni_2] = unique(cols_not,'last');
    rows_second = rows_not(uni_2);    
    
    % merge
    cols = [cols_first;cols_second];
    rows = [rows_first;rows_second];
    
    % sort
    [cols,I] = sort(cols);
    rows = rows(I,:);
    
    
    ma_rows = rows(1:2:end-1); ma_cols = cols(1:2:end-1);
    li_rows = rows(2:2:end);  li_cols = cols(2:2:end);

%         imshow(bwper);hold on
%         plot(li_cols(:),li_rows(:),'c','LineWidth',1.5),hold on,
%         plot(ma_cols(:),ma_rows(:),'y','LineWidth',1.5),hold on,
%         close all

    %% we cut the extrema of the profiles to avoid border effect and make a 
    %% profile predition more similar to the human one
    
    cut_ratio = 33;
    LI_AUTO = [li_cols,li_rows]; MA_AUTO = [ma_cols,ma_rows]; 
    
    cutLI = round(length(LI_AUTO)/cut_ratio); if cutLI == 0; cutLI =1; end
    cutMA = round(length(MA_AUTO)/cut_ratio); if cutMA == 0; cutMA =1; end
    
    LI_AUTO = LI_AUTO(cutLI:length(LI_AUTO)-cutLI,:);
    MA_AUTO = MA_AUTO(cutMA:length(MA_AUTO)-cutMA,:);

    if sum(LI_AUTO(:,2) > MA_AUTO(:,2)) > length(LI_AUTO)/2
        temp = MA_AUTO;
        MA_AUTO = LI_AUTO;
        LI_AUTO = temp;
    end
    
    %% correct for bwperim selection
    LI_AUTO(:,2) = LI_AUTO(:,2)-1;
    MA_AUTO(:,2) = MA_AUTO(:,2)+1;
    
    %% smooth profiles
    LI_AUTO(:,2) = movmean(LI_AUTO(:,2),length(LI_AUTO)/16);
    MA_AUTO(:,2) = movmean(MA_AUTO(:,2),length(MA_AUTO)/16);

%     imshow(out);hold on
%     plot(LI_AUTO(:,1),LI_AUTO(:,2),'c','LineWidth',1),hold on,
%     plot(MA_AUTO(:,1),MA_AUTO(:,2),'y','LineWidth',1),hold on,
%     close all
end

% cs = regionprops(bwper,'all')