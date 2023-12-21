close all
clear
clc

%% script for the organization and interpolation of profiles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% part 1: organize and interpolate manual profiles from multi-op
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% folders

dataDir = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA';

outProfilesDir = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES/GT1-INTERPOLATED';
outCommonProfilesDir = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES/GT1-COMMON';

outProfilesDir2 = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES/GT2-INTERPOLATED';
outCommonProfilesDir2 = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES/GT2-COMMON';

outProfilesDir3 = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES/GT3-INTERPOLATED';
outCommonProfilesDir3 = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES/GT3-COMMON';

outputMerge = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES/GT-CONSENSUS';
outputCommonMerge = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES/GT-CONSENSUS-COMMON';

outputNet = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES/UNET-96x96-v4';
outputCommonNet = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/PROFILES/UNET-96x96-v4-COMMON';

load('/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/interOP.mat')

%% files

interOP = gt2Table{:,1};

%% loop across centers

interpProfiles = struct('name',[],'center',[],'GT1_LI',[],'GT1_MA',[], ...
    'GT2_LI',[],'GT2_MA',[], ...
    'GT3_LI',[],'GT3_MA',[], ...
    'cons_LI',[],'cons_MA',[],...
    'unet_LI',[],'unet_MA',[]);
q=1;
% figure,
imgDir = fullfile(dataDir,"DEVELOPMENT/TEST/IMAGES");
CFdir = fullfile(dataDir,'DEVELOPMENT/TEST/CF');

for ii = 1 : length(interOP)
    fprintf('%d\n',ii);

    filename = interOP{ii,1};
    filename_img = [interOP{ii,1} '.png'];
    img = imread(fullfile(imgDir,filename_img));

    try
        tempCF = load(fullfile(CFdir,[interOP{ii,1} '_CF.txt']));
    catch
        tempCF = 0.06;
    end

    %% load profiles
    try
        GT1_LI = load(fullfile(outProfilesDir,[filename '-GTLI.txt']));
        GT1_MA = load(fullfile(outProfilesDir,[filename '-GTMA.txt']));
        [temp.IMT_GT1,~,~] = PolyDistMethod(GT1_LI,GT1_MA);
        temp.IMT_GT1 = temp.IMT_GT1*tempCF;

    catch
        GT1_MA = nan;
        GT1_LI = nan;
        temp.IMT_GT1 = nan;
    end
    try
        GT2_LI = load(fullfile(outProfilesDir2,[filename '-LI.txt']));
        GT2_MA = load(fullfile(outProfilesDir2,[filename '-MA.txt']));
        [temp.IMT_GT2,~,~] = PolyDistMethod(GT2_LI,GT2_MA);
        temp.IMT_GT2 = temp.IMT_GT2*tempCF;
        
        try
            % common support on interp profiles
            temp_GT1_LI = GT1_LI(:,ismember(round(GT1_LI(1,:)),round(GT2_LI(1,:))));       
            temp_GT2_LI = GT2_LI(:,ismember(round(GT2_LI(1,:)),round(GT1_LI(1,:))));       
            temp_GT1_MA = GT1_MA(:,ismember(round(GT1_MA(1,:)),round(GT2_MA(1,:))));
            temp_GT2_MA = GT2_MA(:,ismember(round(GT2_MA(1,:)),round(GT1_MA(1,:))));

            %% Hausdorff LI
            HM_LI_12 = HaussdorfDistance(temp_GT1_LI,temp_GT2_LI);
            tempHD.HM_LI_12 = HM_LI_12*tempCF;
        
            %% Hausdorff MA
            HM_MA_12 = HaussdorfDistance(temp_GT1_MA,temp_GT2_MA);
            tempHD.HM_MA_12 = HM_MA_12*tempCF;
        catch
            tempHD.HM_LI_12 = nan;
            tempHD.HM_MA_12 = nan;
        end

    catch
        GT2_MA = nan;
        GT2_LI = nan;
        temp.IMT_GT2 = nan;
        tempHD.HM_LI_12 = nan;
        tempHD.HM_MA_12 = nan;
    end
    try
        GT3_LI = load(fullfile(outProfilesDir3,[filename '-LI.txt']));
        GT3_MA = load(fullfile(outProfilesDir3,[filename '-MA.txt']));
        [temp.IMT_GT3,~,~] = PolyDistMethod(GT3_LI,GT3_MA);
        temp.IMT_GT3 = temp.IMT_GT3*tempCF;

        try
            % common support on interp profiles
            temp_GT1_LI = GT1_LI(:,ismember(round(GT1_LI(1,:)),round(GT3_LI(1,:))));      
            temp_GT3_LI = GT3_LI(:,ismember(round(GT3_LI(1,:)),round(GT1_LI(1,:))));      
            temp_GT1_MA = GT1_MA(:,ismember(round(GT1_MA(1,:)),round(GT3_MA(1,:))));
            temp_GT3_MA = GT3_MA(:,ismember(round(GT3_MA(1,:)),round(GT1_MA(1,:))));

            %% Hausdorff LI
            HM_LI_13 = HaussdorfDistance(temp_GT1_LI,temp_GT3_LI);
            tempHD.HM_LI_13 = HM_LI_13*tempCF;
        
            %% Hausdorff MA
            HM_MA_13 = HaussdorfDistance(temp_GT1_MA,temp_GT3_MA);
            tempHD.HM_MA_13 = HM_MA_13*tempCF;
        catch
            tempHD.HM_LI_13 = nan;
            tempHD.HM_MA_13 = nan;
        end

        try
            % common support on interp profiles
            temp_GT2_LI = GT2_LI(:,ismember(round(GT2_LI(1,:)),round(GT3_LI(1,:))));      
            temp_GT3_LI = GT3_LI(:,ismember(round(GT3_LI(1,:)),round(GT2_LI(1,:))));       
            temp_GT2_MA = GT2_MA(:,ismember(round(GT2_MA(1,:)),round(GT3_MA(1,:))));
            temp_GT3_MA = GT3_MA(:,ismember(round(GT3_MA(1,:)),round(GT2_MA(1,:))));

            %% Hausdorff LI
            HM_LI_23 = HaussdorfDistance(temp_GT2_LI,temp_GT3_LI);
            tempHD.HM_LI_23 = HM_LI_23*tempCF;
        
            %% Hausdorff MA
            HM_MA_23 = HaussdorfDistance(temp_GT2_MA,temp_GT3_MA);
            tempHD.HM_MA_23 = HM_MA_23*tempCF;
        catch
            tempHD.HM_LI_23 = nan;
            tempHD.HM_MA_23 = nan;
        end

    catch
        GT3_MA = nan;
        GT3_LI = nan;
        temp.IMT_GT3 = nan;
        tempHD.HM_LI_23 = nan;
        tempHD.HM_MA_23 = nan;
        tempHD.HM_LI_13 = nan;
        tempHD.HM_MA_13 = nan;

    end
    try
        GTC_LI = load(fullfile(outputMerge,[filename '-LI.txt']));
        GTC_MA = load(fullfile(outputMerge,[filename '-MA.txt']));
        [temp.IMT_GTC,~,~] = PolyDistMethod(GTC_LI,GTC_MA);
        temp.IMT_GTC = temp.IMT_GTC*tempCF;
    catch
        GTC_MA = nan;
        GTC_LI = nan;
        temp.IMT_GTC = nan;

    end
    try
        U_LI = load(fullfile(outputNet,[filename '-LI.txt']));
        U_MA = load(fullfile(outputNet,[filename '-MA.txt']));
        U_LI = TurnRow(U_LI);
        U_MA = TurnRow(U_MA);

        [temp.IMT_U,~,~]   = PolyDistMethod(U_LI,U_MA);
        temp.IMT_U   = temp.IMT_U*tempCF;

        try
            % common support on interp profiles
            temp_U_LI   = U_LI(:,ismember(round(U_LI(1,:)),round(GT1_LI(1,:))));       
            temp_GT1_LI = GT1_LI(:,ismember(round(GT1_LI(1,:)),round(U_LI(1,:))));       
            temp_U_MA   = U_MA(:,ismember(round(U_MA(1,:)),round(GT1_MA(1,:))));
            temp_GT1_MA = GT1_MA(:,ismember(round(GT1_MA(1,:)),round(U_MA(1,:))));

            %% Hausdorff LI
            HM_LI_U1 = HaussdorfDistance(temp_U_LI,temp_GT1_LI);
            tempHD.HM_LI_U1 = HM_LI_U1*tempCF;
        
            %% Hausdorff MA
            HM_MA_U1 = HaussdorfDistance(temp_U_MA,temp_GT1_MA);
            tempHD.HM_MA_U1 = HM_MA_U1*tempCF;
        catch
            tempHD.HM_LI_U1 = nan;
            tempHD.HM_MA_U1 = nan;
        end
        try
            % common support on interp profiles
            temp_U_LI = U_LI(:,ismember(round(U_LI(1,:)),round(GT2_LI(1,:))));       
            temp_GT2_LI = GT2_LI(:,ismember(round(GT2_LI(1,:)),round(U_LI(1,:))));       
            temp_U_MA = U_MA(:,ismember(round(U_MA(1,:)),round(GT2_MA(1,:))));
            temp_GT2_MA = GT2_MA(:,ismember(round(GT2_MA(1,:)),round(U_MA(1,:))));

            %% Hausdorff LI
            HM_LI_U2 = HaussdorfDistance(temp_U_LI,temp_GT2_LI);
            tempHD.HM_LI_U2 = HM_LI_U2*tempCF;
        
            %% Hausdorff MA
            HM_MA_U2 = HaussdorfDistance(temp_U_MA,temp_GT2_MA);
            tempHD.HM_MA_U2 = HM_MA_U2*tempCF;
        catch
            tempHD.HM_LI_U2 = nan;
            tempHD.HM_MA_U2 = nan;
        end
        try
            % common support on interp profiles
            temp_U_LI = U_LI(:,ismember(round(U_LI(1,:)),round(GT3_LI(1,:))));       
            temp_GT3_LI = GT3_LI(:,ismember(round(GT3_LI(1,:)),round(U_LI(1,:))));       
            temp_U_MA = U_MA(:,ismember(round(U_MA(1,:)),round(GT3_MA(1,:))));
            temp_GT3_MA = GT3_MA(:,ismember(round(GT3_MA(1,:)),round(U_MA(1,:))));

            %% Hausdorff LI
            HM_LI_U3 = HaussdorfDistance(temp_U_LI,temp_GT3_LI);
            tempHD.HM_LI_U3 = HM_LI_U3*tempCF;
        
            %% Hausdorff MA
            HM_MA_U3 = HaussdorfDistance(temp_U_MA,temp_GT3_MA);
            tempHD.HM_MA_U3 = HM_MA_U3*tempCF;
        catch
            tempHD.HM_LI_U3 = nan;
            tempHD.HM_MA_U3 = nan;
        end

    catch
        U_LI = nan;
        U_MA = nan;
        temp.IMT_U = nan;
        tempHD.HM_LI_U1 = nan;
        tempHD.HM_MA_U1 = nan;
        tempHD.HM_LI_U2 = nan;
        tempHD.HM_MA_U2 = nan;
        tempHD.HM_LI_U3 = nan;
        tempHD.HM_MA_U3 = nan;
    end       

    %% save profiles

    temp.name = filename;

    IMT(q) = temp;
    HD(q) = tempHD;
    clear temp 
    clear GT1 GT1_MA GT1B_LI GT1B_MA GT2_LI GT2_MA GT3_LI GT3_MA
    
    q = q + 1;

end

%% Stats on black box tempCFel
IMTvalues = struct2cell(IMT');
IMTvalues = struct2cell(IMT);
IMTvalues = cell2mat(IMTvalues([1,2,3,5],:))';
IMTtable = struct2table(IMT);

IMTtable = IMTtable(:,[1,2,3,5]);
writetable(IMTtable,fullfile("BA-Fig3/toBA.xlsx"))

countNan = sum(isnan(IMTvalues));
nanRows = sum(IMTvalues,2);
IMTvalues = IMTvalues(~isnan(nanRows),:);

for i = 1 : 5
    for j = 1 : 5
        mBias(i,j) = mean(abs(IMTvalues(:,i)-IMTvalues(:,j)),"all","omitnan");
        sBias(i,j) = std(abs(IMTvalues(:,i)-IMTvalues(:,j)),'omitnan');
        [r(i,j), LB(i,j), UB(i,j), ~, ~, ~, p(i,j)] = ICC([IMTvalues(:,i),IMTvalues(:,j)], '1-1');
%         [rA(i,j), LBA(i,j), UBA(i,j), ~, ~, ~, p(i,j)] = ICC([IMTvalues(:,i),IMTvalues(:,j),IMTvalues(:,3)], 'A-1');
%         [rA(i,j), LBA(i,j), UBA(i,j), ~, ~, ~, p(i,j)] = ICC([mean([IMTvalues(:,i),IMTvalues(:,j),IMTvalues(:,j+1)],2),IMTvalues(:,5)], 'A-1');

        cS(i,j) = min(min(corr(IMTvalues(:,i),IMTvalues(:,j),'type','Spearman')));
    end
end
R = corrcoef(IMTvalues,'Rows','complete');

toComp = [1 2 3 5];
mBias = mBias(toComp,toComp);
sBias = sBias(toComp,toComp);
R = R(toComp,toComp);
r = r(toComp,toComp);
% rA = rA(toComp,toComp);
[ICCops, LBA(i,j), UBA(i,j), ~, ~, ~, p(i,j)] = ICC([IMTvalues(:,1),IMTvalues(:,2),IMTvalues(:,3)], 'A-1');
[ICCU, LBA(i,j), UBA(i,j), ~, ~, ~, p(i,j)] = ICC([mean([IMTvalues(:,1),IMTvalues(:,2),IMTvalues(:,3)],2),IMTvalues(:,4)], 'A-1');


%% Bland Altman between IMT measurments

corrinfo = {'n','SSE','rho (p)','eq'}; % stats to display of correlation scatter plot
BAinfo = {'RPC(%)','IQR'}; % stats to display on Bland-ALtman plot
limits = [0 3 0 3];
colors = 'brg';

names = {'GT1','GT2','GT3','UNET'};
saveDir = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/RESULTS/MATLAB/BA-Fig3';

for i = 1 : 4
    for j = i+1 : 4
        dataA = IMTvalues(:,toComp(i));
        dataB = IMTvalues(:,toComp(j));
        
        tit = []; % figure title
        gnames = {names{i},names{j}};
        label = {names{i},names{j}}; % Names of data sets
        
        [cr, fig, statsStruct] = BlandAltman(dataA, dataB,label,tit,gnames,...
                                            'corrInfo',corrinfo,'baInfo',BAinfo,...
                                            'axesLimits',limits,'colors',colors,...
                                            'data1mode','Compare',...
                                            'forceZeroIntercept','off',...
                                            'showFitCI','on',...
                                            'baStatsmode','Gaussian',...
                                            'legend','off');
        saveas(fig,fullfile(saveDir,sprintf('%s_vs_%s.png',names{i},names{j})));
    end
end

%% Hausdorff

HDvalues = struct2cell(Struct_Empty_To_Nan(HD));
HDvalues = cell2mat(HDvalues)';
nanRows = sum(HDvalues,2);
HDvalues = HDvalues(~isnan(nanRows),:);

namesHD = {'LI_12','MA_12','LI_13','MA_13','LI_23','MA_23','LI_1U','MA_1U','LI_2U','MA_2U','LI_3U','MA_3U'};
mHD = mean(HDvalues);
sHD = std(HDvalues);
p50HD = prctile(HDvalues,50);
p75HD = prctile(HDvalues,75);
p90HD = prctile(HDvalues,90);
p95HD = prctile(HDvalues,95);

tableHD = array2table([mHD;sHD;p50HD;p75HD;p90HD;p95HD]);
tableHD.Properties.VariableNames = namesHD;