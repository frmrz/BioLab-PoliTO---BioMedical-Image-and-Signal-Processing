close all
clear
clc

%% load results

load('/media/francesco/DEV001/PROJECT-CUBS-SPIE/RESULTS/MATLAB/masksEvaluationConsensus-INTER-commonSupport.mat');
load('/media/francesco/DEV001/PROJECT-CUBS-SPIE/RESULTS/MATLAB/profileEvaluationConsensus-INTER-biasProfiles-CommonSupport.mat');

load('/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/intraOP.mat')
load('/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/interOP.mat')

intraOP = gt1bbTable{:,1};
interOP = gt2Table{:,1};

v1 = statsDice;
names = fieldnames(IMTstats);

%% 2) BA plots on IMT Bias

set(0, 'DefaultFigureVisible', 'off')

saveDirintraOP = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/RESULTS/MATLAB/Bland-AltmanConsensus-INTRA1';
saveDirinterOP = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/RESULTS/MATLAB/Bland-AltmanConsensus-INTER1';

if ~isfolder(saveDirintraOP)
    mkdir(saveDirintraOP)
end
if ~isfolder(saveDirinterOP)
    mkdir(saveDirinterOP)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% INTRA-OP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% comparisons between GTconsensus and Operator

for ii = 1 : length(names)
    data = IMTstats.(names{ii});
    names{ii}
    if ~contains(names{ii},'ORIGINAL') && ~isempty(data)
    
        corrinfo = {'n','SSE','rho (p)','eq'}; % stats to display of correlation scatter plot
        BAinfo = {'RPC(%)','IQR'}; % stats to display on Bland-ALtman plot
        limits = [0 4 0 4];
        colors = 'brg';
        
        A = Struct_Empty_To_Nan(data);
        
        %% make intra
        dataA = [];
        dataB = [];
        for k = 1 : length(A)
            tempName = A(k).filename;
            for w = 1 : length(intraOP)
                if strcmp(tempName,intraOP{w,1})
                    dataA = [dataA; A(k).IMTgt];
                    dataB = [dataB; A(k).IMTauto];
                end
            end
        end

%         dataA = [A.IMTgt]'; dataB = [A.IMTauto]';
    
        % make plot intra
        %             tit = [ref ' vs ' strrep(comp,'_','\_') ' - IMT(mm) ' fold ' dataset']; % figure title
        tit = []; % figure title
        gnames = {'GTconsensus',strrep(names{ii},'_','\_')};
        label = {'GTconsensus',strrep(names{ii},'_','\_')}; % Names of data sets
        
        [cr, fig, statsStruct] = BlandAltman(dataA, dataB,label,tit,gnames,...
                                            'corrInfo',corrinfo,'baInfo',BAinfo,...
                                            'axesLimits',limits,'colors',colors,...
                                            'data1Mode','Compare',...
                                            'forceZeroIntercept','off',...
                                            'showFitCI','on',...
                                            'baStatsMode','Gaussian',...
                                            'legend','off');
        saveas(fig,fullfile(saveDirintraOP,sprintf('%s.png',names{ii})));

        %% make inter
        dataA = [];
        dataB = [];
        for k = 1 : length(A)
            tempName = A(k).filename;
            for w = 1 : length(interOP)
                if strcmp(tempName,interOP{w,1})
                    dataA = [dataA; A(k).IMTgt];
                    dataB = [dataB; A(k).IMTauto];
                end
            end
        end

%         dataA = [A.IMTgt]'; dataB = [A.IMTauto]';
    
        % make plot interOP
        %             tit = [ref ' vs ' strrep(comp,'_','\_') ' - IMT(mm) ' fold ' dataset']; % figure title
        tit = []; % figure title
        gnames = {'GTconsensus',strrep(names{ii},'_','\_')};
        label = {'GTconsensus',strrep(names{ii},'_','\_')}; % Names of data sets
        
        [cr, fig, statsStruct] = BlandAltman(dataA, dataB,label,tit,gnames,...
                                            'corrInfo',corrinfo,'baInfo',BAinfo,...
                                            'axesLimits',limits,'colors',colors,...
                                            'data1Mode','Compare',...
                                            'forceZeroIntercept','off',...
                                            'showFitCI','on',...
                                            'baStatsMode','Gaussian',...
                                            'legend','off');
        saveas(fig,fullfile(saveDirinterOP,sprintf('%s.png',names{ii})));
    
    end
end


%% comparison between OP1-OP1b

data1 = IMTstats.GT1_INTERPOLATED;
data1s = IMTstats.GT1B_INTERPOLATED;

A1 = Struct_Empty_To_Nan(data1);
A1s = Struct_Empty_To_Nan(data1s);

dataA = [];
dataB = [];
for k = 1 : length(A1)
    tempName = A1(k).filename;
    for w = 1 : length(intraOP)
        if strcmp(tempName,intraOP{w,1})
            dataA = [dataA; A1(k).IMTauto];
        end
    end
end
for k = 1 : length(A1s)
    tempName = A1s(k).filename;
    for w = 1 : length(intraOP)
        if strcmp(tempName,intraOP{w,1})
            dataB = [dataB; A1s(k).IMTauto];
        end
    end
end

tit = []; % figure title
gnames = {'GT1','GT1s'};
label = {'GT1','GT1s'}; % Names of data sets

[cr, fig, statsStruct] = BlandAltman(dataA, dataB,label,tit,gnames,...
                                    'corrInfo',corrinfo,'baInfo',BAinfo,...
                                    'axesLimits',limits,'colors',colors,...
                                    'data1Mode','Compare',...
                                    'forceZeroIntercept','off',...
                                    'showFitCI','on',...
                                    'baStatsMode','Gaussian',...
                                    'legend','off');
saveas(fig,fullfile(saveDirintraOP,'intraOP.png'));

%% comparison between OP1-UNET

data1 = IMTstats.GT1_INTERPOLATED;
data1s = IMTstats.UNET_96x96_v4;

A1 = Struct_Empty_To_Nan(data1);
A1s = Struct_Empty_To_Nan(data1s);

dataA = [];
dataB = [];
for k = 1 : length(A1)
    tempName = A1(k).filename;
    for w = 1 : length(intraOP)
        if strcmp(tempName,intraOP{w,1})
            dataA = [dataA; A1(k).IMTauto];
        end
    end
end
for k = 1 : length(A1s)
    tempName = A1s(k).filename;
    for w = 1 : length(intraOP)
        if strcmp(tempName,intraOP{w,1})
            dataB = [dataB; A1s(k).IMTauto];
        end
    end
end

tit = []; % figure title
gnames = {'GT1','UNET_96x96_v4'};
label = {'GT1','UNET_96x96_v4'}; % Names of data sets

[cr, fig, statsStruct] = BlandAltman(dataA, dataB,label,tit,gnames,...
                                    'corrInfo',corrinfo,'baInfo',BAinfo,...
                                    'axesLimits',limits,'colors',colors,...
                                    'data1Mode','Compare',...
                                    'forceZeroIntercept','off',...
                                    'showFitCI','on',...
                                    'baStatsMode','Gaussian',...
                                    'legend','off');
saveas(fig,fullfile(saveDirintraOP,'UNET_96x96_v4_vs_GT1.png'));

%% comparison between OP1-UNET

data1 = IMTstats.GT1B_INTERPOLATED;
data1s = IMTstats.UNET_96x96_v4;

A1 = Struct_Empty_To_Nan(data1);
A1s = Struct_Empty_To_Nan(data1s);

dataA = [];
dataB = [];
for k = 1 : length(A1)
    tempName = A1(k).filename;
    for w = 1 : length(intraOP)
        if strcmp(tempName,intraOP{w,1})
            dataA = [dataA; A1(k).IMTauto];
        end
    end
end
for k = 1 : length(A1s)
    tempName = A1s(k).filename;
    for w = 1 : length(intraOP)
        if strcmp(tempName,intraOP{w,1})
            dataB = [dataB; A1s(k).IMTauto];
        end
    end
end

tit = []; % figure title
gnames = {'GT1s','UNET_96x96_v4'};
label = {'GT1s','UNET_96x96_v4'}; % Names of data sets

[cr, fig, statsStruct] = BlandAltman(dataA, dataB,label,tit,gnames,...
                                    'corrInfo',corrinfo,'baInfo',BAinfo,...
                                    'axesLimits',limits,'colors',colors,...
                                    'data1Mode','Compare',...
                                    'forceZeroIntercept','off',...
                                    'showFitCI','on',...
                                    'baStatsMode','Gaussian',...
                                    'legend','off');
saveas(fig,fullfile(saveDirintraOP,'UNET_96x96_v4_vs_GT1s.png'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% INTER-OP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% comparison between OP1-OP2

data1 = IMTstats.GT1_INTERPOLATED;
data1s = IMTstats.GT2_INTERPOLATED;

A1 = Struct_Empty_To_Nan(data1);
A1s = Struct_Empty_To_Nan(data1s);

dataA = [];
dataB = [];
for k = 1 : length(A1)
    tempName = A1(k).filename;
    for w = 1 : length(interOP)
        if strcmp(tempName,interOP{w,1})
            dataA = [dataA; A1(k).IMTauto];
        end
    end
end
for k = 1 : length(A1s)
    tempName = A1s(k).filename;
    for w = 1 : length(interOP)
        if strcmp(tempName,interOP{w,1})
            dataB = [dataB; A1s(k).IMTauto];
        end
    end
end

tit = []; % figure title
gnames = {'GT1','GT2'};
label = {'GT1','GT2'}; % Names of data sets

[cr, fig, statsStruct] = BlandAltman(dataA, dataB,label,tit,gnames,...
                                    'corrInfo',corrinfo,'baInfo',BAinfo,...
                                    'axesLimits',limits,'colors',colors,...
                                    'data1Mode','Compare',...
                                    'forceZeroIntercept','off',...
                                    'showFitCI','on',...
                                    'baStatsMode','Gaussian',...
                                    'legend','off');
saveas(fig,fullfile(saveDirinterOP,'OP1vsOP2.png'));

%% comparison between OP1-OP3

data1 = IMTstats.GT1_INTERPOLATED;
data1s = IMTstats.GT3_INTERPOLATED;

A1 = Struct_Empty_To_Nan(data1);
A1s = Struct_Empty_To_Nan(data1s);

dataA = [];
dataB = [];
for k = 1 : length(A1)
    tempName = A1(k).filename;
    for w = 1 : length(interOP)
        if strcmp(tempName,interOP{w,1})
            dataA = [dataA; A1(k).IMTauto];
        end
    end
end
for k = 1 : length(A1s)
    tempName = A1s(k).filename;
    for w = 1 : length(interOP)
        if strcmp(tempName,interOP{w,1})
            dataB = [dataB; A1s(k).IMTauto];
        end
    end
end

tit = []; % figure title
gnames = {'GT1','GT3'};
label = {'GT1','GT3'}; % Names of data sets

[cr, fig, statsStruct] = BlandAltman(dataA, dataB,label,tit,gnames,...
                                    'corrInfo',corrinfo,'baInfo',BAinfo,...
                                    'axesLimits',limits,'colors',colors,...
                                    'data1Mode','Compare',...
                                    'forceZeroIntercept','off',...
                                    'showFitCI','on',...
                                    'baStatsMode','Gaussian',...
                                    'legend','off');
saveas(fig,fullfile(saveDirinterOP,'OP1vsOP3.png'));

%% comparison between OP2-OP3

data1 = IMTstats.GT2_INTERPOLATED;
data1s = IMTstats.GT3_INTERPOLATED;

A1 = Struct_Empty_To_Nan(data1);
A1s = Struct_Empty_To_Nan(data1s);

dataA = [];
dataB = [];
for k = 1 : length(A1)
    tempName = A1(k).filename;
    for w = 1 : length(interOP)
        if strcmp(tempName,interOP{w,1})
            dataA = [dataA; A1(k).IMTauto];
        end
    end
end
for k = 1 : length(A1s)
    tempName = A1s(k).filename;
    for w = 1 : length(interOP)
        if strcmp(tempName,interOP{w,1})
            dataB = [dataB; A1s(k).IMTauto];
        end
    end
end

tit = []; % figure title
gnames = {'GT2','GT3'};
label = {'GT2','GT3'}; % Names of data sets

[cr, fig, statsStruct] = BlandAltman(dataA, dataB,label,tit,gnames,...
                                    'corrInfo',corrinfo,'baInfo',BAinfo,...
                                    'axesLimits',limits,'colors',colors,...
                                    'data1Mode','Compare',...
                                    'forceZeroIntercept','off',...
                                    'showFitCI','on',...
                                    'baStatsMode','Gaussian',...
                                    'legend','off');
saveas(fig,fullfile(saveDirinterOP,'OP2vsOP3.png'));

%% comparison between OP1-UNET

data1 = IMTstats.GT1_INTERPOLATED;
data1s = IMTstats.UNET_96x96_v4;

A1 = Struct_Empty_To_Nan(data1);
A1s = Struct_Empty_To_Nan(data1s);

dataA = [];
dataB = [];
for k = 1 : length(A1)
    tempName = A1(k).filename;
    for w = 1 : length(interOP)
        if strcmp(tempName,interOP{w,1})
            dataA = [dataA; A1(k).IMTauto];
        end
    end
end
for k = 1 : length(A1s)
    tempName = A1s(k).filename;
    for w = 1 : length(interOP)
        if strcmp(tempName,interOP{w,1})
            dataB = [dataB; A1s(k).IMTauto];
        end
    end
end

tit = []; % figure title
gnames = {'GT1','UNET_96x96_v4'};
label = {'GT1','UNET_96x96_v4'}; % Names of data sets

[cr, fig, statsStruct] = BlandAltman(dataA, dataB,label,tit,gnames,...
                                    'corrInfo',corrinfo,'baInfo',BAinfo,...
                                    'axesLimits',limits,'colors',colors,...
                                    'data1Mode','Compare',...
                                    'forceZeroIntercept','off',...
                                    'showFitCI','on',...
                                    'baStatsMode','Gaussian',...
                                    'legend','off');
saveas(fig,fullfile(saveDirinterOP,'UNET_96x96_v4_vs_GT1.png'));

%% comparison between OP2-UNET

data1 = IMTstats.GT2_INTERPOLATED;
data1s = IMTstats.UNET_96x96_v4;

A1 = Struct_Empty_To_Nan(data1);
A1s = Struct_Empty_To_Nan(data1s);

dataA = [];
dataB = [];
for k = 1 : length(A1)
    tempName = A1(k).filename;
    for w = 1 : length(interOP)
        if strcmp(tempName,interOP{w,1})
            dataA = [dataA; A1(k).IMTauto];
        end
    end
end
for k = 1 : length(A1s)
    tempName = A1s(k).filename;
    for w = 1 : length(interOP)
        if strcmp(tempName,interOP{w,1})
            dataB = [dataB; A1s(k).IMTauto];
        end
    end
end

tit = []; % figure title
gnames = {'GT2','UNET_96x96_v4'};
label = {'GT2','UNET_96x96_v4'}; % Names of data sets

[cr, fig, statsStruct] = BlandAltman(dataA, dataB,label,tit,gnames,...
                                    'corrInfo',corrinfo,'baInfo',BAinfo,...
                                    'axesLimits',limits,'colors',colors,...
                                    'data1Mode','Compare',...
                                    'forceZeroIntercept','off',...
                                    'showFitCI','on',...
                                    'baStatsMode','Gaussian',...
                                    'legend','off');
saveas(fig,fullfile(saveDirinterOP,'UNET_96x96_v4_vs_GT2.png'));

%% comparison between OP3-UNET

data1 = IMTstats.GT3_INTERPOLATED;
data1s = IMTstats.UNET_96x96_v4;

A1 = Struct_Empty_To_Nan(data1);
A1s = Struct_Empty_To_Nan(data1s);

dataA = [];
dataB = [];
for k = 1 : length(A1)
    tempName = A1(k).filename;
    for w = 1 : length(interOP)
        if strcmp(tempName,interOP{w,1})
            dataA = [dataA; A1(k).IMTauto];
        end
    end
end
for k = 1 : length(A1s)
    tempName = A1s(k).filename;
    for w = 1 : length(interOP)
        if strcmp(tempName,interOP{w,1})
            dataB = [dataB; A1s(k).IMTauto];
        end
    end
end

tit = []; % figure title
gnames = {'GT3','UNET_96x96_v4'};
label = {'GT3','UNET_96x96_v4'}; % Names of data sets

[cr, fig, statsStruct] = BlandAltman(dataA, dataB,label,tit,gnames,...
                                    'corrInfo',corrinfo,'baInfo',BAinfo,...
                                    'axesLimits',limits,'colors',colors,...
                                    'data1Mode','Compare',...
                                    'forceZeroIntercept','off',...
                                    'showFitCI','on',...
                                    'baStatsMode','Gaussian',...
                                    'legend','off');
saveas(fig,fullfile(saveDirinterOP,'UNET_96x96_v4_vs_GT3.png'));












%% boxplots HD

saveDir = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/RESULTS/MATLAB/BoxplotsConsensus-v4';

if ~isfolder(saveDir)
    mkdir(saveDir)
end

HM_LI_Vector = [];
HM_MA_Vector = [];
experimentVector = [];

for ii = 1 : length(names)
    data = IMTstats.(names{ii});

    if ~contains(names{ii},'ORIGINAL') && ~isempty(data) && ~contains(names{ii},'64') && ~contains(names{ii},'480')
        HM_LI_Vector = [HM_LI_Vector; [data(:).HM_LI]'];
        HM_MA_Vector = [HM_MA_Vector; [data(:).HM_MA]'];
        category = categorical(cellstr(names(ii)));
        experimentVector = [experimentVector; repmat(category,length([data(:).HM_LI]),1)];
    end
end

h=figure;
h.Position = [ 0 0 1800 900];
h.Visible = 'off';subplot(211),boxchart(experimentVector,HM_LI_Vector),ylabel('HM\_LI'),ylim([0 1]),set(gca,'YTick',0:0.1:1);
subplot(212),boxchart(experimentVector,HM_MA_Vector),ylabel('HM\_MA'),ylim([0 1]),set(gca,'YTick',0:0.1:1);
sgtitle('HAUSDORFF');

saveas(h,fullfile(saveDir,'Hausdorff96.png'));

%% boxplots Dice
statsDice = Struct_Empty_To_Nan(statsDice)
colsDice = struct2table(statsDice);

% DiceMat = table2array(colsDice(:,2:end)); %all
% DiceMat = table2array(colsDice(:,2:6)); %480
% DiceMat = table2array(colsDice(:,7:11)); %64
DiceMat = table2array(colsDice(:,12:end)); %96

DiceVect = reshape(DiceMat,[numel(DiceMat),1]);

exps = {...
%     '480x480-v0','480x480-v1','480x480-v2','480x480-v3','480x480-v4',...
%         '64x64-v0','64x64-v1','64x64-v2','64x64-v3','64x64-v4',...
        '96x96-v0','96x96-v1','96x96-v2','96x96-v3','96x96-v4',
        };
experimentVector = [];

for ii = 1 : length(exps)
    category = categorical(cellstr(exps(ii)));
    experimentVector = [experimentVector; repmat(category,height(colsDice),1)];
end

h=figure;
h.Position = [ 0 0 1800 900];
h.Visible = 'on';
boxchart(experimentVector,DiceVect),ylabel('Dice'),ylim([0 1]),set(gca,'YTick',0:0.1:1),title('DICE');
saveas(h,fullfile(saveDir,'Dice96.png'));

%% Stats on 96x96


    %% DICE
    [statsDice] = Struct_Empty_To_Nan(statsDice);
    colsDice = struct2table(statsDice);
    Dice96 = table2array(colsDice(:,2:6)); %64
%     Dice96 = cell2mat(Dice96); %64
    
    n96 = sum(isnan(Dice96));

    m96 = mean(Dice96,1,"omitnan");
    s96 = std(Dice96,1,"omitnan");
    
    %% test if distributions are normal : 0 -normal 1 -not normal
    hChiD.v0 = chi2gof(Dice96(:,1));
    hChiD.v1 = chi2gof(Dice96(:,2));
    hChiD.v2 = chi2gof(Dice96(:,3));
    hChiD.v3 = chi2gof(Dice96(:,4));
    hChiD.v4 = chi2gof(Dice96(:,5));
    
    %% not normal -> signrank
    
    for i = 1 : width(Dice96)
        for k = i + 1 : width(Dice96)
            pw(i,k) = signrank(Dice96(:,i),Dice96(:,k));
        end
    end

        %% imHM
        a = struct2table(Struct_Empty_To_Nan(IMTstats.UNET_96x96_v4));
        colshd = struct2table(Struct_Empty_To_Nan(statsHD));
        imHAUS96 = table2array(colshd(:,1:5)); %64
        imHAUS96(imHAUS96==inf)=nan;

        m96hd = mean(imHAUS96,1,"omitnan");
        s96hd = std(imHAUS96,1,"omitnan");
        
        %% test if distributions are normal : 0 -normal 1 -not normal
        hChiH.v0 = chi2gof(imHAUS96(:,1));
        hChiH.v1 = chi2gof(imHAUS96(:,2));
        hChiH.v2 = chi2gof(imHAUS96(:,3));
        hChiH.v3 = chi2gof(imHAUS96(:,4));
        hChiH.v4 = chi2gof(imHAUS96(:,5));
        
        %% not normal -> signrank
        
        for i = 1 : width(imHAUS96)
            for k = i + 1 : width(imHAUS96)
                pwGHim(i,k) = signrank(imHAUS96(:,i),imHAUS96(:,k));
            end
        end
    
            %% AbsBiasIMT
        
            temp = Struct_Empty_To_Nan(IMTstats.UNET_96x96_v0);
            AB_IMT96.UNET_96x96_v0 = [temp(:).BiasIMT]';
        
            temp = Struct_Empty_To_Nan(IMTstats.UNET_96x96_v1);
            AB_IMT96.UNET_96x96_v1 = [temp(:).BiasIMT]';
        
            temp = Struct_Empty_To_Nan(IMTstats.UNET_96x96_v2);
            AB_IMT96.UNET_96x96_v2 = [temp(:).BiasIMT]';
        
%             temp = Struct_Empty_To_Nan(IMTstats.UNET_96x96_v3);
%             AB_IMT96.UNET_96x96_v3 = [temp(:).BiasIMT]';
%         
            temp = Struct_Empty_To_Nan(IMTstats.UNET_96x96_v4);
            AB_IMT96.UNET_96x96_v4 = [temp(:).BiasIMT]';
        
            AB_IMT96_mat = struct2array(AB_IMT96);
        
            ab96m = mean(abs(AB_IMT96_mat),1,"omitnan");
            ab96s = std(abs(AB_IMT96_mat),1,"omitnan");
        
            %% test if distributions are normal : 0 -normal 1 -not normal
            hChiB.v0 = chi2gof(AB_IMT96_mat(:,1));
            hChiB.v1 = chi2gof(AB_IMT96_mat(:,2));
            hChiB.v2 = chi2gof(AB_IMT96_mat(:,3));
%             hChiB.v3 = chi2gof(AB_IMT96_mat(:,4));
            hChiB.v4 = chi2gof(AB_IMT96_mat(:,4));
        
            %% not normal -> signrank
            
            for i = 1 : width(AB_IMT96_mat)
                for k = i + 1 : width(AB_IMT96_mat)
                    ABpw(i,k) = signrank(AB_IMT96_mat(:,i),AB_IMT96_mat(:,k));
                end
            end

                %% HM LI
            
                temp = Struct_Empty_To_Nan(IMTstats.UNET_96x96_v0);
                HM_LI96.UNET_96x96_v0 = [temp(:).HM_LI]';
            
                temp = Struct_Empty_To_Nan(IMTstats.UNET_96x96_v1);
                HM_LI96.UNET_96x96_v1 = [temp(:).HM_LI]';
            
                temp = Struct_Empty_To_Nan(IMTstats.UNET_96x96_v2);
                HM_LI96.UNET_96x96_v2 = [temp(:).HM_LI]';
            
%                 temp = Struct_Empty_To_Nan(IMTstats.UNET_96x96_v3);
%                 HM_LI96.UNET_96x96_v3 = [temp(:).HM_LI]';
            
                temp = Struct_Empty_To_Nan(IMTstats.UNET_96x96_v4);
                HM_LI96.UNET_96x96_v4 = [temp(:).HM_LI]';
            
                HM_LI96_mat = struct2array(HM_LI96);
            
                hmli96m = mean(abs(HM_LI96_mat),1,"omitnan");
                hmli96s = std(abs(HM_LI96_mat),1,"omitnan");
            
                %% test if distributions are normal : 0 -normal 1 -not normal
                hChiL.v0 = chi2gof(HM_LI96_mat(:,1));
                hChiL.v1 = chi2gof(HM_LI96_mat(:,2));
                hChiL.v2 = chi2gof(HM_LI96_mat(:,3));
%                 hChiL.v3 = chi2gof(HM_LI96_mat(:,4));
                hChiL.v4 = chi2gof(HM_LI96_mat(:,4));
            
                %% not normal -> signrank
                
                for i = 1 : width(HM_LI96_mat)
                    for k = i + 1 : width(HM_LI96_mat)
                        HM_LIpw(i,k) = signrank(HM_LI96_mat(:,i),HM_LI96_mat(:,k));
                    end
                end

                     %% HM MA
                
                    temp = Struct_Empty_To_Nan(IMTstats.UNET_96x96_v0);
                    HM_MA96.UNET_96x96_v0 = [temp(:).HM_MA]';
                
                    temp = Struct_Empty_To_Nan(IMTstats.UNET_96x96_v1);
                    HM_MA96.UNET_96x96_v1 = [temp(:).HM_MA]';
                
                    temp = Struct_Empty_To_Nan(IMTstats.UNET_96x96_v2);
                    HM_MA96.UNET_96x96_v2 = [temp(:).HM_MA]';
                
%                     temp = Struct_Empty_To_Nan(IMTstats.UNET_96x96_v3);
%                     HM_MA96.UNET_96x96_v3 = [temp(:).HM_MA]';
                
                    temp = Struct_Empty_To_Nan(IMTstats.UNET_96x96_v4);
                    HM_MA96.UNET_96x96_v4 = [temp(:).HM_MA]';
                
                    HM_MA96_mat = struct2array(HM_MA96);
                
                    hmMA96m = mean(abs(HM_MA96_mat),1,"omitnan");
                    hmMA96s = std(abs(HM_MA96_mat),1,"omitnan");
                
                    %% test if distributions are normal : 0 -normal 1 -not normal
                    hChiM.v0 = chi2gof(HM_MA96_mat(:,1));
                    hChiM.v1 = chi2gof(HM_MA96_mat(:,2));
                    hChiM.v2 = chi2gof(HM_MA96_mat(:,3));
%                     hChiM.v3 = chi2gof(HM_MA96_mat(:,4));
                    hChiM.v4 = chi2gof(HM_MA96_mat(:,4));
                
                    %% not normal -> signrank
                    
                    for i = 1 : width(HM_MA96_mat)
                        for k = i + 1 : width(HM_MA96_mat)
                            HM_MApw(i,k) = signrank(HM_MA96_mat(:,i),HM_MA96_mat(:,k));
                        end
                    end

%% excel table

Tab = [m96(1,1),s96(1,1),m96(1,2),s96(1,2),m96(1,3),s96(1,3),m96(1,4),s96(1,4);...
    m96hd(1,1),s96hd(1,1),m96hd(1,2),s96hd(1,2),m96hd(1,3),s96hd(1,3),m96hd(1,4),s96hd(1,4);...
    ab96m(1,1),ab96s(1,1),ab96m(1,2),ab96s(1,2),ab96m(1,3),ab96s(1,3),ab96m(1,4),ab96s(1,4);...
    hmli96m(1,1),hmli96s(1,1),hmli96m(1,2),hmli96s(1,2),hmli96m(1,3),hmli96s(1,3),hmli96m(1,4),hmli96s(1,4);...
    hmMA96m(1,1),hmMA96s(1,1),hmMA96m(1,2),hmMA96s(1,2),hmMA96m(1,3),hmMA96s(1,3),hmMA96m(1,4),hmMA96s(1,4)];

%% order profiles for ICC
