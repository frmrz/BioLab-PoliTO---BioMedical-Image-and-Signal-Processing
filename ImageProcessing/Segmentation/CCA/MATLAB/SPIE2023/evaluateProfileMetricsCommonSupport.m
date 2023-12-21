close all
clear
clc

%% script for the evaluation metrics of profiles on TEST set results
% the profiles are the original or the interpolated GT ones
% for the unet outputs, the profiles are extracted fromthe images resized
% to original size

mainDir = '/media/francesco/DEV001/PROJECT-CUBS-SPIE';

imgDir = fullfile(mainDir,"DATA","DEVELOPMENT","TEST","IMAGES");
profilesDir = fullfile(mainDir,"DATA","DEVELOPMENT","TEST","PROFILES");
cfDir = fullfile(mainDir,"DATA","DEVELOPMENT","TEST","CF");

load('/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/interOP.mat')

%% loop across measurments
mode = 'mm'; % measure distances in mm or pxl

folders = dir(profilesDir);

% IMTstats = struct('GT1B_INTERPOLATED',[],'GT1B_ORIGINAL',[],...
%     'UNET_480X480_v0',[],'UNET_480X480_v1',[],'UNET_480X480_v2',[],'UNET_480X480_v3',[],...
%     'UNET_64x64_v0',[],'UNET_64x64_v1',[],'UNET_64x64_v2',[],'UNET_64x64_v3',[]);

% IMTstats = struct('UNET_480x480_v0',[],'UNET_480x480_v1',[],'UNET_480x480_v2',[],'UNET_480x480_v3',[],'UNET_480x480_v4',[],...
%     'UNET_64x64_v0',[],'UNET_64x64_v1',[],'UNET_64x64_v2',[],'UNET_64x64_v3',[],'UNET_64x64_v4',[],...
%     'UNET_96x96_v0',[],'UNET_96x96_v1',[],'UNET_96x96_v2',[],'UNET_96x96_v3',[],'UNET_96x96_v4',[]);

IMTstats = struct('UNET_96x96_v0',[],'UNET_96x96_v1',[],'UNET_96x96_v2',[],'UNET_96x96_v3',[],'UNET_96x96_v4',[]);

BIASstatsInterp = IMTstats;
BIASstats = IMTstats;
imgs = dir(imgDir);
no_auto_pf = 0;
no_gt_pf = 0;

toCompare = [22,23,24,26];

for k = 1 : length(toCompare)
    ii = toCompare(k);
    tempFolder = fullfile(profilesDir,folders(ii).name);
    
    intGTFolder = fullfile(profilesDir,'GT-CONSENSUS');

    GTfiles = dir(intGTFolder);

    % eval only on Unet folders

    if 1>0 %~contains(tempFolder,'GT')
            
        [~,folderName] = fileparts(tempFolder);
        folderName = strrep(folderName,'-','_');
        fprintf('\n\n ... Evaluating results of %s  %d of %d ... \n\n',folderName,k,length(toCompare))

        files = dir(tempFolder);
        q = 1;

        fxs = table2struct(gt2Table(:,1));
        for i = 1 :length(fxs)
            fxs(i).name = [fxs(i).name '.png'];
        end
        
        fxs = [struct('name','.');struct('name','..');fxs];

        % loop across GT files, NaN if not found
        for jj = 3 : length(fxs)
            
            filename = extractBefore(fxs(jj).name,'.png');
            fprintf('\n\t ... File %s ... exp %d of %d',filename,k,length(toCompare))

            try
                CF = load(fullfile(cfDir,strrep(fxs(jj).name,'.png','_CF.txt')));
            catch
                CF = 0.06;
            end

            try
                %% retrieve GT profiles
%                 LI_GT = TurnColumn(load(fullfile(GTFolder,strrep(fxs(jj).name,'.png','-LI.txt'))));
%                 MA_GT = TurnColumn(load(fullfile(GTFolder,strrep(fxs(jj).name,'.png','-MA.txt'))));
                
                try
                    LI_GT_interp = TurnColumn(load(fullfile(intGTFolder,strrep(fxs(jj).name,'.png','-LI.txt'))));
                    MA_GT_interp = TurnColumn(load(fullfile(intGTFolder,strrep(fxs(jj).name,'.png','-MA.txt'))));
    
                %% compare with current profiles
                    LI = TurnColumn(load(fullfile(tempFolder,strrep(fxs(jj).name,'.png','-LI.txt'))));
                    MA = TurnColumn(load(fullfile(tempFolder,strrep(fxs(jj).name,'.png','-MA.txt'))));
        
                    lenGT = numel(LI)/2*CF;

                    supportCommon = LI(ismember(round(LI(:,1)),round(LI_GT_interp(:,1))),:);
                    percentOverlapPred = numel(supportCommon(:,1))/numel(LI(:,1));

                    %% get general common support
                    t1 = TurnColumn(load(fullfile(fullfile(profilesDir,folders(toCompare(1)).name),strrep(fxs(jj).name,'.png','-LI.txt'))));
                    t2 = TurnColumn(load(fullfile(fullfile(profilesDir,folders(toCompare(2)).name),strrep(fxs(jj).name,'.png','-LI.txt'))));
                    t3 = TurnColumn(load(fullfile(fullfile(profilesDir,folders(toCompare(3)).name),strrep(fxs(jj).name,'.png','-LI.txt'))));
                    t4 = TurnColumn(load(fullfile(fullfile(profilesDir,folders(toCompare(4)).name),strrep(fxs(jj).name,'.png','-LI.txt'))));
                    tGT = LI_GT_interp;
                    
                    supp = 1:1000;
                    common = true(5,1000);

                    common(1,:) = ismember(supp,round(t1(:,1))) & common(1,:);
                    common(2,:) = ismember(supp,round(t2(:,1))) & common(2,:);
                    common(3,:) = ismember(supp,round(t3(:,1))) & common(3,:);
                    common(4,:) = ismember(supp,round(t4(:,1))) & common(4,:);
                    common(5,:) = ismember(supp,round(tGT(:,1))) & common(5,:);
                    
                    common = min(common);%sum(common)
                    indexes = find(common==1);

                    % if profile is interp make both compare else make only bias
                    if ~contains(tempFolder,'ORIGINAL') && lenGT > 1 && percentOverlapPred > 0.1
                        
                        % common support on interp profiles
                        comm.temp_LI = LI(ismember(round(LI(:,1)),indexes),:);               comm.temp_MA = MA(ismember(round(MA(:,1)),indexes),:);
                        comm.gt_LI = LI_GT_interp(ismember(round(LI_GT_interp(:,1)),indexes),:);       comm.gt_MA = MA_GT_interp(ismember(round(MA_GT_interp(:,1)),indexes),:);
                        
                        if length(comm.temp_LI)>2
                            tempIMTstats(q)        = LI_MA_stats_light(comm.gt_LI,comm.gt_MA,comm.temp_LI,comm.temp_MA,CF,filename,mode);
                            tempBIASstatsInterp(q) = LI_MA_stats_not_interp(comm.gt_LI,comm.gt_MA,comm.temp_LI,comm.temp_MA,CF,filename,mode);
                        else
                            tempIMTstats(q) = struct('filename',filename,'CF',CF,'IMTauto',[],'IMTVauto',[],'IMTgt',[],'IMTVgt',[],...
                                'BiasIMT',[],'BiasIMTV',[],'HM_LI',[],'HM_MA',[],'PDM_LI',[],'PDM_MA',[],'commonLenPxl',[],'commonLenMM'...
                                ,'mPointBiasLI',[],'sPointBiasLI',[],'pPointBiasLI',[],'mPointBiasMA',[],'sPointBiasMA',[],'pPointBiasMA',[]);
                            tempBIASstatsInterp(q) = struct('filename',filename,'CF',CF,...
                                'PDMauto',[],'PDMgt',[],'BiasPDM',[],'AbsBiasPDM',[],'EDauto',[],'EDgt',[],'BiasED',[],'AbsBiasED',[]);
                        end
        
                        clear comm
        
%                         % common support on original profiles
%                         [LI,MA,~] = km_CommonSupport(LI,MA);
%                         [LI_GT,MA_GT,~] = km_CommonSupport(LI_GT,MA_GT);
%         
%                         comm.temp_LI = LI(ismember(round(LI(:,1)),round(LI_GT(:,1))),:);   comm.temp_MA = MA(ismember(round(MA(:,1)),round(MA_GT(:,1))),:);
%                         comm.gt_LI = LI_GT(ismember(round(LI_GT(:,1)),round(LI(:,1))),:);  comm.gt_MA = MA_GT(ismember(round(MA_GT(:,1)),round(MA(:,1))),:);
%         
%                         if length(comm.temp_LI)>2
%                             tempBIASstats(q) = LI_MA_stats_not_interp(comm.gt_LI,comm.gt_MA,comm.temp_LI,comm.temp_MA,CF,filename,mode);
%                         else
%                             tempBIASstats(q) = struct('filename',filename,'CF',CF,...
%                                 'PDMauto',[],'PDMgt',[],'BiasPDM',[],'AbsBiasPDM',[],'EDauto',[],'EDgt',[],'BiasED',[],'AbsBiasED',[]);
%                         end
                        
                    else
        
%                         % common support on original profiles
%                         [LI,MA,~] = km_CommonSupport(LI,MA);
%                         [LI_GT,MA_GT,~] = km_CommonSupport(LI_GT,MA_GT);
%         
%                         comm.temp_LI = LI(ismember(round(LI(:,1)),round(LI_GT(:,1))),:);   comm.temp_MA = MA(ismember(round(MA(:,1)),round(MA_GT(:,1))),:);
%                         comm.gt_LI = LI_GT(ismember(round(LI_GT(:,1)),round(LI(:,1))),:);  comm.gt_MA = MA_GT(ismember(round(MA_GT(:,1)),round(MA(:,1))),:);
%         
%                         if length(comm.temp_LI)>2
%                             tempBIASstats(q) = LI_MA_stats_not_interp(LI_GT,MA_GT,LI,MA,CF,filename,mode);
%                         end

                        tempIMTstats(q) = struct('filename',filename,'CF',CF,'IMTauto',[],'IMTVauto',[],'IMTgt',[],'IMTVgt',[],...
                            'BiasIMT',[],'BiasIMTV',[],'HM_LI',[],'HM_MA',[],'PDM_LI',[],'PDM_MA',[],'commonLenPxl',[],'commonLenMM',[]...
                            ,'mPointBiasLI',[],'sPointBiasLI',[],'pPointBiasLI',[],'mPointBiasMA',[],'sPointBiasMA',[],'pPointBiasMA',[]);
                        tempBIASstatsInterp(q) = struct('filename',filename,'CF',CF,...
                            'PDMauto',[],'PDMgt',[],'BiasPDM',[],'AbsBiasPDM',[],'EDauto',[],'EDgt',[],'BiasED',[],'AbsBiasED',[]);

                    end % end contains original for profile comparison choice
                catch
                    fprintf('\tNo profile in file %s',fxs(jj).name)
                    tempIMTstats(q) = struct('filename',filename,'CF',CF,'IMTauto',[],'IMTVauto',[],'IMTgt',[],'IMTVgt',[],...
                        'BiasIMT',[],'BiasIMTV',[],'HM_LI',[],'HM_MA',[],'PDM_LI',[],'PDM_MA',[],'commonLenPxl',[],'commonLenMM',[]...
                        ,'mPointBiasLI',[],'sPointBiasLI',[],'pPointBiasLI',[],'mPointBiasMA',[],'sPointBiasMA',[],'pPointBiasMA',[]);

                    tempBIASstatsInterp(q) = struct('filename',filename,'CF',CF,...
                        'PDMauto',[],'PDMgt',[],'BiasPDM',[],'AbsBiasPDM',[],'EDauto',[],'EDgt',[],'BiasED',[],'AbsBiasED',[]);
                    tempBIASstats(q) = struct('filename',filename,'CF',CF,...
                        'PDMauto',[],'PDMgt',[],'BiasPDM',[],'AbsBiasPDM',[],'EDauto',[],'EDgt',[],'BiasED',[],'AbsBiasED',[]);

                    no_auto_pf = no_auto_pf + 1;

                end
            catch
                fprintf('\tNo GT profile in file %s',fxs(jj).name)

                tempIMTstats(q) = struct('filename',filename,'CF',CF,'IMTauto',[],'IMTVauto',[],'IMTgt',[],'IMTVgt',[],...
                    'BiasIMT',[],'BiasIMTV',[],'HM_LI',[],'HM_MA',[],'PDM_LI',[],'PDM_MA',[],'commonLenPxl',[],'commonLenMM',[]...
                    ,'mPointBiasLI',[],'sPointBiasLI',[],'pPointBiasLI',[],'mPointBiasMA',[],'sPointBiasMA',[],'pPointBiasMA',[]);
                tempBIASstatsInterp(q) = struct('filename',filename,'CF',CF,...
                    'PDMauto',[],'PDMgt',[],'BiasPDM',[],'AbsBiasPDM',[],'EDauto',[],'EDgt',[],'BiasED',[],'AbsBiasED',[]);
                tempBIASstats(q) = struct('filename',filename,'CF',CF,...
                    'PDMauto',[],'PDMgt',[],'BiasPDM',[],'AbsBiasPDM',[],'EDauto',[],'EDgt',[],'BiasED',[],'AbsBiasED',[]);   

                no_gt_pf = no_gt_pf + 1;

            end

           q = q + 1;
        end % end loop across GT files

        if exist("tempIMTstats","var")
            IMTstats.(folderName) = tempIMTstats;
            BIASstatsInterp.(folderName) = tempBIASstatsInterp;
        end
        
        BIASstats.(folderName) = tempBIASstats;

        clear tempIMTstats tempBIASstatsInterp tempBIASstats

    end % end if ~contain GT1

end % end loop across folders

%% performance summary

comparisons = fieldnames(IMTstats);

tableSummaryIMT = table();
tableSummaryBIASinterp = table();
tableSummaryBIAS = table();

for ii = 1 : length(comparisons)

    if ~isempty(IMTstats.(comparisons{ii}))
        tempSummary =  get_dataset_summary(IMTstats.(comparisons{ii}));
    
        tempSummary = struct2table(tempSummary);
        tempSummary = tempSummary(:,3:end-1);
        tempSummary = [cell2table(comparisons(ii)),tempSummary];

        tableSummaryIMT = [tableSummaryIMT;tempSummary];
        
        clear tempSummary
    end

    if ~isempty(BIASstatsInterp.(comparisons{ii}))
        
        tempSummary =  get_dataset_summary(BIASstatsInterp.(comparisons{ii}));
        
        tempSummary = struct2table(tempSummary);
        tempSummary = tempSummary(:,3:end-1);
        tempSummary = [cell2table(comparisons(ii)),tempSummary];

        tableSummaryBIASinterp = [tableSummaryBIASinterp;tempSummary];

        clear tempSummary
    end

%     if ~isempty(BIASstats.(comparisons{ii}))
%     
%         tempSummary =  get_dataset_summary(BIASstats.(comparisons{ii}));
%     
%         tempSummary = struct2table(tempSummary);
%         tempSummary = tempSummary(:,3:end-1);
%         tempSummary = [cell2table(comparisons(ii)),tempSummary];
% 
%         tableSummaryBIAS = [tableSummaryBIAS;tempSummary];
% 
%         clear tempSummary
%     end

end

save('/media/francesco/DEV001/PROJECT-CUBS-SPIE/RESULTS/MATLAB/profileEvaluationConsensus-INTER-biasProfiles-CommonSupport.mat', ...
    "IMTstats","BIASstatsInterp","tableSummaryIMT","tableSummaryBIASinterp");

