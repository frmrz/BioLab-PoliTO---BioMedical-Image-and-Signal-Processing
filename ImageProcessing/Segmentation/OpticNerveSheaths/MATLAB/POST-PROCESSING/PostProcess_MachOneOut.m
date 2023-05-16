close all
clear 

%% Post-processing ensemble network
% Takes as input the folder with N runs and perform weighted ensemble of
% predictions

%% 1. Selecting data

mainFolder =  '/media/francesco/DEV001/PROJECT-GITHUB/Image_Processing/Segmentation/Ultrasound-OpticNerveSheaths';

dataFolder = fullfile(mainFolder,'DATA');
masksFolder = fullfile(dataFolder,'LABELS_256');
imageFolder = fullfile(dataFolder,'IMAGES_256');

resultsFolder = fullfile(mainFolder,'RESULTS','CHECKPOINTS_MACH_OUT');

outputFolder = fullfile(resultsFolder,'PREDICTIONS_MACH_OUT');
outputFinalFolder = '/home/francesco/PROJECT-OPTIC-NERVE/CODE/Segmentation_models/OUTPUT-MachOut/';

% mkdir(outputFolder)
%% 2. Get runs amd loop

runsFodler = dir(mainFolder); % get all folders with N runs for M folds
DiceMat = zeros(464,8);
HausMat = zeros(464,2);
lenMAch = [207 113 93 51];

for Mach = 1 : 4

    tempRun = fullfile(resultsFolder,sprintf('Unet-resnet50-f%s-mach%s',num2str(0),num2str(Mach)),'Test_output');
    filenames = dir(tempRun);

%     figure(),
    for ii = 3 : length(filenames)
%         close all
        fprintf('Machine %s\n', num2str(Mach));

        EnsembleMatrix = uint8(zeros(256,256,5)); % number of folds or repetitions

        GT = imread(fullfile(masksFolder,filenames(ii).name)) == 1;
        I = imread(fullfile(imageFolder,filenames(ii).name));

        q = 1;

%         figure,
        for Fold = 0 : 4

            tempRun = fullfile(resultsFolder,sprintf('Unet-resnet50-f%s-mach%s',num2str(Fold),num2str(Mach)),'Test_output');

            tempRunOutputs = dir(tempRun);
            tempOut= imread(fullfile(tempRunOutputs(ii).folder,tempRunOutputs(ii).name));
            
            %% light morph op on OU
%             figure(1);
%             subplot(221),imshow(tempOut)
            OUT = imopen(tempOut,strel('disk',3));
%             subplot(222),imshow(OUT)
            OUT1 = imclose(OUT,strel('disk',5));
%             subplot(223),imshow(OUT1)
            
            if range(OUT1) > 1; OUT1 = OUT1>0.5; else; OUT1 = OUT1>127.5; end

%             subplot(224),imshow(OUT1)

            DiceMat(ii-2 + sum(lenMAch(1:Mach-1)), Fold+1) = dice(GT,OUT1);
                
            stats=regionprops(logical(OUT1),'Area','Eccentricity','Solidity','PixelIdxList');
%             stats=regionprops(logical(OUT),'all');

            %% cut very small regions

            if length(stats) > 2
                stats = table2struct(sortrows(struct2table(stats), 'Area','descend'));
                meanTwoTop = mean([stats(1).Area,stats(2).Area]);

                if stats(3).Area < 0.25*meanTwoTop
                    for ww = 3 : length(stats)
                        tempOut(stats(ww).PixelIdxList)=0;
                    end
%                     OUT1 = bwareaopen(OUT1,stats(3).Area+1);
                end
            end

%             for kk = 1 : length(stats)
%                 statsRunsTemp = stats(kk);
%                 statsRunsTemp.dice = DiceMat(ii-2 + sum(lenFolds(1:Folds + 1)), Run+1);
%                 statsRunsTemp.run = Run;
% 
%                 statsRuns(q) = statsRunsTemp;
%                 q = q + 1;
%                 subplot(2,4,Run+1), imshow(OUT1);
%             end

            EnsembleMatrix(:,:,Fold+1) = uint8(OUT1)*255;
%             EnsembleMatrix(:,:,Run+1) = uint8(OUT1)*255;
%             subplot(4,3,Run+1),imshow(EnsembleMatrix(:,:,Run+1));
        end

%         mStatsRun = mean(table2array(struct2table(statsRuns)));

%         subplot(2,4,8), imshow(GT);
%         clear statRuns
        close all
        %% create ensemble
        meanEnsemble = uint8(mean(EnsembleMatrix,3));
        meanEnsemblePick = uint8(mean(EnsembleMatrix(:,:,:),3));
        
%         subplot(221),imshow(meanEnsemble);

        OUT = imopen(meanEnsemble,strel('disk',5));
%         subplot(222),imshow(OUT);

        OUT = imclose(OUT,strel('disk',4));
%         subplot(223),imshow(OUT);


        %% post process ensemble
        OUT1 = OUT>0.4*255;

        OUT2 = OUT>0.6*255;
        stats=regionprops(logical(OUT2),'Area','Eccentricity','Solidity','PixelIdxList');

        if length(stats) < 2
%                         imshow(OUT2)

            OUT = imopen(meanEnsemble,strel('disk',3));
            OUT = imclose(OUT,strel('disk',5));
            OUT2 = OUT>0.5*255;
            OUT2 = imopen(OUT2,strel('disk',1));

%                         imshow(OUT2)
        a=0;
        end
        stats=regionprops(logical(OUT2),'Area','Eccentricity','Solidity','PixelIdxList');

        if length(stats) > 2
            stats = table2struct(sortrows(struct2table(stats), 'Area','descend'));
            meanTwoTop = mean([stats(1).Area,stats(2).Area]);
            meanThreeTop = mean([stats(1).Area,stats(2).Area,stats(3).Area]);
            medThreeTop = median([stats(1).Area,stats(2).Area,stats(3).Area]);

%             imshow(OUT2)
            if stats(3).Area < 0.5*meanTwoTop && stats(3).Area < 0.5*medThreeTop
                for ww = 3 : length(stats)
                    OUT2(stats(ww).PixelIdxList)=0;
                end
            else
                BWmerge = false(size(OUT2));
                BWmerge(stats(2).PixelIdxList) = true;
                BWmerge(stats(3).PixelIdxList) = true;
                BWmerge = bwconvhull(BWmerge);
                OUT2(BWmerge==1)=1;
            end
%             imshow(OUT2)
            a=0;
        end

        OUT3 = OUT>0.6*255;
        
        %% dice scores
        DiceMat(ii-2 + sum(lenMAch(1:Mach-1)), 6) = dice(GT,OUT1);
        DiceMat(ii-2 + sum(lenMAch(1:Mach-1)), 7) = dice(GT,OUT2);
        DiceMat(ii-2 + sum(lenMAch(1:Mach-1)), 8) = dice(GT,OUT3);

        %% Hausdorff on OUT2
        gtLabel = uint8(bwlabel(GT));
        outLabel = uint8(bwlabel(OUT2));

        commonColumns = ~((sum(outLabel==1,2) > 1 & sum(outLabel==2,2) > 1) & (sum(gtLabel==1,2) > 1 & sum(gtLabel==2,2) > 1)); %find common support columns

        if sum(commonColumns) < 256
        gtLabel(commonColumns,:)=0; % zero out rows
        outLabel(commonColumns,:)=0;
        
%         subplot(121),imshow(label2rgb(gtLabel));
%         subplot(122),imshow(label2rgb(outLabel));

        HausMat(ii-2 + sum(lenMAch(1:Mach-1)),:)=imhausdorff(gtLabel,outLabel,'euclidean');
 
        else
            HausMat(ii-2 + sum(lenMAch(1:Mach-1)),:) = nan;
        end

% 
% %         OUT = imopen(meanEnsemblePick,strel('disk',3));
% %         OUT = imclose(OUT,strel('disk',5));
% % 
% %         OUT1 = OUT>0.4*255;
% %         OUT2 = OUT>0.5*255;
% %         OUT3 = OUT>0.6*255;
% % 
% %         DiceMat(ii-2 + sum(lenFolds(1:Folds + 1)), 12) = dice(GT,OUT1);
% %         DiceMat(ii-2 + sum(lenFolds(1:Folds + 1)), 13) = dice(GT,OUT2);
% %         DiceMat(ii-2 + sum(lenFolds(1:Folds + 1)), 14) = dice(GT,OUT3);
        debugImage = cat(2,I,uint8(GT)*255);    
        debugImage1 = cat(2,OUT,uint8(OUT2)*255);    
        debugImage2 = cat(1,debugImage,debugImage1);    
        
        %% open sheets to avoid bias
%         OUT2 = imopen(OUT2,strel('disk',5));
%         imshow(OUT2)
% 
%         imwrite(debugImage2,fullfile(outputFolder,filenames(ii).name));
%         imwrite(OUT2,fullfile(outputFinalFolder,filenames(ii).name))

%         subplot(4,3,9),imshow(meanEnsemble);
%         subplot(4,3,10),imshow(OUT2);
%         subplot(4,3,11),imshow(GT);
    end
end

olenMAch = [207 113 93 51];
MachineGroups = [zeros(464,1);ones(lenMAch(1),1);ones(lenMAch(2),1)*2;ones(lenMAch(3),1)*3;ones(lenMAch(4),1)*4];

binEdges = -0.5:1:4.5;
bins = {'Total','Mach 1','Mach 2','Mach 3','Mach 4'};
groupMach = discretize(MachineGroups,binEdges,'categorical',bins);

figure, 
boxchart(DiceMat)

figure, 
subplot(121),boxchart(groupMach,[DiceMat(:,7);DiceMat(:,7)]),ylabel('Dice score');
subplot(122),boxchart(groupMach,[max(HausMat,[],2);max(HausMat,[],2)]),ylabel('Hausdorff distance (pixels)');
sgtitle('Testing on left out machine')

%% paper results

dScores = DiceMat(:,7);
d1Scores = dScores(1:207,:);
d2Scores = dScores(208:320,:);
d3Scores = dScores(321:413,:);
d4Scores = dScores(414:end,:);

avgDiceMach1 = mean(d1Scores) 
stdDiceMach1 = std(d1Scores)

avgDiceMach2 = mean(d2Scores) 
stdDiceMach2 = std(d2Scores)

avgDiceMach3 = mean(d3Scores) 
stdDiceMach3 = std(d3Scores)

avgDiceMach4 = mean(d4Scores) 
stdDiceMach4 = std(d4Scores)
