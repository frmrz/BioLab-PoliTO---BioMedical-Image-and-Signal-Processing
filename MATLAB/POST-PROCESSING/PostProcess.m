close all
clear 

%% Post-processing ensemble network
% Takes as input the folder with N runs and perform weighted ensemble of
% predictions

%% 1. Selecting data

mainFolder =  '/media/francesco/DEV001/PROJECT-GITHUB/Image_Processing/Segmentation/Ultrasound-OpticNerveSheaths';

dataFolder = fullfile(mainFolder,'DATA');
resultsFolder = fullfile(mainFolder,'RESULTS','CHECKPOINTS');

masksFolder = fullfile(dataFolder,'LABELS_256');
imageFolder = fullfile(dataFolder,'IMAGES_256');
outputFolder = fullfile(resultsFolder,'PREDICTIONS_new');

%% 2. Get runs amd loop

runsFodler = dir(mainFolder); % get all folders with N runs for M folds
DiceMat = zeros(464,9);
HausMat = zeros(464,9);
lenFolds = 0;

for Folds = 0:4

    tempRun = fullfile(resultsFolder,sprintf('Unet-resnet50-f%s-r%s',num2str(Folds),num2str(0)),'Test_output'); 

    filenames = dir(tempRun);

%     figure(),
    for ii = 3 : length(filenames)
        
%         close all
        fprintf('Fold %s\n', num2str(Folds));
        EnsembleMatrix = uint8(zeros(256,256,8));

        GT = imread(fullfile(masksFolder,filenames(ii).name)) == 1;
        I = imread(fullfile(imageFolder,filenames(ii).name));

        runWeights = zeros(1,8);
        q = 1;
%         figure,
        for Run = 0 : 7

            tempRun = fullfile(resultsFolder,sprintf('Unet-resnet50-f%s-r%s',num2str(Folds),num2str(Run)),'Test_output');
            tempRunOutputs = dir(tempRun);
            tempOut= imread(fullfile(tempRunOutputs(ii).folder,tempRunOutputs(ii).name));
            
            %% light morph op on OU
%             figure(1);
%             subplot(221),imshow(tempOut)
%             imwrite(tempOut,sprintf('pre%d.png',Run));

            OUT = imopen(tempOut,strel('disk',3));
%             subplot(222),imshow(OUT)
            OUT1 = imclose(OUT,strel('disk',5));
%             subplot(223),imshow(OUT1)
            
            if range(OUT1) > 1; OUT1 = OUT1>0.5; else; OUT1 = OUT1>127.5; end

            stats=regionprops(logical(OUT1),'Area','Eccentricity','Solidity','PixelIdxList');

            %% cut very small regions

            if length(stats) > 2
                stats = table2struct(sortrows(struct2table(stats), 'Area','descend'));
                meanTwoTop = mean([stats(1).Area,stats(2).Area]);

                if stats(3).Area < 0.25*meanTwoTop
                    for ww = 3 : length(stats)
                        OUT1(stats(ww).PixelIdxList)=0;
                    end
                end
            end

%             subplot(224),imshow(OUT1)
%             imwrite(OUT1,sprintf('post%d.png',Run));

            EnsembleMatrix(:,:,Run+1) = uint8(OUT1)*255;

            %% Dice score
            DiceMat(ii-2 + lenFolds, Run+1) = dice(GT,OUT1);

            %% Hausdorff on OUT2
            gtLabel = uint8(bwlabel(GT));
            outLabel = uint8(bwlabel(OUT1));
            commonRows = ~((sum(outLabel==1,2) > 1 & sum(outLabel==2,2) > 1) & (sum(gtLabel==1,2) > 1 & sum(gtLabel==2,2) > 1)); %find common support columns
    
            if sum(commonRows) < 256
                gtLabel(commonRows,:)=0; % zero out rows
                outLabel(commonRows,:)=0;
        
                HausMat(ii-2 + lenFolds,Run+1)=max(imhausdorff(gtLabel,outLabel,'euclidean'));
     
            else
                HausMat(ii-2 + lenFolds,Run+1) = nan;
            end

        end

        %% create ensemble

        meanEnsemble = uint8(mean(EnsembleMatrix,3));
        meanEnsemblePick = uint8(mean(EnsembleMatrix(:,:,:),3));
        
%         subplot(221),imshow(meanEnsemble);
%         imwrite(meanEnsemble,'meanEns.png');

        OUT = imopen(meanEnsemble,strel('disk',5));
%         subplot(222),imshow(OUT);

        OUT = imclose(OUT,strel('disk',4));
%         subplot(223),imshow(OUT);

        %% post process ensemble

        OUT2 = OUT>0.6*255;
        stats=regionprops(logical(OUT2),'Area','Eccentricity','Solidity','PixelIdxList');

        if length(stats) < 2

            OUT = imopen(meanEnsemble,strel('disk',3));
            OUT = imclose(OUT,strel('disk',5));
            OUT2 = OUT>0.5*255;
            OUT2 = imopen(OUT2,strel('disk',1));

        end

        stats=regionprops(logical(OUT2),'Area','Eccentricity','Solidity','PixelIdxList');

        if length(stats) > 2
            stats = table2struct(sortrows(struct2table(stats), 'Area','descend'));
            meanTwoTop = mean([stats(1).Area,stats(2).Area]);
            meanThreeTop = mean([stats(1).Area,stats(2).Area,stats(3).Area]);
            medThreeTop = median([stats(1).Area,stats(2).Area,stats(3).Area]);

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
        end

        OUT3 = OUT>0.6*255;

        %% dice scores
        DiceMat(ii-2 + lenFolds, 9) = dice(GT,OUT2);

        %% Hausdorff on OUT2
        gtLabel = uint8(bwlabel(GT));
        outLabel = uint8(bwlabel(OUT2));

        commonRows = ~((sum(outLabel==1,2) > 1 & sum(outLabel==2,2) > 1) & (sum(gtLabel==1,2) > 1 & sum(gtLabel==2,2) > 1)); %find common support columns

        if sum(commonRows) < 256
            gtLabel(commonRows,:)=0; % zero out rows
            outLabel(commonRows,:)=0;
    
            HausMat(ii-2 + lenFolds,9)=max(imhausdorff(gtLabel,outLabel,'euclidean'));
        else
            HausMat(ii-2 + lenFolds,9) = nan;
        end
        
        %% debug image
%         debugImage = cat(2,I,uint8(GT)*255);    
%         debugImage1 = cat(2,OUT,uint8(OUT2)*255);    
%         debugImage2 = cat(1,debugImage,debugImage1);    

        %% write images
%         imwrite(debugImage2,fullfile(outputFolder,filenames(ii).name));
%         imwrite(OUT2,fullfile(outputFinalFolder,filenames(ii).name))

%         subplot(224),imshow(OUT2);
%                 imwrite(OUT2,'outEns.png');
% 
%         save_image_jpg_eps(filenames(ii).name,outputFolder,outputFolder,0)

    end
    lenFolds = lenFolds + length(filenames) - 2;
end

%% plotting runs comparison

DiceVect = reshape(DiceMat,[],1);
HausVect = reshape(HausMat,[],1);

RunsGroups = reshape(repmat(1:9,464,1),[],1);

binEdges = 0.5:1:9.5;
bins = {'Run 1','Run 2','Run 3','Run 4','Run 5','Run 6','Run 7','Run 8','Ensemble'};
groupMach = discretize(RunsGroups,binEdges,'categorical',bins);

% diceTab = array2table([DiceVect,RunsGroups]);
% writetable(diceTab,'DiceRuns.xlsx');
% 
% HausTab = array2table([HausVect,RunsGroups]);
% writetable(HausTab,'HausRuns.xlsx');

% ll='';

h = figure; h.Position = [100 100 1000 800];

subplot(121),boxchart(groupMach,DiceVect),ylabel('Dice score'),
set(gca,'FontName','Calibri','FontSize',14, 'XColor','k','YColor','k');

grid on
ax = gca; % Get handle to current axes.
ax.XColor = 'k'; % Red
ax.YColor = 'k'; % Blue
ax.GridLineStyle = "-";
ax.GridAlpha = 0.1;  % Make grid lines less transparent.
ax.GridColor = 'k'; % Dark Green.

subplot(122),boxchart(groupMach,HausVect),ylabel('Hausdorff distance (pixels)'),set(gca,'FontSize',12);
set(gca,'FontName','Calibri','FontSize',14, 'XColor','k','YColor','k');

grid on
ax = gca; % Get handle to current axes.
ax.XColor = 'k'; % Red
ax.YColor = 'k'; % Blue
ax.GridLineStyle = "-";
ax.GridAlpha = 0.1;  % Make grid lines less transparent.
ax.GridColor = 'k'; % Dark Green.


% saveas(h,fullfile(ll,'DiceRunsEPS.eps'),'epsc')
% save_image_jpg_eps('DiceRunsEPSv1.eps',ll,ll)
% sgtitle('Testing on full dataset (by runs)')

%% plotting machine comparison

lenMAch = [207 113 93 51];
MachineGroups = [zeros(464,1);ones(lenMAch(1),1);ones(lenMAch(2),1)*2;ones(lenMAch(3),1)*3;ones(lenMAch(4),1)*4];

binEdges = -0.5:1:4.5;
bins = {'Total','Mach 1','Mach 2','Mach 3','Mach 4'};
groupMach = discretize(MachineGroups,binEdges,'categorical',bins);

h = figure; h.Position = [100 100 900 800];

subplot(121),boxchart(groupMach,[DiceMat(:,9);DiceMat(:,9)]),ylabel('Dice score');
set(gca,'FontName','Calibri','FontSize',14, 'XColor','k','YColor','k');

grid on
ax = gca; % Get handle to current axes.
ax.XColor = 'k'; % Red
ax.YColor = 'k'; % Blue
ax.GridLineStyle = "-";
ax.GridAlpha = 0.1;  % Make grid lines less transparent.
ax.GridColor = 'k'; % Dark Green.
subplot(122),boxchart(groupMach,[HausMat(:,9);HausMat(:,9)]),ylabel('Hausdorff distance (pixels)');
set(gca,'FontName','Calibri','FontSize',14, 'XColor','k','YColor','k');

grid on
ax = gca; % Get handle to current axes.
ax.XColor = 'k'; % Red
ax.YColor = 'k'; % Blue
ax.GridLineStyle = "-";
ax.GridAlpha = 0.1;  % Make grid lines less transparent.
ax.GridColor = 'k'; % Dark Green.

% ll='';
% 
% saveas(h,fullfile(ll,'MachRunsEPS.eps'),'epsc')
% sgtitle('Testing on full dataset (by machines)')

%% statistical comparison

DiceNo = sum(DiceMat,2);
DiceMatNum = DiceMat(~isinf(DiceNo) & ~isnan(DiceNo),:);

HausNo = sum(HausMat,2);
HausMatNum = HausMat(~isinf(HausNo) & ~isnan(HausNo),:);

for i = 1 : size(DiceMat,2)
    for j = 1 : size(DiceMat,2)
        tempAd = DiceMat(:,i); tempAh = HausMat(:,i);
        tempBd = DiceMat(:,j); tempBh = HausMat(:,j);

        p.Dice(i,j) = signrank(tempAd,tempBd);
        p.Haus(i,j) = signrank(tempAh,tempBh);        

    end
end

for i = 1 : size(HausMatNum,2)
    for j = 1 : size(HausMatNum,2)
        tempAd = HausMatNum(:,i);
        tempBd = HausMatNum(:,j);

        p.Haus(i,j) = signrank(tempAd,tempBd);

    end
end

%% paper results

% Results : dice score
DiceAvg = mean(DiceMat,1);

min(DiceAvg(1:8))
max(DiceAvg(1:8))
DiceAvg(9)

DiceStd = std(DiceMat,1);

min(DiceStd(1:8))
max(DiceStd(1:8))
DiceStd(9)
DiceAvg = mean(DiceMat,1);

% Results : hausdorff
HausAvg = mean(HausMatNum,1,"omitnan");

min(HausAvg(1:8))
max(HausAvg(1:8))
HausAvg(9)

HausStd = std(HausMatNum,1,"omitnan");

min(HausStd(1:8))
max(HausStd(1:8))
HausStd(9)