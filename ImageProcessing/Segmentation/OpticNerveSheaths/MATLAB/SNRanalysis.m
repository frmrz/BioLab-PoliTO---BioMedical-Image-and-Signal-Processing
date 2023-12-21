% FILEPATH: /media/francesco/DEV001/CENTRAL_PROJECTS_REPOSITORY/BioLab-PoliTO---BioMedical-Image-and-Signal-Processing/ImageProcessing/Segmentation/OpticNerveSheaths/MATLAB/SNRanalysis.m

% This script calculates the Signal-to-Noise Ratio (SNR) in ultrasound (US) images of the optic nerve. It also calculates the Dice scores and Hausdorff distances for the segmented optic nerve regions. The script loops through a set of images and performs the following steps:
% 1. Load the images, labels, and background images.
% 2. Extract the regions of interest (sheets, optic nerve, and background) from the images.
% 3. Calculate the SNR parameters for the left and right regions.
% 4. Calculate the Dice scores and Hausdorff distances for the segmented optic nerve regions.
% 5. Discretize the SNR values and analyze the differences in Dice scores and Hausdorff distances.

% The script uses various functions such as imread, bwconvhull, regionprops, extractLBPFeatures, dice, imhausdorff, and kruskalwallis.

% Note: Some code sections are commented out and can be uncommented to visualize intermediate results or save figures.

% Author: Francesco Marzola

close all
clear

project_root = '/media/francesco/DEV001/PROJECT-GITHUB/Image_Processing/Segmentation/Ultrasound-OpticNerveSheaths';

%% script for the calculation of SNR in US images of optic nerve
mainFolder =  fullfile(project_root, 'DATA');
imgFolder  = fullfile(mainFolder,'IMAGES_256');
labFolder = fullfile(mainFolder,'LABELS_256');
backFolder = fullfile(mainFolder,'BACKGROUND_256');
cfFolder   = fullfile(mainFolder,'CF');

predictionsFolder = fullfile(project_root,'RESULTS','PREDICTIONS');

%% loop throught images and calc SNR


files = dir(imgFolder);
% figure,

for i = 3:length(files)

    img = imread(fullfile(imgFolder,files(i).name));
    lab = imread(fullfile(labFolder,files(i).name));
    bck = imread(fullfile(backFolder,files(i).name));
    cf(i-2) = load(fullfile(cfFolder,strrep(files(i).name,'.png','.txt')));

    bwS = true(size(img));  bwS(lab~=1)=0; 
    bwON = true(size(img)); bwON(lab~=2)=0;
    bwbk = true(size(img)); bwbk(bck~=3 | lab==1)=0;

    S = img;    S(lab~=1)=0;
    ON = img;   ON(lab~=2)=0;
    bk = img;   bk(bck~=3 | lab==1)=0;
 
    Sb = bwconvhull(bwS>0); 
    st = regionprops(Sb);

    % Create a RECT where to evaluate image
    leftLim = round(max(1,st.BoundingBox(1)-st.BoundingBox(3)*0.5));
    rightLim = round(min(256,st.BoundingBox(1)+st.BoundingBox(3)*1.5));

    % Create background reference
    bwbk(sum(bwS,2)==0,:)=0;            % set rows outside the sheets to zero
    bwbk(:,[1:leftLim,rightLim:end])=0; % set cols outside limits to zero
    bk(bwbk==0)=0;                      % set background pixels to zero
    
%     figure,
%     subplot(131),imshow(S);
%     subplot(132),imshow(ON);
%     subplot(133),imshow(bk);

    %% look at sheets plus background
    shbk = img;
    shbk(S==0 & bk==0) = 0;
    

% %     subplot(224),imshow(bwgxd);
% 
%     %% LBP features
%     ImgUint = shbk;
%     features1 = extractLBPFeatures(ImgUint);
% 
%     Sb1 = bwconvhull(shbk>0); 
%     st1 = regionprops(Sb1);
% 
%     st1.BoundingBox(1,3)=st1.BoundingBox(1,3)/2;
% 
%     ImgCropLeft = imcrop(shbk,st1.BoundingBox);
%     ImgCropLeft = double(ImgCropLeft);
%     ImgCropLeft(ImgCropLeft==0)=nan;
% 
%     featuresLeft = extractLBPFeatures(ImgCropLeft,"NumNeighbors",8,"Radius",2,"Upright",false);
% 
%     imshow(uint8(ImgCropLeft))


    %% left 
    leftS = bwS; leftS(:,round(st.Centroid(1):end))=0;    % set right part of the image to zero
    leftON = bwON;
    leftBK = bwbk; leftBK(:,round(st.Centroid(1):end))=0; % set right part of the background to zero

%     subplot(131),imshow(leftS);
%     subplot(132),imshow(leftON);
%     subplot(133),imshow(leftBK);

    ms = single(S(leftS~=0));
    mnIN = single(ON(leftON~=0));
    mnOUT = single(bk(leftBK~=0));
    
    for j = 1 : 10
        temp(j) = calcSNRparams(ms,mnIN,mnOUT);
    end

    temp = mean(table2array(struct2table(temp)));
    temp = array2table(temp(:)');
    temp.Properties.VariableNames= {'SNRin','SNRout','cIN','cOUT'};
    temp = [cell2table(cellstr(files(i).name)),temp];
    temp.Properties.VariableNames= {'name','SNRin','SNRout','cIN','cOUT'};
    temp = table2struct(temp);

    summaryLeft(i-2) = temp;

    clear temp

    %% right 
    rightS = bwS; rightS(:,round(1:st.Centroid(1)))=0; % set left part of the image to zero
    rightON = bwON;
    rightBK = bwbk; rightBK(:,round(1:st.Centroid(1)))=0; % set left part of the background to zero

    ms = single(S(rightS~=0));
    mnIN = single(ON(rightON~=0));
    mnOUT = single(bk(rightBK~=0));
    
    for j = 1 : 20
        temp(j) = calcSNRparams(ms,mnIN,mnOUT);
    end

    temp = mean(table2array(struct2table(temp)));
    temp = array2table(temp(:)');
    temp.Properties.VariableNames= {'SNRin','SNRout','cIN','cOUT'};
    temp = [cell2table(cellstr(files(i).name)),temp];
    temp.Properties.VariableNames= {'name','SNRin','SNRout','cIN','cOUT'};
    temp = table2struct(temp);

    summaryRight(i-2) = temp;

%     subplot(131),imshow(leftS);
%     subplot(132),imshow(leftON);
%     subplot(133),imshow(leftBK);
    clear temp

    %% Draw SNR areas

%     h=figure, imshow(img)
%     hold on
%     visboundaries(bwON,'Color','r','LineWidth',0.5);
%     visboundaries(leftS,'Color','g','LineWidth',0.5);
%     visboundaries(rightS,'Color','g','LineWidth',0.5);
%     save_image_jpg_eps('SNRarea_ONS2ON',outputFolder,outputFolder,0)
% 
%     h=figure, imshow(img)
%     hold on
%     visboundaries(leftS,'Color','g','LineWidth',0.5);
%     visboundaries(rightS,'Color','g','LineWidth',0.5);
%     visboundaries(imopen(rightBK,strel('disk',3)),'Color','r','LineWidth',0.5);
%     visboundaries(imopen(leftBK,strel('disk',3)),'Color','r','LineWidth',0.5);
%     save_image_jpg_eps('SNRarea_ONS2BACK',outputFolder,outputFolder,0)


end

%% get Dice scores

files = dir(predictionsFolder);

for i = 3:length(files)
    out = imread(fullfile(predictionsFolder,files(i).name));
    gt = imread(fullfile(labFolder,files(i).name));
    gt(gt~=1)=0;
    gt=logical(gt);

    Sb = bwconvhull(gt>0); 
    st = regionprops(Sb);

    gtLeft=gt; gtLeft(:,round(st.Centroid(1):end))=0;
    gtRight=gt; gtRight(:,round(1:st.Centroid(1)))=0;

    outLeft=out; outLeft(:,round(st.Centroid(1):end))=0;
    outRight=out; outRight(:,round(1:st.Centroid(1)))=0;

    summaryLeft(i-2).DiceLeft = dice(outLeft,gtLeft);
    summaryRight(i-2).DiceRight = dice(outRight,gtRight);

    diceFromSource(i-2) = dice(out,gt);

    summaryLeft(i-2).HaussLeft = imhausdorff(outLeft,gtLeft,'euclidean');
    summaryRight(i-2).HaussRight = imhausdorff(outRight,gtRight,'euclidean');

    summaryLeft(i-2).cf = cf(i-2);
    summaryRight(i-2).cf = cf(i-2);
end

diceAndSnrMatLeft = [[summaryLeft(:).SNRin]',...
                 [summaryLeft(:).SNRout]',...
                 mean([abs([summaryLeft(:).SNRin]);abs([summaryLeft(:).SNRout])])',...
                 [summaryLeft(:).cIN]',...
                 [summaryLeft(:).cOUT]',...
                 [summaryLeft(:).DiceLeft]',...
                 [summaryLeft(:).HaussLeft]'.*[summaryLeft(:).cf]'];

corrSNRdiceLeft = array2table(corr(diceAndSnrMatLeft));
corrSNRdiceLeft.Properties.VariableNames= {'SNRin','SNRout','SNRmean','cIN','cOUT','Dice','Hauss'};

diceAndSnrMatRight = [[summaryRight(:).SNRin]',...
                 [summaryRight(:).SNRout]',...
                 mean([abs([summaryRight(:).SNRin]);abs([summaryRight(:).SNRout])])',...
                 [summaryRight(:).cIN]',...
                 [summaryRight(:).cOUT]',...
                 [summaryRight(:).DiceRight]',...
                 [summaryRight(:).HaussRight]'.*[summaryLeft(:).cf]'];

corrSNRdiceRight = array2table(corr(diceAndSnrMatRight));
corrSNRdiceRight.Properties.VariableNames= {'SNRin','SNRout','SNRmean','cIN','cOUT','Dice','Hauss'};

%% discretize SNR and see if there are any differences - Left, Right, and Avg - wr Dice
compCol = 1

SNRcategoryLeft = diceAndSnrMatLeft(:,[compCol,6]);
SNRcategoryRight = diceAndSnrMatRight(:,[compCol,6]);

bins = {'Low','Average','High'};

SNRcategoryLeft = discretize(SNRcategoryLeft(:,1), ...
    [-inf,mean(SNRcategoryLeft(:,1))-std(SNRcategoryLeft(:,1)),mean(SNRcategoryLeft(:,1))+std(SNRcategoryLeft(:,1)),+inf],...
    'categorical',bins); %#ok<*NASGU>

SNRcategoryRight = discretize(SNRcategoryRight(:,1), ...
    [-inf,mean(SNRcategoryRight(:,1))-std(SNRcategoryRight(:,1)),mean(SNRcategoryRight(:,1))+std(SNRcategoryRight(:,1)),+inf],...
    'categorical',bins);


% pLeft = kruskalwallis(diceAndSnrMatLeft(:,6),SNRcategoryLeft,'off');
% pRight = kruskalwallis(diceAndSnrMatRight(:,6),SNRcategoryRight,'off');
% 
% figure, 
% h = figure, h.Position = [100 100 1000 800];

% subplot(131),boxchart(SNRcategoryLeft,diceAndSnrMatLeft(:,6)),ylabel(sprintf('Dice score (Left) - p-value (KW) = %.3f',pLeft)),ylim([0 1]),set(gca,'YTick',0:0.1:1);
% set(gca,'FontName','Calibri','FontSize',14, 'XColor','k','YColor','k');
% 
% grid on
% ax = gca; % Get handle to current axes.
% ax.XColor = 'k'; % Red
% ax.YColor = 'k'; % Blue
% ax.GridLineStyle = "-";
% ax.GridAlpha = 0.1;  % Make grid lines less transparent.
% ax.GridColor = 'k'; % Dark Green.
% 
% subplot(132),boxchart(SNRcategoryRight,diceAndSnrMatRight(:,6)),ylabel(sprintf('Dice score (Right) - p-value (KW) = %.3f',pRight)),ylim([0 1]),set(gca,'YTick',0:0.1:1);
% set(gca,'FontName','Calibri','FontSize',14, 'XColor','k','YColor','k');
% 
% grid on
% ax = gca; % Get handle to current axes.
% ax.XColor = 'k'; % Red
% ax.YColor = 'k'; % Blue
% ax.GridLineStyle = "-";
% ax.GridAlpha = 0.1;  % Make grid lines less transparent.
% ax.GridColor = 'k'; % Dark Green.

tryCorr = corr(diceAndSnrMatLeft(:,[2,6]));

meanMatrix = mean([diceAndSnrMatLeft(:,compCol),diceAndSnrMatRight(:,compCol)],2);
meanMatrix = [meanMatrix,mean([diceAndSnrMatLeft(:,6),diceAndSnrMatRight(:,6)],2)];


bins = {'Low','Average','High'};

SNRcategoryMean = discretize(meanMatrix(:,1), ...
    [-inf,mean(meanMatrix(:,1))-std(meanMatrix(:,1)),mean(meanMatrix(:,1))+std(meanMatrix(:,1)),+inf],...
    'categorical',bins);

h = figure, h.Position = [100 100 400 800];

pAvg = kruskalwallis(meanMatrix(:,2),SNRcategoryMean,'off');

subplot(121),boxchart(SNRcategoryMean,diceFromSource'),ylabel(sprintf('Dice score - p-value (KW) = %.3f',pAvg)),ylim([0 1]),set(gca,'YTick',0:0.1:1);
set(gca,'FontName','Calibri','FontSize',14, 'XColor','k','YColor','k');

grid on
ax = gca; % Get handle to current axes.
ax.XColor = 'k'; % Red
ax.YColor = 'k'; % Blue
ax.GridLineStyle = "-";
ax.GridAlpha = 0.1;  % Make grid lines less transparent.
ax.GridColor = 'k'; % Dark Green.

% saveas(h,fullfile(ll,'DiceSnrEPS.eps'),'epsc')
% 
% ll='D:\ProjectOpticNerveLinux\PROJECT-OPTIC-NERVE\PROJECT-OPTIC-NERVE\FIGURES\Paper\Figure7';
% tryCorr = corr(meanMatrix);

%% discretize SNR and see if there are any differences - Left, Right, and Avg - wr HD
compCol = 1

SNRcategoryLeft = diceAndSnrMatLeft(:,[compCol,7]);
SNRcategoryRight = diceAndSnrMatRight(:,[compCol,7]);

bins = {'Low','Average','High'};

SNRcategoryLeft = discretize(SNRcategoryLeft(:,1), ...
    [-inf,mean(SNRcategoryLeft(:,1))-std(SNRcategoryLeft(:,1)),mean(SNRcategoryLeft(:,1))+std(SNRcategoryLeft(:,1)),+inf],...
    'categorical',bins);

SNRcategoryRight = discretize(SNRcategoryRight(:,1), ...
    [-inf,mean(SNRcategoryRight(:,1))-std(SNRcategoryRight(:,1)),mean(SNRcategoryRight(:,1))+std(SNRcategoryRight(:,1)),+inf],...
    'categorical',bins);


pLeft = kruskalwallis(diceAndSnrMatLeft(:,7),SNRcategoryLeft,'off');
pRight = kruskalwallis(diceAndSnrMatRight(:,7),SNRcategoryRight,'off');

% figure, 
% h = figure, h.Position = [100 100 1000 800];

% subplot(131),boxchart(SNRcategoryLeft,diceAndSnrMatLeft(:,7)),ylabel(sprintf('HD (Left) - p-value (KW) = %.3f',pLeft)),ylim([0 7]),set(gca,'YTick',0:0.5:7);
% set(gca,'FontName','Calibri','FontSize',14, 'XColor','k','YColor','k');
% 
% grid on
% ax = gca; % Get handle to current axes.
% ax.XColor = 'k'; % Red
% ax.YColor = 'k'; % Blue
% ax.GridLineStyle = "-";
% ax.GridAlpha = 0.1;  % Make grid lines less transparent.
% ax.GridColor = 'k'; % Dark Green.
% 
% subplot(132),boxchart(SNRcategoryRight,diceAndSnrMatRight(:,7)),ylabel(sprintf('HD (Right) - p-value (KW) = %.3f',pRight)),ylim([0 7]),set(gca,'YTick',0:0.5:7);
% set(gca,'FontName','Calibri','FontSize',14, 'XColor','k','YColor','k');
% 
% grid on
% ax = gca; % Get handle to current axes.
% ax.XColor = 'k'; % Red
% ax.YColor = 'k'; % Blue
% ax.GridLineStyle = "-";
% ax.GridAlpha = 0.1;  % Make grid lines less transparent.
% ax.GridColor = 'k'; % Dark Green.
% 
% tryCorr = corr(diceAndSnrMatLeft(:,[2,7]));
% 
meanMatrix = mean([diceAndSnrMatLeft(:,compCol),diceAndSnrMatRight(:,compCol)],2);
meanMatrix = [meanMatrix,mean([diceAndSnrMatLeft(:,7),diceAndSnrMatRight(:,7)],2)];

bins = {'Low','Average','High'};

SNRcategoryMean = discretize(meanMatrix(:,1), ...
    [-inf,mean(meanMatrix(:,1))-std(meanMatrix(:,1)),mean(meanMatrix(:,1))+std(meanMatrix(:,1)),+inf],...
    'categorical',bins);

pAvg = kruskalwallis(meanMatrix(:,2),SNRcategoryMean,'off');

subplot(122),boxchart(SNRcategoryMean,mean([diceAndSnrMatLeft(:,7),diceAndSnrMatRight(:,7)],2)),ylabel(sprintf('HD - p-value (KW) = %.3f',pAvg)),ylim([0 7]),set(gca,'YTick',0:0.5:7);
set(gca,'FontName','Calibri','FontSize',14, 'XColor','k','YColor','k');

grid on
ax = gca; % Get handle to current axes.
ax.XColor = 'k'; % Red
ax.YColor = 'k'; % Blue
ax.GridLineStyle = "-";
ax.GridAlpha = 0.1;  % Make grid lines less transparent.
ax.GridColor = 'k'; % Dark Green.

% ll=fullfile();
% saveas(h,fullfile(ll,'DiceHdSnrEPS.eps'),'epsc')
% % 
% sgtitle('HD by SNR category')
% 
% tryCorr = corr(meanMatrix);

%% Find SNR categories and evaluate impact of SNR on clinical correctness
diametersErrors = readtable(fullfile(mainFolder,'EXCEL','DiameterMeasurements.xlsx'));

Errors = diametersErrors(:,["unet50OND","unet50ONSD","manualOND","manualONSD"]);

eOND = table2array(Errors(:,"unet50OND"))-table2array(Errors(:,"manualOND"));
eONSD = table2array(Errors(:,"unet50ONSD"))-table2array(Errors(:,"manualONSD"));

h = figure, h.Position = [100 100 400 800];

bins = {'Under','OK','Over'};

DiamONDquality = discretize(eOND,[-inf,-0.2,0.2,+inf],'categorical',bins);
DiamONSDquality = discretize(eONSD,[-inf,-0.2,0.2,+inf],'categorical',bins);

% OND

SNRinMean = mean([diceAndSnrMatLeft(:,1),diceAndSnrMatRight(:,1)],2);
pAvg = kruskalwallis(SNRinMean,DiamONDquality,'off');

subplot(121),boxchart(DiamONDquality,SNRinMean),ylabel(sprintf('SNR (sheets to nerve) dB - p-value (KW) = %.3f',pAvg));%,ylim([0 1]),set(gca,'YTick',0:0.1:1);
set(gca,'FontName','Calibri','FontSize',14, 'XColor','k','YColor','k');

grid on
ax = gca; % Get handle to current axes.
ax.XColor = 'k'; % Red
ax.YColor = 'k'; % Blue
ax.GridLineStyle = "-";
ax.GridAlpha = 0.1;  % Make grid lines less transparent.
ax.GridColor = 'k'; % Dark Green.

SNRoutMean = mean([diceAndSnrMatLeft(:,2),diceAndSnrMatRight(:,2)],2);
pAvg = kruskalwallis(SNRoutMean,DiamONDquality,'off');

subplot(122),boxchart(DiamONDquality,SNRoutMean),ylabel(sprintf('SNR (sheets to exterior) dB - p-value (KW) = %.3f',pAvg));%,ylim([0 1]),set(gca,'YTick',0:0.1:1);
set(gca,'FontName','Calibri','FontSize',14, 'XColor','k','YColor','k');

grid on
ax = gca; % Get handle to current axes.
ax.XColor = 'k'; % Red
ax.YColor = 'k'; % Blue
ax.GridLineStyle = "-";
ax.GridAlpha = 0.1;  % Make grid lines less transparent.
ax.GridColor = 'k'; % Dark Green.
% SNRavgMean = mean([diceAndSnrMatLeft(:,3),diceAndSnrMatRight(:,3)],2);
% pAvg = kruskalwallis(SNRavgMean,DiamONDquality,'off');
% 
% subplot(133),boxchart(DiamONDquality,SNRavgMean),ylabel(sprintf('SNR (averag in/out) dB - p-value (KW) = %.3f',pAvg));%,ylim([0 1]),set(gca,'YTick',0:0.1:1);
% sgtitle('SNR by clinical correctness of OND (thresh = +/- 0.2 mm)')

% ll='';
% saveas(h,fullfile(ll,'DiameterSnrEPSond.eps'),'epsc')

% ONSD

h = figure, h.Position = [100 100 400 800];

SNRinMean = mean([diceAndSnrMatLeft(:,1),diceAndSnrMatRight(:,1)],2);
pAvg = kruskalwallis(SNRinMean,DiamONSDquality,'off');

subplot(121),boxchart(DiamONSDquality,SNRinMean),ylabel(sprintf('SNR (sheets to nerve) dB - p-value (KW) = %.3f',pAvg));%,ylim([0 1]),set(gca,'YTick',0:0.1:1);
set(gca,'FontName','Calibri','FontSize',14, 'XColor','k','YColor','k');

grid on
ax = gca; % Get handle to current axes.
ax.XColor = 'k'; % Red
ax.YColor = 'k'; % Blue
ax.GridLineStyle = "-";
ax.GridAlpha = 0.1;  % Make grid lines less transparent.
ax.GridColor = 'k'; % Dark Green.

SNRoutMean = mean([diceAndSnrMatLeft(:,2),diceAndSnrMatRight(:,2)],2);
pAvg = kruskalwallis(SNRoutMean,DiamONSDquality,'off');

subplot(122),boxchart(DiamONSDquality,SNRoutMean),ylabel(sprintf('SNR (sheets to exterior) dB - p-value (KW) = %.3f',pAvg));%,ylim([0 1]),set(gca,'YTick',0:0.1:1);
set(gca,'FontName','Calibri','FontSize',14, 'XColor','k','YColor','k');

grid on
ax = gca; % Get handle to current axes.
ax.XColor = 'k'; % Red
ax.YColor = 'k'; % Blue
ax.GridLineStyle = "-";
ax.GridAlpha = 0.1;  % Make grid lines less transparent.
ax.GridColor = 'k'; % Dark Green.

% ll='';
% saveas(h,fullfile(ll,'DiameterSnrEPSonsd.eps'),'epsc')


% 
% SNRavgMean = mean([diceAndSnrMatLeft(:,3),diceAndSnrMatRight(:,3)],2);
% pAvg = kruskalwallis(SNRavgMean,DiamONDquality,'off');
% 
% subplot(133),boxchart(DiamONDquality,SNRavgMean),ylabel(sprintf('SNR (averag in/out) dB - p-value (KW) = %.3f',pAvg));%,ylim([0 1]),set(gca,'YTick',0:0.1:1);
% sgtitle('SNR by clinical correctness of ONSD (thresh = +/- 0.2 mm)')

%% make matrix for R analysis

% matR = [diceAndSnrMatLeft(:,1:2),diceAndSnrMatRight(:,1:2),diceAndSnrMatLeft(:,6),diceAndSnrMatRight(:,6),diceFromSource'];
% tabR = array2table(matR);
% tabR.Properties.VariableNames= {'L_SNRin','L_SNRout','R_SNRin','R_SNRout','DiceLeft','DiceRight','Dice'};
% writetable(tabR,'/home/francesco/PROJECT-OPTIC-NERVE/ExcelFiles/SNRandDiceSummary.xl')

function [summary] = calcSNRparams(ms,mnIN,mnOUT)

    if length(ms) > length(mnIN)
        sample = randi(length(ms),1,length(mnIN));
        ms = ms(sample);
        else if length(ms) > length(mnOUT) %#ok<SEPEX,ALIGN>
            sample = randi(length(ms),1,length(mnOUT));
            ms = ms(sample);
        end
    end

    if length(mnOUT) > length(ms)
        sample = randi(length(mnOUT),1,length(ms));
        mnOUT = mnOUT(sample);
    end

    if length(mnIN) > length(ms)
        sample = randi(length(mnIN),1,length(ms));
        mnIN = mnIN(sample);
    end

    %% noise IN
    s.ps = sum(ms.^2); s.pn = sum(mnIN.^2);
    f.pot_log = 20*log10(s.ps/s.pn); 
%     f.snr = (mean(ms)-mean(mnIN))/sqrt(std(ms)^2-std(mnIN)^2);
    
    summary.pot_log_IN = f.pot_log;

    %% noise OUT
    s.ps = sum(ms.^2); s.pn = sum(mnOUT.^2);
    f.pot_log = 20*log10(s.ps/s.pn); f.pot = s.ps/s.pn;

    summary.pot_log_OUT = f.pot_log;

    %% contrast IN
    s.ps = median(ms); s.pn = median(mnIN);
    f.contr = (s.ps-s.pn)/255;

    summary.contrast_IN = f.contr;

    %% contrast IN
    s.ps = median(ms); s.pn = median(mnOUT);
    f.contr = (s.ps-s.pn)/255;

    summary.contrast_OUT = f.contr;
end

