close all
clear

%% script for the evaluation metrics of segmentation on TEST set results
% for now the evaluation is on the RESIZED images DICE/HD
% BiasIMT need to be on the originals

mainDir = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/';
saveDir = fullfile(mainDir,"RESULTS","MATLAB","TEST-OVERLAY");

imgDir = fullfile(mainDir,"DATA","DEVELOPMENT","TEST","IMAGES-RESIZED");
gtDir = fullfile(mainDir,"DATA","DEVELOPMENT","TEST","MASKS-CONSENSUS-RESIZED");
cfDir = fullfile(mainDir,"DATA","DEVELOPMENT","TEST","CF");

load('/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/TEST/interOP.mat')

% v0_dir = fullfile(mainDir,"RESULTS","PYTHON","UNET-480x480-v0","Test_output");
% v1_dir = fullfile(mainDir,"RESULTS","PYTHON","UNET-480x480-v1","Test_output");
% v2_dir = fullfile(mainDir,"RESULTS","PYTHON","UNET-480x480-v2","Test_output");
% v3_dir = fullfile(mainDir,"RESULTS","PYTHON","UNET-480x480-v3","Test_output");
% v4_dir = fullfile(mainDir,"RESULTS","PYTHON","UNET-480x480-v4","Test_output");
% 
% patch_v0_dir = fullfile(mainDir,"RESULTS","PYTHON","UNET-64x64-v0","Test_output");
% patch_v1_dir = fullfile(mainDir,"RESULTS","PYTHON","UNET-64x64-v1","Test_output");
% patch_v2_dir = fullfile(mainDir,"RESULTS","PYTHON","UNET-64x64-v2","Test_output");
% patch_v3_dir = fullfile(mainDir,"RESULTS","PYTHON","UNET-64x64-v3","Test_output");
% patch_v4_dir = fullfile(mainDir,"RESULTS","PYTHON","UNET-64x64-v4","Test_output");

patch96_v0_dir = fullfile(mainDir,"RESULTS","PYTHON","UNET-96x96-v0","Test_output");
patch96_v1_dir = fullfile(mainDir,"RESULTS","PYTHON","UNET-96x96-v1","Test_output");
patch96_v2_dir = fullfile(mainDir,"RESULTS","PYTHON","UNET-96x96-v2","Test_output");
patch96_v3_dir = fullfile(mainDir,"RESULTS","PYTHON","UNET-96x96-v3","Test_output");
patch96_v4_dir = fullfile(mainDir,"RESULTS","PYTHON","UNET-96x96-v4","Test_output");

%% get files
files = table2struct(gt2Table(:,1));
for i = 1 :length(files)
    files(i).name = [files(i).name '.png'];
end

files = [struct('name','.');struct('name','..');files];

imgs = dir(imgDir);
gts = dir(gtDir);

% v0 = dir(v0_dir);
% v1 = dir(v1_dir);
% v2 = dir(v2_dir);
% v3 = dir(v3_dir);
% v4 = dir(v4_dir);
% 
% pv0 = dir(patch_v0_dir);
% pv1 = dir(patch_v1_dir);
% pv2 = dir(patch_v2_dir);
% pv3 = dir(patch_v3_dir);
% pv4 = dir(patch_v4_dir);

p96v0 = dir(patch96_v0_dir);
p96v1 = dir(patch96_v1_dir);
p96v2 = dir(patch96_v2_dir);
p96v3 = dir(patch96_v3_dir);
p96v4 = dir(patch96_v4_dir);

% h=figure;
% h.Position = [ 0 0 1800 900];
% h.Visible = 'off';

for ii = 3 : length(files)
    fprintf('%d\n',ii);

    temp_img = imread(fullfile(imgDir,files(ii).name));
    try
        temp_gt = imread(fullfile(gtDir,files(ii).name));
    catch
        temp_gt = zeros(size(temp_img));
    end

    try
        CF = load(fullfile(cfDir,strrep(files(ii).name,'.png','_CF.txt')));
    catch
        CF = 0.06;
%         write_txt_file(CF,strrep(files(ii).name,'.png','_CF.txt'),'/media/francesco/FMZ_archive/DIDATTICA/CHALLENGE/US-IMT-SEGMENTATION/DOC/CF');
    end

%     temp_v0 = imread(fullfile(v0_dir,v0(ii).name));
%     temp_v1 = imread(fullfile(v1_dir,v1(ii).name));
%     temp_v2 = imread(fullfile(v2_dir,v2(ii).name));
%     temp_v3 = imread(fullfile(v3_dir,v3(ii).name));
%     temp_v4 = imread(fullfile(v4_dir,v4(ii).name));
% 
%     patch_temp_v0 = imread(fullfile(patch_v0_dir,pv0(ii).name));
%     patch_temp_v1 = imread(fullfile(patch_v1_dir,pv1(ii).name));
%     patch_temp_v2 = imread(fullfile(patch_v2_dir,pv2(ii).name));
%     patch_temp_v3 = imread(fullfile(patch_v3_dir,pv3(ii).name));
%     patch_temp_v4 = imread(fullfile(patch_v4_dir,pv4(ii).name));

    patch96_temp_v0 = imread(fullfile(patch96_v0_dir,files(ii).name));
    patch96_temp_v1 = imread(fullfile(patch96_v1_dir,files(ii).name));
    patch96_temp_v2 = imread(fullfile(patch96_v2_dir,files(ii).name));
    patch96_temp_v3 = imread(fullfile(patch96_v3_dir,files(ii).name));
    patch96_temp_v4 = imread(fullfile(patch96_v4_dir,files(ii).name));

    %% process automated masks

%     temp_v0 = temp_v0>255*0.5; st = regionprops(temp_v0,'Area'); if ~isempty(st) ; temp_v0 = imfill(bwareaopen(temp_v0,max([st(:).Area])-1),'holes');end
%     temp_v1 = temp_v1>255*0.5; st = regionprops(temp_v1,'Area'); if ~isempty(st) ; temp_v1 = imfill(bwareaopen(temp_v1,max([st(:).Area])-1),'holes');end
%     temp_v2 = temp_v2>255*0.5; st = regionprops(temp_v2,'Area'); if ~isempty(st) ; temp_v2 = imfill(bwareaopen(temp_v2,max([st(:).Area])-1),'holes');end
%     temp_v3 = temp_v3>255*0.5; st = regionprops(temp_v3,'Area'); if ~isempty(st) ; temp_v3 = imfill(bwareaopen(temp_v3,max([st(:).Area])-1),'holes');end
%     temp_v4 = temp_v4>255*0.5; st = regionprops(temp_v4,'Area'); if ~isempty(st) ; temp_v4 = imfill(bwareaopen(temp_v4,max([st(:).Area])-1),'holes');end
% 
%     patch_temp_v0 = patch_temp_v0>255*0.5; st = regionprops(patch_temp_v0,'Area'); if ~isempty(st) ; patch_temp_v0 = imfill(bwareaopen(patch_temp_v0,max([st(:).Area])-1),'holes');end
%     patch_temp_v1 = patch_temp_v1>255*0.5; st = regionprops(patch_temp_v1,'Area'); if ~isempty(st) ; patch_temp_v1 = imfill(bwareaopen(patch_temp_v1,max([st(:).Area])-1),'holes');end
%     patch_temp_v2 = patch_temp_v2>255*0.5; st = regionprops(patch_temp_v2,'Area'); if ~isempty(st) ; patch_temp_v2 = imfill(bwareaopen(patch_temp_v2,max([st(:).Area])-1),'holes');end
%     patch_temp_v3 = patch_temp_v3>255*0.5; st = regionprops(patch_temp_v3,'Area'); if ~isempty(st) ; patch_temp_v3 = imfill(bwareaopen(patch_temp_v3,max([st(:).Area])-1),'holes');end
%     patch_temp_v4 = patch_temp_v4>255*0.5; st = regionprops(patch_temp_v4,'Area'); if ~isempty(st) ; patch_temp_v4 = imfill(bwareaopen(patch_temp_v4,max([st(:).Area])-1),'holes');end

    patch96_temp_v0 = patch96_temp_v0>255*0.5; st = regionprops(patch96_temp_v0,'Area'); if ~isempty(st) ; patch96_temp_v0 = imfill(bwareaopen(patch96_temp_v0,max([st(:).Area])-1),'holes');end
    patch96_temp_v1 = patch96_temp_v1>255*0.5; st = regionprops(patch96_temp_v1,'Area'); if ~isempty(st) ; patch96_temp_v1 = imfill(bwareaopen(patch96_temp_v1,max([st(:).Area])-1),'holes');end
    patch96_temp_v2 = patch96_temp_v2>255*0.5; st = regionprops(patch96_temp_v2,'Area'); if ~isempty(st) ; patch96_temp_v2 = imfill(bwareaopen(patch96_temp_v2,max([st(:).Area])-1),'holes');end
    patch96_temp_v3 = patch96_temp_v3>255*0.5; st = regionprops(patch96_temp_v3,'Area'); if ~isempty(st) ; patch96_temp_v3 = imfill(bwareaopen(patch96_temp_v3,max([st(:).Area])-1),'holes');end
    patch96_temp_v4 = patch96_temp_v4>255*0.5; st = regionprops(patch96_temp_v4,'Area'); if ~isempty(st) ; patch96_temp_v4 = imfill(bwareaopen(patch96_temp_v4,max([st(:).Area])-1),'holes');end

    %% check length of GT measurement

    %% Dice score
    statsDice(ii-2).name = imgs(ii).name;

    supportGT = sum(temp_gt)>0;
    supportP = sum(patch96_temp_v4)>0;
    supportCommon1 = sum([supportP;supportGT]);

    lenGT = numel(supportGT(supportGT>0))*CF;

    %% common support
    supportGT = sum(temp_gt)>0;
    supportP0 = sum(patch96_temp_v0)>0;
    supportP1 = sum(patch96_temp_v1)>0;
    supportP2 = sum(patch96_temp_v2)>0;
    supportP4 = sum(patch96_temp_v4)>0;
    supportCommon = sum([supportGT;supportP0;supportP1;supportP2;supportP4]);
    
    temp_gt(:,supportCommon~=5)=0;
    patch96_temp_v0(:,supportCommon~=5)=0;
    patch96_temp_v1(:,supportCommon~=5)=0;
    patch96_temp_v2(:,supportCommon~=5)=0;
    patch96_temp_v4(:,supportCommon~=5)=0;

    % percentage of predicted support overlapped to GT, if low probably the
    % predicted profiles are shifted or detection phase gone wrong
    percentOverlapPred = numel(supportCommon(supportCommon1==2))/numel(supportP(supportP>0));

    if lenGT > 1 && percentOverlapPred > 0.1
%         statsDice(ii-2).v0 = diceCommonSupport(logical(temp_gt),logical(temp_v0));
%         statsDice(ii-2).v1 = diceCommonSupport(logical(temp_gt),logical(temp_v1));
%         statsDice(ii-2).v2 = diceCommonSupport(logical(temp_gt),logical(temp_v2));
%         statsDice(ii-2).v3 = diceCommonSupport(logical(temp_gt),logical(temp_v3));
%         statsDice(ii-2).v4 = diceCommonSupport(logical(temp_gt),logical(temp_v4));
%     
%         statsDice(ii-2).pv0 = dice(logical(temp_gt),logical(patch_temp_v0));
%         statsDice(ii-2).pv1 = dice(logical(temp_gt),logical(patch_temp_v1));
%         statsDice(ii-2).pv2 = dice(logical(temp_gt),logical(patch_temp_v2));
%         statsDice(ii-2).pv3 = dice(logical(temp_gt),logical(patch_temp_v3));
%         statsDice(ii-2).pv4 = dice(logical(temp_gt),logical(patch_temp_v4));
    
%         statsDice(ii-2).p96v0 = dice(logical(temp_gt),logical(patch96_temp_v0));
%         statsDice(ii-2).p96v1 = dice(logical(temp_gt),logical(patch96_temp_v1));
%         statsDice(ii-2).p96v2 = dice(logical(temp_gt),logical(patch96_temp_v2));
%         statsDice(ii-2).p96v3 = dice(logical(temp_gt),logical(patch96_temp_v3));
%         statsDice(ii-2).p96v4 = dice(logical(temp_gt),logical(patch96_temp_v4));

        statsDice(ii-2).p96v0 = diceCommonSupport(logical(temp_gt),logical(patch96_temp_v0));
        statsDice(ii-2).p96v1 = diceCommonSupport(logical(temp_gt),logical(patch96_temp_v1));
        statsDice(ii-2).p96v2 = diceCommonSupport(logical(temp_gt),logical(patch96_temp_v2));
%         statsDice(ii-2).p96v3 = diceCommonSupport(logical(temp_gt),logical(patch96_temp_v3));
        statsDice(ii-2).p96v4 = diceCommonSupport(logical(temp_gt),logical(patch96_temp_v4));
    
        %% HD score
%         statsHD(ii-2).name = imgs(ii).name;
%     
%         statsHD(ii-2).v0 = imhausdorffCommonSupport(logical(temp_gt),logical(temp_v0));
%         statsHD(ii-2).v1 = imhausdorffCommonSupport(logical(temp_gt),logical(temp_v1));
%         statsHD(ii-2).v2 = imhausdorffCommonSupport(logical(temp_gt),logical(temp_v2));
%         statsHD(ii-2).v3 = imhausdorffCommonSupport(logical(temp_gt),logical(temp_v3));
%         statsHD(ii-2).v4 = imhausdorffCommonSupport(logical(temp_gt),logical(temp_v4));
%         
%         statsHD(ii-2).pv0 = imhausdorffCommonSupport(logical(temp_gt),logical(patch_temp_v0));
%         statsHD(ii-2).pv1 = imhausdorffCommonSupport(logical(temp_gt),logical(patch_temp_v1));
%         statsHD(ii-2).pv2 = imhausdorffCommonSupport(logical(temp_gt),logical(patch_temp_v2));
%         statsHD(ii-2).pv3 = imhausdorffCommonSupport(logical(temp_gt),logical(patch_temp_v3));
%         statsHD(ii-2).pv4 = imhausdorffCommonSupport(logical(temp_gt),logical(patch_temp_v4));
%     
        statsHD(ii-2).p96v0 = imhausdorffCommonSupport(logical(temp_gt),logical(patch96_temp_v0));
        statsHD(ii-2).p96v1 = imhausdorffCommonSupport(logical(temp_gt),logical(patch96_temp_v1));
        statsHD(ii-2).p96v2 = imhausdorffCommonSupport(logical(temp_gt),logical(patch96_temp_v2));
%         statsHD(ii-2).p96v3 = imhausdorffCommonSupport(logical(temp_gt),logical(patch96_temp_v3));
        statsHD(ii-2).p96v4 = imhausdorffCommonSupport(logical(temp_gt),logical(patch96_temp_v4));
    else 
        fprintf(' In image %s : GT length < 1mm\n',imgs(ii).name);
    end

    %% visualize
%     myMap = [ 0.75 0 0; 1 1 0; 0 0.75 0];
% 
%     ov0 = labeloverlay(temp_img,temp_gt/255 + temp_v0/255*2,'Colormap',myMap, 'Transparency',0.7);
%     ov1 = labeloverlay(temp_img,temp_gt/255 + temp_v1/255*2,'Colormap',myMap, 'Transparency',0.7);
%     ov2 = labeloverlay(temp_img,temp_gt/255 + temp_v2/255*2,'Colormap',myMap, 'Transparency',0.7);
%     ov3 = labeloverlay(temp_img,temp_gt/255 + temp_v3/255*2,'Colormap',myMap, 'Transparency',0.7);
%     
%     subplot(4,8,[1,2,9,10]),imshow(ov0),title(sprintf('v0-Dice-%.2f-HD-%.2f',statsDice(ii-2).v0,statsHD(ii-2).v0));
%     subplot(4,8,[3,4,11,12]),imshow(ov1),title(sprintf('v1-Dice-%.2f-HD-%.2f',statsDice(ii-2).v1,statsHD(ii-2).v1));
%     subplot(4,8,[5,6,13,14]),imshow(ov2),title(sprintf('v2-Dice-%.2f-HD-%.2f',statsDice(ii-2).v2,statsHD(ii-2).v2));
%     subplot(4,8,[7,8,15,16]),imshow(ov3),title(sprintf('v3-Dice-%.2f-HD-%.2f',statsDice(ii-2).v3,statsHD(ii-2).v3));
% 
%     ov0 = labeloverlay(temp_img,temp_gt/255 + patch_temp_v0/255*2,'Colormap',myMap, 'Transparency',0.7);
%     ov1 = labeloverlay(temp_img,temp_gt/255 + patch_temp_v1/255*2,'Colormap',myMap, 'Transparency',0.7);
%     ov2 = labeloverlay(temp_img,temp_gt/255 + patch_temp_v2/255*2,'Colormap',myMap, 'Transparency',0.7);
%     ov3 = labeloverlay(temp_img,temp_gt/255 + patch_temp_v3/255*2,'Colormap',myMap, 'Transparency',0.7);
%     
%     subplot(4,8,[17,18,25,26]),imshow(ov0),title(sprintf('pv0-Dice-%.2f-HD-%.2f',statsDice(ii-2).pv0,statsHD(ii-2).pv0));
%     subplot(4,8,[19,20,27,28]),imshow(ov1),title(sprintf('pv1-Dice-%.2f-HD-%.2f',statsDice(ii-2).pv1,statsHD(ii-2).pv1));
%     subplot(4,8,[21,22,29,30]),imshow(ov2),title(sprintf('pv2-Dice-%.2f-HD-%.2f',statsDice(ii-2).pv2,statsHD(ii-2).pv2));
%     subplot(4,8,[23,24,31,32]),imshow(ov3),title(sprintf('pv3-Dice-%.2f-HD-%.2f',statsDice(ii-2).pv3,statsHD(ii-2).pv3));
% 
%     saveas(h,fullfile(saveDir,gt(ii).name))
end

statsDiceSummary = get_dataset_summary(statsDice);
statsHDSummary = get_dataset_summary(statsHD);

save('/media/francesco/DEV001/PROJECT-CUBS-SPIE/RESULTS/MATLAB/masksEvaluationConsensus-INTER-commonSupport.mat',"statsHDSummary","statsDiceSummary","statsDice","statsHD");

%% stat test dice and hd


function [diceScore] = diceCommonSupport(gt,pred)

    supportGT = sum(gt)>0;
    supportP = sum(pred)>0;

    commonSupp = sum([supportP;supportGT]);

%     initGT = find(supportGT,1,"first");
%     endGT = find(supportGT,1,"last");
% 
%     pred(:,1:initGT) = 0;
%     pred(:,endGT:end) = 0;
    gt(:,commonSupp~=2)=0;
    gt(:,commonSupp~=2)=0;

    pred(:,commonSupp~=2)=0;
    pred(:,commonSupp~=2)=0;

    diceScore = dice(gt,pred);
end

function [diceScore] = imhausdorffCommonSupport(gt,pred)
%     initGT = find(sum(gt)>0,1,"first");
%     endGT = find(sum(gt)>0,1,"last");
% 
%     pred(:,1:initGT) = 0;
%     pred(:,endGT:end) = 0;

    supportGT = sum(gt)>0;
    supportP = sum(pred)>0;

    commonSupp = sum([supportP;supportGT]);

    gt(:,commonSupp~=2)=0;
    gt(:,commonSupp~=2)=0;

    pred(:,commonSupp~=2)=0;
    pred(:,commonSupp~=2)=0;

    diceScore = imhausdorff(gt,pred);
end
