close all
clear

mainDir = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/FULL';

bboxDir = fullfile(mainDir,'BBOX');
imgDir  = fullfile(mainDir,'IMAGES-RESIZED');
A1Dir   = fullfile(mainDir,'MASKS-HYBRID-RESIZED-v1');
% TUMDir  = fullfile(mainDir,'MASKS-TUM-RESIZED');

outImg = fullfile(mainDir,'IMAGES-DETECT');
outA1 = fullfile(mainDir,'MASKS-HYBRID-DETECT-v1');
% outTUM = fullfile(mainDir,'MASKS-TUM-DETECT');

if ~isfolder(outImg)
    mkdir(outImg);
end
if ~isfolder(outA1)
    mkdir(outA1);
end
% if ~isfolder(outTUM)
%     mkdir(outTUM);
% end
% 
files = dir(imgDir);

for i = 3 : length(files)
    try
        fprintf('%d\n',i);
        img = imread(fullfile(imgDir,files(i).name));
        a1 = imread(fullfile(A1Dir,files(i).name));
    %     tum = imread(fullfile(TUMDir,files(i).name));
        bb = imread(fullfile(bboxDir,files(i).name));
    
        box = regionprops(logical(bb),'BoundingBox');
    
    %     imshow(bb)
    
    %     ov = labeloverlay(img,logical(bb));
    %     imshow(ov)
    
        imgDetected = imcrop(img,box.BoundingBox);
        a1Detected  = imcrop(a1,box.BoundingBox);
    %     tumDetected = imcrop(tum,box.BoundingBox);
    
%         imwrite(imgDetected,fullfile(outImg,files(i).name));
        imwrite(a1Detected,fullfile(outA1,files(i).name));
    %     imwrite(tumDetected,fullfile(outTUM,files(i).name));
    catch
        fprintf('Error in %s',files(i).name)
    end
end

imgs = dir(outImg);
for i = 3 : length(imgs)
    imm = imread(fullfile(imgs(i).folder,imgs(i).name));
    [r(i-2),c(i-2)] = size(imm);
end

min(r)
min(c)