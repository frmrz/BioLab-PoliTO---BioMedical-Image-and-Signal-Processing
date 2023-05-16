clear 
close all
clc

%% 1. Selecting data

mainFolder =  '/media/francesco/DEV001/PROJECT-GITHUB/Image_Processing/Segmentation/Ultrasound-OpticNerveSheaths';

dataFolder = fullfile(mainFolder,'DATA');

imageFolder = fullfile(dataFolder,'IMAGES');
outputFolder = fullfile(dataFolder,'IMAGES_256_new');
resultsFolder = fullfile(mainFolder,'RESULTS');

debugFolder = fullfile(resultsFolder,'CROP_DEBUG');
rectFolder = fullfile(resultsFolder,'RECT');

%Operazione di Cropping per ognuna delle immagini
filenames = dir(imageFolder);

for i = 3:length(filenames)
        
    I = imread(fullfile(imageFolder, filenames(i).name));
    
    %convert the input image to a grayscale intensity image. 
    I=im2gray(I);

    %initial black border removal 
    %falg_dim assume valore zero se l'immagine ha un numero di righe
    %inferiore a 700, diversamente assume valore 1
    [Ic1, rect,flag_dim]=BlackBorderCrop(I);
 
    %another crop to avoid really black vertical sides --> basato
    %sull'intensità dei pixel delle colonne
    sumC = sum(Ic1,1);
    ind = find(sumC > max(sumC)*0.15);
    sumCbinary = sumC;
    sumCbinary(ind) = 1;
    sumCbinary(setdiff(1:length(sumC),ind)) = 0;
  
    cc = bwconncomp(logical(sumCbinary));
    stats = regionprops(cc);
    idx = find([stats.Area] == max([stats.Area])); 
    sumCbinaryClean = ismember(labelmatrix(cc), idx); 
    if idx~=1
        for j=1:(idx-1)
            dist1=(stats(idx).Centroid(1,1)-(stats(idx).Area\2))-...
                (stats(idx-j).Centroid(1,1)+(stats(idx-j).Area\2));
            if dist1<22
                sumCbinaryClean = sumCbinaryClean+ismember(labelmatrix(cc), idx-j);
            end
        end
    end
    if idx<size(stats,1)
        for j=1:(size(stats,1)-idx)
            dist2=(stats(idx+j).Centroid(1,1)-(stats(idx+j).Area\2))-...
                (stats(idx).Centroid(1,1)+(stats(idx).Area\2));
            if dist2<22
                sumCbinaryClean = sumCbinaryClean+ismember(labelmatrix(cc), idx+j);
            end
        end
    end
    
    xmin2 = find(sumCbinaryClean == 1,1,'first') - 1; 
    width2 = find(sumCbinaryClean == 1,1,'last')-xmin2;
    
    rect2 = [xmin2 1 width2 size(Ic1,1)];
    Ic = imcrop(Ic1,rect2);
    
   %looking for the best circle on ocular bulb 
    BW = zeros(size(Ic));
    BW(Ic < 20) = 1;
    
    [centers,radii,metric] = imfindcircles(BW,[80 120],'Sensitivity',...
        0.985,'ObjectPolarity','bright');
    
    %Setting parameters according to the Calibration Factor
    % --> in base al valore di CF le immagini vengono considerate più o
    % meno zoommate --> in base al livello di zoom verranno ricercate circonferenze 
    % che approssimano il bulbo oculare con raggi diversi 
    CF=getCalibrationFactor(I);
    if flag_dim==0
        if CF<0.07
            Rmax=round(150-(150*((0.1-CF)*10*0.5)));
            Rmin=round(80-(80*((0.1-CF)*10*0.5)));
            S=0.99;
        else
            Rmax=120;
            Rmin=80;
            S=0.99;
        end
    else
        Rmax=100;
        Rmin=78;
        S=0.985;
    end
    
    l = 0.005;
    while isempty(radii)
        [centers,radii,metric] = imfindcircles(Ic,[Rmin Rmax],'Sensitivity',...
            S + l,'ObjectPolarity','dark');
        
        l = l + 0.005;
    end
    
    for j = 1:length(radii)
        mask(:,:,j) = getCircleMask(Ic,centers(j,:),radii(j));
        maskBW(:,:,j) = mask(:,:,j).*BW;
        maskBWsum(j) = sum(sum(maskBW(:,:,j)));
    end
    
    [~,ind] = max(maskBWsum);
    
    BWfinal = maskBW(:,:,ind);
    
    %Individuazione della circonferenza che meglio approssima il bulbo ocoulare
    clear cc stats
    cc = bwconncomp(logical(BWfinal));
    stats = regionprops(cc,'BoundingBox','Area');
    idx = find([stats.Area] == max([stats.Area])); 
    BWfinalBIG = ismember(labelmatrix(cc), idx); 
    BoundingBoxBIG = stats(idx).BoundingBox;
    
    bottomRow = BoundingBoxBIG(2) +  BoundingBoxBIG(4);
    middleCol = (BoundingBoxBIG(3))\2 + BoundingBoxBIG(1);
    
    %Setting bottomRow according to the Calibration Factor
    %--> in base al valore di CF e quindi di zoom, bottomRow viene traslata verso 
    %il basso il modo da includere l'interesa estensione del nervo nella patch256 finale
    if CF<0.09
        bottomRow=bottomRow*(1+((0.1-CF)*10*0.7));
    end
    
    %Initialize 256x256 image
    Patch256 = uint8(zeros(256,256));
    rect256 = floor([middleCol-128 bottomRow-128]);
    rect256(rect256 <= 0) = 1;
    endR = rect256(2) + 255;
    endC = rect256(1) + 255;
    
    if endR > size(Ic,1)
        endR = size(Ic,1);
        rect256(2)=endR-255;
    end
    if endC > size(Ic,2)
        endC = size(Ic,2);
    end
    
    height = floor(endR - rect256(2)+1);
    width = floor(endC - rect256(1)+1);
    
    %Patch finale centrata sul nervo di dimensioni 256x256
    Patch256(1:height,1:width) = Ic(rect256(2):endR,rect256(1):endC);
    
    %% store RECT

    q = 1;

    for row = 1 : size(I,1)/2
        for col = 1 : size(I,2)/2
            test = I(row:row+255,col:col+255);
            if isequal(test,Patch256)
                flag=1;
                RECT(q,1)=row;
                RECT(q,2)=col;
                RECT(q,3)=256;
                RECT(q,4)=256;
                q=q+1;
            end
        end
    end

    %Creazione immagine di debug in cui si visualizzeranno il risultato
    %del primo crop (Ic), tutte le circonferenze individuate su Ic, la
    %circonferenza che meglio approssima il bulbo oculare e la patch finale 
    warning('off','images:initSize:adjustingMag');
%         fig = figure('units','normalized','outerposition',[0 0 1 1],'visible','on');
%         subplot(2,2,1),imshow(Ic),subplot(2,2,2),imshow(BW),...
%             hold on, viscircles(centers,radii);
%         subplot(2,2,3),imshow(BWfinal), hold on, plot(middleCol,bottomRow,'*r');subplot(2,2,4),imshow(Patch256);
 
%         save_image_jpg_eps([D(i).name(1:end-4) '_bulb'],debug_directory,...
%         debug_directory,0);

    ws = warning('off','all');
    warning(ws);
%     close(fig);

    imwrite(Patch256,fullfile(outputFolder,filenames(i).name));

    fileID = fopen(fullfile(rectFolder,strrep(filenames(i).name, '.png','.txt')),'w');
    fprintf(fileID,'%d %d %d %d',RECT(1),RECT(2),RECT(3),RECT(4));
    fclose(fileID);

    clear mask maskBW maskBWsum


end

% rmpath(fullfile(main_directory, 'code'));
% cd(current_directory)

        