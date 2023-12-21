clear
close all
clc

mainFolder =  '/media/francesco/DEV001/PROJECT-GITHUB/Image_Processing/Segmentation/Ultrasound-OpticNerveSheaths';

% --- DIRECTORIES SETUP ---
dataFolder = fullfile(mainFolder,'DATA');
resultsFolder = fullfile(mainFolder,'RESULTS');

masksFolder = fullfile(dataFolder,'LABELS_256');
imageFolder = fullfile(dataFolder,'IMAGES_256');
outputFolder = fullfile(resultsFolder,'PREDICTIONS');

bulbDistanceFolder = fullfile(dataFolder,'DIAMETERS_DISTANCE_TO_BULB');

cfFolder = fullfile(dataFolder,'CF');

diametersFolder=fullfile(resultsFolder,'DIAMETERS');
debugFolder = fullfile(resultsFolder,'DEBUG');

% Image filenames
filenames = dir(imageFolder);

% --- PROCESSING LOOP ---
for ii = 3:length(filenames)
    fprintf('Processing image %i\n', ii);

    % Reading ground truth, image, and prediction
    GT = imread(fullfile(masksFolder,filenames(ii).name)) == 1;
    I = imread(fullfile(imageFolder,filenames(ii).name));
    pred = imread(fullfile(outputFolder,filenames(ii).name));

    % Load calibration factor and bulb distance
    CF = load(fullfile(cfFolder,strrep(filenames(ii).name,'.png','.txt')));
    BulbDist = load(fullfile(bulbDistanceFolder,strrep(filenames(ii).name,'.png','.txt')));
           
    %% -----NERVE MASK CREATION-----
    %Find connected components in binary image (BW = predicted mak)
    BW=pred;
    cc = bwconncomp(BW);

    stats = regionprops(cc,'Extrema','BoundingBox','Centroid','Solidity','PixelList');

    % Checking object's connectivity
    
    if length(stats) == 1
        if BW(round(stats(1).Centroid(2)),...
                round(stats(1).Centroid(1))) == 0
            rows = round(stats(1).BoundingBox(2):stats(1).Centroid(2));
            rows(rows < 1) = [];

            BW(sub2ind(size(BW),rows,round(stats(1).Centroid(1).*...
                ones(1,length(rows))))) = 0;

            cc = bwconncomp(BW);
            stats = regionprops(cc,'Extrema','BoundingBox','Centroid','Solidity');
        end
    end

    %Get some of the extrema points of the objects:
     % --> top/bottom Right of first object and top/bottom Left of second one
    %Structure of the property "Extrema":
     % --> [top-left top-right right-top right-bottom bottom-right bottom-left left-bottom left-top]
     
    topL = stats(2).Extrema(1,:);
    bottomL = stats(2).Extrema(6,:);
    topR = stats(1).Extrema(2,:);
    bottomR = stats(1).Extrema(5,:);
        
    %locating the lowest point between the highest extreme points of the two objects
    lowestHighPoint = floor(max([stats(1).Extrema(1,2) stats(1).Extrema(2,2) ... 
        stats(2).Extrema(1,2) stats(2).Extrema(2,2)]));
    %locating the highest point between the lowest extreme points of the two objects
    highestLowPoint = floor(min([stats(1).Extrema(5,2) stats(1).Extrema(6,2) ...
        stats(2).Extrema(5,2) stats(2).Extrema(6,2)]));

    %Create mask containing the optic nerve
    BWnerve = roipoly(BW,[topL(1) topR(1) bottomR(1) bottomL(1)],...
        [topL(2) topR(2) bottomR(2) bottomL(2)]);
    BWnerve = BWnerve & ~BW;
    BWnerve = imopen(BWnerve,strel('disk',7));
    BWnerve(1:lowestHighPoint,:) = 0;
    BWnerve(highestLowPoint:end,:) = 0;

    % Start creating debug image --> figure with 4 images --> The first one is the nerve mask
    %     warning('off','images:initSize:adjustingMag');
    %     fig = figure('visible','on');
    %     subplot(2,2,1), imshow(BWnerve)
        
    %% -----CREATING THE CENTERLINE-----

    cont = 1;

    %loop that runs all (and only) the length of the nerve

    for kk = lowestHighPoint:highestLowPoint

        %saving of the extreme left and right points of the nerve, respectively in first and last
        first = find(BWnerve(kk,:) == 1,1,'first');
        last = find(BWnerve(kk,:) == 1,1,'last');

        %calculation of the intermediate point between first and last
         % --> saving the coordinates of the points that make up the centerline

        if ~isempty(first)
           centerline(cont,2) = kk;
           centerline(cont,1) = round((last-first)/2 + first);
%                hold on, plot(centerline(cont,1),centerline(cont,2),'*r');
           cont = cont + 1;
        end

    end
        
    %creating a temporary mask that represents the centerline
    skel_tmp = false(size(BWnerve));
    skel_tmp(sub2ind(size(BWnerve),centerline(:,2),centerline(:,1))) = 1;
        
    %Keep largest object to avoid small disconnected dots 
     % --> saving in the variable "skel" the definitive mask that represents the centerline
    cc_skel = bwconncomp(skel_tmp); 
    stats_skel = regionprops(cc_skel,'Area'); 
    idx_skel = find([stats_skel.Area] == max([stats_skel.Area])); 
    skel = ismember(labelmatrix(cc_skel),idx_skel);  
        
    % Continues the construction of the debug figure --> the secon image is the centerline mask
    %     subplot(2,2,2), imshow(skel+BW)

    % Centerline approximated to a line
    f1 = fit(centerline([1 2:20:end-1 end],2),centerline([1 2:20:end-1 end],1),'poly1');
    % Coefficients of the straight line containing the centerline
    a = f1.p1;
    b = f1.p2;

    yy = 1:0.01:256;
    xx = a.*yy + b;

    xx = round(xx);
    yy = round(yy);

    ind = find(xx>0 & xx<256);

    % Coordinates of the points that make up the straight line containing the centerline
    xx = xx(ind);
    yy = yy(ind);

    %% -----SEARCH FOR THE BOTTOM OF THE EYEBALL-----
    % Get some of the extrema points of the objects:
     % --> top/bottom Left of first object and top/bottom Right of second one
    topL = stats(1).Extrema(1,:);
    bottomL = stats(1).Extrema(6,:);
    topR = stats(2).Extrema(2,:);
    bottomR = stats(2).Extrema(5,:);

    %locating a point 50 pixel higher than the highest extreme points of the two objects
    topMin = min([topR(2)-60, topL(2)-60]);
    if topMin <= 0
        topMin = 1;
    end
        
        
    BWbulb = roipoly(BW,[topL(1,1)-20 topR(1,1)+20 topR(1,1)+20 topL(1,1)-20],...
        [topL(1,2) topR(1,2) topMin topMin]);
    
    %Highlighted in the original image the ocula bulb area that is located inside a box, positioned 
    %immediately above the structure identified by predited mask and slightly wider. Everything else is put in black
    Ibulb = im2double(I).*im2double(BWbulb);
    Ibulb = imadd(Ibulb,im2double(~BWbulb));
    
    %Temporary eyeball mask created from "Ibulb" --> Ocular bulb is put to 1
    IbulbBWtmp = zeros(size(Ibulb));
    IbulbBWtmp(Ibulb < 30/255) = 1;
    
    %Keep largest object to avoid small disconnected areas
    % --> saving in the variable "IbulbBW" the definitive mask that represents the ocular bulb above the nerve
    ccBulb = bwconncomp(IbulbBWtmp); 
    statsBulb = regionprops(ccBulb,'Area'); 
    idxBulb = find([statsBulb.Area] == max([statsBulb.Area])); 
    IbulbBW = ismember(labelmatrix(ccBulb),idxBulb);  

    indBulb = IbulbBW(sub2ind(size(IbulbBW),round(yy),round(xx)));
    lastInd = find(indBulb == 1,1,'last');
    firstInd = find(indBulb == 1,1,'first');
    
    if yy(lastInd) > yy(firstInd)
        xBulb = xx(lastInd);
        yBulb = yy(lastInd);
    else
        xBulb = xx(firstInd); 
        yBulb = yy(firstInd);
    end

    %Continues the construction of the debug figure 
     % --> the third image is the original image in which the lowest point of the eyeball is indicated
%         subplot(2,2,3), imshow(I), hold on, plot(xBulb,yBulb,'*r');

    %% -----COORDINATES OF THE POINTS 3 mm AWAY FROM THE OCULAR BULBUS-----
    %Create mask for optic nerve centerline --> BWline
    BWline = false(size(IbulbBW));
    % remove points above bulb point
    ind = find(yy >= yBulb);
    yylow = yy(ind);
    xxlow = xx(ind);
    BWline(sub2ind(size(BWline),yylow,xxlow)) = 1;
    
    %Coefficient of the line perpendicular to optic nerve centerline
    m_perp = -a;

    %Find spots 3mm below ocular bulb
    distances = bwdistgeodesic(BWline,xBulb,yBulb).*CF;
    
    dmin=BulbDist-0.5;
    dmax=BulbDist+2.5;

    [y3mm,x3mm] = find(distances > dmin & distances < dmax);
    y3mm=round(mean(y3mm));
    x3mm=round(mean(x3mm));
      
        
    %% -----CALCULATION OF OND AND ONSD (for each spots 3mm below the ocular bulb)-----
    %Do some initialization for loop
    BWlinesp = false(size(IbulbBW,1),size(IbulbBW,2),length(x3mm));
    
    x3p = zeros(length(x3mm));
    y3p = zeros(length(y3mm));
    
    xxp = zeros(length(x3mm),length(1:0.25:256));
    yyp = zeros(length(y3mm),length(1:0.25:256));
    
    ONS_L = zeros(length(y3mm),2);
    ON_L = zeros(length(y3mm),2);
    ONS_R = zeros(length(y3mm),2);
    ON_R = zeros(length(y3mm),2);
    
    OND = zeros(length(y3mm),1);
    ONSD = zeros(length(y3mm),1);
        
    L_mask=GT+GT+pred;
    %L_mask overlaid on the original image representing the different labels with different colors
    % --> FP=Blue, FN=Red, TP=Green
%         Im_overlay=labeloverlay(I,L_mask,'Colormap',[0 0 1;1 0 0;0 1 0],'Transparency',0.8);
    
%         Continues the construction of the debug figure
%         --> The fourth image is the original image on which the diameters will be plotted (during the loop)
%         subplot(2,2,4), imshow(Im_overlay);
    %Also represent the line containing the centerline on the original image
%         hold on, plot(xx,yy,'-y')
        
    %% loop for each of the 3mm points found
    for k = 1:length(x3mm)
        x3p(k) = x3mm(k);
        y3p(k) = y3mm(k);
        
        %Get the line perpendicular to optic nerve centerline through one of the points at 3mm
        xxp(k,:) = 1:0.25:256;
        yyp(k,:) = m_perp .* (xxp(k,:) - x3p(k)) + y3p(k);
        
        %Representation of the perpendicular line on the original image
        % --> different colors for each 3mm points
        colors = {'c','g','b','y','c','r','g','b','y','c'};
        hold on, plot(xxp(k,:),yyp(k,:),colors{k});
        
        BWlinestmp = false(size(IbulbBW));
        BWlinestmp(sub2ind(size(BWlinestmp),round(yyp(k,:)),...
            round(xxp(k,:)))) = 1;
        BWlinesp(:,:,k) = BWlinestmp;
        %Mask representing the segments corresponding to the overlap between the perpendicular line and the predicted mask
        BWlinesp(:,:,k) = BWlinesp(:,:,k) & logical(BW);
        
        %Extreme points of the segments
        ccD = bwconncomp(BWlinesp(:,:,k));
        statsD = regionprops(ccD,'Extrema','Area');
        plus = 0;
        if length(statsD)~=2
            plus = 0.5;
            try
                while length(statsD)~=2

                    yyp(k,:) = m_perp .* (xxp(k,:) - x3p(k)) + y3p(k) + plus;

                    BWlinestmp = false(size(IbulbBW));
                    BWlinestmp(sub2ind(size(BWlinestmp),round(yyp(k,:)),round(xxp(k,:)))) = 1;
                    BWlinesp(:,:,k) = BWlinestmp;
                    BWlinesp(:,:,k) = BWlinesp(:,:,k) & logical(BW);

                    %Extreme points of the segments
                    ccD = bwconncomp(BWlinesp(:,:,k));
                    statsD = regionprops(ccD,'Extrema','Area');
%                     imshow(BWlinestmp);
                    plus = plus + 0.5;

                end
            catch
            plus = -0.5;
                 while length(statsD)~=2

                    yyp(k,:) = m_perp .* (xxp(k,:) - x3p(k)) + y3p(k) + plus;

                    BWlinestmp = false(size(IbulbBW));
                    BWlinestmp(sub2ind(size(BWlinestmp),round(yyp(k,:)),round(xxp(k,:)))) = 1;
                    BWlinesp(:,:,k) = BWlinestmp;
                    BWlinesp(:,:,k) = BWlinesp(:,:,k) & logical(BW);

                    %Extreme points of the segments
                    ccD = bwconncomp(BWlinesp(:,:,k));
                    statsD = regionprops(ccD,'Extrema','Area');
%                     imshow(BWlinestmp);
                    plus = plus - 0.5;

                end
            end
                
        end
            
            
        %Get Optic Nerve (ON) and Optic Nerve Sheath (ONS) points
        %Get top/bottom right of first object and top/bottom left of second
        %[top-left top-right right-top right-bottom bottom-right bottom-left left-bottom left-top]
        ONS_L(k,:) = statsD(1).Extrema(7,:);
        ON_L(k,:) = statsD(1).Extrema(3,:);
        ON_R(k,:) = statsD(2).Extrema(7,:);
        ONS_R(k,:) = statsD(2).Extrema(3,:);

        %Calculation of the distances between these extreme points of the segments
        % --> Correspond to OND and ONSD (values in pixel)
        OND(k,1) = pdist([ON_L(k,:); ON_R(k,:)]);
        ONSD(k,1) = pdist([ONS_L(k,:); ONS_R(k,:)]);
        
        %Marked the points that represent the extremes of the diameters
        % --> ONDS = Green and OND = Red
%             hold on, plot(ONS_L(k,1),ONS_L(k,2),'*b');
%             hold on, plot(ON_L(k,1),ON_L(k,2),'*r');
%             hold on, plot(ON_R(k,1),ON_R(k,2),'*r');
%             hold on, plot(ONS_R(k,1),ONS_R(k,2),'*b');
    end

    %Final values in mm
    ONDfinal = mean(OND)*CF;
    ONSDfinal = mean(ONSD)*CF;

    %Saving the debug image to the debug directory
%         sgtitle(['OND = ' num2str(ONDfinal) '; ONSD = ' num2str(ONSDfinal) '; Dist = ' num2str(dmin+0.15+plus*CF)]);
%         print(fullfile(debug_directory,DIm(i).name),'-dpng','-r300');
%         save_image_jpg_eps([DIm(i).name(1:end-4) '_overlay.jpg'],debug_directory,...
%              debug_directory,0);
     
    ws = warning('off','all');
    warning(ws);
%     close(fig);
     
    %Creating the "Diameters" structure where image name and OND and ONSD values are saved
    Diameters(ii).Image = filenames(ii).name; 
    Diameters(ii).OND = ONDfinal;
    Diameters(ii).ONSD = ONSDfinal;
 
%     catch ME
%        % Error message displayed in case of error
%        fprintf([D(i).name ': error!\n']);
%        disp(ME.message)    
%     end

clearvars -except current_directory main_directory CNN_directory outputFolder bulbDistanceFolder ...
    debug_directory imageFolder masksFolder GT_directory images_directory D ONDall ONSDall i ...
    CF_pats diametersFolder cfFolder DIm filenames imageFolder Diameters save_directory DGT BulbDist debug_directory1
    
end

%Converting "Diameters" structure to table and creating an excel file
T=struct2table(Diameters);
T = T(3:end,:);
writetable(T,fullfile(diametersFolder,'Diameters_values.xlsx'));

%Saving "Diameters" structure containing all the diameters values for each image as a matlab file
% save(fullfile(save_directory,'Diameters_values.mat'),'Diameters');     