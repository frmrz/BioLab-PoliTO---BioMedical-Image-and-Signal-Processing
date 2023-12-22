function [Ic,Ic2,minX,minY] = fm_autocrop(I,info)
%function [Ic,minX,minY] = fm_autocrop(I,info)
%
% This function is developed for the automatic cropping of ultrasound
% images in DICOM3 format. The specific fields have been designed on the
% DICOM classes of the ATL HDI5000 scanner, altough they should be present
% and valid in any DICOM format.
%
% The field 'SequenceOfUltrasoundRegions' needs to be defined according to
% the DICOM dictionary in use. If it is not, an experience-based autocropping
% is carried out and a warning message is displayed.
%
% Inputs:  I = input image to be cropped
%          info = infos of the DICOM header (usually provided by dicominfo)
% Outputs: Ic = autocropped image

%--------------------
% Author: F. Molinari
% Date: 16 Dec 2010.

%tic
if size(I,3) > 1
    I2=I(:,:,2);
    I=I(:,:,1);
%     
%     R=I(:,:,1);
%     G=I(:,:,2);
%     sub=imsubtract(R,G);
%     
%     thr=0.4;%graythresh(sub);
%     XX=zeros(size(sub));
%     XX(sub>=thr)=1;
%     se=strel('disk',5);
%     I = immultiply(im2double(rgb2gray(I)),im2double(imcomplement(imerode(XX,se))));%(:,:,1);
%     I=im2uint8(I);
else
    I2=I;
end

% test if the DICOM header is correctly formatted
%------------------------------------------------
if isfield(info,'SequenceOfUltrasoundRegions')
    minX = info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMinX0;
    maxX = info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMaxX1;
    minY = info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMinY0;
    maxY = info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMaxY1;
else   %it is not: using standard space
%     disp('The DICOM header is not fully formatted or it is corrupted: trying cropping anyway ...')
    %minX=1;
    minX=round(size(I,2)*0.08);
    %minX = 70;
    maxX=round(size(I,2)*0.85);%maxX = 717;
    minY=round(size(I,1)*0.08);
    %minY = 1;%round(size(I,1)*0.1);
    maxY=round(size(I,1)*0.9);%
    
end
%keyboard
Ic = imcrop(I,[minX minY maxX-minX maxY-minY]);   %first crop
Ic2=imcrop(I2,[minX minY maxX-minX maxY-minY]);
[r c] = size(Ic);
%keyboard
%fprintf('\nfirst image dimensions = %dx%d',size(Ic,2),size(Ic,1));
% vertical refinement based on Sobel
%-----------------------------------
h = fspecial('sobel');
IFiltered = imfilter(Ic,h);

a = find(IFiltered(:,round(c/2)) >= 20);

%keyboard
%fprintf('\ncc = %d',round(c/2));
if ~isempty(a)
    Ic = Ic(a(1)+5:a(end)-5,:);
    
    Ic2 = Ic2(a(1)+5:a(end)-5,:);
    minY = minY + a(1) + 5;
    
    %figure,imshow(Ic);
    %pause
    
    
    %fprintf('\nsecond image dimensions = %dx%d',size(Ic,2),size(Ic,1));
    %fprintf('\na size = %d\ta(1) = %d\ta(end) = %d',length(a),a(1),a(end));
    % horizontal refinement: using zero valued columns
    %-------------------------------------------------
    mean_ind = find(mean(im2double(Ic')) == max(mean(im2double(Ic'))));
    zero_values = find(Ic(mean_ind,:) <= 1);

    %----------ADDED BY KRISTEN-------------
    diffzero=diff(zero_values);
    ind=find(diffzero>50);
    %fprintf('\nind = %d',ind);
    
    %Take out zero_values that have a diffzero over 50 and are at the
    %front, back, or are isolated:
    %keyboard
    if ~isempty(ind)
        if ind(1) == 1
            zero_values = zero_values(2:end); %Take out first
            diffzero = diffzero(2:end);
            ind = ind(2:end)-1; %-1 because the first element was taken out
        end
    end
    
    if ~isempty(ind)
        if ind(end) == length(diffzero)
            zero_values = zero_values(1:end-1); %Take out last
            ind = ind(1:end-1);
        end
        diffind = diff(ind);
        indclose = find(diffind == 1);
        takeout = ind(indclose + 1);
        zero_values(takeout) = [];
    end
    
    diffzero = diff(zero_values);
    ind = find(diffzero>50);
    
    %keyboard
    if ~isempty(ind) && zero_values(ind(1))<round(c/2)
        firstind=ind(1);
    else
        firstind=length(zero_values);
    end
    if ~isempty(ind) && zero_values(ind(end)+1)>round(c/2)
        lastind=ind(end)+1;
    else
        lastind=1;
    end
    
    %keyboard
    %----------ADDED BY KRISTEN-------------
    c0 = find(zero_values(1:firstind) < round(c/2));

    if isempty(c0) %|| c0(end) %< c*0.2
        c1 = 5;
    else
        c1 = zero_values(c0(end))+5;
    end
    c0 = find(zero_values(lastind:end) > round(c/2))+lastind-1;
    %fprintf('\nfirstind = %d lastind = %d',firstind,lastind);

    if isempty(c0) 
        c2 = c-5;
    else
        c2 = zero_values(c0(1))-5;
    end
    Ic = Ic(:,c1:c2);   %final crop
    Ic2 = Ic2(:,c1:c2);   %final crop
    minX = minX + c1;
    
    %keyboard
    %fprintf('\nmean_ind = %d\nzero_values size = %d\nc1 = %d\tc2 = %d',mean_ind,length(zero_values),c1,c2);
    if isempty(Ic)
        minX=90;
        minY=52;
        Ic=imcrop(I,[minX minY maxX-minX maxY-minY]);
        Ic2=imcrop(I2,[minX minY maxX-minX maxY-minY]);
    end

else
%     disp('\nThis image has been probably zoomed and cropping result could be inaccurate.');
    minX = 0;
    maxX = 0;
    minY = 0;
    maxY = 0;
    [r c] = size(I);
    h = fspecial('sobel');
    IFiltered = imfilter(I,h);
    a = find(IFiltered(:,round(size(I,2)/2)) > 20);
    if ~isempty(a)
        Ic = I(a(1)+5:a(end)-5,:);
        minY = minY + a(1) + 5;

        % horizontal refinement: using zero valued columns
        %-------------------------------------------------
        mean_ind = find(mean(im2double(Ic')) == max(mean(im2double(Ic'))));
        zero_values = find(Ic(mean_ind,:) == 0);
        c0 = find(zero_values < c/2);
        if isempty(c0)
            c1 = 5;
        else
            c1 = zero_values(c0(end))+5;
        end
        c0 = find(zero_values > c/2);
        if isempty(c0)
            c2 = c-5;
        else
            c2 = zero_values(c0(1))-5;
        end
        Ic = Ic(:,c1:c2);   %final crop
        Ic2 = Ic2(:,c1:c2);   %final crop
        minX = minX + c1;
    else
%         disp('\nUnable to autocrop. Output image is same as input and original one. Sorry.');
        Ic = I;
        Ic2=Ic;
        minX = 0;
        minY = 0;
    end
end

%toc