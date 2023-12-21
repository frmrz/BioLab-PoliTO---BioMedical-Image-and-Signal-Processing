function [Ic,rect,flag] = BlackBorderCrop(I)
%
% This function is developed for the automatic cropping of ultrasound.
% Authors: Filippo Molinari, Kristen Meiburger, Francesco Marzola
%
% Inputs:  I = input image to be cropped
%
% Outputs: Ic = autocropped image
%          rect = cropping dimensions
%          flag = its value can be 0 or 1. If the rows of the image I are less then 700, 
%                   flag will be equal to 0, otherwise it will be equal to 1.


minX=round(size(I,2)*0.08);
maxX=round(size(I,2)*0.85);
minY=round(size(I,1)*0.17);
maxY=round(size(I,1)*0.9);
    
%keyboard

%First Crop
Ic = imcrop(I,[minX minY maxX-minX maxY-minY]);  
[r c] = size(Ic);
%keyboard

%another crop to avoid really black vertical sides 
mean_ind = find(mean(im2double(Ic')) == max(mean(im2double(Ic'))));
zero_values = find(Ic(mean_ind,:) <= 1);

%----------ADDED BY KRISTEN-------------
diffzero=diff(zero_values);
ind=find(diffzero>50);

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
minX = minX + c1;
maxX=maxX+1-(c-c2);

%keyboard
if isempty(Ic)
    minX=90;
    minY=52;
    maxX=round(size(I,2)*0.85);
    maxY=round(size(I,1)*0.9);
    Ic=imcrop(I,[minX minY maxX-minX maxY-minY]);
end

%another crop to avoid really black horizontal side
if size(I,1)>700
    if size(I,1)>750
        Rowmax=round(size(Ic,1)*0.25);
    else
        Rowmax=round(size(Ic,1)*0.095);
    end
    sumR = sum(Ic(1:Rowmax,:),2);
    ind = find(sumR > max(sumR)*0.1);
    sumRbinary = sumR;
    sumRbinary(ind) = 1;
    sumRbinary(setdiff(1:length(sumR),ind)) = 0;
 
    minY_temp = find(sumRbinary == 1,1,'last');
    minY=minY + minY_temp;
    height_temp = size(Ic,1)-minY_temp-1;
    
    rect_temp = [1 minY_temp size(Ic,2) height_temp];
    Ic = imcrop(Ic,rect_temp);
    flag=1;
else
    flag=0;
end
rect=[minX minY maxX-minX maxY-minY];
end
