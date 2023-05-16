function CF = getCalibrationFactor(I)
%
% La fuzione calcola in modo automatico il Fattore di Calibrazione espresso in mm/pixel, a partire 
% dall0immagine originale non croppata. 
% Ricerca la scala numerica laterale verticale (in alcune immagini a dx e in altre a sx), individua 
% le tacche principali (la cui distanza corrisponde a 10 mm reali), misura la distanza in pixel tra le tacche 
% principali e ne esegue la media (avarage_vertical_dist). 
% Infine ottiene il valore di CF come: 10/avarage_vertical_dist. 
%
% Inputs:  I = input image 
%
% Outputs: CF = Calibration Factor in mm/pixel

%Border removal
minY=round(size(I,1)*0.156);
maxY=round(size(I,1)*0.95);
minX=round(size(I,2)*0.06);
rect=[minX, minY, size(I,2)-minX, maxY-minY];
Ic=imcrop(I, rect);

BW=zeros(size(Ic));
BW(Ic>135)=1;

sumC=sum(BW,1);
ind=find(sumC~=0);
RowMin=round(size(BW,1)*0.58);
ColMin=round(size(BW,2)*0.8);
if sum(BW(RowMin:size(BW,1),ColMin:size(BW,2)))==0
    %Ricerca da sx
    first_ind=find(sumC>14,1,'first');
    last_ind=first_ind+26;
    widthmin=9;
else
    %Ricerca da dx
    last_ind=ind(end);
    first_ind=last_ind-6;
    widthmin=6;
end

BWcol=BW(:,first_ind:last_ind);

%Calibration Factor Calculation (CF)
cc = bwconncomp(logical(BWcol));
stats = regionprops(cc);
BoundingBox=[stats.BoundingBox];
height=BoundingBox(4:4:length(BoundingBox));
stats(height>9)=[];
BoundingBox=[stats.BoundingBox];
width=BoundingBox(3:4:length(BoundingBox));
stats(width<widthmin | width> 10)=[];

Centroid=[stats.Centroid];
if length(Centroid)>10
    a=2:2:10;
else
    a=2:2:length(Centroid);
end
dist_ind=diff(sort(Centroid(a)));
dist_takeout=find(diff(dist_ind)>8)+1;
dist_ind(dist_takeout)=[];
avarage_vertical_dist=mean(dist_ind); %equal to 10 mm
CF=10/avarage_vertical_dist; %mm per pixel
end

