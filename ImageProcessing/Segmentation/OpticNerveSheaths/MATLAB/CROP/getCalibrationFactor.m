function CF = getCalibrationFactor(I)

    %  The function automatically calculates the Calibration Factor in mm/pixel from the original uncropped image.
    %  It searches for the vertical numerical scale (either on the right or left side, depending on the image),
    %  identifies the main markers (which are spaced 10 mm apart), measures the pixel distance between the main markers,
    %  and calculates the average distance (avarage_vertical_dist).

    %  Finally, it obtains the CF value as: 10/avarage_vertical_dist.

    % Inputs: I = input image

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
        % Search from left
        first_ind=find(sumC>14,1,'first');
        last_ind=first_ind+26;
        widthmin=9;
    else
        % Search from right
        last_ind=ind(end);
        first_ind=last_ind-6;
        widthmin=6;
    end

    BWcol=BW(:,first_ind:last_ind);

    % Calibration Factor Calculation (CF)
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
    CF=10/avarage_vertical_dist; % mm per pixel
end

