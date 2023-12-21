% FILEPATH: /media/francesco/DEV001/CENTRAL_PROJECTS_REPOSITORY/BioLab-PoliTO---BioMedical-Image-and-Signal-Processing/ImageProcessing/Segmentation/OpticNerveSheaths/MATLAB/CROP/getCircleMask.m
% getCircleMask - Generates a circular mask for image segmentation.
%
% Syntax:
%   mask = getCircleMask(I, center, radius)
%
% Input Arguments:
%   - I: Input image.
%   - center: Center coordinates of the circle [x, y].
%   - radius: Radius of the circle.
%
% Output Argument:
%   - mask: Binary mask representing the circular region of interest.
%
% Description:
%   This function generates a circular mask for image segmentation. It
%   creates a circular region of interest (ROI) with the specified center
%   and radius on the input image. The resulting mask is a binary image
%   where the pixels inside the circle are set to 1 and the pixels outside
%   the circle are set to 0.
%
% Example:
%   I = imread('image.jpg');
%   center = [100, 100];
%   radius = 50;
%   mask = getCircleMask(I, center, radius);
%
%   figure;
%   subplot(1, 2, 1);
%   imshow(I);
%   title('Original Image');
%
%   subplot(1, 2, 2);
%   imshow(mask);
%   title('Circular Mask');

function mask = getCircleMask(I,center,radius)



    thetaResolution = 2; 
    theta=(0:thetaResolution:360)'*pi/180;

    x = bsxfun(@times,radius',cos(theta));
    x = bsxfun(@plus,x,(center(:,1))');
    x = cat(1,x,nan(1,length(radius)));
    x = x(:);
    x = x(~isnan(x));

    y = bsxfun(@times,radius',sin(theta));
    y = bsxfun(@plus,y,(center(:,2))');
    y = cat(1,y,nan(1,length(radius)));
    y = y(:);
    y = y(~isnan(y));

    mask = roipoly(I,x,y);