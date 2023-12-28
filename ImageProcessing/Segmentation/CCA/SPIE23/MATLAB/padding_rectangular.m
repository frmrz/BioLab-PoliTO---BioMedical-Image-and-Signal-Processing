% padding_rectangular pads the input image I to make it a square image by adding zeros on the sides.
% The function calculates the aspect ratio of the input image and determines the amount of padding required.
% If the aspect ratio is less than 1, padding is added to the top and bottom of the image.
% If the aspect ratio is equal to 1, no padding is added.
% If the aspect ratio is greater than 1, padding is added to the left and right of the image.
%
% Inputs:
%   - I: Input image
%
% Outputs:
%   - I_padded: Padded image
%   - n: Amount of padding added
%   - dims: Dimensions of the padded image
%
% Example usage:
%   I = imread('image.jpg');
%   [I_padded, n, dims] = padding_rectangular(I);
function [I_padded, n, dims] = padding_rectangular(I)

rows = size(I,1);
cols = size(I,2);
aspect_ratio = rows/cols;

%--------------------- padding ----------------------------------------

if aspect_ratio < 1
    n = ceil((cols-rows)/2);     
    I_padded = padarray(I,[n 0],0,'both'); 
elseif aspect_ratio == 1
    I_padded = I;
    n = 0;
else
    n = ceil((rows-cols)/2);     
    I_padded = padarray(I,[0 n],0,'both'); 
end

dims = size(I_padded);