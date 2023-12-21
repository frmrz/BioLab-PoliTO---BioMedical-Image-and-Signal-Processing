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