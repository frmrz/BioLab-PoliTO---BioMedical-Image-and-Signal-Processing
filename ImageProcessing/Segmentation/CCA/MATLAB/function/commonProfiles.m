function [comm_LI,comm_MA] = commonProfiles(img,GT1_LI,GT1_MA,GT2_LI,GT2_MA,GT3_LI,GT3_MA,GTM_LI,GTM_MA,U_LI,U_MA)
    
    commonSuppIdx = 1:size(img,2);
    commonSupp = zeros(1,size(img,2));

    [a,b,idx1] = intersect(GT1_LI(1,:),commonSuppIdx);


end