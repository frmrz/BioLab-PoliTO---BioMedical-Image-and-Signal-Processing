function [PDM,DvB1SB2,DvB2SB1] = PolyDistMethod(B1,B2)
%
% This function calculates the distance between two profiles B1 and B2 by
% using the Polyline Distance Method (see Suri, JS et al. Pattern Analysis
% & Applications, 2000, 3:39-60). For purposes of CCA segmentation
% validation, here we consider B2 as Gorund Truth (even though inessential
% in practice).
%
% First, the distance of the vertices of profile B1 from the segments edges
% of profile B2 are computed. Profiles are then swapped and distances are
% calulated again for every vertex.
% The final PDM is the sum of the distances B1->B2 and B2->B1 over the
% total number of vertices of B1 and B2.
%
% PLEASE NOTE: with respect to the original formula in the manuscript, here
% we take absolute values of "Lambda" and of "distances", since vertexes of
% the profiles can lead to positive or negative distances.
%
% Input(s):     B1    vector of coordinates of vertexes of B1. First row in
%                     x-coord, second row is y-coord.
%               B2    vector of coordinates of vertexes of B2. First row in
%                     x-coord, second row is y-coord.
% Output(s):    PDM   Polyline Distance between B1 and B2.
%               DvB1SB2   Vector of distances of all the vertices of B1 from the segments of B2
%               DvB2SB1   Vector of distances of all the vertices of B2 from the segments of B1

% Author: F. Molinari
% Date: Sept. 22, 2009
% Version 1.01

%Make sure profiles are in row format:
B1 = TurnRow(B1);
B2 = TurnRow(B2);

S1 = size(B1);
S2 = size(B2);

% Start by computing distance of B1 vertexes from B2 segments (B2 is GT)

if S1(2) > 1 && S2(2) > 1
    for j = 1:S1(2)
        %swapping on all the segments of B2
        for k = 1:S2(2)-1
            Lambda = abs(((B2(2,k+1) - B2(2,k))*(B1(2,j) - B2(2,k)) + (B2(1,k+1) - B2(1,k))*(B1(1,j) - B2(1,k))) / ((B2(1,k+1) - B2(1,k)).^2 + (B2(2,k+1) - B2(2,k)).^2));

            if Lambda >= 0 && Lambda <=1
                Dvs(k) = abs(((B2(2,k+1) - B2(2,k))*(-B1(1,j) + B2(1,k)) + (B2(1,k+1) - B2(1,k))*(B1(2,j) - B2(2,k))) / (sqrt((B2(1,k+1) - B2(1,k)).^2 + (B2(2,k+1) - B2(2,k)).^2)));
            else
                d1 = sqrt((B1(1,j) - B2(1,k))^2 + (B1(2,j) - B2(2,k))^2);
                d2 = sqrt((B1(1,j) - B2(1,k+1))^2 + (B1(2,j) - B2(2,k+1))^2);
                Dvs(k) = min([d1 d2]);
            end
        end
        %here I computed all the distances between vertex j of B1 and segments
        %of B2 - Assigning the minimum distance to the variable DvB1sB2
        DvB1SB2(j) = min(Dvs);
    end

    % Swapping the profiles, in order to compute the distance between the B2
    % vertex and the B1 segments (B2 id GT)
    Temp = B1;
    B1 = B2;
    B2 = Temp;
    S1 = size(B1);
    S2 = size(B2);
    clear Temp Dvs

    for j = 1:S1(2)
        %swapping on all the segments of B2
        for k = 1:S2(2)-1
            Lambda = abs(((B2(2,k+1) - B2(2,k))*(B1(2,j) - B2(2,k)) + (B2(1,k+1) - B2(1,k))*(B1(1,j) - B2(1,k))) / ((B2(1,k+1) - B2(1,k)).^2 + (B2(2,k+1) - B2(2,k)).^2));

            if Lambda >= 0 && Lambda <=1
                Dvs(k) = abs(((B2(2,k+1) - B2(2,k))*(-B1(1,j) + B2(1,k)) + (B2(1,k+1) - B2(1,k))*(B1(2,j) - B2(2,k))) / (sqrt((B2(1,k+1) - B2(1,k)).^2 + (B2(2,k+1) - B2(2,k)).^2)));
            else
                d1 = sqrt((B1(1,j) - B2(1,k))^2 + (B1(2,j) - B2(2,k))^2);
                d2 = sqrt((B1(1,j) - B2(1,k+1))^2 + (B1(2,j) - B2(2,k+1))^2);
                Dvs(k) = min([d1 d2]);
            end
        end
        %here I computed all the distances between vertex j of B1 and segments
        %of B2 - Assigning the minimum distance to the variable DvB2sB1 (since
        %I previously swapped
        DvB2SB1(j) = min(Dvs);
    end

    % DvB1SB2 = sort(DvB1SB2);
    % if length(DvB1SB2) >= 20
    %     DvB1SB2 = DvB1SB2(1:20);
    %     S1 = 20;
    % else
    %     S1 = S1(2);
    % end
    % DvB1SB2 = sort(DvB1SB2);
    % if length(DvB2SB1) >= 20
    %     DvB2SB1 = DvB2SB1(1:20);
    %     S2 = 20;
    % else
    %     S2 = S2(2);
    % end
    PDM = (sum(abs(DvB1SB2)) + sum(abs(DvB2SB1)))/(S1(2) + S2(2));
    
else
    PDM = mean(abs(B1(2,:) - B2(2,:)));
    DvB1SB2 = [];
    DvB2SB1 = [];
end
