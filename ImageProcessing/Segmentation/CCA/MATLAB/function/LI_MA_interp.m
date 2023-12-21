function [LI,MA]= LI_MA_interp(li,ma)%,saveoption, filenames, directory)

% this function requires 2 vectors with li and ma profiles to be interpolated
%% avoid repeated columns
[~, idxli] = unique(li(:,1), 'rows');
[~, idxma] = unique(ma(:,1), 'rows');

li = li(idxli,:);
ma = ma(idxma,:);

%% GT profiles interpolation
LI = [li(1,1),ceil(li(1,1)): floor(li(end,1)), li(end,1)];
MA = [ma(1,1),ceil(ma(1,1)):floor(ma(end,1)), ma(end,1);];

[~,idxLI] = unique(round(LI));
[~,idxMA] = unique(round(MA));

LI = LI(idxLI);
MA = MA(idxMA);

LI = [LI; interp1(li(:,1),li(:,2),LI,'pchip')];    
MA = [MA; interp1(ma(:,1),ma(:,2),MA,'pchip')];    

%% avoid repeated columns
LI = unique(LI', 'rows')';
MA = unique(MA', 'rows')';

%% LI and MA intersection
[~,inter_LI,inter_MA] = intersect(LI(1,:),MA(1,:));

LI = LI(:,inter_LI);
MA = MA(:,inter_MA);

% %% writing profiles
% if saveoption == 1
%      write_txt_file(LI,filenames{1},directory)
%      write_txt_file(MA,filenames{2},directory)
% end