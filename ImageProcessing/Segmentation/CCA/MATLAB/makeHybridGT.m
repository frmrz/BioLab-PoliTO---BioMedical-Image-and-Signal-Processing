close all
clear

%% script for the creation of hybrid GT

mainDir = '/media/francesco/DEV001/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT/FULL';

imgDir = fullfile(mainDir,"IMAGES");

A1Dir = fullfile(mainDir,"MASKS-A1");
A1sDir = fullfile(mainDir,"MASKS-A1s");
TUMDir = fullfile(mainDir,"MASKS-TUM");
% HDir = fullfile(mainDir,"MASKS-HYBRID");
HDir = fullfile(mainDir,"MASKS-HYBRID-v1");

CFDir = fullfile(mainDir,"CF");
ProfilesDir = fullfile(mainDir,"PROFILES");

%% get files in struct

imgs = dir(imgDir);

A1 = dir(A1Dir);
A1s = dir(A1sDir);
TUM = dir(TUMDir);

CF = dir(CFDir);

A1p = dir(fullfile(ProfilesDir,'Manual-A1'));
A1sp = dir(fullfile(ProfilesDir,'Manual-A1s'));
TUMp = dir(fullfile(ProfilesDir,'Computerized-TUM_DE'));

%% loop through files
%     figure,
CVli = zeros(length(A1)-2,3);
CVma = zeros(length(A1)-2,3);

for i = 36 : length(A1)
    fprintf('%d\n',i);

    %% load images and profiles
    tempImg = imread(fullfile(imgDir, imgs(i).name));
    tempA1  = imread(fullfile(A1Dir, A1(i).name));
    tempA1s = imread(fullfile(A1sDir, A1(i).name));
    tempTUM = imread(fullfile(TUMDir, TUM(i).name));

    tempCF = load(fullfile(CFDir, CF(i).name));

    LI_A1 = load(fullfile(A1p(i).folder, strrep(imgs(i).name,'.tiff','-LI.txt')));
    MA_A1 = load(fullfile(A1p(i).folder, strrep(imgs(i).name,'.tiff','-MA.txt')));

    LI_A1s = load(fullfile(A1sp(i).folder, strrep(imgs(i).name,'.tiff','-LI.txt')));
    MA_A1s = load(fullfile(A1sp(i).folder, strrep(imgs(i).name,'.tiff','-MA.txt')));

    LI_TUM = TurnRow(load(fullfile(TUMp(i).folder, strrep(imgs(i).name,'.tiff','-LI.txt'))));
    MA_TUM = TurnRow(load(fullfile(TUMp(i).folder, strrep(imgs(i).name,'.tiff','-MA.txt'))));

    %% store profiles with upper and lower bound to calculate intersection and metrics

    A1pMask = false(size(tempA1));
    idxA1p = sub2ind(size(tempA1s),round(LI_A1(2,:)),round(LI_A1(1,:)));
    A1pMask(idxA1p) = true; % create mask of profile

    A1spMask = false(size(tempA1s));
    idxA1sp = sub2ind(size(tempA1s),round(LI_A1s(2,:)),round(LI_A1s(1,:)));
    A1spMask(idxA1sp) = true; % create mask of profile

    TUMpMask = false(size(tempTUM));
    idxTUMp = sub2ind(size(tempTUM),round(LI_TUM(2,:)),round(LI_TUM(1,:)));
    TUMpMask(idxTUMp) = true; % create mask of profile
    
    %% plot profiles
%     h=figure;
%     imshow(tempImg); hold on,
%     plot(LI_A1(1,:),LI_A1(2,:),"Color",'r', 'LineWidth',5);
%     plot(MA_A1(1,:),MA_A1(2,:),"Color",'r', 'LineWidth',5);
%     plot(LI_A1s(1,:),LI_A1s(2,:),"Color",'b', 'LineWidth',5, 'LineStyle','--');
%     plot(MA_A1s(1,:),MA_A1s(2,:),"Color",'b', 'LineWidth',5, 'LineStyle','--');
%     plot(LI_TUM(1,:),LI_TUM(2,:),"Color",'y', 'LineWidth',5);
%     plot(MA_TUM(1,:),MA_TUM(2,:),"Color",'y', 'LineWidth',5);
%            
%     export_fig(h,'/media/francesco/DEV001/PROJECT-CUBS-SPIE/PAPER/FIGURES/FOV.png','-png','-r300');
%     export_fig(h,'/media/francesco/DEV001/PROJECT-CUBS-SPIE/PAPER/FIGURES/FOV.eps','-eps','-r300');

    %% check IMTV, if high IMTV probably there's a plque, so take only the manual measurments

    [IMTa1,IMTa1Li2Ma,IMTa1Ma2Li] = PolyDistMethod(LI_A1,MA_A1);     
    IMTVa1 = sqrt(((std(IMTa1Li2Ma)*tempCF )^2+(std(IMTa1Ma2Li)*tempCF )^2)/(numel(IMTa1Li2Ma) + numel(IMTa1Ma2Li)));

    [IMTa1s,IMTa1sLi2Ma,IMTa1sMa2Li] = PolyDistMethod(LI_A1s,MA_A1s);     
    IMTVa1s = sqrt(((std(IMTa1sLi2Ma)*tempCF )^2+(std(IMTa1sMa2Li)*tempCF )^2)/(numel(IMTa1sLi2Ma) + numel(IMTa1sMa2Li)));

    %% quality factor A1-A1s

    [~,inter_A1,inter_A1s] = intersect(LI_A1(1,:),LI_A1s(1,:));

    iLI_A1 = LI_A1(:,inter_A1);    iLI_A1s = LI_A1s(:,inter_A1s);
    iMA_A1 = MA_A1(:,inter_A1);   iMA_A1s = MA_A1s(:,inter_A1s);

    biasLI = abs(iLI_A1(2,:)-iLI_A1s(2,:)); corrLI = abs(min(corrcoef(iLI_A1(2,:),iLI_A1s(2,:)),[],'all'));
    biasMA = abs(iMA_A1(2,:)-iMA_A1s(2,:)); corrMA = abs(min(corrcoef(iMA_A1(2,:),iMA_A1s(2,:)),[],'all'));

    cvLI = corrLI/mean(biasLI); % sort of a quality coeff for GT agreement
    cvMA = corrMA/mean(biasMA);

    %% quality factor A1-TUM
    [~,inter_A1,inter_TUM] = intersect(LI_A1(1,:),LI_TUM(1,:));

    iLI_A1 = LI_A1(:,inter_A1);    iLI_TUM = LI_TUM(:,inter_TUM);
    iMA_A1 = MA_A1(:,inter_A1);   iMA_TUM = MA_TUM(:,inter_TUM);

    biasLItum = abs(iLI_A1(2,:)-iLI_TUM(2,:)); corrLItum = abs(min(corrcoef(iLI_A1(2,:),iLI_TUM(2,:)),[],'all'));
    biasMAtum = abs(iMA_A1(2,:)-iMA_TUM(2,:)); corrMAtum = abs(min(corrcoef(iMA_A1(2,:),iMA_TUM(2,:)),[],'all'));

    cvLItum = corrLItum/mean(biasLItum);
    cvMAtum = corrMAtum/mean(biasMAtum);    

    %% quality factor A1s-TUM
    [~,inter_A1s,inter_TUM] = intersect(LI_A1s(1,:),LI_TUM(1,:));

    iLI_A1s = LI_A1s(:,inter_A1s);    iLI_TUM = LI_TUM(:,inter_TUM);
    iMA_A1s = MA_A1s(:,inter_A1s);    iMA_TUM = MA_TUM(:,inter_TUM);

    biasLItum1 = abs(iLI_A1s(2,:)-iLI_TUM(2,:)); corrLItum1 = abs(min(corrcoef(iLI_A1s(2,:),iLI_TUM(2,:)),[],'all'));
    biasMAtum1 = abs(iMA_A1s(2,:)-iMA_TUM(2,:)); corrMAtum1 = abs(min(corrcoef(iMA_A1s(2,:),iMA_TUM(2,:)),[],'all'));

    cvLItum1 = corrLItum1/mean(biasLItum1);
    cvMAtum1 = corrMAtum1/mean(biasMAtum1);

    %% calculate hybrid

    LI_Hybrid = [];
    MA_Hybrid = [];

    
    for jj = 1 : size(tempImg,2)

        %% check profiles defined at column jj
        fA1 = sum(LI_A1(1,:) == jj);
        fA1s = sum(LI_A1s(1,:) == jj);
        fTUM = sum(LI_TUM(1,:) == jj);
        
        qA1 = find(LI_A1(1,:) == jj);
        qA1s = find(LI_A1s(1,:) == jj);
        qTUM = find(LI_TUM(1,:) == jj);

        %% merge profiles
        
        %% zero condition: if high IMTV and low concordance with TUM take only manual

        if IMTVa1 > 0.0125 && IMTVa1s > 0.0125 
            if cvLItum < 0.6 ; irregularLI = 1; else ; irregularLI = 0; end
            if cvLItum1 < 0.6 ; irregularLIs = 1; else ; irregularLIs = 0; end
            if cvMAtum < 0.6 ; irregularMA = 1; else ; irregularMA = 0; end
            if cvMAtum1 < 0.6 ; irregularMAs = 1; else ; irregularMAs = 0; end
        else
            irregularLI = 0;
            irregularLIs = 0;
            irregularMA = 0;
            irregularMAs = 0;
        end

        %% first condition: one profile present -> take that profile

        if fA1 == 1 && fA1s == 0 && fTUM == 0
            LI_Hybrid = [LI_Hybrid, LI_A1(:,qA1)];
            MA_Hybrid = [MA_Hybrid, MA_A1(:,qA1)];
        elseif fA1 == 0 && fA1s == 1 && fTUM == 0
            LI_Hybrid = [LI_Hybrid, LI_A1s(:,qA1s)];
            MA_Hybrid = [MA_Hybrid, MA_A1s(:,qA1s)];
        elseif fA1 == 0 && fA1s == 0 && fTUM == 1
            LI_Hybrid = [LI_Hybrid, LI_TUM(:,qTUM)];
            MA_Hybrid = [MA_Hybrid, MA_TUM(:,qTUM)];
        end

        %% second condition (1): two profile present -> take avg between profiles with
        %% best concordance or single manual with best concordance to TUM when
        %% low concordance between all profiles

        if fA1 == 1 && fA1s == 1 && fTUM == 0
            %% LI
            if cvLI > 0.5 || (irregularLI > 0 && irregularLIs >0) % or plaque condition
                LI_Hybrid = [LI_Hybrid, [LI_A1(1,qA1);mean([LI_A1(2,qA1),LI_A1s(2,qA1s)])]]; % a1-a1s
            else
                if cvLItum > cvLItum1
                    LI_Hybrid = [LI_Hybrid, LI_A1(:,qA1)]; % a1-tum
                else
                    LI_Hybrid = [LI_Hybrid, LI_A1s(:,qA1s)]; % a1s-tum
                end
            end
            %% MA
            if cvMA > 0.5 || (irregularMA > 0 && irregularMAs >0)
                MA_Hybrid = [MA_Hybrid, [MA_A1(1,qA1);mean([MA_A1(2,qA1),MA_A1s(2,qA1s)])]]; % a1-a1s
            else
                if cvMAtum > cvMAtum1
                    MA_Hybrid = [MA_Hybrid, MA_A1(:,qA1)];
                else
                    MA_Hybrid = [MA_Hybrid, MA_A1s(:,qA1s)];
                end
            end
        end

        %% second condition (2): two profile present -> take avg between profiles with
        %% best concordance or single manual with best concordance to TUM when
        %% low concordance between all profiles

        if fA1 == 1 && fA1s == 0 && fTUM == 1
            %% LI
            if cvLItum > 0.5 && irregularLI == 0 && irregularLIs == 0 % and no plaque condition
                LI_Hybrid = [LI_Hybrid, [LI_A1(1,qA1);mean([LI_A1(2,qA1),LI_TUM(2,qTUM)])]]; % a1-a1s
            else
                if cvLItum > cvLItum1 || (irregularLI > 0 && irregularLIs >0) % or plaque condition
                    LI_Hybrid = [LI_Hybrid, LI_A1(:,qA1)]; % a1-tum
                else
                    LI_Hybrid = [LI_Hybrid, LI_TUM(:,qTUM)]; % a1-tum
                end
            end
            %% MA
            if cvMAtum > 0.5 && irregularMA == 0 && irregularMAs == 0
                MA_Hybrid = [MA_Hybrid, [MA_A1(1,qA1);mean([MA_A1(2,qA1),MA_TUM(2,qTUM)])]]; % a1-a1s
            else
                if cvMAtum > cvMAtum1 || (irregularMA > 0 && irregularMAs >0)
                    MA_Hybrid = [MA_Hybrid, MA_A1(:,qA1)];
                else
                    MA_Hybrid = [MA_Hybrid, MA_TUM(:,qTUM)];
                end
            end
        end
        
        %% second condition (3): two profile present -> take avg between profiles with
        %% best concordance or single manual with best concordance to TUM when
        %% low concordance between all profiles

        if fA1 == 0 && fA1s == 1 && fTUM == 1
            %% LI
            if cvLItum1 > 0.5 && irregularLI == 0 && irregularLIs == 0 % and no plaque condition
                LI_Hybrid = [LI_Hybrid, [LI_A1s(1,qA1s);mean([LI_A1s(2,qA1s),LI_TUM(2,qTUM)])]]; % a1-a1s
            else
                if cvLItum1 > cvLItum || (irregularLI > 0 && irregularLIs >0) % or plaque condition
                    LI_Hybrid = [LI_Hybrid, LI_A1s(:,qA1s)]; % a1-tum
                else
                    LI_Hybrid = [LI_Hybrid, LI_TUM(:,qTUM)]; % a1-tum
                end
            end
            %% MA
            if cvMAtum1 > 0.5 && irregularMA == 0 && irregularMAs == 0
                MA_Hybrid = [MA_Hybrid, [MA_A1s(1,qA1s);mean([MA_A1s(2,qA1s),MA_TUM(2,qTUM)])]]; % a1-a1s
            else
                if cvMAtum1 > cvMAtum || (irregularMA > 0 && irregularMAs >0)
                    MA_Hybrid = [MA_Hybrid, MA_A1s(:,qA1s)];
                else
                    MA_Hybrid = [MA_Hybrid, MA_TUM(:,qTUM)];
                end
            end
        end

       
        %% third condition: three profiles present -> take avg between profiles with
        %% best concordance or between manual with best concordance and TUM 

        if fA1 == 1 && fA1s == 1 && fTUM == 1
            %% LI
            if cvLI > 0.5 || (irregularLI > 0 && irregularLIs >0) % or plaque condition
                LI_Hybrid = [LI_Hybrid, [LI_A1(1,qA1);mean([LI_A1(2,qA1),LI_A1s(2,qA1s)])]]; % a1-a1s
            elseif cvLI < 0.50001 && cvLItum > 0.5
                LI_Hybrid = [LI_Hybrid, [LI_A1(1,qA1);mean([LI_A1(2,qA1),LI_TUM(2,qTUM)])]]; % a1-tum
            elseif cvLI < 0.50001 && cvLItum1 > 0.5
                LI_Hybrid = [LI_Hybrid, [LI_A1s(1,qA1s);mean([LI_A1s(2,qA1s),LI_TUM(2,qTUM)])]]; % a1s-tum
            else
                if cvLItum > cvLItum1
                    LI_Hybrid = [LI_Hybrid, [LI_A1(1,qA1);mean([LI_A1(2,qA1),LI_TUM(2,qTUM)])]]; % a1-tum
                else
                    LI_Hybrid = [LI_Hybrid, [LI_A1s(1,qA1s);mean([LI_A1s(2,qA1s),LI_TUM(2,qTUM)])]]; % a1s-tum
                end
            end
            %% MA
            if cvMA > 0.5 || (irregularMA > 0 && irregularMAs >0) % or plaque condition
                MA_Hybrid = [MA_Hybrid, [MA_A1(1,qA1);mean([MA_A1(2,qA1),MA_A1s(2,qA1s)])]]; % a1-a1s
            elseif cvMA < 0.50001 && cvMAtum > 0.5
                MA_Hybrid = [MA_Hybrid, [MA_A1(1,qA1);mean([MA_A1(2,qA1),MA_TUM(2,qTUM)])]]; % a1-tum
            elseif cvMA < 0.50001 && cvMAtum1 > 0.5
                MA_Hybrid = [MA_Hybrid, [MA_A1s(1,qA1s);mean([MA_A1s(2,qA1s),MA_TUM(2,qTUM)])]]; % a1s-tum
            else
                if cvMAtum > cvMAtum1
                    MA_Hybrid = [MA_Hybrid, [MA_A1(1,qA1);mean([MA_A1(2,qA1),MA_TUM(2,qTUM)])]]; % a1-tum
                else
                    MA_Hybrid = [MA_Hybrid, [MA_A1s(1,qA1s);mean([MA_A1s(2,qA1s),MA_TUM(2,qTUM)])]]; % a1s-tum
                end
            end
        end

    end

    CVli(i-2,:) = [cvLI,cvLItum,cvLItum1];
    CVma(i-2,:) = [cvMA,cvMAtum,cvMAtum1];

%     %% moving average to smooth profiles
    LI_Hybrid(2,:) = movmean(LI_Hybrid(2,:),5);
    MA_Hybrid(2,:) = movmean(MA_Hybrid(2,:),5);
% 
%     %% write profiles
% 
%     write_txt_file(LI_Hybrid,strrep(imgs(i).name,'.tiff','-LI.txt'),fullfile(ProfilesDir,"Hybrid-v1"));
%     write_txt_file(MA_Hybrid,strrep(imgs(i).name,'.tiff','-MA.txt'),fullfile(ProfilesDir,"Hybrid-v1"));
% 
%     %% write image         
%     mask = uint8(roipoly(tempImg,[LI_Hybrid(1,:),flip(MA_Hybrid(1,:))],[LI_Hybrid(2,:),flip(MA_Hybrid(2,:))])).*255;
%     imwrite(mask,fullfile(HDir,strrep(imgs(i).name,'.tiff','.png')));

    %% plot profiles
        h=figure;

    hold off
    imshow(tempImg); hold on,
    plot(LI_A1(1,:),LI_A1(2,:),"Color",'r', 'LineWidth',8);
    plot(MA_A1(1,:),MA_A1(2,:),"Color",'r', 'LineWidth',8);
    plot(LI_A1s(1,:),LI_A1s(2,:),"Color",'b', 'LineWidth',8, 'LineStyle','--');
    plot(MA_A1s(1,:),MA_A1s(2,:),"Color",'b', 'LineWidth',8, 'LineStyle','--');
    plot(LI_TUM(1,:),LI_TUM(2,:),"Color",'y', 'LineWidth',8);
    plot(MA_TUM(1,:),MA_TUM(2,:),"Color",'y', 'LineWidth',8);
    plot(LI_Hybrid(1,:),LI_Hybrid(2,:),"Color",'g',  'LineWidth',8, 'LineStyle',':');
    plot(MA_Hybrid(1,:),MA_Hybrid(2,:),"Color",'g',  'LineWidth',8, 'LineStyle',':');

    export_fig(h,'/media/francesco/DEV001/PROJECT-CUBS-SPIE/PAPER/FIGURES/three.png','-png','-r300');
    export_fig(h,'/media/francesco/DEV001/PROJECT-CUBS-SPIE/PAPER/FIGURES/three.eps','-eps','-r300');
end

% figure, 
% subplot(121),histogram(IMTVa1),title('IMTVa1');
% subplot(122),histogram(IMTVa1s),title('IMTVa1s');
% 
% p90a1 = prctile(IMTVa1,90);
% p90a1s = prctile(IMTVa1s,90);


figure, 
subplot(321),histogram(CVli(:,1)),title('CVli');
subplot(322),histogram(CVma(:,1)),title('CVma');
subplot(323),histogram(CVli(:,2)),title('CVliT');
subplot(324),histogram(CVma(:,2)),title('CVmaT');
subplot(325),histogram(CVli(:,3)),title('CVliT1');
subplot(326),histogram(CVma(:,3)),title('CVmaT1');

CVli = reshape(CVli,numel(CVli),1);
CVma = reshape(CVma,numel(CVma),1);

max(CVli(isfinite(CVli)))
min(CVli(isfinite(CVli)))

max(CVma(isfinite(CVma)))
min(CVma(isfinite(CVma)))

figure, 
subplot(121),histogram(CVli),title('CVli');
subplot(122),histogram(CVma),title('CVma');

p90CVli = prctile(CVli,90);
p90CVma = prctile(CVma,90);

p50CVli = prctile(CVli,50);
p50CVma = prctile(CVma,50);

p10CVli = prctile(CVli,10);
p10CVma = prctile(CVma,10);

% 
% figure, 
% subplot(321),histogram(cvLI),title('CV LI');
% subplot(322),histogram(cvMA),title('CV MA');
% 
% subplot(323),histogram(cvLItum),title('CV LI t ');
% subplot(324),histogram(cvMAtum),title('CV MA t');
% 
% subplot(325),histogram(cvLItum1),title('CV LI t1');
% subplot(326),histogram(cvMAtum1),title('CV MA t1');


function output1 = SlideFun(x,WindowLength)

    output1 = zeros(length(x)-WindowLength,1);
%     output2 = zeros(length(x)-WindowLength,1);
    
    for idx = 1:length(x)-WindowLength
        Block = x(idx:idx+WindowLength);
        output1(idx) = std(Block)/mean(Block);
%         output2(idx) = fun(Block);
    end

end

