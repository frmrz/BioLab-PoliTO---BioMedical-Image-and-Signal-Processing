function [ImtStat]= LI_MA_stats(LI_GT,MA_GT,LI_AUTO,MA_AUTO,CF,filename,mode)

    % this function requires the profiles in common support
    ImtStat.filename = filename;
    ImtStat.CF = CF;
    
    if strcmp(mode,'mm')
        mod = CF;
    end
    
    if strcmp(mode,'pxl')
        mod = 1;
    end

    %% AUTO
    [ImtStat.IMTauto,IMTautoLi2Ma,IMTautoMa2Li] = PolyDistMethod(LI_AUTO,MA_AUTO);     
%     ImtStat.IMTVauto = sqrt(((std(IMTautoLi2Ma)*ImtStat.CF )^2+(std(IMTautoMa2Li)*ImtStat.CF )^2)/(numel(IMTautoLi2Ma) + numel(IMTautoMa2Li)));
    
    ImtStat.IMTauto = ImtStat.IMTauto*mod;
    
    %% GT
    [ImtStat.IMTgt,IMTgtLi2Ma,IMTgtMa2Li] = PolyDistMethod(LI_GT,MA_GT);     
%     ImtStat.IMTVgt = sqrt(((std(IMTgtLi2Ma)*ImtStat.CF )^2+(std(IMTgtMa2Li)*ImtStat.CF )^2)/(numel(IMTgtLi2Ma) + numel(IMTgtMa2Li)));
    
    ImtStat.IMTgt = ImtStat.IMTgt*mod;

    %% Bias
    ImtStat.BiasIMT = ImtStat.IMTauto - ImtStat.IMTgt;
%     ImtStat.BiasIMTV = ImtStat.IMTVauto - ImtStat.IMTVgt;

    %% FoM
    ImtStat.FoM  = 100 - abs((ImtStat.IMTauto - ImtStat.IMTgt) / ImtStat.IMTgt) * 100;

    %% Hausdorff LI
    [ImtStat.HM_LI] = HaussdorfDistance(LI_GT,LI_AUTO);
    ImtStat.HM_LI = ImtStat.HM_LI*mod;
    
    %% Hausdorff MA
    [ImtStat.HM_MA] = HaussdorfDistance(MA_GT,MA_AUTO);
    ImtStat.HM_MA = ImtStat.HM_MA*mod;
    
    %% Profile length
    ImtStat.commonLenPxl = height(LI_GT);
    ImtStat.commonLenMM = height(LI_GT)*CF;
    
    %% Point based Bias
    pointBiasLI = LI_GT(:,2)-LI_AUTO(:,2);
    ImtStat.mPointBiasLI = mean(pointBiasLI)*mod;
    ImtStat.sPointBiasLI = std(pointBiasLI)*mod;
    ImtStat.pPointBiasLI = (sum(pointBiasLI>0)/length(LI_GT))*100;
  
    pointBiasMA = MA_AUTO(:,2)-MA_GT(:,2);
    ImtStat.mPointBiasMA = mean(pointBiasMA)*mod;
    ImtStat.sPointBiasMA = std(pointBiasMA)*mod;
    ImtStat.pPointBiasMA = (sum(pointBiasMA>0)/length(MA_GT))*100;

end
