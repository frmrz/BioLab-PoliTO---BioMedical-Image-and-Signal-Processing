function [ImtStat]= LI_MA_stats_not_interp(LI_GT,MA_GT,LI_AUTO,MA_AUTO,CF,filename,mode)

    % this function requires the profiles in common support
    ImtStat.filename = filename;
    ImtStat.CF = CF;
    
    if strcmp(mode,'mm')
        mod = CF;
    end
    
    if strcmp(mode,'pxl')
        mod = 1;
    end
    
    LI_AUTO = TurnColumn(LI_AUTO);
    MA_AUTO = TurnColumn(MA_AUTO);
    LI_GT = TurnColumn(LI_GT);
    MA_GT = TurnColumn(MA_GT);
    
    %% AUTO
    [ImtStat.PDMauto,~,~] = PolyDistMethod(LI_AUTO,MA_AUTO);     
    ImtStat.PDMauto = ImtStat.PDMauto*mod;
    %% GT
    [ImtStat.PDMgt,~,~] = PolyDistMethod(LI_GT,MA_GT);     
    ImtStat.PDMgt = ImtStat.PDMgt*mod;

    %% Bias
    ImtStat.BiasPDM = ImtStat.PDMauto - ImtStat.PDMgt;
    ImtStat.AbsBiasPDM = abs(ImtStat.PDMauto - ImtStat.PDMgt);

    %% AUTO
    ImtStat.EDauto = mean(pdist2(LI_AUTO,MA_AUTO,'euclidean','Smallest',1));     
    ImtStat.EDauto = ImtStat.EDauto*mod;
    
    %% GT
    ImtStat.EDgt = mean(pdist2(LI_GT,MA_GT,'euclidean','Smallest',1));     
    ImtStat.EDgt = ImtStat.EDgt*mod;
    
    %% Bias
    ImtStat.BiasED = ImtStat.EDauto - ImtStat.EDgt;
    ImtStat.AbsBiasED = abs(ImtStat.EDauto - ImtStat.EDgt);
end
