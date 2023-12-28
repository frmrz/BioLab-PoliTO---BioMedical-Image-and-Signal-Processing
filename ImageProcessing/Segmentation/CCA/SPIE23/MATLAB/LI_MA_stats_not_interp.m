
% LI_MA_stats_not_interp calculates various statistics for LI and MA tracing precision assesment when profile are not interpolated
%
% Input:
%   LI_GT: Ground truth LI (Luminance Index) data
%   MA_GT: Ground truth MA (Microaneurysm) data
%   LI_AUTO: Automatically generated LI data
%   MA_AUTO: Automatically generated MA data
%   CF: Conversion factor for unit conversion (mm or pxl)
%   filename: Name of the file being processed
%   mode: Mode of operation ('mm' or 'pxl')
%
% Output:
%   ImtStat: Structure containing the calculated statistics
%       - filename: Name of the file being processed
%       - CF: Conversion factor for unit conversion (mm or pxl)
%       - PDMauto: PolyDistMethod result for the automatically generated data
%       - PDMgt: PolyDistMethod result for the ground truth data
%       - BiasPDM: Bias in PolyDistMethod result (PDMauto - PDMgt)
%       - AbsBiasPDM: Absolute bias in PolyDistMethod result (|PDMauto - PDMgt|)
%       - EDauto: Mean Euclidean distance for the automatically generated data
%       - EDgt: Mean Euclidean distance for the ground truth data
%       - BiasED: Bias in mean Euclidean distance (EDauto - EDgt)
%       - AbsBiasED: Absolute bias in mean Euclidean distance (|EDauto - EDgt|)
%
% Note: This function assumes that the profiles have a common support.
%
% Example usage:
%   LI_GT = [1, 2, 3];
%   MA_GT = [4, 5, 6];
%   LI_AUTO = [7, 8, 9];
%   MA_AUTO = [10, 11, 12];
%   CF = 0.5;
%   filename = 'image.jpg';
%   mode = 'mm';
%   ImtStat = LI_MA_stats_not_interp(LI_GT, MA_GT, LI_AUTO, MA_AUTO, CF, filename, mode);
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
