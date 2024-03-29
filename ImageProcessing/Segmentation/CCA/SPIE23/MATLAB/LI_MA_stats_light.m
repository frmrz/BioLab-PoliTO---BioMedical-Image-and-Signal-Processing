% LI_MA_stats_light calculates various statistics for LI and MA tracing precision assesment.
%
% Input:
%   - LI_GT: Ground truth profile for the LI (Lumen Intima) boundary
%   - MA_GT: Ground truth profile for the MA (Media Adventitia) boundary
%   - LI_AUTO: Automatically segmented profile for the LI region
%   - MA_AUTO: Automatically segmented profile for the MA region
%   - CF: Conversion factor from pixels to millimeters
%   - filename: Name of the file being processed
%   - mode: Mode of operation ('mm' or 'pxl')
%
% Output:
%   - ImtStat: Structure containing the calculated statistics
%
% Note: This function requires the profiles in common support.
%
% Example usage:
%   LI_GT = imread('LI_GT.png');
%   MA_GT = imread('MA_GT.png');
%   LI_AUTO = imread('LI_AUTO.png');
%   MA_AUTO = imread('MA_AUTO.png');
%   CF = 0.1;
%   filename = 'image.png';
%   mode = 'mm';
%   stats = LI_MA_stats_light(LI_GT, MA_GT, LI_AUTO, MA_AUTO, CF, filename, mode);
%
% Author: Francesco Marzola

function [ImtStat]= LI_MA_stats_light(LI_GT,MA_GT,LI_AUTO,MA_AUTO,CF,filename,mode)

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
    ImtStat.IMTVauto = sqrt(((std(IMTautoLi2Ma)*ImtStat.CF )^2+(std(IMTautoMa2Li)*ImtStat.CF )^2)/(numel(IMTautoLi2Ma) + numel(IMTautoMa2Li)));
    
    ImtStat.IMTauto = ImtStat.IMTauto*mod;
    
    %% GT
    [ImtStat.IMTgt,IMTgtLi2Ma,IMTgtMa2Li] = PolyDistMethod(LI_GT,MA_GT);     
    ImtStat.IMTVgt = sqrt(((std(IMTgtLi2Ma)*ImtStat.CF )^2+(std(IMTgtMa2Li)*ImtStat.CF )^2)/(numel(IMTgtLi2Ma) + numel(IMTgtMa2Li)));
    
    ImtStat.IMTgt = ImtStat.IMTgt*mod;

    %% Bias
    ImtStat.BiasIMT = ImtStat.IMTauto - ImtStat.IMTgt;
    ImtStat.BiasIMTV = ImtStat.IMTVauto - ImtStat.IMTVgt;

    %% FoM
%     ImtStat.FoM  = 100 - abs((ImtStat.IMTauto - ImtStat.IMTgt) / ImtStat.IMTgt) * 100;

    %% Hausdorff LI
    [ImtStat.HM_LI] = HaussdorfDistance(LI_GT,LI_AUTO);
    ImtStat.HM_LI = ImtStat.HM_LI*mod;
    
    %% Hausdorff MA
    [ImtStat.HM_MA] = HaussdorfDistance(MA_GT,MA_AUTO);
    ImtStat.HM_MA = ImtStat.HM_MA*mod;
    
    %% MAD IMT
%     ColDistAUTO = abs(LI_AUTO(:,2)-MA_AUTO(:,2));
%     ColDistGT = abs(LI_GT(:,2)-MA_GT(:,2));

%     ImtStat.MAD_IMT = mad(abs(ColDistAUTO-ColDistGT))*mod;
%     ImtStat.MAD_IMT = mean(abs(ColDistAUTO-ColDistGT))*mod;
%     ImtStat.MAD_IMT = mad(ColDistAUTO-ColDistGT)*mod;
    
    %% PDM LI
    [ImtStat.PDM_LI] = PolyDistMethod(LI_GT,LI_AUTO);
%     [ImtStat.PDM_LI] = mad(abs(LI_GT(:,2)-LI_AUTO(:,2)));
    ImtStat.PDM_LI = ImtStat.PDM_LI*mod;
    
    %% PDM MA
    [ImtStat.PDM_MA] = PolyDistMethod(MA_GT,MA_AUTO);
%     [ImtStat.PDM_MA] = mad(abs(MA_GT(:,2)-MA_AUTO(:,2)));
    ImtStat.PDM_MA = ImtStat.PDM_MA*mod;

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
