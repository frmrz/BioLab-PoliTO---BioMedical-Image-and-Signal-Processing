function [consensuLI,consensuMA,cvLI,cvMA,flag_high_bias] = consensusGT(img,LI_A1,MA_A1,LI_A1B,MA_A1B,LI_A2,MA_A2,LI_A3,MA_A3)
    
    support = 1:size(img,2);
    flag_high_bias = 0;

    %% consensus LI
    ConsensusMatrixLI = double(zeros(4,size(img,2)));

    if ~isempty(LI_A1)
        [~,~,idxSuppA1] = intersect(round(LI_A1(1,:)),support);
        ConsensusMatrixLI(1,idxSuppA1(1):idxSuppA1(end)) = LI_A1(2,:);
    end
    if ~isempty(LI_A1B)
        [~,~,idxSuppA1B] = intersect(round(LI_A1B(1,:)),support);
        ConsensusMatrixLI(2,idxSuppA1B(1):idxSuppA1B(end)) = LI_A1B(2,:);
    end
    if ~isempty(LI_A2)
        [~,~,idxSuppA2] = intersect(round(LI_A2(1,:)),support);
        ConsensusMatrixLI(3,idxSuppA2(1):idxSuppA2(end)) = LI_A2(2,:);
    end
    if ~isempty(LI_A3)
        [~,~,idxSuppA3] = intersect(round(LI_A3(1,:)),support);
        ConsensusMatrixLI(4,idxSuppA3(1):idxSuppA3(end)) = LI_A3(2,:);
    end

    ConsensusMatrixLI(ConsensusMatrixLI==0)=nan;

    %% evaluate quality factor LI

    corr_1_1b = abs(min(corrcoef([ConsensusMatrixLI(1,:)',ConsensusMatrixLI(2,:)'],'rows','complete'),[],'all'));
    bias_1_1b = mean(abs(ConsensusMatrixLI(1,:)-ConsensusMatrixLI(2,:)),2,'omitnan');
    cvLI(1,1) = corr_1_1b/bias_1_1b;

    corr_1_2 = abs(min(corrcoef([ConsensusMatrixLI(1,:)',ConsensusMatrixLI(3,:)'],'rows','complete'),[],'all'));
    bias_1_2 = mean(abs(ConsensusMatrixLI(1,:)-ConsensusMatrixLI(3,:)),2,'omitnan');
    cvLI(1,2) = corr_1_2/bias_1_2;

    corr_1_3 = abs(min(corrcoef([ConsensusMatrixLI(1,:)',ConsensusMatrixLI(4,:)'],'rows','complete'),[],'all'));
    bias_1_3 = mean(abs(ConsensusMatrixLI(1,:)-ConsensusMatrixLI(4,:)),2,'omitnan');
    cvLI(1,3) = corr_1_3/bias_1_3;

    corr_1b_2 = abs(min(corrcoef([ConsensusMatrixLI(2,:)',ConsensusMatrixLI(3,:)'],'rows','complete'),[],'all'));
    bias_1b_2 = mean(abs(ConsensusMatrixLI(2,:)-ConsensusMatrixLI(3,:)),2,'omitnan');
    cvLI(1,4) = corr_1b_2/bias_1b_2;

    corr_1b_3 = abs(min(corrcoef([ConsensusMatrixLI(2,:)',ConsensusMatrixLI(4,:)'],'rows','complete'),[],'all'));
    bias_1b_3 = mean(abs(ConsensusMatrixLI(2,:)-ConsensusMatrixLI(4,:)),2,'omitnan');
    cvLI(1,5) = corr_1b_3/bias_1b_3;

    corr_2_3 = abs(min(corrcoef([ConsensusMatrixLI(3,:)',ConsensusMatrixLI(4,:)'],'rows','complete'),[],'all'));
    bias_2_3 = mean(abs(ConsensusMatrixLI(3,:)-ConsensusMatrixLI(4,:)),2,'omitnan');
    cvLI(1,6) = corr_2_3/bias_2_3;
   
    %% flag high bias images
    if bias_1_1b > 17 ||  bias_1_2 > 17 ||  bias_1_3 > 17 || ...
            bias_1b_2 > 17 ||  bias_1b_3 > 17 ||  bias_2_3 > 17
        flag_high_bias = 1;
    end

    %% consensus MA
    ConsensusMatrixMA = double(zeros(4,size(img,2)));

    if ~isempty(MA_A1)
        [~,~,idxSuppA1] = intersect(round(MA_A1(1,:)),support);
        ConsensusMatrixMA(1,idxSuppA1(1):idxSuppA1(end)) = MA_A1(2,:);
    end
    if ~isempty(MA_A1B)
        [~,~,idxSuppA1B] = intersect(round(MA_A1B(1,:)),support);
        ConsensusMatrixMA(2,idxSuppA1B(1):idxSuppA1B(end)) = MA_A1B(2,:);
    end
    if ~isempty(MA_A2)
        [~,~,idxSuppA2] = intersect(round(MA_A2(1,:)),support);
        ConsensusMatrixMA(3,idxSuppA2(1):idxSuppA2(end)) = MA_A2(2,:);
    end
    if ~isempty(MA_A3)
        [~,~,idxSuppA3] = intersect(round(MA_A3(1,:)),support);
        ConsensusMatrixMA(4,idxSuppA3(1):idxSuppA3(end)) = MA_A3(2,:);
    end

    ConsensusMatrixMA(ConsensusMatrixMA==0)=nan;

    %% evaluate quality factor MA

    corr_1_1b = abs(min(corrcoef([ConsensusMatrixMA(1,:)',ConsensusMatrixMA(2,:)'],'rows','complete'),[],'all'));
    bias_1_1b = mean(abs(ConsensusMatrixMA(1,:)-ConsensusMatrixMA(2,:)),2,'omitnan');
    cvMA(1,1) = corr_1_1b/bias_1_1b;

    corr_1_2 = abs(min(corrcoef([ConsensusMatrixMA(1,:)',ConsensusMatrixMA(3,:)'],'rows','complete'),[],'all'));
    bias_1_2 = mean(abs(ConsensusMatrixMA(1,:)-ConsensusMatrixMA(3,:)),2,'omitnan');
    cvMA(1,2) = corr_1_2/bias_1_2;

    corr_1_3 = abs(min(corrcoef([ConsensusMatrixMA(1,:)',ConsensusMatrixMA(4,:)'],'rows','complete'),[],'all'));
    bias_1_3 = mean(abs(ConsensusMatrixMA(1,:)-ConsensusMatrixMA(4,:)),2,'omitnan');
    cvMA(1,3) = corr_1_3/bias_1_3;

    corr_1b_2 = abs(min(corrcoef([ConsensusMatrixMA(2,:)',ConsensusMatrixMA(3,:)'],'rows','complete'),[],'all'));
    bias_1b_2 = mean(abs(ConsensusMatrixMA(2,:)-ConsensusMatrixMA(3,:)),2,'omitnan');
    cvMA(1,4) = corr_1b_2/bias_1b_2;

    corr_1b_3 = abs(min(corrcoef([ConsensusMatrixMA(2,:)',ConsensusMatrixMA(4,:)'],'rows','complete'),[],'all'));
    bias_1b_3 = mean(abs(ConsensusMatrixMA(2,:)-ConsensusMatrixMA(4,:)),2,'omitnan');
    cvMA(1,5) = corr_1b_3/bias_1b_3;

    corr_2_3 = abs(min(corrcoef([ConsensusMatrixMA(3,:)',ConsensusMatrixMA(4,:)'],'rows','complete'),[],'all'));
    bias_2_3 = mean(abs(ConsensusMatrixMA(3,:)-ConsensusMatrixMA(4,:)),2,'omitnan');
    cvMA(1,6) = corr_2_3/bias_2_3;

    %% flag high bias images
    if bias_1_1b > 17 ||  bias_1_2 > 17 ||  bias_1_3 > 17 || ...
            bias_1b_2 > 17 ||  bias_1b_3 > 17 ||  bias_2_3 > 17
        flag_high_bias = 1;
    end

    %% merge profiles
    if sum(ConsensusMatrixLI(1,:)-ConsensusMatrixLI(2,:),"all","omitnan") ~= 0

        LI = mean(ConsensusMatrixLI(1:2,:),1,"omitnan");
        xLI = find(~isnan(LI));
        consensuLI = [xLI;movmean(LI(~isnan(LI)),5)];

        MA = mean(ConsensusMatrixMA(1:2,:),1,"omitnan");
        xMA = find(~isnan(MA));
        consensuMA = [xMA;movmean(MA(~isnan(MA)),5)];

    elseif sum(sum(ConsensusMatrixLI(2:end,:),2,"omitnan")>0) > 1
        LI = mean(ConsensusMatrixLI(2:end,:),1,"omitnan");
        xLI = find(~isnan(LI));
        consensuLI = [xLI;movmean(LI(~isnan(LI)),5)];

        MA = mean(ConsensusMatrixMA(2:end,:),1,"omitnan");
        xMA = find(~isnan(MA));
        consensuMA = [xMA;movmean(MA(~isnan(MA)),5)];

    else
        consensuLI = [];
        consensuMA = [];
    end


%     %% merge MA
%     if sum(ConsensusMatrixMA(1,:)-ConsensusMatrixMA(2,:),"all","omitnan") ~= 0
%         
%         MA = mean(ConsensusMatrixMA(1:2,:),1,"omitnan");
%         xMA = find(~isnan(MA));
%         consensuMA = [xMA;movmean(MA(~isnan(MA)),5)];
% 
%     elseif sum(sum(ConsensusMatrixMA(2:end,:),2,"omitnan")>0) > 1
%         MA = mean(ConsensusMatrixMA(2:end,:),1,"omitnan");
%         xMA = find(~isnan(MA));
%         consensuMA = [xMA;movmean(MA(~isnan(MA)),5)];
% 
%     else
%         consensuMA = [];
%     end

end
