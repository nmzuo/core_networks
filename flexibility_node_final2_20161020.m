%% this is to analynize the the flexibility of the node across different mental state
%% begin with the WM (0bk and 2bk) task and its resting

%% 1. Task the average of the LR/RL, after Fisher-z transform
%% 2. Remove-head and Append-tail to reduce HRF delay, by HeYong-YangYihong PNAS paper
%% 3. Remove the mean activity from the time series
%% Date: Oct. 20, 2016

function  flexibility_node_final2_20161020

outdat='/datc/flex/data_CompCor/';
outfig='/datc/flex/fig_CompCor/';
recpath='/datc/dynNet/code/';
addpath(genpath(recpath));
addpath(genpath('/datc/dynNet/code/GenLouvain2.0'));
addpath(genpath('/datc/software/brat/ccm')); %% compute clustering coeficient
%addpath(genpath('/datc/software/SurfStat')); % to use @term function
addpath(genpath('/datc/dynNet/code/aboxplot'));

% roi_rm = load('roi264_Power_uncertain.txt'); % 28 uncertain ROIs
% roiInd = setdiff([1:264], roi_rm);
% nROI = length(roiInd);
nROI = 264;
%Rest: 1200
%WM:   405
%LANGUAGE: 316
%MOTOR: 284
listfile = load([recpath 'hcp_S453_sex_age_avg.txt']);

nSubj = size(listfile, 1);

%%LR
tname = {'rfMRI_REST1_LR','tfMRI_GAMBLING_LR','tfMRI_MOTOR_LR','tfMRI_SOCIAL_LR',  ...
         'tfMRI_EMOTION_LR',  'tfMRI_LANGUAGE_LR',  'tfMRI_RELATIONAL_LR', 'tfMRI_WM_LR'};
load([outdat 'task_blockSeries_LR453.mat']);  
% %RL
% tname = {'rfMRI_REST1_RL','tfMRI_GAMBLING_RL','tfMRI_MOTOR_RL','tfMRI_SOCIAL_RL',  ...
%         'tfMRI_EMOTION_RL',  'tfMRI_LANGUAGE_RL',  'tfMRI_RELATIONAL_RL', 'tfMRI_WM_RL'};
% load([outdat 'task_blockSeries_RL453.mat']);  

tshort2= {'Rest', 'Gambling', 'Motor', 'Social', ...
          'Emotion', 'Language', 'Relational', 'WM'};
     
corrR_rest_WM = nan(2,nSubj,nROI,nROI); % rest, WM
corrR_WM_block = nan(2,nSubj,nROI,nROI); % bk0, bk2 in WM
logTriu = logical(triu(ones(nROI,nROI),1));


%% get the block-WM (0bk and 2bk correlation matrix)
TR=0.72; % TR=0.72s for HCP fMRI data
%% discard the first 8s (~0.72*10) and append 4s (~0.72*5)
%% concatenate the similar conditions into one longer time series, Stewart H. Mostofsky, Brain 2009

load('../data_CompCor/WM_RT_Avg_Correct_Trials_avg.mat');
bk0(isnan(bk0))=0; bk2(isnan(bk2))=0; %[1,2,3,4]=[body, face, places, tools];
RTime=[mean(bk0,2), mean(bk2,2)];

% % disp('Reading original connectivity matrix ...');
% load([outdat 'revise_avg_corrR_all_remMActivation.mat']);  % corrR_rest_WM = nan(2,nSubj,nROI,nROI);
% corrR_rest_WM = nan(2,nSubj,nROI,nROI); % only for rest and WM
% corrR_rest_WM(1,:,:,:)=corrR_all(1,:,:,:); 
% corrR_rest_WM(2,:,:,:)=corrR_all(8,:,:,:); 
% clear corrR_all;
% load([outdat 'corrR_WM_block_original_avg.mat']); % corrR_WM_block = nan(2,nSubj,nROI,nROI);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% how to measure the distance (for NCD)
corr_dist = 1;
hamming_dist = 0;
%% test different density threshold, thr=15,10,5 .
thr = 5;
if hamming_dist == 1
    preout=['hamming_thr' num2str(thr) '_'];
elseif corr_dist == 1
    preout=['thr' num2str(thr) '_'];
end
%preout=[]; % no threshold, gvcN and NCD_3contrat_1mCORR.mat, pcoeff only for binary
% %
% addpath(genpath('/datc/dynNet/code'));  % network_thr_bin()
% disp('Threshold and binarize ...');
% corrR_rest_WM(1,:,:,:) = network_thr_bin(squeeze(corrR_rest_WM(1,:,:,:)), thr);
% corrR_rest_WM(2,:,:,:) = network_thr_bin(squeeze(corrR_rest_WM(2,:,:,:)), thr);
% corrR_WM_block(1,:,:,:) = network_thr_bin(squeeze(corrR_WM_block(1,:,:,:)), thr);
% corrR_WM_block(2,:,:,:) = network_thr_bin(squeeze(corrR_WM_block(2,:,:,:)), thr);

% Compute the modularity and partitions (by Mucha's method)
if 0
    addpath(genpath('/datc/dynNet/code'));
    totQ=nan(2, nSubj); % [0bk, 2bk]*nSubj
    totS=nan(2, nSubj, nROI);
    for i=[1:nSubj] % 1:nSubj
        fprintf('%d ', i);
        [Q, S] = Mucha_2D(squeeze(corrR_WM_block(1,i,:,:)), 100);
        totQ(1,i)=Q;
        totS(1,i,:)=S;
        clear Q S
        [Q, S] = Mucha_2D(squeeze(corrR_WM_block(2,i,:,:)), 100);
        totQ(2,i)=Q;
        totS(2,i,:)=S;
        clear Q S
    end
    fprintf('\n');
    save([outdat preout 'avg_totQS_WM_block_bin.mat'], 'totQ', 'totS');   
end
load([outdat preout 'avg_totQS_WM_block_bin.mat']); % totQ=nan(2, nSubj);  totS=nan(2, nSubj, nROI);
%rmQS=load([outdat 'revise_avg_totQS_all_bin_thr5.mat']);

nNet=10;
%% predefine the subNet names and 13 net index
netName={'Sensory', 'CON', 'Auditory','DMN', 'VAN', 'Visual', ...
         'FPN', 'Salience', 'DAN', 'Subcort.', 'Memory','Cerebellar', 'Uncertain'};
Sensory=[13:46,255]; % 35
Cingulo=[47:60]; % 14
Auditory=[61:73]; % 13
DMN=[74:83,86:131,137,139]; % 58
Memory=[133:136,221]; % 5
Vatt=[138,235:242]; % 9
Visual=[143:173]; % 31
FPC=[174:181,186:202]; % 25
Salience=[203:220]; % 18
Subcort=[222:234]; % 13
Cerebellar=[243:246]; % 4
Datt=[251:252,256:264]; % 11
Uncertain=[ 1:12, 84:85,132, 140:142, 182:185, 247:250, 253:254]; % 28
powerPart=zeros(nROI,1); %% partition assignment of Power 264
powerPart(Sensory)=1;
powerPart(Cingulo)=2;
powerPart(Auditory)=3;
powerPart(DMN)=4;
powerPart(Vatt)=5;
powerPart(Visual)=6;
powerPart(FPC)=7;
powerPart(Salience)=8;
powerPart(Datt)=9;

powerPart(Subcort)=10;
powerPart(Memory)=11;
powerPart(Cerebellar)=12;
powerPart(Uncertain)=13;

%% to check the correlation with the behaviour performance
s453 = load([recpath 'hcp_S453_sex_age_avg.txt']);
behave=load('S500_Selected_Cognition.txt');
behave453=zeros(size(s453,1),size(behave,2));
for i=1:size(s453,1)
    behave453(i,:) = behave(find(behave(:,1)==s453(i,1)),:);
end
behave453(isnan(behave453))=0; % remove nan;

Colormap=[0.8 0.8 0.8];

%% the following is to examine the GVC(cole) with the behaviroul performance
%% including the out-scanner (18 scales) and in-scanner (RT)
if 0
    disp('Computing GVC(Cole) matrix ...')
    gvcN=nan(nSubj,nROI);
    for i=[1:nSubj ] % 1:nSubj
        if mod(i, 50) == 0
                fprintf('%d ', i);
        end
        matState=nan(nROI, nROI,3);
        matState(:,:,1)=corrR_rest_WM(1,i,:,:);
        matState(:,:,2)=corrR_WM_block(1,i,:,:);
        matState(:,:,3)=corrR_WM_block(2,i,:,:);
        gvcNode=gvc_cole(matState);
        gvcN(i,:)=gvcNode;
    end
    fprintf('\n');
    save([outdat preout 'gvcN.mat'], 'gvcN'); 
    load([outdat preout 'gvcN.mat']);% , 'gvcN'); 

    repeatCole=nan(nSubj,nNet);
    for i=1:nNet
        repeatCole(:,i)=mean(gvcN(:,powerPart==i),2);
    end

    figure('Position',[100,100, 700,500]), aboxplot(repeatCole, 'labels',netName(1:nNet), 'Colormap',Colormap); 
    set(gca,  'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',16);
    xticklabel_rotate([], 45, [], 'fontsize', 16);
    ylabel('Nodal GVC (averaged in networks)', 'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',16);
    %title('GVC across subjects for each ROI (across 3 states) repeatCole');
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperPosition', [100 100 700 500]);
    print(gcf, '-dtiff', '-r300', [outfig preout 'gvc_13net_across3states_repeatCole.tif']);
    %%above results: 13 net,no significant for DMN and FPC; 10 net: DMN is
    %%smallest but not significant for FPC
end

%% the following is to examine the change of flexibility of node aross states,
%% rest, 0bk, 2bk, by participation coefficient (Sporns BCT toolbox)
if 0 % pcoeff based on Power 14-network partition
    disp('Computing pcoeff matrix ...');
    addpath(genpath('/datc/software/BCT/BCT_20150125'));
    pcoeff=nan(nSubj,nROI,3); % 3 states: rest, 0bk, 2bk
    for i= [  ] % 1:nSubj
        if mod(i,50) ==0
            fprintf('%d  ', i);
        end
        % rest
        P=participation_coef(abs(squeeze(corrR_rest_WM(1,i,:,:))), powerPart);
        pcoeff(i,:,1)=P(:);  clear P;
        % 0bk
        P=participation_coef(abs(squeeze(corrR_WM_block(1,i,:,:))), powerPart);
        pcoeff(i,:,2)=P(:);  clear P;
        % 2bk
        P=participation_coef(abs(squeeze(corrR_WM_block(2,i,:,:))), powerPart);
        pcoeff(i,:,3)=P(:);  clear P;
    end
    fprintf('\n');
%    save([outdat preout 'pcoeff.mat'], 'pcoeff');
    load([outdat preout 'pcoeff.mat']); %, 'pcoeff'); % [nSubj,nROI,3]
end
if 1 % pcoeff based on Individual partition stored in totQ/totS
    disp('Computing pcoeff matrix ...');
    addpath(genpath('/datc/software/BCT/BCT_20150125'));
    pcoeff=nan(nSubj,nROI,3); % 3 states: rest, 0bk, 2bk
    for i= [  ] % 1:nSubj
        if mod(i,50) ==0
            fprintf('%d  ', i);
        end
        % rest
        P=participation_coef(abs(squeeze(corrR_rest_WM(1,i,:,:))), squeeze(rmQS.totS(1,i,:)));
        pcoeff(i,:,1)=P(:);  clear P;
        % 0bk
        P=participation_coef(abs(squeeze(corrR_WM_block(1,i,:,:))), squeeze(totS(1,i,:)));
        pcoeff(i,:,2)=P(:);  clear P;
        % 2bk
        P=participation_coef(abs(squeeze(corrR_WM_block(2,i,:,:))), squeeze(totS(2,i,:)));
        pcoeff(i,:,3)=P(:);  clear P;
    end
    fprintf('\n');
%    save([outdat preout 'pcoeff_individualPartition.mat'], 'pcoeff');
    load([outdat preout 'pcoeff_individualPartition.mat']); %, 'pcoeff'); % [nSubj,nROI,3]
end


%% the following is to test nodal connection distribution (NCD)
%% this is similar to GVC by cole, but SHOULD be more precise
%% diversity of the node across different mental states
if 1
    disp('Computing nodal connection distribution (NCD) ...');
    ncd12=nan(nSubj,nROI); % rest vs. 0bk % connectivity are all positive
    ncd13=nan(nSubj,nROI); % rest vs. 2bk
    ncd23=nan(nSubj,nROI); % 0bk vs. 2bk
    for i=[   ] % 1:nSubj
        if mod(i,50) ==0
            fprintf('%d  ', i);
        end
        
        if corr_dist == 1
            for j=1:nROI
                [rho, pval]=corr(reshape(corrR_rest_WM(1,i,j,:), nROI,1), reshape(corrR_WM_block(1,i,j,:), nROI,1));
                %%%%%%%%%%%%%%%  corrR_all(1,4,78,:) are all 0
                ncd12(i,j)=1-abs(rho);
                [rho, pval]=corr(reshape(corrR_rest_WM(1,i,j,:), nROI,1), reshape(corrR_WM_block(2,i,j,:), nROI,1));
                ncd13(i,j)=1-abs(rho);
                [rho, pval]=corr(reshape(corrR_WM_block(1,i,j,:), nROI,1), reshape(corrR_WM_block(2,i,j,:), nROI,1));
                ncd23(i,j)=1-abs(rho);
            end
        elseif hamming_dist == 1
            for j=1:nROI
                warning off; % pdist2/hamming is stupid to accept double
                D=pdist2(reshape(corrR_rest_WM(1,i,j,:), 1,nROI), reshape(corrR_WM_block(1,i,j,:), 1,nROI), 'hamming');
                %%%%%%%%%%%%%%%  corrR_all(1,4,78,:) are all 0
                ncd12(i,j)=D;
                D=pdist2(reshape(corrR_rest_WM(1,i,j,:), 1,nROI), reshape(corrR_WM_block(2,i,j,:), 1,nROI), 'hamming');
                ncd13(i,j)=D;
                D=pdist2(reshape(corrR_WM_block(1,i,j,:), 1,nROI), reshape(corrR_WM_block(2,i,j,:), 1,nROI), 'hamming');
                ncd23(i,j)=D;  
                warning on;
            end
        end
    end
    ncd12(isnan(ncd12))=1; ncd13(isnan(ncd13))=1; ncd23(isnan(ncd23))=1;
    fprintf('\n');
%    save([outdat preout 'NCD_3contrat_1mCORR.mat'], 'ncd12', 'ncd23', 'ncd13');
    load([outdat preout 'NCD_3contrat_1mCORR.mat']); %%, 'ncd12', 'ncd13', 'ncd23');
end

%% the following is to repeat the figure 4b (mean GVC) and figure 5a (mean pcoeff) 
%% in two ways in Cole paper (Nat Neur. 2013) stem bars for the 13 subnetworks
% pcoeff and ncd
conLab={'resting vs. 0back', '0back vs. 2back', 'resting vs. 2back'};
pcstem=nan(nSubj,nNet,3); % 13 subnetwork
ncstem=nan(nSubj,nNet,3); % 13 subnetwork
tmppc=cat(3, squeeze(pcoeff(:,:,1)-pcoeff(:,:,2)),squeeze(pcoeff(:,:,2)-pcoeff(:,:,3)), squeeze(pcoeff(:,:,1)-pcoeff(:,:,3)));
tmppc=abs(tmppc);
tmpnc=cat(3,ncd12,ncd23,ncd13);
for i=1:nNet %% DMN is the biggest
    pcstem(:,i,:)=nanmean(tmppc(:,powerPart==i,:), 2);
    ncstem(:,i,:)=nanmean(tmpnc(:,powerPart==i,:), 2);
end
%save([outdat preout 'NCD_pcoeff_13net.mat'], 'pcstem', 'ncstem');
if 1
    figure('Position',[100,100, 500,750]);
    for i=1:2 %% FPC is the smallest
        subplot(2,1,i);aboxplot(squeeze(pcstem(:,1:nNet,i)), 'labels', netName(1:nNet), 'Colormap',Colormap); 
        set(gca,  'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',16);
        title(['PCD: ' conLab{i}]);
        %xticklabel_rotate([], 45, [], 'fontsize', 14);   
        rotateXLabels(gca, 45);
        if i==2 || i==1
            ylabel('PCD (averaged in networks)', 'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',16);
        end
    end
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperPosition', [100 100 500 750]);
    print(gcf, '-dtiff', '-r300', [outfig preout 'repeatCole_13net_across3states_pcoeff_individualPartition.tif']);
end
[h,p,c,s]=ttest(repmat(squeeze(pcstem(:,4,2)), [1,10]), squeeze(pcstem(:,:,2))); dispn(sign(s.tstat).*p,6);
figure('Position',[100,100, 500,750]);
for i=1:2
    subplot(2,1,i);aboxplot(squeeze(ncstem(:,1:nNet,i)), 'labels', netName(1:nNet), 'Colormap',Colormap); 
    set(gca,  'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',16);
    title(['NCD: ' conLab{i}]);
    %xticklabel_rotate([], 45, [], 'fontsize', 14);   
    rotateXLabels(gca, 45);
    if i==2 | i==1
        ylabel('NCD (averaged in networks)', 'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',16);
    end
end   
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperPosition', [100 100 500 750]);
print(gcf, '-dtiff', '-r300', [outfig preout 'repeatCole_13net_across3states_NCD.tif']);

% test whether DMN and FPC is significant in each group
disp('test whether DMN and FPC is significant in each group');
p1=nan(nNet,3); p2=nan(nNet,3);
for i=1:nNet
    [h,p,c,s]=ttest(squeeze(pcstem(:,4,1)), squeeze(pcstem(:,i,1)));
    p1(i,1)=sign(s.tstat)*p;
    [h,p,c,s]=ttest(squeeze(pcstem(:,4,2)), squeeze(pcstem(:,i,2)));
    p1(i,2)=sign(s.tstat)*p;
    [h,p,c,s]=ttest(squeeze(pcstem(:,4,3)), squeeze(pcstem(:,i,3)));
    p1(i,3)=sign(s.tstat)*p;
    [h,p,c,s]=ttest(squeeze(ncstem(:,7,1)), squeeze(ncstem(:,i,1)));
    p2(i,1)=sign(s.tstat)*p;
    [h,p,c,s]=ttest(squeeze(ncstem(:,7,2)), squeeze(ncstem(:,i,2)));
    p2(i,2)=sign(s.tstat)*p;
    [h,p,c,s]=ttest(squeeze(ncstem(:,7,3)), squeeze(ncstem(:,i,3)));
    p2(i,3)=sign(s.tstat)*p;
end
disp('p1: DMN'); dispn(p1, 6);
disp('p2: FPC'); dispn(p2, 6);
% the above output is:  %% ALL ARE p<e-6 (thr=15)
% p1'=
%     0.0000    0.0000    0.0000       NaN    0.0000    0.0000    0.0000    0.0000    0.0000   -0.0000    0.0000   -0.0000   -0.0000
%     0.0000    0.0000    0.0000       NaN    0.0000    0.0000    0.0000    0.0000    0.0000   -0.0000    0.0000   -0.0000   -0.0000
%     0.0000    0.0000    0.0000       NaN    0.0000    0.0000    0.0000    0.0000    0.0000   -0.2127    0.0000   -0.0000   -0.0000
% p2'=
%    -0.0000   -0.0000   -0.0000   -0.0000   -0.0000   -0.0000       NaN   -0.0000   -0.0000   -0.0000   -0.0000   -0.0000   -0.0000
%    -0.0000   -0.0000   -0.0000   -0.0000   -0.0000   -0.0040       NaN   -0.0000   -0.0000   -0.0000   -0.0000   -0.0000   -0.0000
%    -0.0000   -0.0000   -0.0000   -0.0000   -0.0000   -0.0289       NaN   -0.0000   -0.0000   -0.0000   -0.0000   -0.0000   -0.0000

%% to print out the pval
%% DMN and FPC, correlation to behave
disp('DMN and FPC, correlation to behave');
[r1,p1]=corr(pcstem(:,4,1), RTime(:,1)); disp(p1*sign(r1)); % p1=0.0418; DMN  % thr=15
[r2,p2]=corr(pcstem(:,4,2), RTime(:,2)); disp(p2*sign(r2)); % p2=-0.1944; DMN
[r3,p3]=corr(ncstem(:,7,1), RTime(:,1)); disp(p3*sign(r3)); % p3=6.5e-6; FPC
[r4,p4]=corr(ncstem(:,7,2), RTime(:,2)); disp(p4*sign(r4)); % p4=2.9e-4; FPC
disp('Test whether the above correlation also occure for other 12 networks');
if 1 % test whether the above correlation also occure for other 12 networks
    %results: pcstem, no significant corr.
    %ncstem for <rest bk0>: DMN: 0.001, Visual: 0.005, FPC: 0.001, Salience: 0.015, Datt: 0.011 
    %ncstem for <rest bk2>: DMN: 0.003, FPC: 0.000, Salience: 0.003
    disp('*************** pcstem: 0bk **************'); % no significant
    for i=1:nNet
        [r,p]=corr(pcstem(:,i,1), RTime(:,1)); fprintf('%s: %0.4f \n', netName{i}, p*sign(r));
    end
    disp('*************** pcstem: 2bk **************'); % no significant
    for i=1:nNet
        [r,p]=corr(pcstem(:,i,3), RTime(:,2)); fprintf('%s: %0.4f \n', netName{i}, p*sign(r));
    end
    disp('*************** ncstem: 0bk **************'); % significant, DMN and FPC are most
    for i=1:nNet
        [r,p]=corr(ncstem(:,i,1), RTime(:,1)); fprintf('%s: %0.4f \n', netName{i}, p*sign(r));
    end
    disp('*************** ncstem: 2bk **************'); % significant, DMN and FPC are most
    for i=1:nNet
        [r,p]=corr(ncstem(:,i,3), RTime(:,2)); fprintf('%s: %0.4f \n', netName{i}, p*sign(r));
    end
end

%% So far, we found: for global connection change (dynamics), DMN is big but FPC is small,
%% the above is for change, now we examine the connectivity strength (node strength)
%% The strength was normalized by node number
disp('Reading original connectivity matrix ...');
load([outdat 'corrR_rest_WM_original_avg.mat']);  % corrR_rest_WM = nan(2,nSubj,nROI,nROI);
load([outdat 'corrR_WM_block_original_avg.mat']); % corrR_WM_block = nan(2,nSubj,nROI,nROI);

disp('Threshold and binarize ...');
corrR_rest_WM(1,:,:,:) = network_thr_bin(squeeze(corrR_rest_WM(1,:,:,:)), thr);
corrR_rest_WM(2,:,:,:) = network_thr_bin(squeeze(corrR_rest_WM(2,:,:,:)), thr);
corrR_WM_block(1,:,:,:) = network_thr_bin(squeeze(corrR_WM_block(1,:,:,:)), thr);
corrR_WM_block(2,:,:,:) = network_thr_bin(squeeze(corrR_WM_block(2,:,:,:)), thr);

restStren=nan(nSubj,nNet); bk0Stren=nan(nSubj,nNet); bk2Stren=nan(nSubj,nNet);
conType=[];
for i=1:nNet
    if 0 % To global
        conType='toGlobal_';
        tmp=sum(squeeze(corrR_rest_WM(1,:,:,:)),3);
        restStren(:,i)=mean(tmp(:,powerPart==i),2);
        tmp=sum(corrR_WM_block,4);
        tmp=mean(tmp(:,:,powerPart==i),3);
        bk0Stren(:,i)=tmp(1,:);
        bk2Stren(:,i)=tmp(2,:);
    else % To others
        conType='toOthers_';
        extrM=squeeze(corrR_rest_WM(1,:,:,:));
        extrM(:,:,powerPart==i)=[];
        tmp=sum(extrM,3);
        restStren(:,i)=mean(tmp(:,powerPart==i),2);
        extrM=corrR_WM_block;
        extrM(:,:,:,powerPart==i)=[];
        tmp=sum(extrM,4);
        tmp=mean(tmp(:,:,powerPart==i),3);
        bk0Stren(:,i)=tmp(1,:);
        bk2Stren(:,i)=tmp(2,:);
    end
end
figure('Position',[100,100, 500,1000]); 
subplot(3,1,1); aboxplot(restStren, 'labels', netName(1:nNet), 'Colormap',Colormap); 
    set(gca,  'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',16);
    title(['Connectivity strength: resting']);
    %xticklabel_rotate([], 45, [], 'fontsize', 14);
    rotateXLabels(gca, 45);
subplot(3,1,2); aboxplot(bk0Stren, 'labels', netName(1:nNet), 'Colormap',Colormap);   
    set(gca,  'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',16);
    title('Connectivity strength: 0back');
    %xticklabel_rotate([], 45, [], 'fontsize', 14);  
    rotateXLabels(gca, 45);
    ylabel('Nodal NCD (averaged in networks)', 'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',16);
subplot(3,1,3); aboxplot(bk2Stren, 'labels', netName(1:nNet), 'Colormap',Colormap); 
    set(gca,  'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',16);
    title('Connectivity strength: 2back');
    %xticklabel_rotate([], 45, [], 'fontsize', 14); 
    rotateXLabels(gca, 45);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperPosition', [100 100 500 1000]);
print(gcf, '-dtiff', '-r300', [outfig preout conType 'NodeStrength_13net_3state.tif']);
%The above output: resting, DMN and FPC is small; bk0 and bk2, FPC is
%biggest and DMN is the smallest (no significance)

%% The following is to test, during bk0 and bk2, since FPC NodeStrength is strong,
%% FPC connect mostly to which network (13) ?
%% Results: FPC connect to Salience and Datt
confpc=nan(nSubj, nNet,3);
for i=1:nNet
    normSize=length(find(powerPart==7))*length(find(powerPart==i));
    tmp=squeeze(corrR_rest_WM(1,:,:,:));
    tmp=squeeze(tmp(:,powerPart==7, powerPart==i)); % 7 is the FPC
    confpc(:,i,1)=sum(sum(tmp,3),2)/normSize;
    tmp=squeeze(corrR_WM_block(1,:,:,:));
    tmp=squeeze(tmp(:,powerPart==7, powerPart==i));
    confpc(:,i,2)=sum(sum(tmp,3),2)/normSize;
    tmp=squeeze(corrR_WM_block(2,:,:,:));
    tmp=squeeze(tmp(:,powerPart==7, powerPart==i));
    confpc(:,i,3)=sum(sum(tmp,3),2)/normSize;
end
figure('Position',[100,100, 500,1000]); 
subplot(3,1,1); aboxplot(squeeze(confpc(:,1:nNet,1)), 'labels', netName(1:nNet), 'Colormap',Colormap);     
    set(gca,  'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',16);    
    %xticklabel_rotate([], 45, [], 'fontsize', 14); 
    rotateXLabels(gca, 45);
    title('Connectivity strength (FPC to others): resting'); 
subplot(3,1,2); aboxplot(squeeze(confpc(:,1:nNet,2)), 'labels', netName(1:nNet), 'Colormap',Colormap);     
    set(gca,  'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',16);   
    %xticklabel_rotate([], 45, [], 'fontsize', 14); 
    rotateXLabels(gca, 45);
    title('Connectivity strength (FPC to others): 0bk');
    ylabel('Connectivity strength', 'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',16);
subplot(3,1,3); aboxplot(squeeze(confpc(:,1:nNet,3)), 'labels', netName(1:nNet), 'Colormap',Colormap);     
    set(gca,  'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',16);    
    %xticklabel_rotate([], 45, [], 'fontsize', 14); 
    rotateXLabels(gca, 45);
    title('Connectivity strength (FPC to others): 2bk');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperPosition', [100 100 500 1000]);
print(gcf, '-dtiff', '-r300', [outfig preout 'conn_NodeStrength_FPCto13net_3state.tif']);
disp('FPC connect mostly to Salience and Datt: whether they are significant ? ');
for J=[7,8,9]
    disp(netName{J});
    for i=1:nNet
        tmp=squeeze(confpc(:,:,1));
        [h1,p1,c1,s1]=ttest(tmp(:,J), tmp(:,i));
        tmp=squeeze(confpc(:,:,2));
        [h2,p2,c2,s2]=ttest(tmp(:,J), tmp(:,i));
        tmp=squeeze(confpc(:,:,3));
        [h3,p3,c3,s3]=ttest(tmp(:,J), tmp(:,i));
        fprintf('%d:%d\t%1.6f\t%1.6f\t%1.6f\n', J,i, p1*sign(s1.tstat), p2*sign(s2.tstat), p3*sign(s3.tstat));
        %% thr=15, all p<e-6
    end
end
fprintf('\n');
% the above output are all p<e-6 (for thr=15); COOL!
% so sort the mean value

disp('print pval for the NodeStrength compare DMN&FPC with others');
disp('To examine whether the DMN has similar significance to FPC ?');
% details of the above results
fprintf('%s\t%s\t\t%s\t\t%s\n', '#', 'resting', '0bk','2bk');
for i=1:nNet
    [h1,p1,c1,s1]=ttest(restStren(:,4), restStren(:,i));
    [h2,p2,c2,s2]=ttest(restStren(:,7), restStren(:,i));
    [h3,p3,c3,s3]=ttest(bk0Stren(:,4), bk0Stren(:,i));
    [h4,p4,c4,s4]=ttest(bk0Stren(:,7), bk0Stren(:,i));
    [h5,p5,c5,s5]=ttest(bk2Stren(:,4), bk2Stren(:,i));
    [h6,p6,c6,s6]=ttest(bk2Stren(:,7), bk2Stren(:,i));
    fprintf('%d\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n', i, p1*sign(s1.tstat), p2*sign(s2.tstat), ...
           p3*sign(s3.tstat), p4*sign(s4.tstat), p5*sign(s5.tstat), p6*sign(s6.tstat));
end
fprintf('\n');
%% The above output: FPC strength is all significant, but DMN is not always
% #	resting         0bk             2bk
% 1	-0.000	-0.000	0.000	0.000	-0.003	0.000
% 2	-0.000	-0.000	-0.196	0.000	-0.149	0.000
% 3	-0.000	-0.000	0.000	0.000	0.905	0.000
% 4	NaN	0.000	NaN	0.000	NaN	0.000
% 5	-0.000	-0.000	-0.000	0.000	-0.000	0.000
% 6	-0.000	-0.000	-0.000	0.000	-0.000	0.000
% 7	-0.000	NaN	-0.000	NaN	-0.000	NaN
% 8	-0.000	-0.001	-0.000	0.000	-0.000	0.000
% 9	-0.000	-0.000	-0.000	0.000	-0.000	0.000
% 10	0.000	0.000	0.000	0.000	0.000	0.000


%% test the correlation between NodeStrength and ReactionTime (WM)
disp('Test the correlation between NodeStrength and ReactionTime (WM)');
[r,p]=corr(restStren, RTime);
disp('rest vs. RTime'); (p.*sign(r))'
[r,p]=corr(bk0Stren, RTime);
disp('0bk vs. RTime'); (p.*sign(r))'
[r,p]=corr(bk2Stren, RTime);
disp('2bk vs. RTime'); (p.*sign(r))'
%% The above output: only FPC/0-2bk, Sensory/0-2bk are all significant;
%% but FPC/Salience/Datt all show negative significant trend, or no correlation.


%% the following is to investigate the partition difference (measured by 
%% zscore of Rand Coefficient) is correltated to performance
if 0
    disp('The following is to investigate the partition difference (measured by ');
    disp('zscore of Rand Coefficient) is correltated to performance');
    QSpath='/datc/flex/data_CompCor/';
    load([QSpath 'revise_avg_totQS_all_bin_thr' num2str(thr) '.mat']);
    restS=squeeze(totS(1,:,:)); % [463,236]
    wmS=squeeze(totS(8,:,:));
    load(['../data_CompCor/thr' num2str(thr) '_avg_totQS_WM_block_bin.mat']);
    bkS=totS;
    clear totQ totS;
    zsc=nan(nSubj,3); % 3: [rest:bk0, rest:bk2, rest:WM]
    addpath(genpath('/datc/dynNet/code/NCT_Bassett'));
    for i=[ 1:nSubj ] %[1:nSubj]
        tmp1=zrand(restS(i,:), squeeze(bkS(1,i,:)));
        tmp2=zrand(restS(i,:), squeeze(bkS(2,i,:)));
        tmp3=zrand(restS(i,:), wmS(i,:));
        zsc(i,:)=[tmp1; tmp2; tmp3];
    end
    save([outdat preout 'zscoreRC_rest_0bk_2bk_WM.mat'], 'zsc');
end
load([outdat preout 'zscoreRC_rest_0bk_2bk_WM.mat']); %, 'zsc'

%zsc(287,:)=[]; RTime(287,:)=[]; % 3 sigma outlier
zsc(255,:)=[]; RTime(255,:)=[]; % 4 sigma outlier

[r1,p1]=corr(zsc(:,1), RTime(:,1)); 
[r2,p2]=corr(zsc(:,2), RTime(:,2)); 
[r3,p3]=corr(zsc(:,3), mean(RTime,2)); 
fprintf('%f; %f; %f;\n', sign(r1)*p1, sign(r2)*p2,sign(r3)*p3);
%% output: -0.000725; -0.106602; -0.000015; % for thr=15
%% this means connection similarity is positively correlated to the performance.
stat3={'WM vs. resting', '0-bk vs. resting', '2-bk vs. resting'};
if thr==15
    figure('Position',[100,100, 1200,400]);
    %% rest vs. WM
    subplot(1,3,1), plot_regress(zsc(:,3), mean(RTime,2))
    ylabel('Reaction Time (ms)', 'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    xlabel({'Network Similarity','resting vs. WM'}, 'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    set(gca,  'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    xlim([0,70]); ylim([500, 1400]); set(gca, 'xtick',[0 20 40 60], 'ytick',[400, 800,1200]);
    logtext={sprintf('R = %0.2f',r3); sprintf('P = %0.2e', p3)};    
    ax=xlim; ay=ylim;
    text(ax(1)+0.55*(ax(2)-ax(1)), ay(1)+0.85*(ay(2)-ay(1)), logtext, ...
        'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);   
    %% rest vs. 0bk
    subplot(1,3,2), plot_regress(zsc(:,1), RTime(:,1))
    ylabel('Reaction Time (ms)', 'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    xlabel({'Network Similarity','resting vs. 0-back'}, 'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    set(gca,  'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    xlim([0,70]); ylim([400, 1400]); set(gca, 'xtick',[0 20 40 60], 'ytick',[400, 800,1200]);
    logtext={sprintf('R = %0.2f',r1); sprintf('P = %0.2e', p1)};    
    ax=xlim; ay=ylim;
    text(ax(1)+0.55*(ax(2)-ax(1)), ay(1)+0.85*(ay(2)-ay(1)), logtext, ...
        'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);   
    %% rest vs. 2bk
    subplot(1,3,3), plot_regress(zsc(:,2), RTime(:,2))
    ylabel('Reaction Time (ms)', 'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    xlabel({'Network Similarity','resting vs. 2-back'}, 'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    set(gca,  'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    xlim([0,60]); ylim([600, 1500]); set(gca, 'xtick',[0 20 40 60], 'ytick',[600, 900,1200]);
    logtext={sprintf('R = %0.2f',r2); sprintf('P = %0.2e', p2)};    
    ax=xlim; ay=ylim;
    text(ax(1)+0.55*(ax(2)-ax(1)), ay(1)+0.85*(ay(2)-ay(1)), logtext, ...
        'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);   
     
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperPosition', [0 0 1200 400]);
    print(gcf, '-dtiff', '-r300', [outfig preout 'plot_NetSim_RT_full_rest_WM_0bk_2bk.tif']);    
elseif thr==10
    figure('Position',[100,100, 1200,400]);
    %% rest vs. WM
    subplot(1,3,1), plot_regress(zsc(:,3), mean(RTime,2))
    ylabel('Reaction Time (ms)', 'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    xlabel({'Network Similarity','resting vs. WM'}, 'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    set(gca,  'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    xlim([0,60]); ylim([500, 1400]); set(gca, 'xtick',[0 20 40 60], 'ytick',[400, 800,1200]);
    logtext={sprintf('R = %0.2f',r3); sprintf('P = %0.2e', p3)};    
    ax=xlim; ay=ylim;
    text(ax(1)+0.55*(ax(2)-ax(1)), ay(1)+0.85*(ay(2)-ay(1)), logtext, ...
        'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);   
    %% rest vs. 0bk
    subplot(1,3,2), plot_regress(zsc(:,1), RTime(:,1))
    ylabel('Reaction Time (ms)', 'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    xlabel({'Network Similarity','resting vs. 0-back'}, 'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    set(gca,  'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    xlim([0,70]); ylim([400, 1400]); set(gca, 'xtick',[0 20 40 60], 'ytick',[400, 800,1200]);
    logtext={sprintf('R = %0.2f',r1); sprintf('P = %0.2e', p1)};    
    ax=xlim; ay=ylim;
    text(ax(1)+0.55*(ax(2)-ax(1)), ay(1)+0.85*(ay(2)-ay(1)), logtext, ...
        'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);   
    %% rest vs. 2bk
    subplot(1,3,3), plot_regress(zsc(:,2), RTime(:,2))
    ylabel('Reaction Time (ms)', 'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    xlabel({'Network Similarity','resting vs. 2-back'}, 'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    set(gca,  'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    xlim([0,45]); ylim([600, 1500]); set(gca, 'xtick',[0 20 40], 'ytick',[600, 900,1200]);
    logtext={sprintf('R = %0.2f',r2); sprintf('P = %0.2e', p2)};    
    ax=xlim; ay=ylim;
    text(ax(1)+0.55*(ax(2)-ax(1)), ay(1)+0.85*(ay(2)-ay(1)), logtext, ...
        'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);   

    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperPosition', [0 0 1200 400]);
    print(gcf, '-dtiff', '-r300', [outfig preout 'plot_NetSim_RT_full_rest_WM_0bk_2bk.tif']);     
elseif thr==5
    figure('Position',[100,100, 1200,400]);     
    %% rest vs. WM
    subplot(1,3,1), plot_regress(zsc(:,3), mean(RTime,2))
    ylabel('Reaction Time (ms)', 'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    xlabel({'Network Similarity','resting vs. WM'}, 'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    set(gca,  'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    xlim([0,70]); ylim([500, 1400]); set(gca, 'xtick',[0 20 40 60], 'ytick',[400, 800,1200]);
    logtext={sprintf('R = %0.2f',r3); sprintf('P = %0.2e', p3)};    
    ax=xlim; ay=ylim;
    text(ax(1)+0.55*(ax(2)-ax(1)), ay(1)+0.85*(ay(2)-ay(1)), logtext, ...
        'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);   
    %% rest vs. 0bk
    subplot(1,3,2), plot_regress(zsc(:,1), RTime(:,1))
    ylabel('Reaction Time (ms)', 'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    xlabel({'Network Similarity','resting vs. 0-back'}, 'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    set(gca,  'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    xlim([0,60]); ylim([400, 1400]); set(gca, 'xtick',[0 20 40 60], 'ytick',[400, 800,1200]);
    logtext={sprintf('R = %0.2f',r1); sprintf('P = %0.2e', p1)};    
    ax=xlim; ay=ylim;
    text(ax(1)+0.55*(ax(2)-ax(1)), ay(1)+0.85*(ay(2)-ay(1)), logtext, ...
        'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);   
    %% rest vs. 2bk
    subplot(1,3,3), plot_regress(zsc(:,2), RTime(:,2))
    ylabel('Reaction Time (ms)', 'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    xlabel({'Network Similarity','resting vs. 2-back'}, 'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    set(gca,  'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);
    xlim([0,45]); ylim([600, 1500]); set(gca, 'xtick',[0 20 40], 'ytick',[600, 900,1200]);
    logtext={sprintf('R = %0.2f',r2); sprintf('P = %0.2e', p2)};    
    ax=xlim; ay=ylim;
    text(ax(1)+0.55*(ax(2)-ax(1)), ay(1)+0.85*(ay(2)-ay(1)), logtext, ...
        'FontName', 'Courier', 'FontWeight', 'bold', 'Fontsize',16);   

    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperPosition', [0 0 1200 400]);
    print(gcf, '-dtiff', '-r300', [outfig preout 'plot_NetSim_RT_full_rest_WM_0bk_2bk.tif']);   
end



end %% end for ALL


function [TK, TK2]= setTicks(xx, numTick)
    x1=xx(1);  x2=xx(2);
    if ~exist('numTick', 'var')
        numTick=2;
    end
    xdist=abs(x2-x1);
    xspan=linspace(x1,x2,numTick+1);
    TK = [];
    if xdist>numTick
        TK=fix(xspan(1)):1:xspan(end);
    elseif xdist/0.5>numTick
        TK=fix(xspan(1)/0.5)*0.5:0.5:xspan(end);
    elseif xdist/0.4>numTick
        TK=fix(xspan(1)/0.4)*0.4:0.4:xspan(end);
    elseif xdist/0.3>numTick
        TK=fix(xspan(1)/0.3)*0.3:0.3:xspan(end);
    elseif xdist/0.2>numTick
        TK=fix(xspan(1)/0.2)*0.2:0.2:xspan(end);
    elseif xdist/0.1>numTick
        TK=fix(xspan(1)/0.1)*0.1:0.1:xspan(end);
    elseif xdist/0.05>numTick
        TK=fix(xspan(1)/0.05)*0.05:0.05:xspan(end);
    end
    
    TK=unique(TK);
  
    TK2=arrayfun( @(x) sprintf('%g',x), TK , 'uniformoutput', false);

    for i=1:length(TK2)
        if abs(str2double(TK2{i}))<1.e-6
            TK2{i}='0';
        end
    end

end
