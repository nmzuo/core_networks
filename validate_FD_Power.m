%% this is to analynize the the flexibility of the node across different mental state
%% begin with the WM (0bk and 2bk) task and its resting

%% 1. Task the average of the LR/RL, based on Reviewer #3, after Fisher-z transform
%% 2. Remove-head and Append-tail to reduce HRF delay, by HeYong-YangYihong PNAS paper
%% 3. Remove the mean activity from the time series
%% Date: Oct. 20, 2016

function  validate_FD_Power

outdat='/datc/flex/data_CompCor/';
outfig='/datc/flex/fig_CompCor/';
recpath='/datc/dynNet/code/';
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

% %%LR
% tname = {'rfMRI_REST1_LR','tfMRI_GAMBLING_LR','tfMRI_MOTOR_LR','tfMRI_SOCIAL_LR',  ...
%          'tfMRI_EMOTION_LR',  'tfMRI_LANGUAGE_LR',  'tfMRI_RELATIONAL_LR', 'tfMRI_WM_LR'};
% load([outdat 'task_blockSeries_LR453.mat']);  
%RL
tname = {'rfMRI_REST1_RL','tfMRI_GAMBLING_RL','tfMRI_MOTOR_RL','tfMRI_SOCIAL_RL',  ...
        'tfMRI_EMOTION_RL',  'tfMRI_LANGUAGE_RL',  'tfMRI_RELATIONAL_RL', 'tfMRI_WM_RL'};
load([outdat 'task_blockSeries_RL453.mat']);  

tshort2= {'Rest', 'Gambling', 'Motor', 'Social', ...
          'Emotion', 'Language', 'Relational', 'WM'};



%% get the resting and entire-WM correlation-matrix
for j=[1:8]  %% j=[1,8], for rest and 7 tasks %% only for resting and WM
    tsp1 = '/data/hcp_S500_defil/';
    tsp2 = ['/MNINonLinear/Results/' tname{j} '/'];
    disp(tname{j});
    fdCount=[];
    for i=1:nSubj
        if mod(i, 50) == 0
            fprintf('%d ', i);
        end
        clear ts_rois FD;
        load([tsp1 num2str(listfile(i)) tsp2  num2str(listfile(i))  '_ts_CompCor.mat']);
        fdCount=[fdCount, FD];
    end
    fprintf('\nFD remove ratio: %0.5f\n', length(find(fdCount>0.5))/length(fdCount(:)));
end %%j=[1,8]

%% output: 
% rfMRI_REST1_RL
% FD remove ratio: 0.01571
% tfMRI_GAMBLING_RL
% FD remove ratio: 0.01099
% tfMRI_MOTOR_RL
% FD remove ratio: 0.02434
% tfMRI_SOCIAL_RL
% FD remove ratio: 0.01562
% tfMRI_EMOTION_RL
% FD remove ratio: 0.01263
% tfMRI_LANGUAGE_RL
% FD remove ratio: 0.01528
% tfMRI_RELATIONAL_RL
% FD remove ratio: 0.02091
% tfMRI_WM_RL
% FD remove ratio: 0.01471


end

