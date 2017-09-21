function validate_modularity_Mucha_gamma
    outdat='/datc/flex/data_CompCor/';
    outfig='/datc/flex/fig_CompCor/';
    recpath='/datc/dynNet/code/';
    addpath(genpath('/datc/dynNet/code/GenLouvain2.0'));
    addpath(genpath('/datc/software/brat/ccm')); %% compute clustering coeficient
    %addpath(genpath('/datc/software/SurfStat')); % to use @term function
    addpath(genpath('/datc/dynNet/code/aboxplot'));
    nROI = 264;
    listfile = load([recpath 'hcp_S453_sex_age_avg.txt']);

    nSubj = size(listfile, 1);
    tname = {'rfMRI_REST1_LR','tfMRI_GAMBLING_LR','tfMRI_MOTOR_LR','tfMRI_SOCIAL_LR',  ...
             'tfMRI_EMOTION_LR',  'tfMRI_LANGUAGE_LR',  'tfMRI_RELATIONAL_LR', 'tfMRI_WM_LR'};
    %tname = {'rfMRI_REST1_RL','tfMRI_GAMBLING_RL','tfMRI_MOTOR_RL','tfMRI_SOCIAL_RL',  ...
    %         'tfMRI_EMOTION_RL',  'tfMRI_LANGUAGE_RL',  'tfMRI_RELATIONAL_RL', 'tfMRI_WM_RL'};

    % for the global network reconfiguration
    load('../data_CompCor/WM_RT_Avg_Correct_Trials_avg.mat');
    bk0(isnan(bk0))=0; bk2(isnan(bk2))=0; %[1,2,3,4]=[body, face, places, tools];
    RTime=[mean(bk0,2), mean(bk2,2)];
    %% the following is to investigate the partition difference (measured by 
    %% zscore of Rand Coefficient) is correltated to performance
    preout='validate_globalReconf_zscore_';
    nGa=7;
    if 0
        disp('The following is to investigate the partition difference (measured by ');
        disp('zscore of Rand Coefficient) is correltated to performance');

        zsc=nan(nGa, nSubj,3); % 3: [rest:bk0, rest:bk2, rest:WM]
        Q = nan(nGa, nSubj, 4); % 4: [rest, WM, 0bk, 2bk]
        addpath(genpath('/datc/dynNet/code/NCT_Bassett'));
        for j=1:nGa
            strGa=sprintf('%1.1f', 0.4+(j-1)*0.2);
            disp(strGa);
            QSfile=[outdat 'avg_totQS_rest_WM_block_thr15_gamma' strGa '.mat'];
            load(QSfile); % totQ, totS
            if ~isempty(find(isnan(totQ))) || ~isempty(find(isnan(totS)))
                ttt=0;
            end
            for i=[ 1:nSubj ] %[1:nSubj]
                tmp1=zrand(squeeze(totS(1,i,:)), squeeze(totS(3,i,:)));
                tmp2=zrand(squeeze(totS(1,i,:)), squeeze(totS(4,i,:)));
                tmp3=zrand(squeeze(totS(1,i,:)), squeeze(totS(2,i,:)));
                zsc(j,i,:)=[tmp1; tmp2; tmp3];
                Q(j,i,:)=totQ(:,i);
                if ~isempty(find(isnan(zsc))) || ~isempty(find(isnan(zsc)))
                    ttt=0;
                end
            end
        end
        save([outdat preout 'rest_0bk_2bk_WM.mat'], 'zsc', 'Q');
    end
    load([outdat preout 'rest_0bk_2bk_WM.mat']); %, 'zsc'
    for i=1:nGa
        if ~isempty(find(isnan(zsc(i,:,:)))) || ~isempty(find(isnan(Q(i,:,:))))
            disp(i); %% only gamma=0.4, there is nan in totS
        end
    end

    zsc(isnan(zsc))=0;
    Q(isnan(Q))=0;

    
    %% aboxplot the Q in 4 states [rest, WM, 0bk, 2bk]
    xlab={'rest','WM','0bk','2bk'};
    figure('Position',[100,100, 700,700]); % decrease across the WM loads
    for i=1:3
        for j=1:3
            II=(i-1)*3+j;
            if II<=nGa
                subplot(3,3,II), aboxplot(squeeze(Q(II,:,:))); 
                title(['\gamma' sprintf('=%1.1f', (II-1)*0.2+0.4)],'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',16);
                set(gca,  'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',16);
                set(gca,'xticklabel', xlab, 'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',16);
                rotateXLabels(gca, 45);
            end
        end
    end
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperPosition', [100 100 700 700]);
    print(gcf, '-dtiff', '-r300', [outfig preout 'Q4states.tif']);
    %zsc(287,:)=[]; RTime(287,:)=[]; % 3 sigma outlier
    zsc(:,255,:)=[]; RTime(255,:)=[]; % 4 sigma outlier    
    pMat=nan(nGa,3);
    for i=1:nGa
       [r,p] =corr(squeeze(zsc(i,:,3))', mean(RTime,2)); pMat(i,1)=sign(r)*p;
       [r,p] =corr(squeeze(zsc(i,:,1))', RTime(:,1)); pMat(i,2)=sign(r)*p;
       [r,p] =corr(squeeze(zsc(i,:,2))', RTime(:,2)); pMat(i,3)=sign(r)*p;
    end
    figure;
    plot(pMat(:,1), 'r*'); hold on; plot(pMat(:,2), 'go'); hold on;  plot(pMat(:,3), 'b^'); hold on;
    plot([0:1:8], -ones(1,nGa+2)*0.05, 'c--');
    xlim([0,8]); set(gca, 'xtick',[1:7], 'xticklabel',[0.4:0.2:1.6]);
    xlabel('\gamma in modularity computing','FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',16);
    ylabel('Signed {\it{P}}-value of the correlations','FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',16);
    legend({'rest vs. WM', 'rest vs. 0back', 'rest vs. 2back'}, 'Location','East');
    set(gca,  'FontName', 'Arial', 'FontWeight', 'bold', 'Fontsize',16);
    
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperPosition', [100 100 700 700]);
    print(gcf, '-dtiff', '-r300', [outfig preout 'pVal_corr_behave.tif']);    
end


