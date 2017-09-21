function modularity_compute
%% compute the modularity [264, 264] of the WM_block data

outdat='/datc/flex/code/data/';
load('data/corrR_rest_WM_signed.mat'); % corrR_rest_WM: [2   463   264   264]

nSubj=463; 
nROI = 264;
totQ = nan(2, nSubj); totS = nan(2,nSubj,nROI);

for j=1:2  %% j=1:8, for rest and 7 tasks
    for i=1:nSubj
        if mod(i, 50) == 0
            fprintf('%d ', i);
        end
        corrR=squeeze(corrR_rest_WM(j,i,:,:));
        [Q, S] = Mucha_2D(corrR,100);
%        [Q, S] = Mucha_2D(abs(corrR),100);
        totQ(j,i) = Q;
        totS(j,i,:) = S(:);
    end
end
fprintf('\n');

save([outdat 'totQS_FDR_rest_WM_signed.mat'], 'totQ', 'totS');

end