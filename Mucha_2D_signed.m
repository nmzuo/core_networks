%% For the Mucha'code: Ordered Multislice Matrix
%% http://netwiki.amath.unc.edu/GenLouvain/GenLouvain
%% For the A+/- rewritting, please refer Bassett et al, PNAS, 2011,
%% "Dynamic reconfiguration of human brain networks during learning"

function [Q, S] = Mucha_2D_signed(A, T, gamma)
    addpath(genpath('/DATA/239/nmzuo/dynNet/code/NCT_Bassett'));
    
    if ~exist('gamma', 'var')
        gamma=1.0;
    end
    
    if ~exist('T', 'var')
        T=1;
        [Q, S] = Mucha_2D_signed_onestep(A);
    else
        numNode = size(A, 1);
        Qt=nan(T, 1);
        St=nan(numNode, T);
        for i=1:T
            [Qtmp, Stmp]=Mucha_2D_signed_onestep(A, gamma);
            Qt(i)=Qtmp;
            St(:,i)=Stmp(:);
        end
%         % to generate the nodal association matrix
%         % Ref: Bassett_Chaos2013_Robust detection of dynamic community structure in networks
%         [S2, Q2] = consensus_iterative(St');
%         S = consensus_similarity(S2);
%         if all(all(A>=0)) || all(all(A<=0))
%             disp('WARNING: validate_Mod in a unsigned matrix!');
%             Q = validate_Mod(A, S);
%         else
%             Q = validate_Mod_signed(A, S);
%         end
        %% based on Vatansever, JN 2015, Default Mode Dynamics for Global Functional Integration
        %% August 31, 2016
        [Q, iQ]=max(Qt);
        S=St(:,iQ);
    end
end

function [Q, S] = Mucha_2D_signed_onestep(A, gamma)
    if ~exist('gamma', 'var')
        gamma=1.0;
    end
%     k = full(sum(A));
%     twom = sum(k);
%     B = full(A - gamma*k'*k/twom);
    if all(all(A>=0)) || all(all(A<=0))
        [Q,S] = Mucha_2D(abs(A));
    else
        Ap = zeros(size(A)); An = zeros(size(A));
        Ap(A>0) = A(A>0); An(A<0) = -A(A<0); % the negative - is necessary!
        kp = full(sum(Ap)); kn = full(sum(An));
        twomp = sum(kp); twomn = sum(kn);
        B = full(A - (gamma*kp'*kp/twomp - gamma*kn'*kn/twomn));
        
        [S,Q] = genlouvain(B, 10000, 0); % no verbose output
        Q = Q/(twomp+twomn);
    end
end

function vMod = validate_Mod_signed(M, S, gamma) 
    if ~exist('gamma', 'var')
        gamma=1.0;
    end
    
    nNode=length(M);
    vMod = 0.0;
    Mp = M .* (M>0);    Mn = - M.*(M<0);  % the negative - is necessary!
    twomp = sum(Mp(:)); twomn = sum(Mn(:)); %% Count all the stubs within A or B and between A&B
    
    for i=1:nNode
        for j=1:nNode
            if S(i) == S(j)
                vMod = vMod + M(i,j)- (gamma*sum(Mp(i,:))*sum(Mp(j,:))/twomp - gamma*sum(Mn(i,:))*sum(Mn(j,:))/twomn);
            end
        end
    end
    vMod = vMod/(twomp+twomn);
end

function vMod = validate_Mod(M, S) 
    nNode=length(M);
    gamma = 1.0;
    vMod = 0.0;
    twom = sum(M(:)); 
    
    for i=1:nNode
        for j=1:nNode
            if S(i) == S(j)
                vMod = vMod + M(i,j)- gamma*sum(M(i,:))*sum(M(j,:))/twom;
            end
        end
    end
    vMod = vMod/twom;
end

