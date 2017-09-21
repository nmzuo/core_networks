%% For the Mucha'code: Ordered Multislice Matrix
%% http://netwiki.amath.unc.edu/GenLouvain/GenLouvain
%
% by nmzuo, Aug 2015

function [Q, S] = Mucha_2D(A, T,gamma)
    addpath(genpath('/DATA/239/nmzuo/dynNet/code/NCT_Bassett')); % only for GenLouvain
    
    if all(all(A>=0)) || all(all(A<=0))
        A=abs(A);
    else
        error('validate_Mod: input A cotains both positive and negative.');
    end
    
    if ~exist('gamma', 'var')
        gamma=1.0;
    end
    
    if ~exist('T', 'var')
        T=1;
        [Q, S] = Mucha_2D_onestep(A);
    else
        numNode = size(A, 1);
        Qt=nan(T, 1);
        St=nan(numNode, T);
        for i=1:T
            [Qtmp, Stmp]=Mucha_2D_onestep(A, gamma);
            Qt(i)=Qtmp;
            St(:,i)=Stmp(:);
        end
%         % to generate the nodal association matrix
%         % Ref: Bassett_Chaos2013_Robust detection of dynamic community structure in networks
%         [S2, Q2] = consensus_iterative(St');
%         S = consensus_similarity(S2);
%         Q = validate_Mod(A, S);
        %% based on Vatansever, JN 2015, Default Mode Dynamics for Global Functional Integration
        %% August 31, 2016
        [Q, iQ]=max(Qt);
        S=St(:,iQ);
    end
end

function [Q, S] = Mucha_2D_onestep(A, gamma)
    if ~exist('gamma', 'var')
        gamma=1.0;
    end

    k = full(sum(A));
    twom = sum(k);
    B = full(A - gamma*k'*k/twom);
    [S,Q] = genlouvain(B, 10000, 0); % no verbose outpu
    Q = Q/twom;
end

function vMod = validate_Mod(M, S, gamma) 
    if ~exist('gamma', 'var')
        gamma=1.0;
    end
    
    nNode=length(M);
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

