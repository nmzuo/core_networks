function pval = permutation_test_paired(x,y, nPerm, Verb)
% x = [m,n]; y = [m,1]
% Compare the mean value between each column of x and y
% The x could be only one colum (e.g. n=1) if you only want to compare two
% vectors
% Calling invprctile(...)
% by nmzuo, Jan. 13, 2016
    [m n] = size(x);
    x = single(x); y = single(y);
    
    if ~exist('nPerm', 'var')
        nPerm = 10000;
    end
    if ~exist('Verb', 'var')
        Verb=0;
    end
    diffRecord=zeros(nPerm,n);
    for i=1:nPerm
       if (mod(i,100) == 0 && Verb ~=0)
           fprintf('%d ', i);
           if (mod(i,1000) == 0 )
               fprintf('\n');
           end
       end
       
       %rInd = randperm(mm);
       randswap=rand(m,1);
       randswap(randswap<0.5)=0; randswap=logical(randswap);
       for j=1:n
           x1=x(:,j); y1=y;
           x1(randswap)=y(randswap); y1(randswap)=x(randswap,j);
           diffRecord(i,j)=mean(x1) - mean(y1);
       end
    end
    
    orgDiff = abs(mean(x,1)-mean(y));
    
    pval=zeros(1,n);
    for i=1:n
        pval(i) = invprctile(diffRecord(:,i),orgDiff(i));
    end
    
    pval = 1.0-pval/100.0;

end
