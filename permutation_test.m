function pval = permutation_test(x,y)
% x = [m1,n]; y = [m2,1]
% Compare the mean value between each column of x and y
% The x could be only one colum (e.g. n=1) if you only want to compare two
% vectors
% Calling invprctile(...)
% by nmzuo, Jan. 13, 2016
    [m1 n] = size(x);
    x = single(x); y = single(y);
    xy = [x; repmat(y,[1,n])]; % so xy=[m1+m2, n]
    clear x y;
    mm = size(xy, 1); % mm=m1+m2
    
    nPerm = 5000;
    diffRecord=zeros(nPerm,n);
    for i=1:nPerm
       if (mod(i,100) == 0)
           fprintf('%d ', i);
           if (mod(i,1000) == 0)
               fprintf('\n');
           end
       end
       
       rInd = randperm(mm);
       %xyRand=xy(rInd, :);
       diffRecord(i,:)=mean(xy(rInd(1:m1),:),1)-mean(xy(rInd(m1+1:mm),:),1);
    end
    
    orgDiff = abs(mean(xy(1:m1,:),1) - mean(xy(m1+1:mm,:),1));
    
    pval=zeros(1,n);
    for i=1:n
        pval(i) = invprctile(diffRecord(:,i),orgDiff(i));
    end
    
    pval = 1.0-pval/100.0;

end