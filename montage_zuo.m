function montage_zuo(X)
    if length(size(X)) == 3
       [nx,ny,nt]=size(X);
       sx=ceil(sqrt(nt));
       sy=sx;
       for i=1:sx
           for j=1:sy
               sCount=(i-1)*sy+j;
               if sCount<=nt
                   subplot(sx,sy,sCount); imshow(squeeze(X(:,:,sCount)),[]); 
                   title(num2str(sCount));
               end
           end
       end
    end

end