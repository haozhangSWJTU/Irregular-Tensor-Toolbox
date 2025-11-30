function [error,mpsnr] = my_metrics(X_rec,Xtrue)

 [XTi,oo2] =Irunfold(Xtrue);
 [Xreci,oo2] =Irunfold(X_rec);
 error = norm(Xreci(:)-XTi(:))/norm(XTi(:));

     
 mpsnr  = PSNR(Xreci,XTi,max(XTi(:)));

  
end