clear all
close all
clc
addpath(genpath(cd));

rng(0);
name{1} = 'superpixel1';
load('superpixel1.mat')

[M,N,B]=size(Xtrue);

idx = isnan(Xtrue);
figureTrue = Xtrue;
figureTrue (idx) = 1000;
ratio = [0.3];
XTi =Irunfold(Xtrue);
%% add sparse noise
Yi=imnoise(XTi,'salt & pepper',ratio);
Xb =Irfold(Yi,Xtrue);
[error,psnr] = my_metrics(Xb,Xtrue);
idx = isnan(Xb);
figureObs =   Xb;
figureObs (idx) = 1000;
imname=[num2str(name{1}),'_ratio_',num2str(ratio),'_result_Obs_error_',num2str(error,'%.4f'),'_psnr_',num2str(psnr,'%.2f'),'.mat'];
save(imname,'Xb','Xtrue','figureObs','figureTrue'); 

%load('superpixel1_ratio_0.3_result_Obs_error_0.5636_psnr_14.88.mat');
                    
%% NMF(unfold the spatial-irregular tensor to a matrix)
L = size(XTi,2);
for mu1 = [1] % 0.01,0.05,0.1,0.5,1,5,10,50,100
for lambda1 = [0.03] % 0.01 0.03 0.05 0.07 0.09 0.1 0.3 0.5 0.7
for rank = [1] % [1:L] L = size(XTi,2);
 MaxIt1 = 500;
 mu = mu1;
 Tol1 = 1e-6;
 lambda = lambda1;
 idx = isnan(Xb);
 Yi =Irunfold(Xb);

 tic;
 [Xi,iter,chgrec] = non_LRMF(Yi,lambda,mu,Tol1,MaxIt1,rank);
 X_rec = Irfold(Xi,Xtrue);
 Time = toc;
 
 [error,psnr] = my_metrics(X_rec,Xtrue);
 
 disp(['Final error :  ',num2str(error)])
 disp(['Final psnr :  ',num2str(psnr)])
 disp(['Time:  ',num2str(Time)])
 disp('  ')
 figureLRMF =   X_rec;
 figureLRMF(idx) = 1000;
 imname=[num2str(name{1}),'_ratio_',num2str(ratio),'_result_non_LRMF_error_',num2str(error,'%.4f'),'_psnr_',num2str(psnr,'%.2f'),'_Time_',num2str(Time,'%.4f'),'_mu_',num2str(mu),'_lambda_',num2str(lambda),'_rank_',num2str(rank),'.mat'];
 save(imname,'X_rec','Time','error','psnr','mu','lambda','figureLRMF'); 
end
end
end 
   
%% TNN (zeropadding)
for mu1 = [5] % 0.01,0.05,0.1,0.5,1,5,10,50,100       
for lambda1 = [0.03] % 0.01,0.03,0.05,0.07,0.1,0.3,0.5,0.7
 MaxIt1 = 500;
 mu = mu1;
 Tol1 = 1e-6;
 lambda = lambda1;
 idx = isnan(Xb);
 Xpadd = Xb;
 Xpadd(idx) = 0;
 %lambda = 1/sqrt(max(size(Xpadd,1),size(Xpadd,2))*size(Xpadd,3));
 tic;
 X_rec = trpca_tnn(Xpadd,lambda,mu,MaxIt1,Tol1);
 Time = toc;
 X_rec(idx)= NaN;
 [error,psnr] = my_metrics(X_rec,Xtrue);
 disp(['Final error :  ',num2str(error)])
 disp(['Final psnr :  ',num2str(psnr)])
 disp(['Time:  ',num2str(Time)])
 disp('  ')
 figureTNNZ =   X_rec;
 figureTNNZ(idx) = 1000;
 imname=[num2str(name{1}),'_ratio_',num2str(ratio),'_result_TNN(zeropadding)_error_',num2str(error,'%.4f'),'_psnr_',num2str(psnr,'%.2f'),'_Time_',num2str(Time,'%.4f'),'_mu_',num2str(mu),'_lambda_',num2str(lambda),'.mat'];
 save(imname,'X_rec','Time','error','psnr','mu','lambda','figureTNNZ');      
end
end
      
%% TNN (interpolation)
for mu1 = [0.1] % 0.01,0.05,0.1,0.5,1,5,10,50,100  
for lambda1 = [0.1]% 0.01,0.03,0.05,0.07,0.1,0.3,0.5,0.7
 MaxIt1 = 500;
 mu = mu1;
 Tol1 = 1e-6;
 XI =  interpolate_nans(Xb); 
 lambda = lambda1;
 %lambda = 1/sqrt(max(size(XI,1),size(XI,2))*size(XI,3));
 tic;
 X_rec = trpca_tnn(XI,lambda,mu,MaxIt1,Tol1);
 Time = toc;
 X_rec(idx)= NaN;
 [error,psnr] = my_metrics(X_rec,Xtrue);
 disp(['Final error :  ',num2str(error)])
 disp(['Final psnr :  ',num2str(psnr)])
 disp(['Time:  ',num2str(Time)])
 disp('  ')
 figureTNNI =   X_rec;
 figureTNNI(idx) = 1000;
 imname=[num2str(name{1}),'_ratio_',num2str(ratio),'_result_TNN(interpolation)_error_',num2str(error,'%.4f'),'_psnr_',num2str(psnr,'%.2f'),'_Time_',num2str(Time,'%.4f'),'_mu_',num2str(mu),'_lambda_',num2str(lambda),'.mat'];
 save(imname,'X_rec','Time','error','psnr','mu','lambda','figureTNNI'); 
end
end
      
%% TNN (Weighted)
for mu1 = [5] % 0.01,0.05,0.1,0.5,1,5,10,50,100  
for lambda1 = [0.03]% 0.01,0.03,0.05,0.07,0.1,0.3,0.5,0.7
 MaxIt1 = 500;
 mu = mu1;
 Tol1 = 1e-6;
 idx = isnan(Xb);
 idx2 = ~isnan(Xb);
 Weighted_Tens = zeros(size(Xb));
 Weighted_Tens(idx2) = 1;
 %lambda = 1/sqrt(max(size(Xb,1),size(Xb,2))*size(Xb,3));
 lambda = lambda1;
 tic;
 X_rec = trpca_tnn_w(Xb,Weighted_Tens,lambda,mu,MaxIt1,Tol1);
 Time = toc;
 X_rec(idx)= NaN;
 [error,psnr] = my_metrics(X_rec,Xtrue);
 disp(['Final error :  ',num2str(error)])
 disp(['Final psnr :  ',num2str(psnr)])
 disp(['Time:  ',num2str(Time)])
 disp('  ')
 figureTNNW =   X_rec;
 figureTNNW(idx) = 1000;
 imname=[num2str(name{1}),'_ratio_',num2str(ratio),'_result_TNN(Weighted)_error_',num2str(error,'%.4f'),'_psnr_',num2str(psnr,'%.2f'),'_Time_',num2str(Time,'%.4f'),'_mu_',num2str(mu),'_lambda_',num2str(lambda),'.mat'];
 save(imname,'X_rec','Time','error','psnr','mu','lambda','figureTNNW');  
end
end
       
%% SIR-TRPCA 
% initialzation    
 MaxIt1 = 500;
 tau = 5;
 lambda = 0.1;
 mu = 0.1;
 Tol1 = 1e-5;
 prosize(1)= 30;
 prosize(2) = 30;
 % using SIR-TRPCA with fixed transform
 [X_rec_i,~,~,Xi] = SIR_TRPCA_PAM_ini(Xb,lambda,tau,mu,Tol1,MaxIt1,prosize);
 [error,psnr] = my_metrics(X_rec_i,Xtrue);
    
 Xini = Xi;
for pro = 30
for tau = [0.5] %0.07,0.09,0.1,0.3,0.5,0.7,0.9,1,3,5,7,9
for lambda = [0.0003]%0.0001,0.0003,0.0005,0.0007,0.001
for mu = [1] % 0.01,0.1,1,10,100
 prosize(1)= pro;
 prosize(2)= 30;
 tic;
 [X_rec,iter,chgrec,~] = SIR_TRPCA_PAM(Xb,Xini,lambda,tau, mu,Tol1,MaxIt1,prosize);
 Time = toc;
 [error,psnr] = my_metrics(X_rec,Xtrue);
 disp(['Final error :  ',num2str(error)])
 disp(['Final psnr :  ',num2str(psnr)])
 disp(['Time:  ',num2str(Time)])
 disp('  ')
 idx = isnan(Xb);
 figureSIR =   X_rec;
 figureSIR(idx) = 1000;
 imname=[num2str(name{1}),'_ratio_',num2str(ratio),'_result_SIR-TRPCA_error_',num2str(error,'%.4f'),'_psnr_',num2str(psnr,'%.2f'),'_Time_',num2str(Time,'%.4f'),'_tau_',num2str(tau),'_lambda_',num2str(lambda),'_mu_',num2str(mu),'_bn1_',num2str(prosize(1)),'_bn2_',num2str(prosize(2)),'.mat'];
 save(imname,'X_rec','Time','error','error','chgrec','figureSIR');
end
end
end
end

% show visual results
imtool([figureObs,figureTrue,figureLRMF,figureTNNZ,figureTNNI,figureTNNW,figureSIR]);
