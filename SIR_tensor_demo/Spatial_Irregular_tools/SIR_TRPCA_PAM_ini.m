    function  [Xrec,iter,chgrec,Xi] = SIR_TRPCA_PAM_ini(Xb,lambda,tau,mu,Tol1,MaxIt1,prosize)

    Yi = Irunfold(Xb);
    rho = 0.001;
    gamma = 1/tau-rho;
   tic;
   % size of latent tensor
   h = prosize(1);
   w = prosize(2);
   b= h*w;
   Si = zeros(size(Yi));
   B = size(Xb,3);
   Xi = Yi-Si; 
   K = size(Xi,1); 
   
   % initial transform
   if K >=b
   Ditemp = dctmtx(K);
   Di =  Ditemp(:,1:b)'; 
   disp(['is incomplet dictionary ']) 
   else
   Ditemp = dctmtx(b);
   Di =  Ditemp(:,1:K); 
   disp(['is overcomplet dictionary '])   
   end
   Gi = Di;
   Zi = Di*Xi;
   start_iter=1000;
   start_iter2=0;
   max_mu = 1e10;
   max_gamma = 1e10;
   iter = 0;
   stop = 0;
   while stop==0
   oldXi =    Xi;
   oldSi =    Si;
   oldZi =    Zi;
   oldDi = Di;
   % update Zi
   [Zi,~] = prox_tnn(reshape((gamma*(Di*Xi)+rho*oldZi)/(gamma+rho),[h w B]),1/(gamma+rho));
   Zi = reshape(Zi,[h*w B]);
   % update Di
   if iter >= start_iter
   Di = (gamma*Zi*Xi'+rho*oldDi)/(gamma*(Xi*Xi')+rho*eye(size(Xi,1)));
   Di =    Di./sqrt(sum(Di.^2,2));
   end
   % update Xi
   %Xi = (gamma*Di'*(Zi)+mu*(Yi-Si)+rho*oldXi)/(mu+rho+gamma);
   if iter >= start_iter2
   %Xi = pinv((mu+rho)*eye(size(Di,2))+gamma*(Di'*Di))*(gamma*Di'*(Zi)+mu*(Yi-Si)+rho*oldXi);
   Xi = ((mu+rho)*eye(size(Di,2))+gamma*(Di'*Di))\(gamma*Di'*(Zi)+mu*(Yi-Si)+rho*oldXi);
   end
   % update Si
   Si = prox_l1((mu*(Yi-Xi)+rho*oldSi)/(mu+rho),lambda/(mu+rho));
   iter = iter+1;
   % update mu gamma 
   mu = min(1.05*mu,max_mu);    
   gamma = min(1.05*gamma,max_gamma); 

   % Calculate the new fit to the model and decide to stop or not
     chgX = norm(Xi(:)-oldXi(:))/norm(oldXi(:));
     chgS = norm(Si(:)-oldSi(:))/norm(oldSi(:));
     chg = max([chgX chgS]);
     disp(['at iteration ',num2str(iter), '\ chg ',num2str(chg)]) 
     chgrec(iter) = chgX;
     if iter >10
     if chg<Tol1 || (iter==MaxIt1) 
        stop=1;
     end
     end
   end
      disp(['finish initialization']) 
   Time = toc;
   Xrec = Irfold(Xi,Xb);
end

