    function  [Xrec,iter,chgrec,Xi] = SIR_TRPCA(Xb, Xini,lambda,tau,mu,Tol1,MaxIt1,prosize)

    Yi = Irunfold(Xb);
    gamma = 1/tau;
tic;
 beta = 100;
 
% size of latent tensor
   h = prosize(1);
   w = prosize(2);
   b= h*w;
   Si = zeros(size(Yi));
   B = size(Xb,3);
if Xini == 0
    Xi = Yi-Si; 
else
    Xi = Xini;
    disp(['using initial ']) 
end 
   K = size(Xi,1); 
   E1 = 0;
   E2 = 0;
   E3 = 0;
  start_iter=0;
  start_iter2 = 3;
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
    % 3.0 ini
%      Ditemp = dctmtx(b);
%      Di =  Ditemp(:,1:K)';    
  

   max_mu = 1e10;
   max_gamma = 1e10;
   max_beta = 1e10;
   iter = 0;
   stop = 0;
   while stop==0
   oldXi =    Xi;
   oldSi =    Si;
   oldZi =    Zi;
   
 
   % update Zi
  %  [Zi,~] = prox_tnn(reshape(Di*Xi-E2/gamma,[h w B]),tau/gamma);
  [Zi,~] = prox_tnn(reshape(Di*Xi-E2/gamma,[h w B]),1/gamma);
   Zi = reshape(Zi,[h*w B]);
   if iter > start_iter
   % update Di
   Di = (gamma*(Zi+E2/gamma)*Xi'+beta*(Gi+E3/beta))/(gamma*(Xi*Xi')+beta*eye(size(Xi,1)));
   % update Gi
   Gi = (Di-E3/beta);
   Gi =    Gi./sqrt(sum(Gi.^2,2));
   end


%    
%      if iter > start_iter
%    % update Di
%    Gi = ((Zi+E2/gamma)*Xi')/((Xi*Xi'));
%    % update Gi
%    Di =    Gi./sqrt(sum(Gi.^2,2));
%    end
   
      if iter > start_iter2
      % update Xi
      Xi = (mu*eye(size(Di,2))+gamma*(Di'*Di))\(gamma*Di'*(Zi+E2/gamma)+mu*(Yi-Si-E1/mu));
      end
      
    % update Si
    Si = prox_l1(Yi-Xi-E1/mu,lambda/mu);
    iter = iter+1;
    % update E1 E2
  
 
   
   E1 = E1 + mu*(Xi+Si-Yi);
   E2 = E2 + gamma*(Zi-Di*Xi);
   E3 = E3 + beta*(Gi-Di);
   % update mu gamma 
     mu = min(1.1*mu,max_mu);    
     gamma = min(1.1*gamma,max_gamma); 
     beta = min(1.1*beta,max_beta); 
      
    
       % Calculate the new fit to the model and decide to stop or not
       
        chgX = norm(Xi(:)-oldXi(:))/norm(oldXi(:));
        chgS = norm(Si(:)-oldSi(:))/norm(oldSi(:));
        chgZ = norm(Zi(:)-oldZi(:))/norm(oldZi(:));
        chg = max([chgX chgS]);
        disp(['at iteration ',num2str(iter), '\ chg ',num2str(chg)]) 
     
       chgrec(iter) = chgX;
       if iter >10
     if chg<Tol1 || (iter==MaxIt1) 
            stop=1;
     end
       end
   end
   Time = toc;
   
   Xrec = Irfold(Xi,Xb);
end

