function  [Xi,iter,chgrec] = non_LRMF(Yi,lambda,mu,Tol1,MaxIt1,rank)
   

    
   % initial dictionary 
%    h = ceil(r/ave);
%    w = ceil(c/ave);

  [I,J]=size(Yi);
  Xi = Yi;
  Si = zeros(size(Yi));
   max_mu = 1e10;
   iter = 0;
   stop=0;
   while stop==0
   oldXi =    Xi;
   oldSi =    Si;
   % update Xi
    options.verbose = 1;
    % MU
    options.alg = 'mu';
    [M, ~] = fro_mu_nmf(Yi-Si, rank, options);
    Xi = M.W*M.H;
   % update Si 
      Si = prox_l1(Yi-Xi,lambda/mu);
    
      iter = iter+1;
  
    
       % Calculate the new fit to the model and decide to stop or not
       
        chgX = norm(Xi(:)-oldXi(:))/norm(oldXi(:));
         chgS = norm(Si(:)-oldSi(:))/norm(oldSi(:));
         chg = max([chgX chgS]);
       chgrec(iter) = chgX;   
       % disp(['at iteration ',num2str(iter), '\ chg ',num2str(chg)]) 
        mu = min(1.1*mu,max_mu);     
     if chg<Tol1 || (iter==MaxIt1) 
            stop=1;
     end
     
   end
   Time = toc;
 end

