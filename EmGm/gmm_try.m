% GMM code

function varargout = gmm_try(X, K_or_centroids)

      % input X:N-by-D data matrix
      % input K_or_centroids: K-by-D centroids
      
      % set the threshold of convergence
      threshold =  1e-6;
      % get the size of X
     [N, D] = size(X);
       % judge if the input is centroids or number of cluster
       if isscalar(K_or_centroids)
           % choose k centorids randomly
           K = K_or_centroids;
           rnpm = randperm(N); 
           centroids = X(rnpm(1:K),:);
       else   % set the centroids
           K = size(K_or_centroids,1);
           centroids = K_or_centroids;
       end
       
       % set the initial parameters of EM GMM
       [pMiu pPi pSigma] = init_params();
       
       Lprev = -inf; maxiter = 50; iter = 1;
       while iter<=maxiter
           % E-step
           % Px: N-by-K
           
           Px = calc_prob();
           
           % pPi: -by-K     pGamma:N-by-K
           pGamma = Px./repmat(pPi,N,1);
           %
           pGamma = pGamma./repmat(sum(pGamma,2),1,K);
           
           % M-step
           % 
           Nk = sum(pGamma,1);
           pMiu = diag(1./Nk)*pGamma'*X;
           pPi = (1/K) * ones(1,K) ;%Nk/N;
           for kk =  1:K
              Xshift = X - repmat(pMiu(kk,:),N,1);
              pSigma(:,:,kk) = (Xshift'*(diag(pGamma(:,kk))*Xshift))/Nk(kk);
           end
           
           % to see if L is converged
           L = sum(log(Px*pPi'));
           if L-Lprev < threshold
               break;
           end
           Lprev = L;
           iter = iter+1;
           
       end
       
       % output set 
       if nargout == 1 
           varargout = {Px};
       else
           model = [];
           model.Miu = pMiu;
           model.Sigma = pSigma;
           model.Pi = pPi;
           varargout = {Px, model};
       end
       
       function [pMiu pPi pSigma] = init_params()
          pMiu = centroids; % centroids of each cluster
          pPi = zeros(1,K); % probability
          pSigma = zeros(D,D,K); % covariances D-by-D
          
          % (X - pMiu)^2  = X^2 + pMiu^2  -  2*X*pMiu
          %distmat = repmat(sum(X.*X,2),1,K) + repmat(sum(pMiu.*pMiu,2)',N,1) -  2*X*pMiu';
          distmat = zeros(N,K);
          for ii = 1:N
              for jj = 1:K
                distmat(ii,jj) = norm(X(ii,:)-pMiu(jj,:));
              end
          end
          [dummy labels] = min(distmat,[],2); % find min value of each row, and its cluster
          %maybe something wrong here so some empty clusters
          AA = unique(labels);
          for k= 1:K 
              if ismember(k,AA) == 1
              Xk = X(labels == k,:);
              else
                  Xk = X(randperm(N, 31),:);
              end
              pPi(k) = 1/K; %size(Xk,1)/N;
              pSigma(:, :, k) = cov(Xk);
          end 
          
          % add something to make sure that pSigma will not be NaN;
       end
   
       % calculate probability
       function Px = calc_prob()
           Px = zeros(N,K); 
           A = (1e-60) * rand(N,1) .* ones(N,1);
           for k= 1:K
               Xshift = X - repmat(pMiu(k,:),N,1);
               inv_pSigma = inv(pSigma(:,:,k)+diag(repmat(threshold,1,size(pSigma(:,:,k),1))));
               tmp = sum((Xshift*inv_pSigma).*Xshift,2);
               coef = (2*pi)^(-D/2)*sqrt(det(inv_pSigma));
               Px(:,k) = (coef * exp(-0.5*tmp)) + A;
              
           end
       end
         
   
   end

