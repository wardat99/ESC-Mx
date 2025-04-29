function [Plabel,Timecost] = ESC_Mx(A,mu,lambda,numClust,pr)
% Enhanced Community Detection in Multiplex Networks via Tensor Weighted
% Schatten p-norm code
% Definition:
%     [Plabel,Timecost]=ExMSC(A,mu,lambda,numClust,pr)
%
% Inputs:
% A               3rd mode tensor [n*n*V], non-negative adjacency tensor.
% mu              scalar, regularization parameter that controls the L_{2,1} norm.
% lambda          scalar, regularization parameter that controls the subspace clustering term.
% numClust        scalar, the number of clusters.
% pr              scalar, the probability of the norm
%
% Outputs:
% Plabel          vector [1*n], the final clustering labels for each node in the netwrok.
% Timecost        scalar, the time cost for each attempt.
%
%
% required external files:
%  [1] wshrinkObj, "On unifying multi-view self-representations for clustering by tensor multi-rank minimization."
%  [2] Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

%%%%
%
% Copyright (C)   <Esraa Al-Sharoa & Mohammad Al-Wardat>
%  
%%%%


[n1,n2,n3]=size(A);
[beta] = getBeta(n3);
I=eye(n1,n2);
Z=eps*ones(size(A));
E=Z; Y1=Z; Y2=Z;  
dY1=Z;    
iter = 1;
tol=0.05;
gamma1=0.01;
gamma2=0.01;
rho = 2;           
Uopt = zeros(n1,numClust);
w_v = eps*ones(n3,1);
max_gamma = 1e10; 
max_iter=50;


tic;
for i=1:max_iter
                                      
        
        % update Q
        Q = reshape(wshrinkObj_weight_lp(Z+Y2/gamma2,beta./gamma2,[n1 n2 n3],0,3,pr),[n1 n2 n3]);
        
        % update E
        
        for v=1:n3
            Emod(:,:,v)=A(:,:,v)-A(:,:,v)*Z(:,:,v)+Y1(:,:,v)/gamma1;
        end
        
        E = l21_operator(Emod,mu,gamma1);
        
        % update Z
        
        
            for v=1:1:n3
            t(:,:,v)=gamma2*I+gamma1*A(:,:,v)'*A(:,:,v);
            B(:,:,v)=A(:,:,v)-E(:,:,v)+Y1(:,:,v)/gamma1;
            Di(:,:,v)=1./getD(Z(:,:,v));
            Di(Di==Inf)=0;
            Z(:,:,v)=inv(t(:,:,v))*(gamma2*Q(:,:,v)-Y2(:,:,v)+ gamma1*A(:,:,v)'*B(:,:,v)...
                +lambda*w_v(v)*Di(:,:,v)*(Uopt*Uopt'));
            end 
        
        
        % make sure Z is symmetric and non-negative
        for v=1:n3
            
            Z(:,:,v)=(Z(:,:,v)+Z(:,:,v)')/2;
        end
        Z(Z<0)=0;
        
        % Compute Z_{rw}
        for v=1:n3
           Dinv(:,:,v)=1./getD(Z(:,:,v)); 
           Dinv(Dinv==Inf)=0;
           Z_rw(:,:,v)=Dinv(:,:,v)*Z(:,:,v);
        end
        
        % Tucker decomposition 
        C = [numClust numClust  1];
        [Uf]= tucker_als(Z_rw,C);  % HOOI
        Uopt = Uf{1};   % Common subspace
        w_v = Uf{3};
        
        % update Y1
        for j=1:n3
           dY1(:,:,j)=A(:,:,j)-A(:,:,j)*Z(:,:,j)-E(:,:,j); 
        end
        Y1 = Y1 + gamma1*dY1;
        
        % update Y2
        dY2=Z-Q;
        Y2 = Y2+gamma2*dY2;
        
        
        err1=norm(dY1(:),'Inf'); err2=norm(dY2(:),'Inf');  % stopping criteria 
        error1(iter)=err1; error2(iter)=err2; 
        err = max(err1,err2);
        
        gamma1 = min(rho*gamma1,max_gamma);
        gamma2 = min(rho*gamma2,max_gamma);
        
        %if mod(iter,10) == 0
        disp(['iter ' num2str(iter) ', err=' num2str(err)])
        %end
    if err < tol
       break;
    end
    iter = iter + 1;
end
Timecost=toc;
      % Compute normalized common subspace
norm_mat = repmat(sqrt(sum(Uopt.*Uopt,2)),1,size(Uopt,2));  
for i=1:size(norm_mat,1)
   if (norm_mat(i,1)==0)
        norm_mat(i,:) = 1;
    end
end
Uf1 = Uopt./norm_mat;
          % Evaluation metrics 
        Plabel = kmeans(Uf1(:,1:numClust),numClust,'replicates',100,'emptyAction','singleton');   
end

