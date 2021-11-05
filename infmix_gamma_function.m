function [quantile_matrix_x,GLD_matrix_x,LD_matrix_x,den_pm_store,M_store,alpha_store,mean_store,gini_store]=infmix_gamma_function(income_data)


%load('data_dominance_paper.mat');

Ngrid = 1000;
%y_space = linspace(1,max(income_data),Ngrid);           % discretized domain for density
y_space=logspace(0,4,Ngrid);  
v_alpha = 10;
s_alpha = 10;

y=income_data;
T=length(y);

mu_alpha=2.2;
mu_beta=100;

theta=0.1;
%m0 = 0;
%h0 = 1;
%v0 = 4.5 ;
%s0 = (v0-2);

K = 3;
a = ones(1, K);
S = ones(T, 1);
U = ones(1, T); %initial value does not make sense, need to simulate it first.
M = 1; %only one active regime at beginning

m_gamma=1./random('gam',mu_alpha,1./mu_beta,1,K);
v_gamma=exprnd(1/theta,1,K);
w_gamma = sy_dir(a); % need to write your own random variable generator from Dirichlet distribution.

alpha = 1; %concentration
K = K-1; % the last value is the residual

burn = 1000;               %burn-in sample
nit = 5000;               % total number of iterations
s = burn+nit;

%den_pm_store = zeros(nit, Ngrid);   % density drawn from posterior
%M_store = zeros(nit, 1);
%alpha_store = zeros(nit, 1);
x_store = zeros(nit, 1);
r_gammaprop=100;
rng(123)
num_rep=100000;

p_prop = (0.001:0.001:0.999)';
n_prop = size(p_prop,1);


ii=1;
for i = 1:s
    
    % show iterations
%     if (mod(i, 1000)==0 )
%         fprintf('iteration %d\n', i)
%     end
    
    % step 1: update w and U
    % step 1.1: update w
    % find the n vector to count obs in regimes
    
    Sind = repmat(S, 1, M) == repmat(1:M, T, 1); % T by M matrix
    n = sum(Sind); %n is a row vector
    w_gamma(1:M+1) = sy_dir([n, alpha]); % w is a row vector
    % step 1.2: update U
    U = rand(1, T) .* w_gamma(S); % U is a row vector
    % step 1.3: expand K until umin > residual w_{K'}
    umin = min(U);
    wres = w_gamma(K+1);
    while (umin < wres)
        %stick breaking
        v = betarnd(1, alpha);
        w_gamma(K+1:K+2) = w_gamma(K+1) * [v, 1-v];
        %generate associated theta
        m_gamma(K+1) = 1 ./ random('gam', mu_alpha, 1./ mu_beta);
        v_gamma(K+1) = exprnd(1/theta);
        %update residual value
        wres = w_gamma(K+1); 
        %increase K by 1
        K = K + 1;
    end
    % residual weight's value
    m_gamma(K+2) = 1 ./ random('gam', mu_alpha, 1./ mu_beta);
    v_gamma(K+2) = exprnd(1/theta);

   
    % step 2: update theta: (mu and sigma2)
    Y = (repmat(S, 1, M) == repmat(1:M, T, 1)) .* repmat(y, 1, M);
    logY = (repmat(S, 1, M) == repmat(1:M, T, 1)) .* repmat(log(y), 1, M);
    
    for j=1:M
        %sampling m_gamma
        m_gamma(j)=1/(gamrnd(mu_alpha+n(j)*v_gamma(j),1/(mu_beta+sum(Y(:,j))*v_gamma(j))));
        
        %sampling v_gamma
        
        A_rand=rand(1);
        v_prop=gamrnd(r_gammaprop,v_gamma(j)/r_gammaprop);
        
        logvv1=v_gamma(j)*n(j)*log(v_gamma(j))-n(j)*gammaln(v_gamma(j))-v_gamma(j)*(theta+(sum(Y(:,j))/m_gamma(j))+n(j)*log(m_gamma(j))-sum(logY(:,j)));
        logvc1=v_prop*n(j)*log(v_prop)-n(j)*gammaln(v_prop)-v_prop*(theta+(sum(Y(:,j))/m_gamma(j))+n(j)*log(m_gamma(j))-sum(logY(:,j)));
        logV_v1=log(gampdf(v_gamma(j),r_gammaprop,v_prop/r_gammaprop));
        logv_V1=log(gampdf(v_prop,r_gammaprop,v_gamma(j)/r_gammaprop));
        MH1=exp(logvc1+logV_v1-logvv1-logv_V1);
        C1=min(1,MH1);
        if A_rand<=C1
           v_gamma(j)=v_prop;
        else
           v_gamma(j)=v_gamma(j);
        end
    end
    

    
    % step 3: update S
    % step 3.1: S
    prob = gampdf(repmat(y, 1, K), repmat(v_gamma(1:K), T, 1), repmat(m_gamma(1:K), T, 1)./repmat(v_gamma(1:K), T, 1));
    prob = prob .* (repmat(U', 1, K) < repmat(w_gamma(1:K), T, 1));
    prob = cumsum(prob, 2); % T by K matrix
    S = 1 + sum(repmat(rand(T,1) .* prob(:,K), 1, K) > prob, 2); % S is a column vector
    
    % step 3.2: update active regime number M and relabel S, mu, sigma2, w
    % and K
    ind = repmat(S, 1, K) == repmat(1:K, T, 1); %T by K matrix
    ind = sum(ind) > 0; % logical. 1 means active
    M = sum(ind);
    m_gamma(1:M) = m_gamma(ind);
    v_gamma(1:M) = v_gamma(ind);
    w_gamma(1:M) = w_gamma(ind);
    w_gamma(M+1) = 1 - sum(w_gamma(1:M));
    perm = (1:K) .* ind;
    perm(ind) = (1:M);
    S = (repmat(S, 1, K) == repmat(1:K, T, 1)) * perm';
    K = M;
    
    m_gamma(K+1) = 1 ./ random('gam', mu_alpha, 1./ mu_beta);
    v_gamma(K+1) = exprnd(1/theta);

    %sigma2(K+1) = 1 ./ random('gam', v0/2, 2./ s0);
    %mu(K+1) = normrnd(m0, sqrt(sigma2(K+1) / h0));
    
    % step 4: update alpha
    % step 4.1: auxillliary variable x
    x = betarnd(alpha, T);
    % step 4.2: alpha
    alpha = gamrnd(v_alpha+M, 1 / (s_alpha - log(x)));
    
   
    % results
    if (i>burn) & mod(i,100)==0
        den_pm_store(ii,:) = (sum(repmat(w_gamma(1:K+1), Ngrid, 1) .* gampdf(repmat(y_space', 1, K+1), repmat(v_gamma(1:K+1), Ngrid, 1),repmat(m_gamma(1:K+1), Ngrid, 1)./repmat(v_gamma(1:K+1), Ngrid, 1)), 2))';
        
        M_store(ii,1) = M;
        alpha_store(ii,1) = alpha;
        %x(ii,1) = x;
        
        u1=rand(num_rep,1);
        obs_x=[];
        for j=1:K+1
            %id=u>sum(probs(1:(i-1),:),1) & u<=sum(probs(1:i,:),1);
            z_x=(u1>sum(w_gamma(1:(j-1)))) & (u1<=sum(w_gamma(1:j)));
            n_x=sum(z_x);
            obs_temp=gamrnd(v_gamma(j),m_gamma(j)/v_gamma(j),n_x,1);
            obs_x=[obs_x;obs_temp];
        end
        obs_x=sort(obs_x,1);
        cdf_x=sum(repmat(w_gamma(1:K+1), num_rep, 1) .* gamcdf(repmat(obs_x, 1, K+1), repmat(v_gamma(1:K+1), num_rep, 1),repmat(m_gamma(1:K+1), num_rep, 1)./repmat(v_gamma(1:K+1), num_rep, 1)), 2);
        cdf_x(end,1)=1;
        table_x=[obs_x,cdf_x];
        
        for t=1:n_prop
            indx_x=find(p_prop(t,1)<=table_x(:,2),1,'first');
            if indx_x == 1
               quan_x(t,1) = table_x(indx_x(1,1),1);
            else
               quan_x(t,1) = unifrnd(table_x(indx_x(1,1)-1,1),table_x(indx_x(1,1),1));
            end
        end    
        quantile_matrix_x(ii,:)=quan_x';
        
        mean_x=sum(w_gamma(1:K+1).*m_gamma(1:K+1))-...
            sum(w_gamma(1:K+1).*m_gamma(1:K+1).* gamcdf(repmat(quan_x(1,1), 1, K+1), v_gamma(1:K+1)+1,m_gamma(1:K+1)./v_gamma(1:K+1)));
        
        GLD_x=sum(repmat(w_gamma(1:K+1), n_prop, 1).*repmat(m_gamma(1:K+1), n_prop, 1).* gamcdf(repmat(quan_x, 1, K+1), repmat(v_gamma(1:K+1)+1, n_prop, 1),repmat(m_gamma(1:K+1), n_prop, 1)./repmat(v_gamma(1:K+1), n_prop, 1)), 2)-...
            sum(w_gamma(1:K+1).*m_gamma(1:K+1).* gamcdf(repmat(quan_x(1,1), 1, K+1), v_gamma(1:K+1)+1,m_gamma(1:K+1)./v_gamma(1:K+1)));
        LD_x=GLD_x./mean_x;
        GLD_matrix_x(ii,:) = (GLD_x)';
        LD_matrix_x(ii,:) = (LD_x)';
        
        mean_store(ii,1) = sum(w_gamma(1:K+1).*m_gamma(1:K+1));
        gini_store(ii,1) = samplegini(obs_x);
        ii=ii+1;
        

    end
end


end
%save('infmix_gam1999.mat','Post');
% subplot(2, 1, 1)
% plot(y_space,Post.den_pm,'b--',y_space,f0,'r'); 
% xlim([-7,5])
% 
% subplot(2, 1, 2)
% hist(y,30); xlim([-7,5]);
% 
% 
% 
% plot(Post.M)
% legend('\alpha\sim G(100, 10)')
% 
% plot(Post.alpha)
% plot(Post.x)
% 
% 
% 
