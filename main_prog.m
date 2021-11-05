num_FPBB=200;
load('income_overall_2001_used.mat');
load('income_overall_2006_used.mat');

income_x=income_overall_2001_CPI;
samplingweight_x=weight_overall_2001;
n_x=length(income_x);

income_y=income_overall_2006_CPI;
samplingweight_y=weight_overall_2006;
n_y=length(income_y);

%parpool(28)
for j=1:num_FPBB
    
    [obs_FPBB_x]=generate_pseudo_representative(income_x,samplingweight_x,n_x);
    [quantile_matrix_x_temp(:,:,j),GLD_matrix_x_temp(:,:,j),LD_matrix_x_temp(:,:,j),den_pm_x_temp(:,:,j),M_x_temp(:,:,j),alpha_x_temp(:,:,j),mean_x_temp(:,:,j),...
        gini_x_temp(:,:,j)]=infmix_gamma_function(obs_FPBB_x);
    
    [obs_FPBB_y]=generate_pseudo_representative(income_y,samplingweight_y,n_y);
    [quantile_matrix_y_temp(:,:,j),GLD_matrix_y_temp(:,:,j),LD_matrix_y_temp(:,:,j),den_pm_y_temp(:,:,j),M_y_temp(:,:,j),alpha_y_temp(:,:,j),...
        mean_y_temp(:,:,j),gini_y_temp(:,:,j)]=infmix_gamma_function(obs_FPBB_y);
    
end
den_pm_x=[];
den_pm_y=[];
M_x=[];
M_y=[];
alpha_x=[];
alpha_y=[];
mean_x=[];
mean_y=[];
gini_x=[];
gini_y=[];
quantile_matrix_x=[];
quantile_matrix_y=[];
GLD_matrix_x=[];
GLD_matrix_y=[];
LD_matrix_x=[];
LD_matrix_y=[];

for j=1:num_FPBB

    den_pm_x=[den_pm_x;den_pm_x_temp(:,:,j)];
    den_pm_y=[den_pm_y;den_pm_y_temp(:,:,j)];

    M_x=[M_x;M_x_temp(:,:,j)];
    M_y=[M_y;M_y_temp(:,:,j)];

    alpha_x=[alpha_x;alpha_x_temp(:,:,j)];
    alpha_y=[alpha_y;alpha_y_temp(:,:,j)];

    mean_x=[mean_x;mean_x_temp(:,:,j)];
    mean_y=[mean_y;mean_y_temp(:,:,j)];
    
    gini_x=[gini_x;gini_x_temp(:,:,j)];
    gini_y=[gini_y;gini_y_temp(:,:,j)];
    
    quantile_matrix_x=[quantile_matrix_x;quantile_matrix_x_temp(:,:,j)];
    quantile_matrix_y=[quantile_matrix_y;quantile_matrix_y_temp(:,:,j)];
    
    GLD_matrix_x=[GLD_matrix_x;GLD_matrix_x_temp(:,:,j)];
    GLD_matrix_y=[GLD_matrix_y;GLD_matrix_y_temp(:,:,j)];
    
    LD_matrix_x=[LD_matrix_x;LD_matrix_x_temp(:,:,j)];
    LD_matrix_y=[LD_matrix_y;LD_matrix_y_temp(:,:,j)];
end
m = size(mean_x,1);
for j = 1:1000
    
    R = (randperm(m))';
    for i = 1:m
         GLD_matrix_y_perm(i,:) = GLD_matrix_y(R(i,1),:);
         LD_matrix_y_perm(i,:) = LD_matrix_y(R(i,1),:);
         quantile_matrix_y_perm(i,:) = quantile_matrix_y(R(i,1),:);
    end
    
    xGLDy = GLD_matrix_x>=GLD_matrix_y_perm;
    yGLDx = GLD_matrix_y_perm>=GLD_matrix_x;
    xGLDy = double(xGLDy);
    yGLDx = double(yGLDx);
    prop_xGLDy(j,:) = mean(xGLDy);
    prop_yGLDx(j,:) = mean(yGLDx);
    overall_xGLDy(j,1) = mean(prod(xGLDy,2));
    overall_yGLDx(j,1) = mean(prod(yGLDx,2));
    overall_xGLDy_20lowest(j,1) = mean(prod(xGLDy(:,1:200),2));
    overall_yGLDx_20lowest(j,1) = mean(prod(yGLDx(:,1:200),2));
    overall_xGLDy_10lowest(j,1) = mean(prod(xGLDy(:,1:100),2));
    overall_yGLDx_10lowest(j,1) = mean(prod(yGLDx(:,1:100),2));
    
    xLDy = LD_matrix_x>=LD_matrix_y_perm;
    yLDx = LD_matrix_y_perm>=LD_matrix_x;
    xLDy = double(xLDy);
    yLDx = double(yLDx);
    prop_xLDy(j,:) = mean(xLDy);
    prop_yLDx(j,:) = mean(yLDx);
    overall_xLDy(j,1) = mean(prod(xLDy,2));
    overall_yLDx(j,1) = mean(prod(yLDx,2));
    overall_xLDy_20lowest(j,1) = mean(prod(xLDy(:,1:200),2));
    overall_yLDx_20lowest(j,1) = mean(prod(yLDx(:,1:200),2));
    overall_xLDy_10lowest(j,1) = mean(prod(xLDy(:,1:100),2));
    overall_yLDx_10lowest(j,1) = mean(prod(yLDx(:,1:100),2));
    
    xFSDy = quantile_matrix_x>=quantile_matrix_y_perm;
    yFSDx = quantile_matrix_y_perm>=quantile_matrix_x;
    xFSDy = double(xFSDy);
    yFSDx = double(yFSDx);
    prop_xFSDy(j,:) = mean(xFSDy);
    prop_yFSDx(j,:) = mean(yFSDx);
    overall_xFSDy(j,1) = mean(prod(xFSDy,2));
    overall_yFSDx(j,1) = mean(prod(yFSDx,2));
    overall_xFSDy_20lowest(j,1) = mean(prod(xFSDy(:,1:200),2));
    overall_yFSDx_20lowest(j,1) = mean(prod(yFSDx(:,1:200),2));
    overall_xFSDy_10lowest(j,1) = mean(prod(xFSDy(:,1:100),2));
    overall_yFSDx_10lowest(j,1) = mean(prod(yFSDx(:,1:100),2));
    
end

mean_prop_xGLDy = mean(prop_xGLDy);
mean_prop_yGLDx = mean(prop_yGLDx);
mean_overall_xGLDy = mean(overall_xGLDy);
max_overall_xGLDy = max(overall_xGLDy);
min_overall_xGLDy = min(overall_xGLDy);
mean_overall_yGLDx = mean(overall_yGLDx);
max_overall_yGLDx = max(overall_yGLDx);
min_overall_yGLDx = min(overall_yGLDx);
mean_overall_xGLDy_20lowest = mean(overall_xGLDy_20lowest);
max_overall_xGLDy_20lowest = max(overall_xGLDy_20lowest);
min_overall_xGLDy_20lowest = min(overall_xGLDy_20lowest);
mean_overall_yGLDx_20lowest = mean(overall_yGLDx_20lowest);
max_overall_yGLDx_20lowest = max(overall_yGLDx_20lowest);
min_overall_yGLDx_20lowest = min(overall_yGLDx_20lowest);
mean_overall_xGLDy_10lowest = mean(overall_xGLDy_10lowest);
max_overall_xGLDy_10lowest = max(overall_xGLDy_10lowest);
min_overall_xGLDy_10lowest = min(overall_xGLDy_10lowest);
mean_overall_yGLDx_10lowest = mean(overall_yGLDx_10lowest);
max_overall_yGLDx_10lowest = max(overall_yGLDx_10lowest);
min_overall_yGLDx_10lowest = min(overall_yGLDx_10lowest);



mean_prop_xLDy = mean(prop_xLDy);
mean_prop_yLDx = mean(prop_yLDx);
mean_overall_xLDy = mean(overall_xLDy);
max_overall_xLDy = max(overall_xLDy);
min_overall_xLDy = min(overall_xLDy);
mean_overall_yLDx = mean(overall_yLDx);
max_overall_yLDx = max(overall_yLDx);
min_overall_yLDx = min(overall_yLDx);
mean_overall_xLDy_20lowest = mean(overall_xLDy_20lowest);
max_overall_xLDy_20lowest = max(overall_xLDy_20lowest);
min_overall_xLDy_20lowest = min(overall_xLDy_20lowest);
mean_overall_yLDx_20lowest = mean(overall_yLDx_20lowest);
max_overall_yLDx_20lowest = max(overall_yLDx_20lowest);
min_overall_yLDx_20lowest = min(overall_yLDx_20lowest);
mean_overall_xLDy_10lowest = mean(overall_xLDy_10lowest);
max_overall_xLDy_10lowest = max(overall_xLDy_10lowest);
min_overall_xLDy_10lowest = min(overall_xLDy_10lowest);
mean_overall_yLDx_10lowest = mean(overall_yLDx_10lowest);
max_overall_yLDx_10lowest = max(overall_yLDx_10lowest);
min_overall_yLDx_10lowest = min(overall_yLDx_10lowest);

mean_prop_xFSDy = mean(prop_xFSDy);
mean_prop_yFSDx = mean(prop_yFSDx);
mean_overall_xFSDy = mean(overall_xFSDy);
max_overall_xFSDy = max(overall_xFSDy);
min_overall_xFSDy = min(overall_xFSDy);
mean_overall_yFSDx = mean(overall_yFSDx);
max_overall_yFSDx = max(overall_yFSDx);
min_overall_yFSDx = min(overall_yFSDx);
mean_overall_xFSDy_20lowest = mean(overall_xFSDy_20lowest);
max_overall_xFSDy_20lowest = max(overall_xFSDy_20lowest);
min_overall_xFSDy_20lowest = min(overall_xFSDy_20lowest);
mean_overall_yFSDx_20lowest = mean(overall_yFSDx_20lowest);
max_overall_yFSDx_20lowest = max(overall_yFSDx_20lowest);
min_overall_yFSDx_20lowest = min(overall_yFSDx_20lowest);
mean_overall_xFSDy_10lowest = mean(overall_xFSDy_10lowest);
max_overall_xFSDy_10lowest = max(overall_xFSDy_10lowest);
min_overall_xFSDy_10lowest = min(overall_xFSDy_10lowest);
mean_overall_yFSDx_10lowest = mean(overall_yFSDx_10lowest);
max_overall_yFSDx_10lowest = max(overall_yFSDx_10lowest);
min_overall_yFSDx_10lowest = min(overall_yFSDx_10lowest);  

quantile_est_x=mean(quantile_matrix_x);
quantile_est_y=mean(quantile_matrix_y);

GLD_est_x=mean(GLD_matrix_x);
GLD_est_y=mean(GLD_matrix_y);

LD_est_x=mean(LD_matrix_x);
LD_est_y=mean(LD_matrix_y);

den_est_x=mean(den_pm_x);
den_est_y=mean(den_pm_y);

save('/short/jz21/dg2271/overall_2001_2006.mat','M_x','M_y','alpha_x','alpha_y','mean_x','mean_y','gini_x','gini_y',...
     'quantile_est_x','quantile_est_y','GLD_est_x','GLD_est_y','LD_est_x','LD_est_y','den_est_x','den_est_y',...
     'prop_xGLDy','prop_yGLDx','prop_xLDy','prop_yLDx','prop_xFSDy','prop_yFSDx',...
     'mean_overall_xGLDy','mean_overall_yGLDx','mean_overall_xGLDy_20lowest','mean_overall_yGLDx_20lowest',...
     'mean_overall_xGLDy_10lowest','mean_overall_yGLDx_10lowest','mean_overall_xLDy','mean_overall_yLDx','mean_overall_xLDy_20lowest',...
     'mean_overall_yLDx_20lowest','mean_overall_xLDy_10lowest','mean_overall_yLDx_10lowest','mean_overall_xFSDy','mean_overall_yFSDx',...
     'mean_overall_xFSDy_20lowest','mean_overall_yFSDx_20lowest','mean_overall_xFSDy_10lowest','mean_overall_yFSDx_10lowest',...
     'max_overall_xGLDy','max_overall_yGLDx','max_overall_xGLDy_20lowest','max_overall_yGLDx_20lowest',...
     'max_overall_xGLDy_10lowest','max_overall_yGLDx_10lowest','max_overall_xLDy','max_overall_yLDx','max_overall_xLDy_20lowest',...
     'max_overall_yLDx_20lowest','max_overall_xLDy_10lowest','max_overall_yLDx_10lowest','max_overall_xFSDy','max_overall_yFSDx',...
     'max_overall_xFSDy_20lowest','max_overall_yFSDx_20lowest','max_overall_xFSDy_10lowest','max_overall_yFSDx_10lowest',...
     'min_overall_xGLDy','min_overall_yGLDx','min_overall_xGLDy_20lowest','min_overall_yGLDx_20lowest',...
     'min_overall_xGLDy_10lowest','min_overall_yGLDx_10lowest','min_overall_xLDy','min_overall_yLDx','min_overall_xLDy_20lowest',...
     'min_overall_yLDx_20lowest','min_overall_xLDy_10lowest','min_overall_yLDx_10lowest','min_overall_xFSDy','min_overall_yFSDx',...
     'min_overall_xFSDy_20lowest','min_overall_yFSDx_20lowest','min_overall_xFSDy_10lowest','min_overall_yFSDx_10lowest');

%save('overall.mat','overall_xGLDy','overall_yGLDx','overall_xGLDy_20lowest','overall_yGLDx_20lowest',...
%     'overall_xGLDy_10lowest','overall_yGLDx_10lowest','overall_xLDy','overall_yLDx','overall_xLDy_20lowest','overall_yLDx_20lowest',...
%     'overall_xLDy_10lowest','overall_yLDx_10lowest','overall_xFSDy','overall_yFSDx','overall_xFSDy_20lowest','overall_yFSDx_20lowest',...
%     'overall_xFSDy_10lowest','overall_yFSDx_10lowest')