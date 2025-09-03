% Copyright: Panagiotis Delis
% pandelis@unipi.gr

clc
clear
load('data_final.mat')
path_folder_images='C:\Users\Eurhope\Desktop\Research projects\ECB_backup_29062021\FULL_MOST_UPDATED_ONE\1st_round\Personal_folders\Research\paper_Oil_volatility_Uncertainty\results_25042024\';
models={'HAR-OVX';'HAR-OVX-GPR';'DMA ALL'};
est_date='2014-01-15';
h=[1 5 10 15 22];
maxh=max(h);
i=h(1);
p=[1;5;22];
est_sample=1000;
dma_method='SP';
lambda=0.99;
alpha=0.99;
method='direct';
logorlevel='log';
SUM11=[];
p_VaR=0.05;


MCS_final=cell(max(h),1);
mPost_Omega_ALL_h=cell(max(h),1);
idx_ALL1_h=cell(max(h),1);
idx_ALL2_h=cell(max(h),1);
idx_ALL3_h=cell(max(h),1);
ERRORS=cell(size(RV_ann,2),1);
FORECASTS=cell(size(RV_ann,2),1);
MSPE_WTI=cell(size(RV_ann,2),1);
MAE_WTI=cell(size(RV_ann,2),1);
DoC_WTI=cell(size(RV_ann,2),1);
prob_OIL1_DMA=zeros(1408,size(RV_ann,2),1);
prob_OIL2_DMA=zeros(1408,size(RV_ann,2),1);
prob_OIL3_DMA=zeros(1408,size(RV_ann,2),1);
prob_MACRO1_DMA=zeros(1408,size(RV_ann,2),1);
prob_MACRO2_DMA=zeros(1408,size(RV_ann,2),1);
prob_MACRO3_DMA=zeros(1408,size(RV_ann,2),1);
prob_UNC1_DMA=zeros(1408,size(RV_ann,2),1);
prob_UNC2_DMA=zeros(1408,size(RV_ann,2),1);
prob_UNC3_DMA=zeros(1408,size(RV_ann,2),1);
prob_UNC4_DMA=zeros(1408,size(RV_ann,2),1);
prob_ALL1_DMA=zeros(1408,size(RV_ann,2),1);
prob_ALL2_DMA=zeros(1408,size(RV_ann,2),1);
prob_ALL3_DMA=zeros(1408,size(RV_ann,2),1);
CR_all=[];
RV_CR_OIL=[];
RV_CR_MACRO=[];
RV_CR_UNC=[];

X_GPR=[GPRD_OLD GPRD_THREAT_OLD GPRD_ACT_OLD];
X_OTH=[EPU ADS US_FSI];   
X_IV=[log(VIX) log(VXD) log(VXN)]; 

ERR=zeros(size(RV_ann,1)-max(p),size(models,1),max(h));
FOR=zeros(size(RV_ann,1)-max(p),size(models,1),max(h));

USO_eval=USO(max(p)+1:end,1);

%HAR-OVX
[~,~,err_AR3_OVX,RV_for_OVX,Y,MSPE_AR3_OVX,MAE_AR3_OVX,~,DoC_AR3_OVX,IP_for_all_OVX,ssize2,~]=HAR_RV(...
RV_ann,[log(OVX)],0,1,p,logorlevel,est_date,est_sample,i,method,maxh,DATE,USO);

%HAR-OVX-GPR
[~,~,err_AR3_OIL,RV_forX,Y,MSPE_AR3_OIL,MAE_AR3_OIL,~,DoC_AR3_OIL,IP_for_all1,ssize2,~]=HAR_RV(...
RV_ann,[log(OVX) X_GPR(:,1)],0,1,p,logorlevel,est_date,est_sample,i,method,maxh,DATE,USO);    

%DMA-ALL
[err_DMA_ALL,mComb_ALL,mPost_Omega_ALL,y_t_pr_ALL,mPredLik_ALL,BestModel_ALL,...
    MSPE_DMA_ALL,MSPE_DMS_ALL,MSPE_MCF_ALL,MSPE_COMB_ALL,Y_C_for_ALL,...
    MAE_DMA_ALL,MAE_DMS_ALL,MAE_MCF_ALL,MAE_COMB_ALL,RV_for_DMA_ALL,ssize4,DoC_DMA_ALL,Y_USO] = DMA_logRV_dir(...
    RV_ann,[log(OVX) X_GPR X_OTH X_IV],p,lambda,0.95,alpha,est_date,est_sample,i,method,DATE,logorlevel,dma_method,maxh,[],'',USO);
 
idx_ALL1 = cellfun(@(x)any(ismember(x,[5:7])),mComb_ALL(:,3))';
idx_ALL2 = cellfun(@(x)any(ismember(x,[8:10])),mComb_ALL(:,3))';
idx_ALL3 = cellfun(@(x)any(ismember(x,[11:13])),mComb_ALL(:,3))';

prob_ALL1_DMA(:,1,1)=sum(mPost_Omega_ALL(ssize4+2*max(h):end,idx_ALL1),2);
prob_ALL2_DMA(:,1,1)=sum(mPost_Omega_ALL(ssize4+2*max(h):end,idx_ALL2),2);
prob_ALL3_DMA(:,1,1)=sum(mPost_Omega_ALL(ssize4+2*max(h):end,idx_ALL3),2);

FOR(:,:,i)=[RV_for_OVX RV_forX RV_for_DMA_ALL];
ERR(:,:,i)=[err_AR3_OVX err_AR3_OIL err_DMA_ALL];

MSPE=mean(squeeze(ERR(ssize4+2*max(h):end,:,i)).^2,1)';