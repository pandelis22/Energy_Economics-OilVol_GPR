% Copyright: Panagiotis Delis
% pandelis@unipi.gr

function [ERR_DMA,mComb,mPost_Omega,y_t_pr,mPredLik,BestModel,MSPE_DMA,MSPE_DMS,MSPE_MCF,MSPE_COMB,...
    Y_C_for,MAE_DMA,MAE_DMS,MAE_MCF,MAE_COMB,vMean_DMA,ssize_general,DoC,Y_USO] = DMA_logRV_dir(...
    y,x,p1,lambda,vKappa,dAlpha,est_date,est_sample,h,method,date1,logorlevel,meth_forg_fact_self_pert,maxh,FOR_COND,SSHOCKS_VAR,USO)

if strcmpi(logorlevel,'log')
    y=log(y);
else
end  
if ~isempty(p1)
    if size(p1,1)==1 && size(p1,2)>=1 
        p=[p1' p1'];
    elseif min(size(p1))==1 
        p=[ones(size(p1)) p1];
    else
    end
end
numP = size(p,1);
maxP = max(max(p));    
[Y,X] = lagmatrix(y,maxP,0);
[Y_USO,X_x_USO] = lagmatrix(USO,maxP,0);
T = length(Y);
newX_y = zeros(T,numP);
newX_x=[];
newX_x_IND=[];
Y_x=[];
for i=1:numP
    newX_y(:,i) = mean(X(:,p(i,1):p(i,2)),2);
end
for k=1:size(x,2)
    [Y_x1,X_x] = lagmatrix(x(:,k),maxP,0);
    new_X_x = zeros(T,numP);
    for i=1:numP
        if isempty(FOR_COND)
            new_X_x(:,i) = mean(X_x(:,p(i,1):p(i,2)),2);
        else
            new_X_x(:,i) = Y_x1;
        end
    end 
    Y_x=[Y_x Y_x1];
    newX_x=[newX_x new_X_x(:,1)];% lags as predictors
    newX_x_IND=[newX_x_IND new_X_x(:,1:end)];% lags as predictors
end
newX_IND = [newX_y newX_x_IND];
newX = [newX_y newX_x];

date=date1(maxP+1:end);
date_num=datenum(date);
date_ssize=datenum(est_date);
[ssize_general,~]=find(date_ssize==date_num);
ssize=ssize_general;
ssize1=ssize-est_sample;
    [cT, cK] = size(newX);
    if isempty(SSHOCKS_VAR)
        N = cK;
    else
    end
    comb = cell(N,1);
    for nn = 1:N
        comb{nn, 1} = nchoosek(1:N, nn);
    end
    K = 2.^N; 
    index = cell(K, 1);
    dim = zeros(1, N+1);
    for nn = 1:N
        dim(:, nn+1) = size(comb{nn,1},1);
        for jj = 1:size(comb{nn,1},1)
            index{jj + sum(dim(:,1:nn)), 1} = comb{nn, 1}(jj, :);
        end
    end
    lambda = num2cell(lambda);
    vKappaC = num2cell(vKappa);
    mComb = allcomb(lambda, vKappaC, index);
    iM = size(mComb, 1);
    omega_predict = ones(1, iM)./iM;
    mpmt = zeros(cT, iM);
    mPredLik = zeros(cT, iM);
    mMeans = zeros(cT, iM);
    merrors = zeros(cT, iM);
    err = zeros(cT, iM);
    DA = zeros(cT, iM); 
    y_t_pr = zeros(cT, iM);
    mVariances = zeros(cT, iM);
    vOutputBeta = zeros(cT, 1);
    vOutputKappa = zeros(cT, 1);
    vtopMmpt = zeros(cT, 1);
    vOutputBestModel = zeros(cT, 1);
    mPost_Omega = zeros(cT, iM);
    MSPE_all=[];
    MAE_all=[];
    DoC_all=[];  
    
    for z=1:size(mComb, 1)
        disp(z)
        if ~isempty(SSHOCKS_VAR)
        else
            index = mComb{z, size(mComb,2)};          
            SSHOCKS1=[];
        end      
        XX_all=[];
           Xnew = [ones(cT,1), newX(:,index)]; 
        idx=find(ismember(mComb{z, size(mComb,2)},4));
        if isempty(idx)||isempty(FOR_COND) 
        RV_for_cond=[{FOR_COND} {[]}];    
        else
        RV_for_cond=[{FOR_COND} {idx}];        
        end      
        ERR_all=[];
        dbeta=[1e-10];
        [~,RV_for_TVP,errors_TVP,~,MSPE,MAE,~,~,VPredLik,Variance,DoC,P_DoC,~]=TVP_HAR_RV_DMA(...
        y,[],0,1,p1,est_date,est_sample,h,mComb{z},dbeta(1),logorlevel,method,meth_forg_fact_self_pert,maxh,date1,Xnew,RV_for_cond,SSHOCKS1);            
    
    MSPE_all=[MSPE_all;MSPE];
    MAE_all=[MAE_all;MAE];
    DoC_all=[DoC_all;DoC];
    
    DA(:,z)=P_DoC;
    err(:,z)=errors_TVP;
    y_t_pr(:,z)=RV_for_TVP;
    mPredLik(:,z)=VPredLik;
    mVariances(:,z)=Variance;
    end  
    min_MSPE=mComb{find(MSPE_all==min(MSPE_all)),3};
    min_MAE=mComb{find(MAE_all==min(MAE_all)),3};
    max_DoC=mComb{find(DoC_all==max(DoC_all)),3};
    
    merrors=err(ssize+h:end,:);
    mPredLik=mPredLik(ssize+h:end,:);
    mMeans=y_t_pr(ssize+h:end,:);
    mVariances=mVariances(ssize+h:end,:);
    for i = 1:cT-ssize-2*(h-1) 
        %disp(['----> h=',num2str(h)])
        if i == 1
           omega_predict = ones(1, iM)./iM;
        else
           dsumprob = sum(vOmega_Update.^dAlpha);
           omega_predict = (vOmega_Update.^dAlpha + 1e-200)./(dsumprob + 1e-200);
        end
         w_t = (mPredLik(i, :).* omega_predict);
         vSum_W = sum(w_t, 2);
         vOmega_Update = (w_t + 1e-200)./(vSum_W + 1e-200);
         mPost_Omega(ssize+i+h-1, :) = vOmega_Update;
         vmodel = 1:1:iM;
         vmodel = [vmodel', vOmega_Update'];
         vselect = (vmodel == max(vmodel(:, 2)));
         vselect = vmodel(:, 1).*vselect(:, 2);
         dBestModel = vselect(vselect~=0, :);
         if size(dBestModel, 1) > 1
            dBestModel = dBestModel(end, 1);
         end
         vtopMmpt(ssize+i) = max(vOmega_Update); 
         vOutputBestModel(ssize+i+h-1) = dBestModel(end, 1);
         vOutputBeta(i) = mComb{dBestModel, 1};
         vOutputKappa(i) = mComb{dBestModel, 2};
    end
    vPredlik_DMA = zeros(rows(vOutputBestModel),1);
    vMean_DMA = zeros(rows(vOutputBestModel),1);
    vVariance_DMA = zeros(rows(vOutputBestModel),1);
    vPredlik_DMS = zeros(rows(vOutputBestModel),1);
    vMean_DMS = zeros(rows(vOutputBestModel),1);
    vVariance_DMS = zeros(rows(vOutputBestModel),1);
    for i=1:rows(vOutputBestModel)-ssize-2*h+1
        % ======================= DMS part ============================== %
        vPredlik_DMS(ssize+i+h) = mPredLik(i+1,vOutputBestModel(ssize+i+h-1));
        vMean_DMS(ssize+i+2*h-1) = mMeans(i+h,vOutputBestModel(ssize+i+h-1));
        vVariance_DMS(ssize+i+h) = mVariances(i+1,vOutputBestModel(ssize+i+h-1));
        % ======================= DMA part ============================== %
        vPredlik_DMA(ssize+i+h) = mPost_Omega(ssize+i+h-1, :) * mPredLik(i+1, :)';
        vMean_DMA(ssize+i+2*h-1) = mPost_Omega(ssize+i+h-1, :) *  mMeans(i+h, :)';
        vVariance_DMA(ssize+i+h) = mPost_Omega(ssize+i+h-1, :) *  mVariances(i+1,:)';  
        % ======================= MCF part ============================== %
        if Y(ssize+i+2*h-1,1)>Y(ssize+i+h-1,1) && vMean_DMA(ssize+i+2*h-1,1)>Y(ssize+i+h-1,1)
            P(ssize+i+2*h-1,1)=1;
        elseif Y(ssize+i+2*h-1,1)<Y(ssize+i+h-1,1) && vMean_DMA(ssize+i+2*h-1,1)<Y(ssize+i+h-1,1)
            P(ssize+i+2*h-1,1)=1;
        else
            P(ssize+i+2*h-1,1)=0;
        end
    end
    DoC=mean(P(ssize+maxh:end,1))*100;
    vMean_MCF = mean(y_t_pr,2);
    vMedian_MCF = median(y_t_pr,2);
    if strcmpi(logorlevel,'log')
        Y=exp(Y);
    else
    end
    Y_oos=Y;
    Y_oos(find(vMean_DMA==0))=0;
    ERR_DMA=vMean_DMA-Y_oos;
    ERR_DMS=vMean_DMS-Y_oos;
    ERR_MCF=vMean_MCF-Y_oos;
    ERR_median_MCF=vMedian_MCF-Y_oos;    
    MSPE_MEDIAN=mean(ERR_median_MCF(ssize+maxh:end,1).^2);    
    MSPE_DMA=mean(ERR_DMA(ssize+maxh:end,1).^2);
    MSPE_DMS=mean(ERR_DMS(ssize+maxh:end,1).^2);
    MSPE_MCF=mean(ERR_MCF(ssize+maxh:end,1).^2);
    MAE_DMA=mean(abs(ERR_DMA(ssize+maxh:end,1)));
    MAE_DMS=mean(abs(ERR_DMS(ssize+maxh:end,1)));
    MAE_MCF=mean(abs(ERR_MCF(ssize+maxh:end,1)));
    BestModel=mComb(vOutputBestModel(ssize+maxh:end-h+1,:),3);
  
MSPE_COMB=[];
MAE_COMB=[];
Y_C_for=[];
end

% FUNCTIONS
function A = allcomb(varargin)
NC = nargin ;
ii = NC:-1:1 ;
isCellInput = cellfun(@iscell,varargin);
if any(isCellInput)
    ix = cellfun(@(c) 1:numel(c), varargin,'un',0) ;
    [ix{ii}] = ndgrid(ix{ii});
    A = cell(numel(ix{1}),NC);
    for k=1:NC,
        A(:,k) = reshape(varargin{k}(ix{k}),[],1) ;
    end
else
    [A{ii}] = ndgrid(varargin{ii});
    A = reshape(cat(NC+1,A{:}),[],NC);
end
end

function r = rows(x)
  [r,c] = size(x);
end

function [y,x]=lagmatrix(x,lags,const)

[T,M]=size(x);
if T<M
    x=x';
end
if lags>0
    lags=lags+1;
    newX=[x;zeros(lags,1)];
    lagmatrix=repmat(newX,lags,1);
    lagmatrix=reshape(lagmatrix(1:size(lagmatrix,1)-lags),T+lags-1,lags);
    lagmatrix=lagmatrix(lags:T,:);
    y=lagmatrix(:,1);
    x=lagmatrix(:,2:lags);
    if const==1
        x=[ones(size(x,1),1) x];
    end
else
    if const==1
        y=x;
        x=ones(T,1);
    else
        y=x;
        x=[];
    end
end
end

function [beta,s] = ols(y,x,const)
N = size(y,1);
M=size(x,2);
if const==1
    x=[ones(N,1) x];
    M=M+1;
end
beta=x\y;
s=(y-x*beta)'*(y-x*beta)/(N-M);
end
