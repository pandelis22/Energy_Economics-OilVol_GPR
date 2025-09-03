% Copyright: Panagiotis Delis
% pandelis@unipi.gr

function [y_t_pred,RV_for,errors,newX,MSPE,MAE,theta_update,R_t,vPredLik,Variance,DoC,P,ssize_general]=TVP_HAR_RV_DMA(y,x,xlags_only,constant,p,ssize,ssize1,h,lambda,dbeta,logorlevel,method,meth_forg_fact_self_pert,maxh,date,Xnew,RV_for_cond,sshocks)

if strcmpi(logorlevel,'level')
else
end
if strcmpi(meth_forg_fact_self_pert,'SP')
    lambda=1;
else
end
if ~isempty(p)
    if size(p,1)==1 && size(p,2)>=1 
        p=[p' p'];
    elseif min(size(p))==1 
        p=[ones(size(p)) p];
    else
    end
end
p = unique(p,'rows');

numP = size(p,1);
maxP = max(max(p));
[Y,X] = lagmatrix(y,maxP,0);
T = length(Y);
if isempty(Xnew)
newX_y = zeros(T,numP);
newX_x=[];
Y_x=[];
for i=1:numP
    newX_y(:,i) = mean(X(:,p(i,1):p(i,2)),2);
end
if isempty(x)
    if constant
        newX = [ones(size(newX_y,1),1) newX_y];
    end  
    Y_x=zeros(T,size(x,2));
    new_X_x=zeros(T,numP);
    newX_x=zeros(T,numP);   
else
    for k=1:size(x,2)
    [Y_x1,X_x] = lagmatrix(x(:,k),maxP,0);
    new_X_x = zeros(T,numP);
    for i=1:numP
        if isempty(RV_for_cond{1,2})
            new_X_x(:,i) = mean(X_x(:,p(i,1):p(i,2)),2);
        else
            new_X_x(:,i) = Y_x1;
        end
    end 
    Y_x=[Y_x Y_x1];
    newX_x=[newX_x new_X_x(:,1)];
    end
    if constant        
        if xlags_only==1
        newX = [ones(size(newX_y,1),1) newX_x];
        else
        newX = [ones(size(newX_y,1),1) newX_y newX_x];    
        end
    else
        newX = [newX_y newX_x];
    end
end
else
    if isempty(sshocks)
        newX=Xnew; 
    else
        newX=[Xnew zeros(size(Xnew,1),size(sshocks{1},2))]; 
    end
end
numP = size(p,1);
maxP = max(max(p));

date=date(maxP+1:end);
date_num=datenum(date);
date_ssize=datenum(ssize);
[ssize_general,~]=find(date_ssize==date_num);
ssize=ssize_general;
ssize1=ssize-ssize1;

T = length(Y);
if isempty(sshocks)
    theta_pred=zeros(size(newX,2),T);
    R_t=zeros(size(newX,2),size(newX,2),T);
    theta_update=zeros(size(newX,2),T);
    S_t=zeros(size(newX,2),size(newX,2),T);    
    m=size(newX,2);
else
    theta_pred=zeros(size(Xnew,2)+size(sshocks{1},2),T);
    R_t=zeros(size(Xnew,2)+size(sshocks{1},2),size(Xnew,2)+size(sshocks{1},2),T);
    theta_update=zeros(size(Xnew,2)+size(sshocks{1},2),T);
    S_t=zeros(size(Xnew,2)+size(sshocks{1},2),size(Xnew,2)+size(sshocks{1},2),T);
    m=size(Xnew,2)+size(sshocks{1},2);
end

y_t_pred=zeros(T,1);
e_t=zeros(1,T);
e=zeros(1,T);
errors=zeros(T,1);
V_t=zeros(1,T);
vPredLik=zeros(T,1);
dx1=[];
for i=1:T-h+1
x_t(i+h-1,:)=newX(i,:);
end
if isempty(sshocks)
else
   newX(ssize1+1:ssize,:)=[newX(ssize1+1:ssize,1:size(Xnew,2)) sshocks{1}(1:end-1,:)];
end
theta_OLS = (Y(ssize1+1:ssize,:)'/newX(ssize1+1:ssize,:)')';
theta_0_prmean = theta_OLS;
temp = (Y(ssize1+1:ssize,:)'/newX(ssize1+1:ssize,:)')';
s_02 = temp(1,1).^2 + var(Y(ssize1+1:ssize,:));
theta_0_prvar = 0.1*diag([s_02./100;(var(Y(ssize1+1:ssize,:))./var(newX(ssize1+1:ssize,2:end)))'./100]);
s2=((Y(ssize1+1:ssize,:)-newX(ssize1+1:ssize,:)*theta_OLS)'*(Y(ssize1+1:ssize,:)-newX(ssize1+1:ssize,:)*theta_OLS)/(ssize-ssize1-m));
theta_0_prvar = s2*inv(newX(ssize1+1:ssize,:)'*newX(ssize1+1:ssize,:));

RV_for=zeros(T,1);
V_0 = s2;
V_t(:,ssize)=V_0;
theta_update(:,ssize)=theta_0_prmean;
S_t(:,:,ssize)=theta_0_prvar;
Variance=zeros(T,1);
dLogF = 0;
dSumSquares = 0;
if h==1||strcmpi(method,'direct')
    newX=x_t;
    for i=1:T-ssize-h+1
        theta_pred(:,ssize+i) = theta_update(:,ssize+i-1);  
        R_t(:,:,ssize+i) = (1./lambda)*squeeze(S_t(:,:,ssize+i-1));    
        if strcmpi(logorlevel,'log')
        y_t_pred(ssize+i+h-1,:) = exp(newX(ssize+i+h-1,:)*theta_pred(:,ssize+i));    
        errors(ssize+i+h-1,:)=y_t_pred(ssize+i+h-1,:)-exp(Y(ssize+i+h-1,:));
        else
            if isempty(RV_for_cond{1,2})
        y_t_pred(ssize+i+h-1,:) = newX(ssize+i+h-1,:)*theta_pred(:,ssize+i); 
            else
                if size(RV_for_cond{1,1},2)==2
                    COND = [RV_for_cond{1,1}(i,1) RV_for_cond{1,1}(i,2)];
                elseif size(RV_for_cond{1,1},2)==1
                    COND = RV_for_cond{1,1}(i,1);
                end                
        y_t_pred(ssize+i+h-1,:) = [1 newX(ssize+i+h-1,2:RV_for_cond{1,2}) COND]*theta_pred(:,ssize+i);        
            end
        errors(ssize+i+h-1,:)=y_t_pred(ssize+i+h-1,1)-Y(ssize+i+h-1,:);         
        end
        RV_for(ssize+i+h-1,1)=y_t_pred(ssize+i+h-1,:);
        Variance(ssize+i+h-1,1)=V_t(:,ssize+i-1) + newX(ssize+i, :) * R_t(:,:,ssize+i) * newX(ssize+i, :)';   
        e_t(:,ssize+i) = Y(ssize+i,:) - newX(ssize+i,:)*theta_pred(:,ssize+i); 
            A_t = ((i-1)/i)*V_t(:,ssize+i-1)+(1/i)*((e_t(:,ssize+i)).^2 - newX(ssize+i,:)*squeeze(R_t(:,:,ssize+i))*newX(ssize+i,:)');
            if A_t>0
                V_t(:,ssize+i) = A_t;
            else
                V_t(:,ssize+i) = V_t(:,ssize+i-1);
            end
        %else
        %    A_t = (e_t(:,ssize+i-1)).^2;  
        %   V_t(:,ssize+i) = 0.94*V_t(:,ssize+i-1) + (1-0.94)*A_t;
        %end
        if strcmpi(meth_forg_fact_self_pert,'SP')
           dx = ((e_t(:,ssize+i).^2)/V_t(:,ssize+i) - 1);
           dx1=[dx1 dx];
           mterm = dbeta * max(0, floor(dx)) * eye(size(theta_pred, 1));             
        else
           mterm=0; 
        end
           dF = V_t(:,ssize+i) + newX(ssize+i, :) * R_t(:,:,ssize+i) * newX(ssize+i, :)';
           dLogF = dLogF + log(dF);
           dSumSquares = dSumSquares + (e_t(:,ssize+i)^2)/dF;
           vPredLik(ssize+i+h-1,1) = mvnpdf(Y(ssize+i,:), newX(ssize+i,:) * theta_pred(:,ssize+i), V_t(:,ssize+i) + newX(ssize+i, :) * R_t(:,:,ssize+i) * newX(ssize+i, :)');    
          theta_update(:,ssize+i) = theta_pred(:,ssize+i) + squeeze(R_t(:,:,ssize+i))*newX(ssize+i,:)'*inv(V_t(:,ssize+i) + newX(ssize+i,:)*squeeze(R_t(:,:,ssize+i))*newX(ssize+i,:)')*e_t(:,ssize+i);     
          S_t(:,:,ssize+i) = squeeze(R_t(:,:,ssize+i)) - squeeze(R_t(:,:,ssize+i))*newX(ssize+i,:)'*inv(V_t(:,ssize+i) + newX(ssize+i,:)*squeeze(R_t(:,:,ssize+i))*newX(ssize+i,:)')*newX(ssize+i,:)*squeeze(R_t(:,:,ssize+i))+mterm;     
    end
dLogLik = 0.5*((T)*log(2 * pi)+ dLogF + dSumSquares);
else
end
if strcmpi(logorlevel,'log')
    Y=exp(Y);
else
end
for i=ssize+h:T
    if Y(i,1)>Y(i-h,1) && RV_for(i,1)>Y(i-h,1)
        P(i,1)=1;
    elseif Y(i,1)<Y(i-h,1) && RV_for(i,1)<Y(i-h,1)
        P(i,1)=1;
    else
        P(i,1)=0;
    end
end
DoC=mean(P(ssize+maxh:end,1))*100;
MSPE=mean(errors(ssize+maxh:end,:).^2);
MAE=mean(abs(errors(ssize+maxh:end,:)));
end

% FUNCTIONS
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
