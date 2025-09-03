% Copyright: Panagiotis Delis
% pandelis@unipi.gr

function [newX,param,errors,RV_for,Y,MSPE_HAR,MAE_HAR,mdl,DoC,RV_for_all,ssize_general,Y_USO]=HAR_RV(y,x,xlags_only,constant,p,logorlevel,ssize,ssize1,h,method,maxh,date,USO)

RV_for_all=[];
if strcmpi(logorlevel,'level')
else
    y=log(y);
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
[Y_USO,X_x_USO] = lagmatrix(USO,maxP,0);
T = length(Y);
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
        new_X_x(:,i) = mean(X_x(:,p(i,1):p(i,2)),2);
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

date=date(maxP+1:end);
date_num=datenum(date);
date_ssize=datenum(ssize);
[ssize_general,~]=find(date_ssize==date_num);
ssize=ssize_general;
ssize1=ssize-ssize1;

T=size(Y,1);
RV_for=zeros(T,1);
param=zeros(size(newX,2),T);
errors=zeros(T,1);
[B,S] = ols(Y(ssize1+1:ssize,1),newX(ssize1+1:ssize,:),0);
param(:,ssize)=B;

if h==1||strcmpi(method,'direct')
    for i=h:T-ssize
[B,S] = ols(Y(ssize1+i-h+1:ssize+i-h,1),newX(ssize1+i-2*h+2:ssize+i-2*h+1,:),0);
param(:,ssize+i-h)=B;
if strcmpi(logorlevel,'log')
RV_for(ssize+i,1)=exp(newX(ssize+i-h+1,:)*param(:,ssize+i-h));
errors(ssize+i,1) = RV_for(ssize+i,1)-exp(Y(ssize+i,1)); 
else
RV_for(ssize+i,1)=newX(ssize+i-h+1,:)*param(:,ssize+i-h);
RV_for_all(ssize+i,1)=RV_for(ssize+i,1);
errors(ssize+i,1) = RV_for(ssize+i,1)-Y(ssize+i,1);     
end
    end
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
errors1=errors(ssize+maxh:end,1);
MSPE_HAR=mean(errors1.^2,1);
MAE_HAR=mean(abs(errors1));
mdl = fitlm(newX(:,2:end),log(Y(1:end,1)),'linear');
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
