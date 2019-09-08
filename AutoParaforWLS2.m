function [map,map2] = AutoParaforWLS2(T3mat)

% The code was written by Yexian Ren (2018-3-30 V1) 
% Please refer to this paper for a more detailed description of the algorithm.
% Yexian Ren, Jie Yang, Lingli Zhao, Pingxiang Li, Zhiqu Liu, Lei Shi. 
% A Global Weighted Least Squares Optimization Framework for Speckle Filtering 
% of PolSAR Imagery. IEEE Transactions on Geoscience And Remote Sensing. 2018 (99): 1-13. doi: 10.1109/TGRS.2018.2865507

% If you have any comment, or question, please contact Yexian Ren at
% renyexian@foxmail.com  or QQ: 2538715345

%% Function£ºAutomatic estimation of filter parameters
disp('Automatic estimation of filter parameters...')

r = 3;
Pmat = T3mat(:,:,1)+T3mat(:,:,6)+T3mat(:,:,9);
output2 = medfilt2(Pmat,[(2*r+1),(2*r+1)]);

[m,n,ch] = size(T3mat);
map = 3*ones(m,n);
map2 = zeros(m,n);
N = m*n;

HomoIndex = TraceENL(T3mat,r);
look = mean2(HomoIndex);

data0 = HomoIndex(:);

[mu,variance,p,label] = PDFfit3(data0,4);

Llabel = reshape(label,m,n);

map(find(Llabel==1)) = 0.25*sqrt(1/look);
map(find(Llabel==2)) = 0.5*sqrt(1/look);
map(find(Llabel==3)) = 1.5*sqrt(1/look);
map(find(Llabel==4)) = 3*sqrt(1/look);

map2(find(Llabel==1)) = 1.2;
map2(find(Llabel==2)) = 1.2;
map2(find(Llabel==3)) = 2;
map2(find(Llabel==4)) = 2;

th = 15;
a = (Pmat./output2)>th;
map(find(a==1)) = 0;

end

function ENLmap = TraceENL(T3mat,r)
[m,n,channel] = size(T3mat);
N = boxfilterV2(ones(m, n, 9), r);
M = boxfilterV2(T3mat, r)./N;
trM2 = (M(:,:,1) + M(:,:,6) + M(:,:,9)).^2;
trTT = TrCC(T3mat);
M_trTT = boxfilterV2(trTT, r)./N(:,:,1);
trMM = TrCC(M);
ENLmap = trM2./(M_trTT - trMM);
ENLmap(isnan(ENLmap)==1)=0;
end

function imDst = boxfilterV2(imSrc, r)
[hei,wid,channel] = size(imSrc);
imDst = zeros(size(imSrc));
imCum = cumsum(imSrc, 1);
imDst(1:r+1, :, :) = imCum(1+r:2*r+1, :, :);
imDst(r+2:hei-r, :, :) = imCum(2*r+2:hei, :, :) - imCum(1:hei-2*r-1, :, :);
imDst(hei-r+1:hei, :, :) = repmat(imCum(hei, :, :), [r, 1]) - imCum(hei-2*r:hei-r-1, :, :);
imCum = cumsum(imDst, 2);
imDst(:, 1:r+1, :) = imCum(:, 1+r:2*r+1, :);
imDst(:, r+2:wid-r, :) = imCum(:, 2*r+2:wid, :) - imCum(:, 1:wid-2*r-1, :);
imDst(:, wid-r+1:wid, :) = repmat(imCum(:, wid, :), [1, r]) - imCum(:, wid-2*r:wid-r-1, :);
end

function TmDet = computeMDet(Tm)
T11 = Tm(:,:,1);
T12 = complex(Tm(:,:,2),Tm(:,:,3));
T21 = complex(Tm(:,:,2),-Tm(:,:,3));
T13 = complex(Tm(:,:,4),Tm(:,:,5));
T31 = complex(Tm(:,:,4),-Tm(:,:,5));
T22 = Tm(:,:,6);
T23 = complex(Tm(:,:,7),Tm(:,:,8));
T32 = complex(Tm(:,:,7),-Tm(:,:,8));
T33 = Tm(:,:,9);
TmDet=T11.*T22.*T33-T11.*T23.*T32-T12.*T21.*T33+T12.*T23.*T31+T13.*T21.*T32-T13.*T22.*T31;
TmDet = abs(TmDet);
end

function trcc = TrCC(C)
trcc = C(:,:,1).^2 + C(:,:,6).^2 + C(:,:,9).^2 + 2*(C(:,:,2).^2 + C(:,:,3).^2 + C(:,:,4).^2 + C(:,:,5).^2 + C(:,:,7).^2 + C(:,:,8).^2);
end

function [mu,variance,p,label] = PDFfit3(data,K)
data0 = sort(data(:));
N = numel(data0);
muInit = ones(K,1);

for i = 1:K
    ceil((i-1)*N/K)
    t1 = data0(floor((i-1)*N/K)+1);
    t2 = data0(floor(i*N/K));
    muInit(i) = mean(data0(find(data0>=t1&data0<=t2)));
end

% [mu,variance,p,iterNum,label] = GaussEM2(data,K,50,0.01,muInit);
 [mu,variance,p,iterNum,label] = GaussEM2(data,K,10,0.01,muInit);

end

function [mu,variance,p,iterNum,label] = GaussEM2(X,numOfComponent,maxIter,tol,muInit,varianceInit,pInit)
X = X(:);
N = numel(X); 
K = numOfComponent; 

if nargin==4
    temp = randperm(N);
    mu = X(temp(1:K));mu = mu(:);
    clear temp;
    variance = var(X) * ones(K,1);
    p = (1/K) * ones(K,1);
end
if nargin==5
    mu = muInit(:);
    variance = var(X) * ones(K,1);
    p = (1/K) * ones(K,1);
end
if nargin==6
    mu = muInit(:);
    variance = varianceInit(:);
    p = (1/K) * ones(K,1);
end
if nargin==7
    mu = muInit(:);
    variance = varianceInit(:);
    p = pInit(:);
end

for i=1:maxIter
    G = repmat(p',N,1) ./ sqrt(repmat(variance',N,1)) .* exp(-((repmat(X,1,K)-repmat(mu',N,1)).^2) ./ (2*repmat(variance',N,1)));
    g = sum(G,2);g=g(:);
    % Dividing operation is prone to abnormal calculation
    g(find(g<1e-12)) = 1e-6;
    
    G = G./ repmat(g,1,K);
    M = sum(G,1);
    M = M(:);
    M = M + 1e-6;
    M = M/sum(M) * N;
    muNEW = (G'* X) ./ M;
    varianceNEW = sum(G.*((repmat(X,1,K)-repmat(muNEW',N,1)).^2),1)'./M + 1e-6;
    pNEW = M/N;
    if (norm(mu-muNEW,1)/norm(mu,1)<tol && norm(varianceNEW-variance,1)/norm(variance,1)<tol && norm(pNEW-p,1)/norm(p,1)<tol)
        break;
    end
    mu = muNEW;
    variance = varianceNEW;
    p = pNEW;
end
iterNum = i;
mu = muNEW
variance = varianceNEW;
p = pNEW;

[Y,label] = max(G');
label = label(:);
end



