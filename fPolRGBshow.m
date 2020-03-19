function RGBdata=fPolRGBshow(varargin)

% This code was adapted from PolSARpro
% The code was written by Yexian Ren £¨renyexian@foxmail.com£© 
% and used in the paper for the PolSAR Image Visualization
% Ren et al. "SIRV-Based High-Resolution PolSAR Image Speckle Suppression via Dual-Domain Filtering" ,IEEE Trans. Geosci. Remote Sens.
% If you have any comment, or question, please contact Yexian Ren at
% renyexian@foxmail.com  or QQ: 2538715345

%% Function: PolSAR Image Display

if nargin == 1
    T3 = varargin{1};
    R=double(T3(:,:,6));
    G=double(T3(:,:,9));
    B=double(T3(:,:,1));
    k = 3;
elseif nargin == 2    
    T3 = varargin{1};
    R=double(T3(:,:,6));
    G=double(T3(:,:,9));
    B=double(T3(:,:,1));
    k = varargin{2};
elseif nargin == 3    
    R=double(varargin{1});
    G=double(varargin{2});
    B=double(varargin{3});
    k = 3;
else
    R=double(varargin{1});
    G=double(varargin{2});
    B=double(varargin{3});
    k = varargin{4};
end

% deal with NaN
R(find(isnan(R))) = eps;
G(find(isnan(G))) = eps;
B(find(isnan(B))) = eps;

% Dynamic range adjustment
R=abs(R);B=abs(B);G=abs(G);
R(R<eps)=eps;B(B<eps)=eps;G(G<eps)=eps;
R=10*log10(R);
G=10*log10(G);
B=10*log10(B);

[rMin,rMax]=GetMinMax(R,k);
[gMin,gMax]=GetMinMax(G,k);
[bMin,bMax]=GetMinMax(B,k);

R(R<rMin)=rMin;R(R>rMax)=rMax;
G(G<gMin)=gMin;G(G>gMax)=gMax;
B(B<bMin)=bMin;B(B>bMax)=bMax;

R=(R-rMin)/(rMax-rMin);
G=(G-gMin)/(gMax-gMin);
B=(B-bMin)/(bMax-bMin);

R(R>1.0)=1;R(R<0)=0;
G(G>1.0)=1;G(G<0)=0;
B(B>1.0)=1;B(B<0)=0;

R=uint8(floor(R*255));
G=uint8(floor(G*255));
B=uint8(floor(B*255));

RGBdata(:,:,1)=R;
RGBdata(:,:,2)=G;
RGBdata(:,:,3)=B;

figure;imshow(RGBdata);

end

function [xMin,xMax]=GetMinMax(S,k)

if nargin<2
    k=3;
end

s=S(S~=10*log10(eps));
% xMin=min(s(:));
% xMax=max(s(:));
med0=median(s(:));
med1=med0;
med2=med0;
% xMin=min(S(:));
% xMax=max(S(:));
for m=1:k
    temp1=s(s<med1);
    if isempty(temp1)
        break
    end
    med1=median(temp1);
end
for m=1:k
    temp2=s(s>med2);
    if isempty(temp2)
        break;
    end
    med2=median(temp2);
end
xMin=med1;
xMax=med2;
end


