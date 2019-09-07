
%%

% The code was written by Yexian Ren (2018-3-30 V1) 
% Please refer to this paper for a more detailed description of the algorithm.
% Yexian Ren, Jie Yang, Lingli Zhao, Pingxiang Li, Zhiqu Liu, Lei Shi. 
% A Global Weighted Least Squares Optimization Framework for Speckle Filtering 
% of PolSAR Imagery. IEEE Transactions on Geoscience And Remote Sensing. 2018 (99): 1-13. doi: 10.1109/TGRS.2018.2865507

% If you have any comment, or question, please contact Yexian Ren at
% renyexian@foxmail.com  or QQ: 2538715345

%% The GWLS experiment demo

% Description: For very large images, it is recommended to do a multi-look processing first.

clear;clc;

load T3mat.mat
% T3mat = cat(3,T11,T12_real,T12_imag,T13_real,T13_imag,T22,T23_real,T23_imag,T33)
% 4-look PolSAR data

% show the original PolSAR image
z0 = fPolRGBshow(T3mat,3);

tic
% adaptive parameter based on ENL
[lambda,alpha] = AutoParaforWLS2(T3mat);
% global weighted least squares fitting
data1 = spanWLSfilterV2(T3mat,lambda,alpha);
toc

% show the filtered PolSAR image
z0 = fPolRGBshow(data1,3);

