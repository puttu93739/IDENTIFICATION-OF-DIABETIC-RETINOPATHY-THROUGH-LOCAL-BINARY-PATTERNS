%  DOA Estimation  
% This demo implemants  a fast Direction of arrival estimation of Uniform 
% Linear array as described in the following paper : 
% 
% " Fast Algorithm for DOA Estimation with Partial Covariance Matrix and 
%   without Eigendecomposition"   
%   Jianfeng Chen, Yuntao Wu,Hui Cao,Hai Wang  
% 
%  (c) 2012,KHMOU Youssef  
 
clear, 
N=10; % Uniform linear array ,  Number of sensors 
P=2; % Number of sources 
Theta=[10 45]; % Directions of arrival 
 
 
L=100; %  time samples  
Fc=1000000000; % signak frequency Fc=1Ghz 
c=299792458 ;% light speed 
lambda=c/Fc; % wavelength in meters 
d=lambda/2; % space between ULA elements 
v=2.9; % Noise variance 
s=[10 10] ; 
f=Partial_doa(N,P,L,Theta,s,v,Fc); 