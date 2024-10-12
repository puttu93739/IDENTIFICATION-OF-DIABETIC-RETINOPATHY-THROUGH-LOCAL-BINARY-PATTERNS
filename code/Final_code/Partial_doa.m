function  f=Partial_doa(N,P,L,Theta,s,v,Fc) 
  
 %NOTE   : this code is subject of modification/correction. 
 % 
 % Reference : 
 %==================================================================== 
 % "Fast algorithm for DOA Estimation with Partial Covariance Matrix   | 
 % and Without EigenDecomposition"                                     | 
 % Jianfeng Chen, Yuntao Wu, Hui Cao, Hai Wang                         | 
 % Journal of Signal and information Processing , 2011,2,266-269       | 
 %==================================================================== 
 % 
 % Usage : 
 % N      : number of sensors . 
 % P      : number of narrowband signals with condition N>2P. 
 % L      : number of snapshots ( samples),  
 % v      : AWGN variance.  
 % theta  : true angles of arrival, theta is vector of length P. 
 % lambda : wavelength with condition d<=lambda/2. 
 % d      : space between ULA elements . 
 % Fc     : Frequency of the incoming signals . 
 %   
 % 
 % Demo : N=10;P=2;L=100; Theta=[40 55]; v=2; s=[10 10];Fc=1e+9 .  
 %        f=Partial_doa(N,P,L,Theta,s,v); 
 % 
 % 
 %        (c) KHMOU Youssef, 2012  
  
  
 A=zeros(N,P); %Array manifold matrix 
 %Fc=1000000000; % signal frequency Fc=1Ghz 
 c=299792458 ;% light's speed 
 lambda=c/Fc; % wavelength in meters 
 d=lambda/2; % space between ULA elements 
 
 for iter=1:length(Theta) 
    phi=2*pi*d*sin(pi/180*Theta(iter))/lambda; 
    A(:,iter)=exp(j*phi*(0:N-1)); 
 end 
Y=A*sqrt(s'); 
I=eye(N-2*P); 
Q=zeros(N-2*P,N); 
for t=1:L 
Noise=sqrt(v/2)*(randn(N,1)+j*randn(N,1)); 
X=Y+Noise; 
R12=X(1:P)*conj(X(P+1:2*P)'); 
R31=X((2*P)+1:N)*conj(X(1:P)'); 
R32=X((2*P)+1:N)*conj(X(P+1:(2*P))'); 
R21=X(P+1:(2*P))*conj(X(1:P)'); 
T=[R32*(pinv(R12)) R31*(pinv(R21)) -2*I]; 
Q=Q+T; 
end 
clear T; 
Q=Q./(L); % Estimated matrix 
% End of data acquition 
 
% Processing : Partitioning the data array under the condtion that  N>=2P 
%A1=A(1:P,:); 
%A2=A(P+1:2*P,:); 
%A3=A(2*P+1:N,:); 
%R12=A1*Rss*conj(A2'); 
%R31=A3*Rss*conj(A1'); 
%R32=A3*Rss*conj(A2'); 
%R21=A2*Rss*conj(A1'); 
 
% Testing : 
theta=-90:0.1:90; 
phi=(2*pi*d*sin(pi/180*theta))./lambda; 
a=exp(j*(0:N-1)'*phi); 
F=1./((abs(Q*a)).^2); 
f=mean(F); 
f=10*log10(f); 
%**************************PLOT RESULTS*********************************** 
figure,plot(theta,f,'r'), %axis([-60 60 0 1]) 
xlabel(' DOA in degrees') 
ylabel('Spectrum in dB') 
grid on, 
title('Doa estimation, for ULA ') 
hold on 
AOA=zeros(1,length(theta)); 
for i=1:length(Theta) 
AOA(Theta(i))=max(f); 
end 
plot(AOA,'-b'), legend('Estimated doas','true doas'),axis([theta(1) theta(end) 0 max(f)+max(f)/10]),hold off 
 
 