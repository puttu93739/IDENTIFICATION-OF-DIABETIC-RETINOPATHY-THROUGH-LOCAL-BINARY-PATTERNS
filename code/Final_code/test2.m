%%%%%%%%%%% 
 
% number of lifetime 
Total_lifetime=30; 
Lifetime=[0:1:Total_lifetime]; 
 
% number of sensor 
Total_N=length(Lifetime); 
N=[1:1:Total_N]; 
 
% constant of exp distribution 
C=1; 
% length of Lifetime 
LL=length(Lifetime); 
% length of sensor 
LN=length(N); 
 
% reliability of system 
R=zeros(LN,LL); 
P=zeros(LN,LL); 
 
for j=1:LL 
for i=1:LN 
    for k=1:(N(i)) 
        R(i,j)=R(i,j)+(C*Lifetime(j))^(N(i)-k)/factorial(N(i)-k); 
    end 
    R(i,j)=R(i,j)*exp(-C*Lifetime(j)); 
end 
end 
%P=1-R; 
 
 
Lifetime1=[5,10,15,20]; 
Rth=[0.01,0.8,0.99]; 
% length of Lifetime 
LL=length(Lifetime1); 
% length of sensor 
LR=length(Rth); 
 
% reliability of system 
N_sensor=zeros(LL,LR); 
 
for i=1:LL 
    for j=1:LR 
        N1=1; 
        N2=100; 
        while(N2-N1>1) 
            N_midium=floor((N1+N2)/2); 
            RR=0; 
            for k=1:N_midium 
                RR=RR+(C*Lifetime1(i))^(N_midium-k)/factorial(N_midium-k); 
            end 
            RR=RR*exp(-C*Lifetime1(i)); 
            if(RR>Rth(j))  
                N2=N_midium; 
            else 
                N1=N_midium; 
            end 
        end 
        N_sensor(i,j)=N1; 
    end 
end 
 
 
 
 
 
%figure(1); 
%mesh(Lifetime,N,R); 
 
% figure(2); 
% surf(Lifetime,N,R); 
 
 
%figure(3); 
%surf(N,Lifetime,R); 
 
figure(3); 
plot(N,R(:,1),'b-*'); 
hold on; 
plot(N,R(:,3),'r-o'); 
plot(N,R(:,5),'g-s'); 
plot(N,R(:,7),'r-v'); 
plot(N,R(:,9),'m-d'); 
plot(N,R(:,11),'c-<'); 
 
 
 
 
 
figure(5); 
plot(Lifetime,R(1,:),'b-*'); 
hold on; 
plot(Lifetime,R(3,:),'r-o'); 
plot(Lifetime,R(5,:),'g-s'); 
plot(Lifetime,R(7,:),'k-.'); 
plot(Lifetime,R(9,:),'m-d'); 
plot(Lifetime,R(11,:),'c-<'); 
 
 
 
% figure(7); 
% Lifetime111=[Lifetime1(1),Lifetime1(1),Lifetime1(1)]; 
% plot(N_sensor(1,:),Lifetime111,'b-s'); 
% hold on; 
% plot(N_sensor(1,2),Lifetime111(2),'r-*'); 
% Lifetime111=[Lifetime1(2),Lifetime1(2),Lifetime1(2)]; 
% plot(N_sensor(2,:),Lifetime111,'b-<'); 
% plot(N_sensor(2,2),Lifetime111(2),'m-*'); 
% Lifetime111=[Lifetime1(3),Lifetime1(3),Lifetime1(3)]; 
% plot(N_sensor(3,:),Lifetime111,'b-d'); 
% plot(N_sensor(3,2),Lifetime111(2),'g-*'); 
% Lifetime111=[Lifetime1(4),Lifetime1(4),Lifetime1(4)]; 
% plot(N_sensor(4,:),Lifetime111,'b-o'); 
% plot(N_sensor(4,2),Lifetime111(2),'c-*'); 
% 
%  