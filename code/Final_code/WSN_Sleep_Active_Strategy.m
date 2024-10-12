function varargout = WSN_Sleep_Active_Strategy(varargin)
% WSN_SLEEP_ACTIVE_STRATEGY MATLAB code for WSN_Sleep_Active_Strategy.fig
%      WSN_SLEEP_ACTIVE_STRATEGY, by itself, creates a new WSN_SLEEP_ACTIVE_STRATEGY or raises the existing
%      singleton*.
%
%      H = WSN_SLEEP_ACTIVE_STRATEGY returns the handle to a new WSN_SLEEP_ACTIVE_STRATEGY or the handle to
%      the existing singleton*.
%
%      WSN_SLEEP_ACTIVE_STRATEGY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WSN_SLEEP_ACTIVE_STRATEGY.M with the given input arguments.
%
%      WSN_SLEEP_ACTIVE_STRATEGY('Property','Value',...) creates a new WSN_SLEEP_ACTIVE_STRATEGY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before WSN_Sleep_Active_Strategy_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to WSN_Sleep_Active_Strategy_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help WSN_Sleep_Active_Strategy

% Last Modified by GUIDE v2.5 02-Apr-2019 21:56:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @WSN_Sleep_Active_Strategy_OpeningFcn, ...
                   'gui_OutputFcn',  @WSN_Sleep_Active_Strategy_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before WSN_Sleep_Active_Strategy is made visible.
function WSN_Sleep_Active_Strategy_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to WSN_Sleep_Active_Strategy (see VARARGIN)

% Choose default command line output for WSN_Sleep_Active_Strategy
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes WSN_Sleep_Active_Strategy wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = WSN_Sleep_Active_Strategy_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global net SVM T S n
WB=waitbar(0.8,'Generating testing set... ');
pd1 = makedist('Normal','mu',T,'sigma',0.035*T);
pd2 = makedist('Normal','mu',0.4*T,'sigma',0.1*T);
Temp = [random(pd1,fix(n/3),1);random(pd2,fix(n/20),1);
        random(pd1,n-fix(n/3)-fix(n/20),1)];
pd1 = makedist('Normal','mu',S,'sigma',0.035*S);
pd2 = makedist('Normal','mu',0.4*S,'sigma',0.1*S);
backet_size = [random(pd1,fix(n/3),1);random(pd2,fix(n/20),1);
        random(pd1,n-fix(n/3)-fix(n/20),1)];
% Separating good nodes from attacked nodes
Status = (Temp>0.9*T)+(Temp<1.1*T)+(backet_size>0.9*S)+(backet_size<1.1*S);
idx = Status==4; Status=-1*ones(n,1); Status(idx)=1;
% Testing the SVM model
Pr = predict(SVM,[Temp,backet_size]);
error = Status-Pr; 
error=error~=0; error=sum(error);
disp(['Testing error = ' num2str(error) ' Nodes']);
close(WB);
axes(handles.axes1); hold on;
attack = find(Status==-1); attack_pr = find(Pr==-1);
plot(net(2,attack),net(3,attack),'ro',...
    'MarkerSize',5,'MarkerFaceColor','r');
plot(net(2,attack_pr),net(3,attack_pr),'ro',...
    'MarkerSize',10,'LineWidth',2);
hold off


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global net SVM T S n
n = str2double(get(handles.edit1,'String'));    %Number of nodes
w = str2double(get(handles.edit2,'String'));    %length of the network
h = str2double(get(handles.edit3,'String'));    %width of the network
T = str2double(get(handles.edit4,'String'));    %Temperature
S = str2double(get(handles.edit5,'String'));    %Size of data backet

R = 0.25*sqrt(w*h);      %radio range
% Creating the network
net = [1:n;rand([1,n])*w;rand([1,n])*h];
% Ploting the network
% plot(net(2,:),net(3,:),'ko','MarkerSize',5,'MarkerFaceColor','k');
% title('Distributed nodes');
% xlabel('\it x \rm [m] \rightarrow');
% ylabel('\it y \rm [m] \rightarrow');
hold on;
axes(handles.axes1)
for i = 1:numel(net(1,:))
    for j = 1:numel(net(1,:))
        X1 = net(2,i);
        Y1 = net(3,i);
        X2 = net(2,j);
        Y2 = net(3,j);
        xSide = abs(X2-X1);
        ySide = abs(Y2-Y1);
        d = sqrt(xSide^2+ySide^2);       
        if (d<R)&&(i~=j)
            vertice1 = [X1,X2];
            vertice2 = [Y1,Y2];
            plot(vertice1,vertice2,'--b','LineWidth',0.1);
            hold on;
        end
    end
end

v = net(1,:)';
s = int2str(v);
text(net(2,:)+1,net(3,:)+1,s,'FontSize',8,'VerticalAlignment','Baseline');
hold off
WB=waitbar(0.8,'Generating Training set... ');

pd1 = makedist('Normal','mu',T,'sigma',0.05*T);
pd2 = makedist('Normal','mu',0.4*T,'sigma',0.1*T);
Temp = [random(pd2,fix(n/3),1);random(pd1,fix(n/3),1);random(pd2,n-2*fix(n/3),1)];
pd1 = makedist('Normal','mu',S,'sigma',0.05*S);
pd2 = makedist('Normal','mu',0.4*S,'sigma',0.1*S);
backet_size = [random(pd2,fix(n/3),1);random(pd1,fix(n/3),1);random(pd2,n-2*fix(n/3),1)];
clc; %clear command window
% Separating good nodes from attacked nodes
Status = (Temp>0.9*T)+(Temp<1.1*T)+(backet_size>0.9*S)+(backet_size<1.1*S);
idx = Status==4; Status=-1*ones(n,1); Status(idx)=1;
Sensor_read = Temp;
training_set=table(Sensor_read,backet_size,Status) %Show training set generated
% Training the SVM model
SVM = fitcsvm([Temp,backet_size],Status);
error = Status-predict(SVM,[Temp,backet_size]); 
error=error~=0; error=100*sum(error)/n;
disp(['Training error = ' num2str(error) '%']);
%%%%%%%%%%%%%%%%%%Own design
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
 
 
Lifetime1=n; 
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
 
axes(handles.axes2); hold on; 
plot(N,R(:,1),'b-*'); 
hold on; 
plot(N,R(:,3),'r-o'); 
plot(N,R(:,5),'g-s'); 
plot(N,R(:,7),'r-v'); 
plot(N,R(:,9),'m-d'); 
plot(N,R(:,11),'c-<'); 
xlabel('Frequency');
ylabel('Throughput');
 
 
 
 
 
axes(handles.axes3); hold on; 
plot(Lifetime,R(1,:),'b-*'); 
% hold on; 
plot(Lifetime,R(3,:),'r-o'); 
plot(Lifetime,R(5,:),'g-s'); 
plot(Lifetime,R(7,:),'k-.'); 
plot(Lifetime,R(9,:),'m-d'); 
plot(Lifetime,R(11,:),'c-<'); 
xlabel('No of Nodes');
ylabel('Energy Consumption in J');
%%%%%%%%%%%%end
%%%%%%%%%%%%%no of dead nodes
%-------------------------------------------------------------------------
%Protocol Name :  WBAN and Modified Local Emergency Detection (LED) with Adaptive Sampling.

%-------------------------------------------------------------------------
% 

% clc
% close all
% clear all
% clear

close(WB);
set(handles.pushbutton1,'Enable','on')



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
flag_first_deadl=0;
flag_teenth_deadl=0;
flag_all_deadl=0;
p = 0.1 ;
pl=0.1;
Prop_delay_SUM = 0 ;
path_loss_SUM = 0 ;
packet_to_BS_SUM = 0; 
packet_drop_SUM = 0;
 dead_SUM = 0 ;
 E_SUM = 0; 
 allive_SUM = 0;
 packet_rcvd_AVG = 0;
nl = 8;
 P_loss_SUM = 0 ;
PACKETS_DROPD_SUM=0;
PKTS_TO_BS_SUM=0;
DEAD_SUM=0;
E_TOTAL_SUM=0;
Pro_delay_SUM = 0 ;

 
% packet_to_BS_SUM = 0; 
% packet_drop_SUM = 0;
%  dead_SUM = 0 ;
%  E_SUM = 0; 
%  allive_SUM = 0;
%  packet_rcvd_AVG = 0;
for i = 1:1:5
    
xylabel=20;
legendsize=18;

sink.x = .25;
sink.y = 1 ;
n = 8;
Eo = 0.5 ;
% energy parameters
ETX=16.7*0.000000001;
ERX=36.1*0.000000001;
Emp=1.97*0.000000001;
EDA=5*0.000000001;
do = 0.1 ;

lambda = .125 ;% f = 2.4 GHz
% speed = 299792458 ;
speed = 2997924587877576 ;
flag_first_dead=0;         
flag_teenth_dead=0;
flag_all_dead=0;

dead=0;
first_dead=0;
teenth_dead=0;
all_dead=0;


flag_first_dead1=0;         
flag_teenth_dead1=0;
flag_all_dead1=0;

dead1=0;
first_dead1=0;
teenth_dead1=0;
all_dead1=0;

dead2 = 0 ;
thr = 0.1 ;
C = 299792458 ;

packet_to_BS = 0;
packet_to_CH = 0;
Packet_to_BS_total = 0 ;

% path_loss_round = 0 ;
% Prop_delay_round = 0;
% total_delay_round = 0;


packet_drop = 0;

rmax = 8000;
allive  = n;
%-------------------------------------------------------------------------
%                   ATTEMP Data 
%-------------------------------------------------------------------------

xm=0.8;   %2.5 feet
ym=1.8;  %6 feet

 %x and y Coordinates of the Sink
    sink1.x=0.4;
    sink1.y=0.9;
    %Number of Nodes in the field
node = 8 ;
alliveA  = node;
    %x and y Coordinates of the Sink
 sink1.x=0.4;
 sink1.y=0.9;
 th=0.7;
 P_loss_round = 0 ;
 P_delay_round = 0 ;
 P_loss = 0 ;
    th_temp=98;
    rng=0.8062;
   % Network establishment
    SA(1).xd=0.2;
    SA(1).yd=1.2;
    SA(2).xd=0.6;
    SA(2).yd=1.1;
    SA(3).xd=0.7;
    SA(3).yd=0.8;
    SA(4).xd=0.5;
    SA(4).yd=0.6;
    SA(5).xd=0.1;
    SA(5).yd=0.8;
    SA(6).xd=0.3;
    SA(6).yd=0.5;
    SA(7).xd=0.5;
    SA(7).yd=0.3;
    SA(8).xd=0.3;
    SA(8).yd=0.1;
    for i = 1:1:node
%        plot(SA(i).xd, SA(i).yd,'*r');
%        hold on
    end
%     plot(sink.x ,sink.y, 'bo' );
%   hold on ;
%   grid on ;
%   figure(1);
for i = 1:1:node
    distanceA(i)=sqrt((SA(i).xd-(sink1.x) )^2 + (SA(i).yd-(sink1.y) )^2 );
end
    % figure(1)
    %% Energy eqquipment and type identification
    for i=1:1:node
        SA(i).type='N';
        SA(i).E=Eo;
    end
    SA(node+1).xd=sink1.x;
    SA(node+1).yd=sink1.y;
    countCHs=0;         %the number of Stateflow objects in the current context.
    cluster=1;              %first cluster is selected
    dead=0;
    allive=node;
    %counter for bit transmitted to Bases Station and to Cluster Heads
    packets_TO_BS=0;
    packets_TO_CH=0;
    Paskets_TO_BS_total=0;
    pd=0;
    d=0;
    s=0;


%----------------------------------------------------------
% node deployment

    S(1).xd=0.3;  % calf
    S(1).yd=0.1;
    S(1).P = 1;
    
    S(2).xd=0.5;     % knee
    S(2).yd=0.3;
    S(2).P = 1;
          
   S(3).xd=0.3;    % lactic acid(thigh) 2.5
   S(3).yd=0.55;
   S(3).P = 1;
   
   S(4).xd= .5; % temp(thigh)  2.5
   S(4).yd= .55;
   S(4).P = 1;
   
   S(7).xd= .37; % glucose
   S(7).yd= .75;
   S(7).P = 2;
   
   S(8).xd= .45; % EEG
   S(8).yd= .9;
   S(8).P = 2;
    
   S(5).xd= .7; % left PALM 4.3 feet 
   S(5).yd= .8;
   S(5).P = 1;
   
   S(6).xd= .1; %  rite palm   5.4
   S(6).yd= .8;
   S(6).P = 1;
  x2 = [];    
   %--------------------------------------------------------------------------
   
for i = 1:1:n
    S(i).E = Eo ;
    S(i).id = i;
    S(i).g = 0 ;
    axis on ;
    axis ([0 0.8 0 2]) ;
%   plot(S(i).xd ,S(i).yd, 'r*' );
%   hold on ;
%   grid on ;
end
%--------------------------------------------------------------------------
% plot(sink.x ,sink.y, 'bo' );
    a = 1; 
    b = 1; 
x0 = [];
x1 = [];

            for i = 1:1:n 
               if S(i).P ==1 
                  x0{a} = S(i).id; 
                  a = a+1 ;
               end
            end
            
   %%%%%modified         
      for i=1:1:n
            for j=1:1:n
           distance1(i,j)=sqrt((S(i).xd-(S(j).xd))^2 + (S(i).yd-(S(j).yd) )^2 );
            distance(i)=sqrt((S(i).xd-(sink.x) )^2 + (S(i).yd-(sink.y) )^2 );
            end
      end
      %%%%%%%%%%
        pd1 = 0;
        Packets_to_BS_total = 0;
        total_delay_round = 0;
         

                
        for r=1:1:rmax
            r

path_loss = 0;
P_loss = 0;
PL = 0 ;
delay = 0;
dead = 0;
E =0 ;
E1 = 0 ;

for i=1:1:n

   if (S(i).E<=0)
       dead=dead+1;
       if (dead==1)
          if(flag_first_dead==0)
             first_dead=r;
             flag_first_dead=1;
          end
       end
             if(dead==n)
                 if(flag_all_dead==0)
                    all_dead=r;
                    flag_all_dead=1 ;
                 end
             end
       
   end
   if S(i).E>0
       S(i).g = 0;
       E = E+S(i).E ;
   end
        
end
%E1 = E/(Eo*8) *100 ;

Prop_delay =0 ;
packet_to_BS_per_round = 0;
packet_to_BS = 0 ;
Dead(r+1)=dead;
Allive(r+1)=n-dead;
energy(r+1)=E;


%-------------------------------------------------------------------------
% cost function calculation to select forwarder
cost_function = 0;
%cost_function1 = 0;
for i=1:1:length(x0)
             if(S(i).E > 0 )
             cost_function(i) = distance(i)/(S(i).E) ;
             end
end
%-------------------------------------------------------------------------         
for i=1:1:length(x0)
 % node with minimum cost function elected  as forwarder
      [min_node,I] = min(cost_function);
      node_num = I ;
      node_sel(r) = node_num ;
end
% if energy of node is greater then threshold energy
if (S(node_num).E>thr )
         S(node_num).g = 1 ;  % forwarder
          packet_to_BS = packet_to_BS+1 ;
          distanceCH=sqrt( (S(node_num).xd-(sink.x) )^2 + (S(node_num).yd-(sink.y) )^2 );
          S(node_num).E = S(node_num).E - ( (ETX+ERX+EDA)*(4000) + Emp*3.38*4000*(distanceCH^3.38));
          PL(node_num)=10*log(((4*pi*do)/lambda)^2+10*4*log(distanceCH/do))+4.1 ;
          delay(node_num) = distanceCH/C ;
         
end

% node forwading data to neighboure faorwarder
for i=1:1:length(x0)
            if ( S(i).g == 0) && (S(i).E>thr && S(i).P ==1 )
                
            temp = sqrt((S(i).xd - S(node_num).xd)^2 + (S(i).yd-S(node_num).yd)^2 ) ;
             temp2 = sqrt((S(i).xd - sink.x)^2 + (S(i).yd-sink.y)^2 ) ;
            
                if (temp <temp2)
                S(i).E = S(i).E - ( (ETX)*(4000) + Emp*3.38*4000*(temp^3.38));  % 3.38 => human body path loss exponent
                packet_to_CH = packet_to_CH+1 ;
                PL(i)=10*log(((4*pi*do)/lambda)^2+10*4*log(temp/do))+4.1 ;
                path_loss = path_loss + PL(i);
                
                delay(i) = temp/C ;
                Prop_delay = Prop_delay + delay(i) ;

                else
                S(i).E = S(i).E - ( (ETX)*(4000) + Emp*3.38*4000*(temp2^3.38));  % 3.38 => human body path loss exponent
                packet_to_BS = packet_to_BS+1; 
                PL(i)=10*log(((4*pi*do)/lambda)^2+10*4*log(temp2/do))+4.1 ;
               path_loss = path_loss + PL(i);
               
               delay(i) = temp2/C ;
                Prop_delay = Prop_delay + delay(i) ;

                end
            end
end
 %-----------------------------------------------------------------------
            for i = 1:1:n
            if (S(i).E > 0 && S(i).P == 2 )
            d1 = sqrt((S(i).xd - sink.x)^2 + (S(i).yd-sink.y)^2 ) ;
            S(i).E = S(i).E - ( (ETX)*(4000) + Emp*3.38*4000*(d1^3.38));  % 3.38 => human body path loss exponent
            packet_to_BS = packet_to_BS+1 ;
            PL(i)=10*log(((4*pi*do)/lambda)^2+10*4*log(d1/do))+4.1 ;
            path_loss = path_loss+PL(i) ;
            delay(i) = d1 /C ;
             Prop_delay = Prop_delay + delay(i) ;
            end
            end
 
%            direct transmission if energy of nodes decreses below threshold 
             for i = 1:1:length(x0)
            if (S(i).E < thr && S(i).E > 0 )
            d = sqrt((S(i).xd - sink.x)^2 + (S(i).yd-sink.y)^2 ) ;
            S(i).E = S(i).E - ( (ETX)*(4000) + Emp*3.38*4000*(d^3.38));  % 3.38 => human body path loss exponent
            packet_to_BS = packet_to_BS+1 ;
          packet_to_BS = packet_to_BS+1 ;
            PL(i)=10*log(((4*pi*do)/lambda)^2+10*4*log(d/do))+4.1 ;
            delay(i) = d/C ;
            Prop_delay = Prop_delay + delay(i) ;
            
            end
             end
        %--------------------------------------------------
        % Modified LED Protocol 
              path_loss_round(r+1) = path_loss ;
             Prop_delay_round(r+1) = Prop_delay ;

         packet_to_BS_per_round(r+1)=packet_to_BS;
        Packet_to_BS_total=Packet_to_BS_total+packet_to_BS;
        packet_to_BS1(r+1)=Packet_to_BS_total;
        
          
        pkts_rcvd1=0;
        pkts_drpd1=0;
        P_opt1=0.3;      %Optimal probability to determine link status
        for j=1:1: packet_to_BS_per_round(r+1)
            pr=rand(1,1);
            if(pr>=P_opt1)
                lnk_status_flag1='good';
                pkts_rcvd1=pkts_rcvd1+1;
            else
                lnk_status_flag1='bad';
                pkts_drpd1=pkts_drpd1+1;
            end
        end
         pd1 = pd1+pkts_drpd1;
        packet_drop(r+1) =pd1;
        %----------------------------------------------------------------
%        Random uniformed model of packet drops ends here
 
%--------------------------------------------------------------------------
%                           MODIFIED LED
%--------------------------------------------------------------------------
deadA=0;
        E_total=0;
        for i=1:1:node
            
            if (SA(i).E<=0)
                deadA=deadA+1;
            end
            if(SA(i).E>0)
                SA(i).type='N';
                E_total=E_total+SA(i).E;
            end
            
        end
%        E_ATTM = E_total/(Eo*8) *100 ;

        E_TOTAL(r+1)=E_total;
        DEADA(r+1)=deadA;
        ALLIVEA(r+1)=alliveA-deadA;
        packets_TO_BS_per_round=0;
        packets_TO_BSA=0;
        
        delay_ATT = 0 ;
        total_delay = 0;
         
        for i=1:1:node
            ct=90+(rand*10);   % ct is between 90 to 100
            if (SA(i).E>0 ) % && ct<th_temp
                cv=rand;
                if (cv>th)
                    distA=sqrt( (SA(i).xd-(SA(node+1).xd) )^2 + (SA(i).yd-(SA(node+1).yd) )^2 );
                      SA(i).E=SA(i).E- ( (ETX)*(4000) + Emp*3.38*4000*(distA^3.38));
                      if(distA > 0)
                      PL_ATT(i)=10*log(((4*pi*do)/lambda)^2+10*4*log(distA/do))+4.1;
                       P_loss = P_loss + PL_ATT(i);
                       delay_ATT(i) = distA/C ;
                       total_delay = total_delay+delay_ATT(i) ;

                      end
         
                    %S(i).E=S(i).E- ( (ETX)*(4000) + Emp*4000*(dist*dist));
                    d=i;
                    packets_TO_BSA=packets_TO_BSA+1;
                end
                if(cv<=th)
                    min_dis=Inf;
                    min_dis_cluster=0;
                    for c=1:1:node
                        if i~= c
                        tempA=min(min_dis,sqrt( (SA(i).xd-SA(c).xd)^2 + (SA(i).yd-SA(c).yd)^2 ) );
                        if ( tempA<min_dis )
                            min_dis=tempA;
                            min_dis_cluster=c;
                        end
                        end
                    end
                    disA=sqrt( (SA(i).xd-(SA(node+1).xd) )^2 + (SA(i).yd-(SA(node+1).yd) )^2 );
                    if (min_dis<disA)
                        disA=sqrt( (SA(i).xd-(SA(min_dis_cluster).xd) )^2 + (SA(i).yd-(SA(min_dis_cluster).yd) )^2 );
                     SA(i).E=SA(i).E- ( (ETX)*(4000) + Emp*3.38*4000*(disA^3.38));
                     
                     PL_ATT(i)=10*log(((4*pi*do)/lambda)^2+10*4*log(disA/do))+4.1 ;
                      P_loss = P_loss + PL_ATT(i);
                      delay_ATT(i) = disA/C ;
                      total_delay = total_delay+delay_ATT(i) ;
                     
         
                    %    S(i).E=S(i).E- ( (ETX)*(4000) + Emp*4000*(dis*dis));
                        s=i;
                    end
                    
                    if(min_dis>disA)
                        disA=sqrt( (SA(i).xd-(SA(n+1).xd) )^2 + (SA(i).yd-(SA(n+1).yd) )^2 );
                        
                        SA(i).E=SA(i).E- ( (ETX)*(4000) + Emp*3.38*4000*(disA^3.38));
                        
                        PL_ATT(i)=10*log(((4*pi*do)/lambda)^2+10*4*log(disA/do))+4.1 ;
                        P_loss = P_loss + PL_ATT(i);
                        delay_ATT(i) = disA/C ;
                        total_delay = total_delay+delay_ATT(i) ;
                        
         
         
         %               S(i).E=S(i).E- ( (ETX)*(4000) + Emp*4000*(dis*dis ));
                        packets_TO_BSA=packets_TO_BSA+1;
                    end
                end
            end
        end
        if (SA(4).E>0)
            SA(4).E=SA(4).E- ( (ERX+EDA)*(4000) );
        end
        if (SA(5).E>0)
            SA(5).E=SA(5).E- ( (ERX+EDA)*(4000) );
        end
        if (SA(7).E>0)
            SA(7).E=SA(7).E- ( (ERX+EDA)*(4000) );
        end
           P_loss_round(r+1) = P_loss ;
           total_delay_round(r+1)= total_delay ;
             %Prop_delay_round(r+1) = Prop_delay ;
        packets_TO_BS_per_round(r+1)=packets_TO_BSA;
        Paskets_TO_BS_total=Paskets_TO_BS_total+packets_TO_BSA;
        PKTS_TO_BS(r+1)=Paskets_TO_BS_total;
        %%   Random uniformed model of packet drops starts from here
        pkts_rcvd=0;
        pkts_drpd=0;
        P_opt=0.3;      %Optimal probability to determine link status
        for j=1:1:packets_TO_BS_per_round(r+1)
            p_r=rand(1,1);
            if(p_r>=P_opt)
                lnk_status_flag='good';
                pkts_rcvd=pkts_rcvd+1;
            else
                lnk_status_flag='bad';
                pkts_drpd=pkts_drpd+1;
            end
        end
        pd=pd+pkts_drpd;
        PACKETS_DROPPED(r+1)= pd;
        % Random uniformed model of packet drops ends here

        end
        
           %Adding five times simulation resuts(Proposed)
           Prop_delay_SUM = Prop_delay_round+Prop_delay_SUM ;
           path_loss_SUM =  path_loss_round+ path_loss_SUM ;
          packet_to_BS_SUM=packet_to_BS1+packet_to_BS_SUM;
          packet_drop_SUM=packet_drop+packet_drop_SUM;
          dead_SUM=Dead+dead_SUM;
           E_SUM=E_SUM+energy;
          allive_SUM = allive_SUM +Allive ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
% 5 times execution

     ddd(i,r) = total_delay_round(r) ;
   
     
            Pro_delay_SUM = total_delay_round+Pro_delay_SUM ;
           P_loss_SUM =  P_loss_round+ P_loss_SUM ;

    PKTS_TO_BS_SUM=PKTS_TO_BS+PKTS_TO_BS_SUM;
    PACKETS_DROPD_SUM=PACKETS_DROPPED+PACKETS_DROPD_SUM;
    DEAD_SUM=DEADA+DEAD_SUM;
    E_TOTAL_SUM=E_TOTAL_SUM+E_TOTAL;
    %ALIVE_SUM = ALIVE_SUM + 
end
%%  Calculating average values(PROPOSED)


path_loss_AVG = path_loss_SUM/5 ;
Pro_delay_AVG = Pro_delay_SUM/5 ;
 packet_to_BS_AVG = packet_to_BS_SUM/5;
 packet_drop_AVG = packet_drop_SUM/5;
 packet_rcvd_AVG = packet_to_BS_AVG-packet_drop_AVG;
 dead_AVG = dead_SUM/5;
E_AVG = E_SUM/5;
Allive_AVG = allive_SUM/5 ;
% pr = 1:300:rmax ;
% E_AVG = E_AVG(pr);

%%  Calculating average values (ATTEMPT)
Prop_delay_AVG = Prop_delay_SUM/5 ;
P_loss_AVG = P_loss_SUM/5 ;
PKTS_TO_BS_AVG=PKTS_TO_BS_SUM/5;
PACKETS_DROPD_AVG=PACKETS_DROPD_SUM/5;
PACKETS_RCVD_AVG=PKTS_TO_BS_AVG-PACKETS_DROPD_AVG;
DEAD_AVG=DEAD_SUM/5;
E_TOTAL_AVG=E_TOTAL_SUM/5;
 pat = 1:305:rmax ;
 E_TOTAL_AVG = E_TOTAL_AVG(pat);

%% plotting
STATISTICS.ENERGY(r+1)= pd ;

r=0:rmax;

 axes(handles.axes4); hold on; 
% plot(Lifetime,R(1,:),'b-*'); 
plot(r,dead_AVG,'r-',r,DEAD_AVG,'b-','linewidth',2);
legend('SIMPLE','ATTEMPT');
xlabel('Rounds (r)','FontSize',xylabel,'FontName','Arial');
ylabel('No. of dead nodes','FontSize',xylabel,'FontName','Arial');
axis([0 8000 0 12]);
grid on

x = r;
 [E_AVG]= hist(x);
%  bar(E_AVG)
Y = [pr,E_AVG] ;
bar(Y);

 for r= 0:1:rmax
    r
     if(mod(r, round(1/p) )==0) %remainder
         for i = 1:1:length(x2)
             S2(i).type = 0;
             S3(i).type = 0;
         end
     end
end
 if(mod(r, round(1/pl) )==0) %remainder
        for i=1:1:nl
            Sl(i).G=0;            % it will assign to the nodes that have not been cluster head .
        end
    end
Etl=0;
E=0;
    deadl=0;
    for i=1:1:nl
        
         if (SA(i).E<=0)
            deadl=deadl+1;
            
            if (deadl==1)
                if(flag_first_deadl==0)
                    first_deadl=r;
                    flag_first_deadl=1;
                end
            end
            
            if(deadl==0.1*nl)
                if(flag_teenth_deadl==0)
                    teenth_deadl=r;
                    flag_teenth_deadl=1;
                end
            end
            if(deadl==n)
                if(flag_all_deadl==0)
                    all_deadl=r;
                    flag_all_deadl=1;
                end
            end
        end
        if SA(i).E>0
            Sl(i).type='N';                      
            Etl = Etl+Sl(i).E ;
            
        end
     end

STATISTICS.El(r+1)=Etl;
    STATISTICS.DEAD(r+1)=dead;
    STATISTICS.ALLIVE(r+1)=allive-dead;
    alive1l = allive-dead ;
    received_packsl =0;
 
 
 %-------------------------------------------------------------------------
 %                    MLED Ends here
 %-------------------------------------------------------------------------
% STATISTICS.recieved(r+1)
%STATISTICS.packet_delivery_ratio(r+1)
 STATISTICS.El(r+1)
STATISTICS.DEAD(r+1)
STATISTICS.ALLIVE(r+1)
% STATISTICS.PACKETS_TO_BSl(r+1)

%------------------------
%STATISTICS.packet_delivery_ratiol(r+1)
% STATISTICS.received(r+1)
%  STATISTICS.PACKETS_TO_BS(r+1)
 STATISTICS.DEAD(r+1)
 STATISTICS.ALLIVE(r+1)
%   STATISTICS.E(r+1)
  STATISTICS.ENERGY(r+1)
 %------------------------------------------------------------
 r=0:rmax;
% plot(r,STAT.P,'r');
 
%  figure(10)
%  plot(r,STATISTICS.DEAD,'r',r,DEADA,'b');
%  xlabel('Rounds','FontSize',xylabel,'FontName','Arial')
%  ylabel('Dead','FontSize',xylabel,'FontName','Arial')
%  legend1=legend('Proposed','MLED');
%  set(legend1,'FontSize',legendsize)

%  axes(handles.axes5); hold on; 
%  plot(Lifetime,R(1,:),'b-*');
axes(handles.axes5); hold on; 
 plot(r,STATISTICS.ALLIVE,'r',r,ALLIVEA,'b');
 xlabel('rounds','FontSize',xylabel,'FontName','Arial')
 ylabel('Percent of Allive Nodes','FontSize',xylabel,'FontName','Arial')
legend1=legend('Proposed','MLED');
set(legend1,'FontSize',legendsize)

% figure(12)
%  plot(r,STATISTICS.E,'r',r,STATISTICS.El,'b','linewidth',2)
%  xlabel('rounds','FontSize',xylabel,'FontName','Arial')
%  ylabel('Energy','FontSize',xylabel,'FontName','Arial')
% legend1=legend('Proposed','LEACH');
% set(legend1,'FontSize',legendsize)


% figure(12)
%  plot(r,STATISTICS.ENERGY,'r',r,E_TOTAL,'b');
%  xlabel('Rounds','FontSize',xylabel,'FontName','Arial')
%  ylabel('Remaining Energy','FontSize',xylabel,'FontName','Arial')
%  legend1=legend('Proposed','MLED');
%  set(legend1,'FontSize',legendsize)
%%%%%%%%%end