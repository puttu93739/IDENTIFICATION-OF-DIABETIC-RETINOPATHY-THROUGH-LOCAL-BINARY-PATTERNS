function varargout = WSN_Sleep(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @WSN_Sleep_OpeningFcn, ...
                   'gui_OutputFcn',  @WSN_Sleep_OutputFcn, ...
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


% --- Executes just before WSN_Sleep is made visible.
function WSN_Sleep_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for WSN_Sleep
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = WSN_Sleep_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function edit1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit3_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% --- Executes on button press in pushbutton1. (Detect attacked)
function pushbutton1_Callback(hObject, eventdata, handles)
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
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% --- Executes on button press in pushbutton2 (Generate network)
function pushbutton2_Callback(hObject, eventdata, handles)
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
plot(net(2,:),net(3,:),'ko','MarkerSize',5,'MarkerFaceColor','k');
title('Distributed nodes');
xlabel('\it x \rm [m] \rightarrow');
ylabel('\it y \rm [m] \rightarrow');
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

close(WB);
set(handles.pushbutton1,'Enable','on')
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



function edit4_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit5_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
