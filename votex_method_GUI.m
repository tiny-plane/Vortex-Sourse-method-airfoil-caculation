function varargout = votex_method_GUI(varargin)
% VOTEX_METHOD_GUI MATLAB code for votex_method_GUI.fig
%      VOTEX_METHOD_GUI, by itself, creates a new VOTEX_METHOD_GUI or raises the existing
%      singleton*.
%
%      H = VOTEX_METHOD_GUI returns the handle to a new VOTEX_METHOD_GUI or the handle to
%      the existing singleton*.
%
%      VOTEX_METHOD_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VOTEX_METHOD_GUI.M with the given input arguments.
%
%      VOTEX_METHOD_GUI('Property','Value',...) creates a new VOTEX_METHOD_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before votex_method_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to votex_method_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help votex_method_GUI

% Last Modified by GUIDE v2.5 24-Oct-2016 13:06:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @votex_method_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @votex_method_GUI_OutputFcn, ...
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
end

% --- Executes just before votex_method_GUI is made visible.
function votex_method_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to votex_method_GUI (see VARARGIN)

% Choose default command line output for votex_method_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes votex_method_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = votex_method_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end


function M_Callback(hObject, eventdata, handles)
% hObject    handle to M (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of M as text
%        str2double(get(hObject,'String')) returns contents of M as a double
end

% --- Executes during object creation, after setting all properties.
function M_CreateFcn(hObject, eventdata, handles)
% hObject    handle to M (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function P_Callback(hObject, eventdata, handles)
% hObject    handle to P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of P as text
%        str2double(get(hObject,'String')) returns contents of P as a double
end

% --- Executes during object creation, after setting all properties.
function P_CreateFcn(hObject, eventdata, handles)
% hObject    handle to P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function T_Callback(hObject, eventdata, handles)
% hObject    handle to T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T as text
%        str2double(get(hObject,'String')) returns contents of T as a double
end

% --- Executes during object creation, after setting all properties.
function T_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in caculate.
function caculate_Callback(hObject, eventdata, handles)
% hObject    handle to caculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global M p T
M=0;
p=0;
T=0;
M=str2double(get(handles.M,'string'))/100.0;
p=str2double(get(handles.P,'string'))/10.0;
T=str2double(get(handles.T,'string'))/100.0;
if(isnan(M))
    M=0;
end
if(isnan(p))
    p=0;
end
if(~isnan(T))
    plott(hObject, eventdata, handles);
else
    set(handles.conclusion,'string','这样不行哦');
    set(handles.conclusion,'visible','on');
end
end

function plott(hObject, eventdata, handles)
global M p T
LONGS = 1;
t='NACA';
t(5)=' ';
t(6)=num2str(int8(100*M));
t(7)=num2str(int8(10*p));
t(8:9)=num2str(int8(100*T));
set(handles.conclusion,'string',t);
set(handles.conclusion,'visible','on');
%------initialize-----
n=floor(get(handles.N,'Value'));
%--------done--------
if(M==0)
    p=0.5;
end
a=[0.2969 -0.126 -0.3516 0.2843 -0.1036];
m=M;
t=T;
N=n;

x=0:1/N:p;
yt=t./0.2.*(a(1).*power(x,0.5)+a(2).*x+a(3).*x.*x+a(4).*power(x,3)+a(5).*power(x,4));
z1=m/p/p.*(2.*p.*x-x.*x);
dz1=2*m/p/p.*(p-x);
td1=atan(dz1);
x1 = x - yt.*sin(td1);
x4 = x + yt.*sin(td1);
y1 = z1 + yt.*cos(td1);
y4 = z1 - yt.*cos(td1);

x = p:1/N:1;
yt = t./0.2.*(a(1).*power(x,0.5)+a(2).*x+a(3).*x.*x+a(4).*power(x,3)+a(5).*power(x,4));
z2 = m/(1-p)/(1-p).*(1-2*p+2*p.*x-x.*x);
dz2 = 2*m/(1-p)/(1-p).*(p-x);
td2 = atan(dz2);
x2 = x - yt.*sin(td2);
x3 = x + yt.*sin(td2);
y2 = z2 + yt.*cos(td2);
y3 = z2 - yt.*cos(td2);

x0 = [(x1+x4)/2 (x2+x3)/2]*LONGS;     %中弧线
y0 = [(y1+y4)/2 (y2+y3)/2]*LONGS;

x1(end)=[];y1(end)=[];   %舍去重复部分，共4个点,顺时针
x2(end)=[];y2(end)=[];
x3(1)=[];y3(1)=[];
x4(1)=[];y4(1)=[];     %该点与第一点（原点）重合，舍弃令图形不封闭，见line 50

%标准2N个点整合表示
x = [x1 x2 fliplr(x3) fliplr(x4)]*LONGS;      %整合排列，共2N项，顺时针
y = [y1 y2 fliplr(y3) fliplr(y4)]*LONGS;
 alfa=get(handles.aoa,'Value')*pi/180;
    X = zeros(1,2*N);
    Y = zeros(1,2*N);
    LX = zeros(1,2*N);
    LY = zeros(1,2*N);
    l = zeros(1,2*N);
    cosbta = zeros(1,2*N);
    

% for alf=1:Q
    %变量
    V = 1;        %自由来流，无关
%     alfa(alf) = (alf-1)/Q*25/180*pi;        %最大25度
    V_inf = [cos(alfa);sin(alfa)];


    %控制点，板块方向向量，法向量，板块长度
    for i=1:(2*N-1)
        X(i) = (x(i+1)+x(i))/2;
        Y(i) = (y(i+1)+y(i))/2;
        LX(i) = x(i+1)-x(i);
        LY(i) = y(i+1)-y(i);
        l(i) = (LX(i)^2+LY(i)^2)^(1/2);
    end
        X(2*N) = x(2*N)/2;
        Y(2*N) = y(2*N)/2;
        LX(2*N) = -x(2*N);
        LY(2*N) = -y(2*N);
        l(2*N) = (LX(2*N)^2+LY(2*N)^2)^(1/2);
    point = [X;Y];          %控制点坐标:2×2N矩阵
    L = [LX;LY];            %板块方向向量
    n = [-LY;LX];       %板块法向量，指向外:L X n > 

    %右端常数项，自由来流数为V
    for i=1:2*N
        cosbta(i) = V_inf' * n(:,i) / l(i);     %cos(beta_i)
    end
    B = (cosbta)'*2*pi;


    %系数项J[ij]，大小为2N方阵
    J = zeros(2*N);
    k = J;
    for i=1:2*N
        for j=1:2*N
            if i~=j
                if X(i)~=X(j)
                    k(i,j) = (Y(i)-Y(j))/(X(i)-X(j));
                    %积分项
                    J(i,j) = (-1/(1+k(i,j)*k(i,j))*k(i,j)/(X(i)-X(j))...
                        *n(1,i) + 1/(1+k(i,j)*k(i,j))/(X(i)-X(j))...
                        *n(2,i)) / l(i) * l(j) ;
                else
                    J(i,j) = (-1/(Y(i)-Y(j))*n(1,i)) / l(i) * l(j) ;
                end
            end
        end
    end

    J(N+1,:) = 0;           %选取后缘一点，代换为库塔条件
    J(N+1,N) = 1;
    J(N+1,N+1) = 1;
    B(N+1) = 0;

    %J(2*N,:) = 0;          %选取前缘一点，代换为后缘库塔条件
    %J(2*N,N) = 1;
    %J(2*N,N+1) = 1;
    %B(2*N) = 0;

    %解
    GAMA = J\B;         %涡分布，2N*1 向量
    C_Lift = 2/V/LONGS * (l * GAMA);

cla
axes(handles.airfoil);
plot([x 0],[y 0]);
%plot(x,y,'markersize',200)
hold on;
%plot(x0,y0);
axis equal;
%---------plot----------
t=num2str(C_Lift);
set(handles.solution,'string',t);
set(handles.solution,'visible','on');
set(handles.cl,'visible','on');
end


% --- Executes on slider movement.
function aoa_Callback(hObject, eventdata, handles)
set(handles.angle,'string',num2str(get(handles.aoa,'Value')));
% hObject    handle to aoa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
end

% --- Executes during object creation, after setting all properties.
function aoa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to aoa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

% --- Executes on slider movement.
function N_Callback(hObject, eventdata, handles)
set(handles.Nn,'string',num2str(floor(get(handles.N,'Value'))));
% hObject    handle to N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
end

% --- Executes during object creation, after setting all properties.
function N_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end
