function varargout = vortex_panel_method(varargin)
% SOURCE_METHOD_GUI MATLAB code for source_method_GUI.fig
%      SOURCE_METHOD_GUI, by itself, creates a new SOURCE_METHOD_GUI or raises the existing
%      singleton
%
%      H = SOURCE_METHOD_GUI returns the handle to a new SOURCE_METHOD_GUI or the handle to
%      the existing singleton*.
%
%      SOURCE_METHOD_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SOURCE_METHOD_GUI.M with the given input arguments.
%
%      SOURCE_METHOD_GUI('Property','Value',...) creates a new SOURCE_METHOD_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before source_method_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to source_method_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help source_method_GUI

% Last Modified by GUIDE v2.5 30-Oct-2016 18:13:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @source_method_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @source_method_GUI_OutputFcn, ...
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

% --- Executes just before source_method_GUI is made visible.
function source_method_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to source_method_GUI (see VARARGIN)

% Choose default command line output for source_method_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes source_method_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = source_method_GUI_OutputFcn(hObject, eventdata, handles) 
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
t='NACA';
t(5)=' ';
t(6)=num2str(int8(100*M));
t(7)=num2str(int8(10*p));
if (T<0.1)
    t(8)=num2str(0);
    t(9)=num2str(int8(100*T));
else
t(8:9)=num2str(int8(100*T));
end
set(handles.conclusion,'string',t);
set(handles.conclusion,'visible','on');
%------initialize-----
a0=0.2969;
a1=-0.126;
a2=-0.3516;
a3=0.2843;
a4=-0.1036;
n=floor(get(handles.N,'Value'));
%--------done--------
if(M==0)
    p=0.5;
end
P=int32(p.*n);
nn=int32(n/10);
nnn=double(nn);
x(1:nn)=linspace(0,0.05,nn);
x(nn+1:P)=linspace(x(nn)+1.0/(n-nnn),double(P-1)/(n-1),P-nnn);
ycf=M./(p.^2).*(2.*p.*x-x.^2);
graycf=(M*(2*p - 2*x))/p^2;
yt=5*T*(a0.*sqrt(x)+a1.*x+a2.*x.^2+a3.*x.^3+a4.*x.^4);
thetaf=atan(graycf);
xuf=x-yt.*sin(thetaf);
yuf=ycf+yt.*cos(thetaf);
xlf=x+yt.*sin(thetaf);
ylf=ycf-yt.*cos(thetaf);


x=linspace(x(P)+1/double(n-nn),1,n-P);
ycb=M/((1-p)^2)*(1-2*p+2*p.*x-x.^2);
graycb=(M*(2*p - 2*x))/(p - 1)^2;
yt=5*T*(a0.*sqrt(x)+a1.*x+a2.*x.^2+a3.*x.^3+a4.*x.^4);
thetab=atan(graycb);
xub=x-yt.*sin(thetab);
yub=ycb+yt.*cos(thetab);
xlb=x+yt.*sin(thetab);
ylb=ycb-yt.*cos(thetab);
%--------done with calculation-----
xu=zeros(1,n);
yu=zeros(1,n);
xl=zeros(1,n);
yl=zeros(1,n);
for i=1:P
    xu(i)=xuf(i);
    yu(i)=yuf(i);
    xl(i)=xlf(i);
    yl(i)=ylf(i);
end
for i=P+1:n-1
    j=i-P;
    xu(i)=xub(j);
    yu(i)=yub(j);
    xl(i)=xlb(j);
    yl(i)=ylb(j);
end
for i=1:n-1
    if xu(i)==xu(i+1)
        xu(i)=xu(i-1)./2+xu(i+1)./2;
        xl(i)=xl(i-1)./2+xl(i+1)./2;
        yu(i)=yu(i-1)./2+yu(i+1)./2;
        yl(i)=yl(i-1)./2+yl(i+1)./2;
    end
end
xu(n)=1;
xl(n)=1;
yu(n)=0;
yl(n)=0;
%set(gca,'unit','normalized');
%  plot(xu,yu,xl,yl,'linewidth',1);
%  hold on;
%  axis equal;
%----------done with generator-------
file=fopen('DATA1.txt','w');
fwrite(file,'NACA');
M=int8(100*M);
p=int8(10*p);
if M==0
    p=0;
end
T=int8(100*T);
fprintf(file,'%d%d%d',M,p,T);
fprintf(file,'\r\nTotal points=%d\r\nNumber                  xu                  yu                  xl                  yl\r\n',n);
fclose(file);
%---------done with file intro part----
file=fopen('DATA1.txt','a');
for i=1:n
    fprintf(file,'\r\n%3d---%20d%20d%20d%20d',i,xu(i),yu(i),xl(i),yl(i));
end
fclose(file); 
%--------end-------
caculate_fun(hObject, eventdata, handles);
end

function caculate_fun(hObject, eventdata, handles)
pi=3.14159265;
LONGS = 1;
%------initialize-----
file=fopen('DATA1.txt','r');
%file=fopen('C:\Users\asus-pc\Documents\Visual Studio 2013\Projects\Console7\NACA.txt','r');
fscanf(file,'%c',4);
M=fscanf(file,'%1d',1)/100.0;
p=fscanf(file,'%1d',1)/10.0;
if p==0
    p=0.05;
end
T=fscanf(file,'%2d',1)/100.0;
fscanf(file,'%c',15);
n=fscanf(file,'%d',1);
fscanf(file,'%s',5);

xu=zeros(1,n);
xl=zeros(1,n);
yu=zeros(1,n);
yl=zeros(1,n);
A=[0;0;0;0];
for i=1:n
    fscanf(file,'%s',1);
    A=fscanf(file,'%f',4);
    A=A';  
    xu(i)=A(1);
    yu(i)=A(2);
    xl(i)=A(3);
    yl(i)=A(4);
end
fclose(file);
a=[0.2969 -0.126 -0.3516 0.2843 -0.1036];
m=M;
t=T;
N=n-1;

x(1:N)=xu(1:n-1);
x(N+1:2*N)=xl(n:-1:2);
y(1:N)=yu(1:n-1);
y(N+1:2*N)=yl(n:-1:2);

 alfa=get(handles.aoa,'Value')*pi/180;
    X = zeros(1,2*N);
    Y = zeros(1,2*N);
    LX = zeros(1,2*N);
    LY = zeros(1,2*N);
    l = zeros(1,2*N);
    cosbta = zeros(1,2*N);

    V_inf = [cos(alfa);sin(alfa)];

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

    

    GAMA = J\B  ;      %涡分布，2N*1 向量
    C_Lift = 2/LONGS * (l * GAMA);
    
    n=N+1;
alpha=get(handles.aoa,'Value');
vinfinity=1;
%-----done with initialize-----
% plot(xu,yu,xl,yl,'linewidth',2);
% hold on;
% axis equal;
%-----this part only for test-----
% gama=zeros(2,n-1);
conpoi=zeros(2,2*n-2);
length=zeros(1,2*n-2);
for i=1:n-1
    conpoi(1,i)=(xu(i)+xu(i+1))/2.0;
    conpoi(2,i)=(yu(i)+yu(i+1))/2.0;
    length(i)=sqrt((yu(i+1)-yu(i))^2+(xu(i+1)-xu(i))^2);
    length(i+n-1)=sqrt((yl(i+1)-yl(i))^2+(xl(i+1)-xl(i))^2);
    conpoi(1,i+n-1)=(xl(i)+xl(i+1))/2.0;
    conpoi(2,n+i-1)=(yl(i)+yl(i+1))/2.0;
end
pandir=zeros(1,2*n-2);
for i=1:n-1
    pandir(i)=(yu(i+1)-yu(i))./(xu(i+1)-xu(i));
    pandir(i+n-1)=(yl(i+1)-yl(i))./(xl(i+1)-xl(i));
end
for i=1:2*n-2
    if(isnan(pandir(i)))
        pandir(i)=(pandir(i-1)+pandir(i+1))./2.0;
    end
end
ni=zeros(1,2*n-2);
for i=1:n-1
    ni(i)=pandir(i)+pi/2;
    ni(n-1+i)=pandir(n-1+i)-pi/2;
end
for i=1:2*n-2
    ni(i)=tan(ni(i));
end
freestreamdirection=alpha/180*pi;
beta=zeros(1,2*n-2);
for i=1:2*n-2                  %beta(i)表示自由来流和板块切向夹角，即来流实际对板块的作用角度
    pandir(i)=atan(pandir(i));
    beta(i)=pandir(i)-freestreamdirection;
end
freestreamspeedn=zeros(1,2*n-2);
for i=1:2*n-2          %-----------流入为负-------
    freestreamspeedn(i)=-vinfinity.*sin(beta(i));
end
thetaij=zeros(2*n-2,2*n-2);
len=zeros(2*n-2,2*n-2);
for i=1:2*n-2
    for j=1:2*n-2
        if(i~=j)
        thetaij(i,j)=(conpoi(2,j)-conpoi(2,i))./(conpoi(1,j)-conpoi(1,i)); 
        thetaij(i,j)=atan(thetaij(i,j));
        len(i,j)=(conpoi(2,j)-conpoi(2,i)).^2+(conpoi(1,j)-conpoi(1,i))^2;
        len(i,j)=sqrt(len(i,j));
        end
    end
end
J=zeros(2*n-2,2*n-2);
for i=1:2*n-2
    for j=1:2*n-2
        if i~=j
            J(i,j)=(conpoi(1,i)-conpoi(1,j).*cos(ni(i))+(conpoi(2,i)-conpoi(2,j)).*sin(ni(i)))./len(i,j).*length(j).^2./2./pi;
        else
            J(i,j)=0.5;
        end
     end
end
% J(2*n-2,:)=0;
% J(:,2*n-2)=0;
% J(2*n-2,n-1)=1;
% J(2*n-2,2*n-2)=1;
% freestreamspeedn(2*n-2)=0;
%----j induce to i-----
% gama(1:N)=GAMA(1:N);
% gama(N+1:2*N)=GAMA(2*N:-1:N+1);
gama=J\freestreamspeedn';
corgama=zeros(2*n-2,3);
maxg=max(gama);
ming=min(gama);
pol=maxg-ming;
numg=(gama-ming)/pol;
lengama=zeros(1,2*n-2);
times=0.07;
for i=1:2*n-2
    if (numg(i)<=0.5)
        corgama(i,1:3)=[(1-2*numg(i))^0.25,(2*numg(i))^3,0];
    else
        corgama(i,1:3)=[0,(2-2*numg(i))^3,(2*numg(i)-1)^0.25];
    end
    if (gama(i)>=0)
        lengama(i)=times*sqrt(gama(i))./sqrt(ni(i)^2+1);
    else
        lengama(i)=-times*sqrt(-gama(i))./sqrt(ni(i)^2+1);
    end
end
cla;
axes(handles.airfoil);
plot(xu,yu,xl,yl,'linewidth',2);
hold on;
%scatter(conpoi(1,:),conpoi(2,:));
axis equal;
if(n<=300)
for i=1:2*n-2
    if ni(i)>0
    line([conpoi(1,i),conpoi(1,i)-lengama(i)],[conpoi(2,i),conpoi(2,i)-ni(i).*lengama(i)],'color',corgama(i,1:3),'linewidth',3);
    else
    line([conpoi(1,i),conpoi(1,i)+lengama(i)],[conpoi(2,i),conpoi(2,i)+ni(i).*lengama(i)],'color',corgama(i,1:3),'linewidth',3);
    end
     %  line([xl(i-n+2),xl(i-n+2)+0.1./sqrt(tan(beta(i))^2+1)],[yl(i-n+2),yl(i-n+2)+0.1*tan(beta(i))./sqrt(tan(beta(i))^2+1)],'linestyle',':','linewidth',3,'color','red');
end
% for i=1:n-1
%     if ni(i)<0
%     line([conpoi(1,i),conpoi(1,i)-lengama(i)],[conpoi(2,i),conpoi(2,i)-lengama(i).*ni(i)],'color',corgama(i,1:3),'linewidth',3);
%     else
%     line([conpoi(1,i),conpoi(1,i)+lengama(i)],[conpoi(2,i),conpoi(2,i)+lengama(i).*ni(i)],'color',corgama(i,1:3),'linewidth',3);
%     end
%  %   line([xu(i+1),xu(i+1)+0.1./sqrt(tan(beta(i))^2+1)],[yu(i+1),yu(i+1)+0.1*tan(beta(i))./sqrt(tan(beta(i))^2+1)],'linestyle',':','linewidth',3,'color','red');
% end
end
%---------plot----------

%-----
tao=0;
for i=1:2*n-2
    tao=tao+gama(i).*length(i);
end
lift=C_Lift;
t=num2str(lift);
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
