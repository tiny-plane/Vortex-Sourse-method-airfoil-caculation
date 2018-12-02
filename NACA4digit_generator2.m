a0=0.2969;
a1=-0.126;
a2=-0.3516;
a3=0.2843;
a4=-0.1036;
M=0.02;
p=0.2;
T=0.12;
n=100;
P=int32(p.*n);
%--------done--------
if(M==0)
    p=0.5;
end
x=linspace(0,p,P);
ycf=M./(p.^2).*(2.*p.*x-x.^2);
graycf=(M*(2*p - 2*x))/p^2;
yt=5*T*(a0.*sqrt(x)+a1.*x+a2.*x.^2+a3.*x.^3+a4.*x.^4);
thetaf=atan(graycf);
xuf=x-yt.*sin(thetaf);
yuf=ycf+yt.*cos(thetaf);
xlf=x+yt.*sin(thetaf);
ylf=ycf-yt.*cos(thetaf);

x=linspace(p+x(1),1,n-P);
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
%  plot(xu,yu,xl,yl,'linewidth',1);
%  hold on;
%  axis equal;
%----------done with generator-------
file=fopen('DATA1.txt','w');
fwrite(file,'NACA');
M=int8(100*M);
p=int8(10*p);
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