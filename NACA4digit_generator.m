clear
clc
syms p x M T a0 a1 a2 a3 a4 xc yc
%ycf(x) ycb(x) graycf(x) graycb(x) yt(x) thetaf(x) thetab(x)  xuf(x) xub(x) yuf(x) yuf(x) xlf(x) xlb(x) ylf(x) ylb(x)
%-----DATA---------
a0=0.2969;
a1=-0.126;
a2=-0.3516;
a3=0.2843;
a4=-0.1036;
%------------------

%-------done with functions-------
xc=0;
M=0.02;
p=0.4;
T=0.12;
%---------done with 4 digit-------
n=100;
p=int8(p*n);
ycb=zeros(1,p);
xcb=zeros(1,p);
thetab=zeros(1,p);
ylb=zeros(1,p);
xlb=zeros(1,p);
yub=zeros(1,p);
xub=zeros(1,p);
ytb=zeros(1,p);



ycb=M/((1-p)^2)*(1-2*p+2*p*x-x^2);
graycb=diff(ycb);

x=0:1/n:p;
ytb=5*T*(a0*sqrt(x)+a1*x+a2*x^2+a3*x^3+a4*x^4);
thetab=atan(graycb);
xub=xc-yt*sin(thetab);
yub=yc+yt*cos(thetab);
xlb=xc+yt*sin(thetab);
ylb=yc-yt*cos(thetab);

ycf=M/(p^2)*(2*p*x-x^2);
graycf=diff(ycf);
x=x:1/n:n;
thetaf=atan(graycf);
xuf=xc-yt*sin(thetaf);
yuf=yc+yt*cos(thetaf);
xlf=xc+yt*sin(thetaf);
ylf=yc-yt*cos(thetaf);

xu=zeros(1,n);
xl=zeros(1,n);
yu=zeros(1,n);
yl=zeros(1,n);

for j=1:p
    x=j/n;

    xu(j)=subs(xuf,x);
    xl(j)=subs(xlf,x);
    yu(j)=subs(yuf,x);
    yl(j)=subs(yuf,x);
end
for j=(p+1):n
    x=j/n;
    xu(j)=subs(xuf);
    xl(j)=subs(xlf);
    yu(j)=subs(yuf);
    yl(j)=subs(yuf);
end
plot(yu,xu);