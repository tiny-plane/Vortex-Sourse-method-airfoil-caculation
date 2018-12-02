
clear
clc
pi=3.14159265;
%------initialize-----
file=fopen('DATA1.txt','r');
fscanf(file,'%c',4);
M=fscanf(file,'%1d',1)/100.0;
p=fscanf(file,'%1d',1)/10.0;
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
alpha=0;
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
% for iop=1:1
% cla;
% plot(xu,yu,xl,yl,'linewidth',2);
% hold on;
% axis equal;
% for i=n:2*n-2
%     if ni(i)>0
%     line([conpoi(1,i),conpoi(1,i)-0.1./sqrt(ni(i)^2+1)],[conpoi(2,i),conpoi(2,i)-0.1*ni(i)./sqrt(ni(i)^2+1)],'color','black');
%     else  if ni(i)<=0
%     line([conpoi(1,i),conpoi(1,i)+0.1./sqrt(ni(i)^2+1)],[conpoi(2,i),conpoi(2,i)+0.1*ni(i)./sqrt(ni(i)^2+1)],'color','black');
%         end
%     end
%        line([xl(i-n+2),xl(i-n+2)+0.1./sqrt(tan(beta(i))^2+1)],[yl(i-n+2),yl(i-n+2)+0.1*tan(beta(i))./sqrt(tan(beta(i))^2+1)],'linestyle',':','linewidth',3,'color','red');
% end
% for i=1:n-1
%     if ni(i)<0
%     line([conpoi(1,i),conpoi(1,i)-0.1./sqrt(ni(i)^2+1)],[conpoi(2,i),conpoi(2,i)-0.1*ni(i)./sqrt(ni(i)^2+1)],'color','black');
%     else  if ni(i)>=0
%     line([conpoi(1,i),conpoi(1,i)+0.1./sqrt(ni(i)^2+1)],[conpoi(2,i),conpoi(2,i)+0.1*ni(i)./sqrt(ni(i)^2+1)],'color','black');
%         end
%     end
%     line([xu(i+1),xu(i+1)+0.1./sqrt(tan(beta(i))^2+1)],[yu(i+1),yu(i+1)+0.1*tan(beta(i))./sqrt(tan(beta(i))^2+1)],'linestyle',':','linewidth',3,'color','red');
% end
% end
%---------plot----------
gama=J\freestreamspeedn';
total=0;
for i=1:2*n-2
    total=total+gama(i)*length(i);
end
total
%-----
tao=0;
for i=1:2*n-2
    tao=tao+gama(i).*length(i);
end
-tao*2/vinfinity