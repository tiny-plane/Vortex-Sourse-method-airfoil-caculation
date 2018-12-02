lift=0.123456789;
j=5;
liftn=zeros(1,j);
for i=1:j
    liftn(i)=int8(floor(lift*10^i));
    lift=lift-liftn(i)/10^i;
end