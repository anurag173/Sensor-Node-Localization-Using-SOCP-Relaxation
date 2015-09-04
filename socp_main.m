addpath(genpath('C:\Users\Anurag\Dropbox\UGP 2014\results'))
sdpsettings('solver','SEDUMI')
n=500;k=1;
for i=.1:.05:.35
  %  tmp(k)=localization_socp(n,.15,.15,0,i);
[err1(k) err2(k) err3(k)]=localization_socp(n,.15,.15,0,i);
%localization_socp(500,.15,.15);
k=k+1;
end
t=.1:.05:.35;
h=figure;
err1
err2
err3

%plot(t,err1)