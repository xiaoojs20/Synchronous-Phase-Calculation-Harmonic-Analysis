%用HHT求信号瞬时频率及幅值
clc;
clear all;
signal = readmatrix('3_9.csv'); 
len=length(signal);
fs=10000;
dt=1/fs;
n=0:len-1;
t=n*dt;

A=signal(1:len,2);
% 做希尔伯特变换
Ah=hilbert(A); 
abs_Ah=abs(Ah);
aangle=unwrap(angle(Ah));     

subplot(2,1,1)
plot(t, abs_Ah)
title('瞬时幅值绝对值')
xlabel('时间/s');
ylabel('瞬时幅值/V')

subplot(2,1,2)
plot(t, aangle)
title('瞬时相位')
xlabel('时间/s');
ylabel('相位/°')



