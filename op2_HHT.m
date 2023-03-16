%��HHT���ź�˲ʱƵ�ʼ���ֵ
clc;
clear all;
signal = readmatrix('3_9.csv'); 
len=length(signal);
fs=10000;
dt=1/fs;
n=0:len-1;
t=n*dt;

A=signal(1:len,2);
% ��ϣ�����ر任
Ah=hilbert(A); 
abs_Ah=abs(Ah);
aangle=unwrap(angle(Ah));     

subplot(2,1,1)
plot(t, abs_Ah)
title('˲ʱ��ֵ����ֵ')
xlabel('ʱ��/s');
ylabel('˲ʱ��ֵ/V')

subplot(2,1,2)
plot(t, aangle)
title('˲ʱ��λ')
xlabel('ʱ��/s');
ylabel('��λ/��')



