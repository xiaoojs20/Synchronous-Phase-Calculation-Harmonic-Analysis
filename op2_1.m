%双谱线插值法第二版，求了A、F、P瞬时值，用于加分2
clc;
clear all;
wave = xlsread('3_9.csv'); %load 数据
A1=zeros(8799,1);
F1=zeros(8799,1);
P1=zeros(8799,1);
P=zeros(8799,1);

%从第401个点（0.0401s）开始到第9199（0.9199s）结束
for m1=401:9199
    s=wave(m1-400:m1+400,2);
    fs=10000; %采样率
    N=length(s);%采样点数
    n=0:N-1;
    w=0.5-0.5*cos(2*pi*(n)/N);%汉宁窗
    r=s.*w';%对原信号加窗
    v=fft(r,N);%进行FFT
    vz=abs(v)/N*2*2;%求幅值并修正
    u=abs(v);
%     stem(abs(v)/N*2*2);%画出FFT后信号的茎状图，用于判断k0、k1、k2
    %以下步骤为双峰谱线插值修正算法
    y1=u(43);%事先通过前面的茎状图发现基波位于u(6)
    y2=u(44);
    y3=u(42);
    max=y2;
    k1=43;
    k2=44;
    if y3>y2
        max=y3;
    end
    if max==y3
        t=y1;
        y1=max;
        y2=t;
        k1=42;
        k2=43;
    end
    b=(y2-y1)/(y2+y1);%相当于参考文献中的参数β
    a=1.5*b;%相当于参考文献中的参数α
    k0=k1-1+a+0.5;%按理来说不用减1，但不知为啥-1后似乎更准
    t=(m1-400)/fs;
    A1(m1-400)=(y1+y2)*(2.35619403+1.15543628*a^2+0.32607873*a^4+0.07891461*a^6)/N;
    F1(m1-400)=k0*fs/N;
    P1(m1-400)=(angle(v(42))+pi/2-pi*(a-(-1)*0.5))-pi/2-2*pi*50*t;
end

T=401:9199;
t=T/fs;
for index=1:8192
    P(index)=2*pi*(F1(index)-50)*t(index)+P1(1);
end
% plot(t,F1)
subplot(2,1,1)
plot(t,A1)
title('幅值随时间变化曲线')
xlabel('时间/s');
ylabel('幅值/V');
subplot(2,1,2)
plot(t,P)
title('相位随时间变化曲线')
xlabel('时间/s');
ylabel('相位/rad');