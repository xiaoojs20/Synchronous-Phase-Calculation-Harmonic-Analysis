%画出茎叶图
%双谱线插值，用于求必做基波及谐波的平均频率、幅值、相位

clc;
clear all;
signal=readmatrix('1_9.csv'); %load 数据
A_original=signal(1:2048,2);%由于在3900点之后会有幅值、频率阶跃，故不能直接对全部信号进行FFT，先截取前2048个点分析
%s=wave(4353:6400,2);%阶跃后取4353-6400
fs=10000; 
N=length(A_original);% 采样点数
n=0:N-1;
w=0.5-0.5*cos(2*pi*(n)/N);% 这里加了汉宁窗，后续加汉明窗
w_signal=A_original.*w';

%a.*b表示矩阵a中的元素与矩阵b中的元素按位置依次相乘，得到的结果将作为新矩阵中相同位置的元素
fft_sig=fft(w_signal,N);%进行FFT，返回N点的DFT
fftsig_z=(abs(fft_sig)/N)*2*2;% 修正系数为2
u=abs(fft_sig);
stem(fftsig_z);% 画出FFT后信号的茎状图，用于判断k0、k1、k2

A=zeros(1,30);% 谐波幅值
F=zeros(1,30);% 谐波频率
P=zeros(1,30);% 谐波相位
freq=zeros(1,30);% 谐波次数

%以下步骤为双峰谱线插值修正算法
M=11;%事先通过前面的茎状图发现基波位于u(11)
for I=0:29 %I+1对应谐波系数，由于只含有小于30次的谐波，故在29截止即可
    %选择最大、次最大的点
    if(u(M-1+(M-1)*I) > u(M+1+(M-1)*I))
        k1=M-1+(M-1)*I;
        k2=M+(M-1)*I;
    else 
        k1=M+(M-1)*I;
        k2=M+1+(M-1)*I;
    end
    y1 = u(k1);
    y2 = u(k2);
    % 插值
    b=(y2-y1)/(y2+y1);% beta
    a=1.5*b;% alpha
    k0=k1-1+a+0.5;
    A(I+1)=(y1+y2)*(2.35619403+1.15543628*a^2+0.32607873*a^4+0.07891461*a^6)/N;
    F(I+1)=k0*fs/N;
    P(I+1)=(angle(fft_sig(11+10*I))+pi/2-pi*(a-(-1)*0.5))/pi*180;
end

for i=1:29
    for j=1:29
         if (F(i)/F(1)>=j*0.99) && (F(i)/F(1)<=j*1.01)% 只有误差在基波频率的j倍的98%-102%以内才算谐波
             freq(i)=j;
         end
    end
end
for i=1:29
    if(freq(i)==0)
        freq(i)=0;
        A(i)=0;
        F(i)=0;
        P(i)=0;
    end
    if(A(i)<0.0005)
        freq(i)=0;
        A(i)=0;
        F(i)=0;
        P(i)=0;
    end
end

main_freq=find(freq>0.5);

A=A(main_freq);
F=F(main_freq);
P=P(main_freq);

total=[main_freq;A;F;P];
display(total);

xlswrite('必做谐波分析.csv',total');