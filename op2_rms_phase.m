signal=readmatrix("3_9.csv");
len=length(signal);
fs=10000;
% 加窗长度400，最终分析范围只有200-9399
window=400; 
range=1:(len-window);
f=zeros(len,1);
A=zeros(len,1);
phase=zeros(len,1);

% 调用必做的计算函数prjt1_fund_fun
% function [f0,A,phi] = prjt1_fund_fun(fs,x,n,draw)
% 输入序列和采样频率，输出该序列对应的基波频率、幅值、相位

for k=range
    [f(k+window/2),A(k+window/2),phase(k)]=prjt1_fund_fun(fs,signal(k:k+window,2),signal(k:k+window,1),2);
    phase(k)=mod(phase(k)-50*2*pi*signal(k,1),2*pi);
end

% 基波前199个点分析不到，取第200个点63.80646344，后200个取57.49611806，但最后并不呈现出来
A(1:window/2)=A(window/2+1);
A(len-window/2+1:end)=A(len-window/2);
% 同步相量相角结果应为：对于每一时刻计算出来的绝对相位2πft+φ与参考相位2π50t的差值，相角范围[-pi, pi];
% 相角的延拓
for kk=len-window+1:len
    phase(kk)=mod(2*pi*f(len-window)*signal(kk-len+window,1)+phase(len-window)-2*pi*50*signal(kk-len+window,1),2*pi);
end

% 绘图
subplot(2,1,1)
plot(signal(window/2:len-window/2,1),A(window/2:len-window/2));
title('RMS')
hold on;
subplot(2,1,2)
plot(signal(1:len-window,1),phase(1:len-window));
title('phase')

% 写入3—_9_solution.csv
writematrix([signal(:,1) A phase],'3_9_solution.csv')

