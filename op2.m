M=readmatrix("3_9.csv");
fs=10000;
len=length(M);
window=400;  % 可以修改 window值 观察效果
k_range=1:(len-window);
f=zeros(len,1);
A=zeros(len,1);
phase=zeros(len,1);
for k=k_range
    [f((2*k+window)/2),A((2*k+window)/2),phase(k)]=myCal_FreFundamental(M(k:k+window,2),M(k:k+window,1),fs);
    phase(k)=mod(phase(k)-50*2*pi*M(k,1),2*pi);
end
%% 该部分代码是为补全 A和phase 并不必要
A(1:window/2)=A(window/2+1);
A(len-window/2+1:end)=A(len-window/2);
for j=len-window+1:len
    phase(j)=mod(2*pi*f(len-window)*M(j-len+window,1)+phase(len-window)-2*pi*50*M(j-len+window,1),2*pi);
end

%% 绘图
subplot(2,1,1)
plot(M(window/2:len-window/2,1),A(window/2:len-window/2));
title('RMS')
hold on;
subplot(2,1,2)
plot(M(1:len-window,1),phase(1:len-window));
title('phase')
%% write csv
writematrix([M(:,1) A phase],'3_9_solution.csv')

