%输入序列和采样频率，输出该序列对应的基波频率、幅值、相位
function [f0,A,phi] = prjt1_fund_fun(fs,x,n,draw)

%构建hamming窗，计算加窗FFT
N = length(n);
win= hamming(N);
X = fft(x.*win);

%采样不同步，双谱线插值计算基波频率，先求Δμ与α的对应关系
delta = 500;%Δμ精度，参数可调
miu_fun = 0:(1/delta):0.5;
W = czt(win,delta+1,exp(-2*pi*j/(delta*N)),1);%窗函数关于任意Δμ的dtft
i = 0:1:0.5*delta;%遍历计数用
alpha_fun = abs(W(end-i))./abs(W(i+1));

%取幅值，找最大级次大级谱线
Xnorm = abs(X);
Xm = max(Xnorm);
k0 = find(Xnorm==Xm,1);%??
if k0 == 1
    k1 = 2;
elseif k0 == N
    k1 = N-1;
elseif Xnorm(k0+1)>Xnorm(k0-1)
    k1 = k0+1;
else 
    k1 = k0-1;   
end

%计算基频
alpha = Xnorm(k1)/Xnorm(k0);
miu = interp1(alpha_fun,miu_fun,alpha);
sign = k1-k0;
f0 = fs*(k0+sign*miu-1)/N;

%计算幅值和相位
t = 0:N-1;
W1 = sum(win'.*exp(2*pi*j*sign*miu*t/N));
A = abs(sqrt(2)*X(k0)/W1);
phi = angle(X(k0))-angle(W1);

%绘幅频特性
if draw == 1
    k = [0:length(X)-1]/length(X);
    plot(k(1:500)*fs,pow2db(Xnorm(1:500)))
    xlabel('frequency');
    ylabel('Rms/dB');
end

    

