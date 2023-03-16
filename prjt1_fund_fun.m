%�������кͲ���Ƶ�ʣ���������ж�Ӧ�Ļ���Ƶ�ʡ���ֵ����λ
function [f0,A,phi] = prjt1_fund_fun(fs,x,n,draw)

%����hamming��������Ӵ�FFT
N = length(n);
win= hamming(N);
X = fft(x.*win);

%������ͬ����˫���߲�ֵ�������Ƶ�ʣ����󦤦�����Ķ�Ӧ��ϵ
delta = 500;%���̾��ȣ������ɵ�
miu_fun = 0:(1/delta):0.5;
W = czt(win,delta+1,exp(-2*pi*j/(delta*N)),1);%�������������⦤�̵�dtft
i = 0:1:0.5*delta;%����������
alpha_fun = abs(W(end-i))./abs(W(i+1));

%ȡ��ֵ������󼶴δ�����
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

%�����Ƶ
alpha = Xnorm(k1)/Xnorm(k0);
miu = interp1(alpha_fun,miu_fun,alpha);
sign = k1-k0;
f0 = fs*(k0+sign*miu-1)/N;

%�����ֵ����λ
t = 0:N-1;
W1 = sum(win'.*exp(2*pi*j*sign*miu*t/N));
A = abs(sqrt(2)*X(k0)/W1);
phi = angle(X(k0))-angle(W1);

%���Ƶ����
if draw == 1
    k = [0:length(X)-1]/length(X);
    plot(k(1:500)*fs,pow2db(Xnorm(1:500)))
    xlabel('frequency');
    ylabel('Rms/dB');
end

    

