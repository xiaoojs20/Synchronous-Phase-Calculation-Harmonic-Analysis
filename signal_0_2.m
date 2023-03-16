clc;clear;
close all;

signal=readmatrix('1_9.csv');
len=length(signal);
T=signal(1:end,1);
A_original=signal(1:end,2);

%totally 9600 lines
figure;
plot(T,A_original);
xlabel('时间/s');
ylabel('幅值');
title('original signal');

%Amplitude, Frequency, Harmonic Order, Sampling Frequency
A=zeros(len,1);
F = zeros(len,1);
order = zeros(len,1);
fs=10000;
L=400;
N=2200;
D=4010;
n=0:N-1;
haar=0.5-0.5*cos(2*pi*n/N);

%3954-4054
%双峰插值计算基波频率、幅值、初相
S = A_original(1:N, 1);
n = 0:N-1;
SS = S .* haar';
FT = fft(SS,N);
FT_A = abs(FT);

%第一组
for M = 1:N-1
    if(FT_A(M+1,1) > FT_A(M,1) && FT_A(M+1,1) > FT_A(M+2,1))
        M = M + 1;
        break;
    end
end
if(FT_A(M-1) > FT_A(M+1))
    k1 = M - 1;
    k2 = M;
else 
    k1 = M;
    k2 = M + 1;
end
y1 = FT_A(k1);
y2 = FT_A(k2);
belta = (y2 - y1) / (y1 + y2);
alpha = 1.5 * belta;
k0 = k1 + alpha + 0.5 - 1;
A(1) = (y1 + y2) * (2.35619403 + 1.15543628 * alpha^2 + 0.32607873 * alpha^4 + 0.07891461 * alpha^6) / N;
F(1) = k0 * fs / N;
order(1) = (angle(FT(M)) + pi/2 - pi * (alpha + (-1) * 0.5)) - pi/2;
while(order(1) < -pi)
    order(1) = order(1) + 2 * pi;
end

stem(FT_A);
display(M);

theta = 2 * pi * F(1) / fs;
tmp = 0;
for ii = 1:N
    tmp = tmp + A_original(ii) * (cos((ii-1)*theta) - 1i * sin((ii-1)*theta));
end
ang = angle(tmp);

display(ang);


%计算相位与幅值随时间的变化
%401到3600个点
for ii = L + 1:D - N
    S = A_original(ii - L + 1:ii - L + N, 1);
    SS = S .* haar';
    FT = fft(SS,N);
    FT_A = abs(FT);
    y1 = FT_A(k1);
    y2 = FT_A(k2);
    belta = (y2 - y1) / (y1 + y2);
    alpha = 1.5 * belta;
    k0 = k1 + alpha + 0.5 - 1;
    A(ii) = (y1 + y2) * (2.35619403 + 1.15543628 * alpha^2 + 0.32607873 * alpha^4 + 0.07891461 * alpha^6) / N;
    F(ii) = k0 * fs / N;
    if(ii == L + 1)
        order(ii) = (angle(FT(M)) + pi/2 - pi * (alpha + (-1) * 0.5)) - pi/2 + 2 * pi * F(ii) * (ii - L) / fs - 2 * pi * 50 * (ii - 1) / fs;
        while(order(ii) < -pi)
            order(ii) = order(ii) + 2 * pi;
        end
    else
        order(ii) = order(ii-1) - 2 * pi * (50 - F(ii)) / fs;
    end
end
%3601到4800个点
for ii = D - N + 1:D
    S = A_original(D - N + 1:D, 1);
    SS = S .* haar';
    FT = fft(SS,N);
    FT_A = abs(FT);
    y1 = FT_A(k1);
    y2 = FT_A(k2);
    belta = (y2 - y1) / (y1 + y2);
    alpha = 1.5 * belta;
    k0 = k1 + alpha + 0.5 - 1;
    A(ii) = (y1 + y2) * (2.35619403 + 1.15543628 * alpha^2 + 0.32607873 * alpha^4 + 0.07891461 * alpha^6) / N;
    F(ii) = k0 * fs / N;
    order(ii) = order(ii-1) - 2 * pi * (50 - F(ii)) / fs;
end
%4801到6000个点
for ii = D + 1:D + N
    S = A_original(D + 1:D + N, 1);
    SS = S .* haar';
    FT = fft(SS,N);
    FT_A = abs(FT);
    y1 = FT_A(k1);
    y2 = FT_A(k2);
    belta = (y2 - y1) / (y1 + y2);
    alpha = 1.5 * belta;
    k0 = k1 + alpha + 0.5 - 1;
    A(ii) = (y1 + y2) * (2.35619403 + 1.15543628 * alpha^2 + 0.32607873 * alpha^4 + 0.07891461 * alpha^6) / N;
    F(ii) = k0 * fs / N;
    order(ii) = order(ii-1) - 2 * pi * (50 - F(ii)) / fs;
end
%6001到9200个点
for ii = D + N + 1:9600 - L
    S = A_original(ii + L - N + 1:ii + L, 1);
    SS = S .* haar';
    FT = fft(SS,N);
    FT_A = abs(FT);
    y1 = FT_A(k1);
    y2 = FT_A(k2);
    belta = (y2 - y1) / (y1 + y2);
    alpha = 1.5 * belta;
    k0 = k1 + alpha + 0.5 - 1;
    A(ii) = (y1 + y2) * (2.35619403 + 1.15543628 * alpha^2 + 0.32607873 * alpha^4 + 0.07891461 * alpha^6) / N;
    F(ii) = k0 * fs / N;
    order(ii) = order(ii-1) - 2 * pi * (50 - F(ii)) / fs;
end
T = (L + 1:9600 - L).' / fs;

xlswrite('out.csv',[T,A(L + 1:9600 - L),order(L + 1:9600 - L)]);

figure;
subplot(2,1,1)
plot(T,A(L + 1:9600 - L))
title('幅值随时间变化曲线')
xlabel('时间/s');
ylabel('幅值/V');
subplot(2,1,2)
plot(T,order(L + 1:9600 - L))
title('相位随时间变化曲线')
xlabel('时间/s');
ylabel('相位/rad');
