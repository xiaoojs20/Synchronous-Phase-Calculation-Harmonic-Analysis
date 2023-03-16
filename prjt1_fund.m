%读取数据和准备参数
Data = readmatrix("1_9.csv");
N = length(Data);
fs = 10000;
draw = 0;
f = zeros(N,1);
A = zeros(N,1);
phi = zeros(N,1);

%设置截断，进行计算，根据HHT结果确定跳变时间为0.400s
%则0.4s前，取该点前两个周期计算；0.4s后，取后两个周期计算
win_len = 400;%取样长度
for i = 401:(N-win_len)
    if i > 4000
        [f(i),A(i),phi(i)] = prjt1_fund_fun(fs,Data(i:i+win_len,2),Data(i:i+win_len,1),draw);
        phi(i) = mod(phi(i)-100*pi*Data(i,1),2*pi);
        if phi(i) > pi
            phi(i) = phi(i)-2*pi;
        end
    else
        [f(i),A(i),phi(i)] = prjt1_fund_fun(fs,Data(i-win_len:i,2),Data(i-win_len:i,1),draw);
        phi(i) = mod(phi(i)+0.08*pi*f(i)-100*pi*Data(i,1),2*pi);
        if phi(i) > pi
            phi(i) = phi(i)-2*pi;
        end        
    end
end

%绘图与结果输出
T = Data(401:N-win_len,1);
phi_out = phi(401:N-win_len);
A_out = A(401:N-win_len);

subplot(2,1,1)
plot(T,phi_out);
title('phase')
hold on;
subplot(2,1,2)
plot(T,A_out);
title('RMS')
writematrix([T A_out phi_out],'1_9_solution.csv')
