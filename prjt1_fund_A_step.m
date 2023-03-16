Data = readmatrix('1_9.csv');
fs = 10000;
[imf,~,~] = emd(Data(:,2));
hht(imf,fs);