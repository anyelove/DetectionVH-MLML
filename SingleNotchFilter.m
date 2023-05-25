clc;
close all;
clear all;

%% Single Notch Filter

Fc = 700;
Fs = 8192;
theta_c = 2*pi*(Fc/Fs); % 0 < Fc < Fs/2
B = 100;
beta = exp(-pi*B/Fs);
TimeInSample = length(signal_1);

H_snf = zeros(1,TimeInSample);

for nn = 3:1:TimeInSample
H_snf(nn) = (signal_1(nn) - (2*beta*cos(theta_c))*(signal_1(nn-1)) + signal_1(nn-2))/(2*(1- beta*cos(theta_c)));
end

b = [1 -(2*beta*cos(theta_c)) 1];
a = [sum(b) 0 0];

y = filter(a,b,signal_1);

[b,a] = butter(4,2000/(Fs/2));
y = filtfilt(b,a,y);
figure
hold on
%plot(n,signal_1)
plot(n,y)
hold off
grid on
box on

w = [y y y y y y y y y y y y y y y y y y y];
m = zeros(1,length(w));
for j = 1:1:length(w)
    m(j) = (j-1) * Ts;
end

figure
hold on
plot (m, w)
hold off
grid on
box on

w = w./(max(w));

figure
hold on
plot (m, w)
hold off
grid on
box on

figure
hold on
plot (m(1:end-1), diff(w))
title('')
hold off
grid on
box on

sound(diff(w),Fs);
%sound(w,Fs);
%clear sound;