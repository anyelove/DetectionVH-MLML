close all;
clear all;
clear vars;
clc;

%% Modelo de Rosenberg

Fs = 8192; %Frecuencia de muestreo modelo

Ts = 1/Fs; %Tiempo de cada muestra

a = 1; %Amplitud


% % Tiempo apertura
% 
% T_milis = 5000; % tiempo en milisegundos
% T_muestras = T_milis*Fs/1000 + 1; % tiempo en samples
% Tp = 0.40 * (T_muestras-1); % porcentaje del tiempo de apertura con pendiente positiva(0.40 sugerido)
% Tn = 0.16 * (T_muestras-1); % porcentaje del tiempo de apertura con pendiente negativa(0.16 sugerido)

%% Glottal Pulse

%  GlottalPulse(N1,N2,N3)
%
% N1 : opening slope in samples
% N2 : closing slope in samples
% N3 : pulse shape index in Rosenberg Paper
%
% N1=25; % TP in Rosenberg paper
% N2=10; % TN-TP in Rosenberg paper
% N3= 4; % Pulse shape index in Rosenberg Paper a=1,b=2,c=3,d=4,e=5,f=6

RandNumber1 = (2*rand)-1;   % Variacion aleatoria para tiempos de apertura y cierre
RandNumber2 = (2*rand)-1;
factor = 0.05;      % Porcentaje maximo de variacion

Tp = 30-(floor(30 * factor * RandNumber1));
Tn = 10-(floor(10 * factor * RandNumber2));

x = GlottalPulse(Tp,Tn,4);

%% Pulse Train 

% Pulse Train Generator
% N : Size Vector in samples
% P : Period of Pulse in samples
% jitter : in percentage of period
% shimmer : in percentage of amplitude
% PulseTrain(N,P,jitter,shimmer)

TotalPulseTrainTime = 3;   % total time in seconds
PulseTrainRep = 0.015;  % between 150 y 300 (1/150 - 1/300)% time of repetitions in seconds, > Tp+Tn+1
Ts_train = floor(Fs * TotalPulseTrainTime);   % total time of pulse train in samples
Ts_rep = floor(Fs * PulseTrainRep); % time between repetitions in samples
jj = PulseTrain(Ts_train,Ts_rep,0.01,0);

signal_1 = conv(x,jj);
signal_2 = conv(diff(x),jj);
%% Vector tiempo

n = zeros(1,length(signal_1));
for j = 1:1:length(signal_1)
    n(j) = (j-1) * Ts;
end

figure('Name','Convolution Glottal Pulse and Pulse Train','NumberTitle','off');
hold on
plot(n,signal_1,'m-')
xlabel('Time [s]')
ylabel('Amplitude')
hold off

figure('Name','Convolution Flow Accel of Glottal Pulse and Pulse Train','NumberTitle','off');
hold on
plot(n(1:end-1),signal_2,'m-')
xlabel('Time [s]')
ylabel('Amplitude')
hold off

%sound(signal_1,Fs);
%sound(signal_2,Fs);

%% Single Notch Filter

Fc = 300;
Fs = 8192;
B = 50;

[voice_signal_2, b1, a1] = SNF_Inv(signal_2, Fc, Fs, B, 1)

figure
freqz(b1,a1);


%sound(voice_signal_1,Fs);
sound(voice_signal_2,Fs);

Y_1 = fft(voice_signal_2);
%Y_2 = fft(signal_2);
L_1 = length(Y_1);
%L_2 = length(Y_2);
P2_1 = abs(Y_1/L_1);
%P2_2 = abs(Y_2/L_2);
P1_1 = P2_1(1:L_1/2+1);
%P1_2 = P2_2(1:L_2/2+1);
P1_1(2:end-1) = 2*P1_1(2:end-1);
%P1_2(2:end-1) = 2*P1_2(2:end-1);

f_1 = Fs*(0:(L_1/2))/L_1;
%f_2 = Fs*(0:(L_2/2))/L_2;

figure
plot(f_1,P1_1) 
title("Single-Sided Amplitude Spectrum of X(t)")
xlabel("f_1 (Hz)")
ylabel("|P1_1(f)|")

%figure
%plot(f_2,P1_2) 
%title("Single-Sided Amplitude Spectrum of X(t)")
%xlabel("f_2 (Hz)")
%ylabel("|P1_2(f)|")