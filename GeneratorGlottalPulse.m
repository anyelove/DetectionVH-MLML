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
