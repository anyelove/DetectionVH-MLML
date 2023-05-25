clc;
clear all;
close all;

% precisar datos de alguna referencia formant chart
% repetir proceso con signal_2
% graficar (?)
%% Modelo de Rosenberg

Fs = 8192; %Frecuencia de muestreo modelo

Ts = 1/Fs; %Tiempo de cada muestra


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

TotalPulseTrainTime = 1;   % total time in seconds
PulseTrainRep = 0.015;  % between 150 y 300 (1/150 - 1/300)% time of repetitions in seconds, > Tp+Tn+1
Ts_train = floor(Fs * TotalPulseTrainTime);   % total time of pulse train in samples
Ts_rep = floor(Fs * PulseTrainRep); % time between repetitions in samples
jj = PulseTrain(Ts_train,Ts_rep,0.01,0);

signal_1 = conv(x,jj);
signal_2 = conv(diff(x),jj);


numsamples = 400; 


%% Vowel /i/

%Parametros para SNF respecto a formant chart

% First formant f0
meanf1_i = 390; % en Hz
meanbw1_i = 100;
varf1_i = 100; % en Hz
varbw1_i = 25;

% Second formant f1
meanf2_i = 2300;
meanbw2_i = 100;
varf2_i = 100;
varbw2_i = 50;

% Crear matrices de 2 x length(numsamples) en el que se muestren las
% frecuencias resonantes de SNF con su respectivo BW

for mm = 1:numsamples
    f_bw_1_i(1,mm) = randn(1,1) * varf1_i + meanf1_i; 
    f_bw_1_i(2,mm) = randn(1,1) * varbw1_i + meanbw1_i;
    f_bw_2_i(1,mm) = randn(1,1) * varf2_i + meanf2_i;
    f_bw_2_i(2,mm) = randn(1,1) * varbw2_i + meanbw2_i;
end


%% SNF para f0 y f1, /i/

signal_1_ivoice = zeros(numsamples, length(signal_1));


for nn = numsamples:numsamples
    [signal_1_i b1 a1] = SNF_Inv(signal_1, f_bw_1_i(1,nn), Fs, f_bw_1_i(2,nn), 1); % Primer Filtro SNF para f0
    figure
    freqz(b1,a1,512,Fs)
    [signal_1_i b1 a1] = SNF_Inv(signal_1_i, f_bw_2_i(1,nn), Fs, f_bw_2_i(2,nn), 1); %Segundo Filtro SNF para f1
    for aa = 1:length(signal_1)
        signal_1_ivoice(nn,aa) = signal_1_i(1,aa); %Matriz de numsamples x length(signal_voice) para guardar las señales correspondientes a la síntesis de voz
    end
    
    
end

%sound(signal_1_ivoice(randi(numsamples),1:end), Fs);

figure
freqz(b1,a1)

%% Vowel /a/

% Parametros para SNF respecto a formant chart

% Los parámetros a utilizar para la vocal /a/ según "blabla (2019)" son:
% para f0, meanf1_a = x , frecuencia promedio de la primera formante,
% varf1_a = varianza para producir la base de datos
% meanbw1_a = ancho de banda promedio del filtro para la primera formante.
% varbw1_a = varianza para producir la base de datos

% Para f1, proceso equivalente reemplazando los indices por 2


% First formant f0
varf1_a = 100; % en Hz
varbw1_a = 25;
meanf1_a = 300; % en Hz
meanbw1_a = 70;

% Second formant f1
varf2_a = 700;
varbw2_a = 25;
meanf2_a = 2800;
meanbw2_a = 70;

for mm = 1:numsamples
    f_bw_1_a(1,mm) = randn(1,1) * varf1_a + meanf1_a;
    f_bw_1_a(2,mm) = randn(1,1) * varbw1_a + meanbw1_a;
    f_bw_2_a(1,mm) = randn(1,1) * varf2_a + meanf2_a;
    f_bw_2_a(2,mm) = randn(1,1) * varbw2_a + meanbw2_a;
end

%% SNF para f0 y f1, /a/

signal_1_avoice = zeros(numsamples, length(signal_1));

for nn = 1:numsamples
    [signal_1_a b1 a1] = SNF_Inv(signal_1, f_bw_1_a(1,nn), Fs, f_bw_1_a(2,nn), 1); % SNF_Inv(signal_2, Fc, Fs, B, 1)
    [signal_1_a b1 a1] = SNF_Inv(signal_1_a, f_bw_2_a(1,nn), Fs, f_bw_2_a(2,nn), 1);
    for aa = 1:length(signal_1)
        signal_1_avoice(nn,aa) = signal_1_a(1,aa);
    end
end

%sound(signal_1_avoice(randi(numsamples),1:end), Fs);
figure
freqz(b1,a1)

%% Vowel /e/

%Parametros para SNF respecto a formant chart

% First formant f0
varf1_e = 100; % en Hz
varbw1_e = 25;
meanf1_e = 300; % en Hz
meanbw1_e = 70;

% Second formant f1
varf2_e = 700;
varbw2_e = 25;
meanf2_e = 2800;
meanbw2_e = 70;


for mm = 1:numsamples
    f_bw_1_e(1,mm) = randn(1,1) * varf1_e + meanf1_e;
    f_bw_1_e(2,mm) = randn(1,1) * varbw1_e + meanbw1_e;
    f_bw_2_e(1,mm) = randn(1,1) * varf2_e + meanf2_e;
    f_bw_2_e(2,mm) = randn(1,1) * varbw2_e + meanbw2_e;
end

%% SNF para f0 y f1, /e/

signal_1_evoice = zeros(numsamples, length(signal_1));

for nn = 1:numsamples
    [signal_1_e b1 a1] = SNF_Inv(signal_1, f_bw_1_e(1,nn), Fs, f_bw_1_e(2,nn), 1); % SNF_Inv(signal_2, Fc, Fs, B, 1)
    [signal_1_e b1 a1] = SNF_Inv(signal_1_e, f_bw_2_e(1,nn), Fs, f_bw_2_e(2,nn), 1);
    for aa = 1:length(signal_1)
        signal_1_evoice(nn,aa) = signal_1_e(1,aa);
    end
end

%sound(signal_1_evoice(randi(numsamples),1:end), Fs);
figure
freqz(b1,a1)


%% Vowel /o/

%Parametros para SNF respecto a formant chart

% First formant f0
varf1_o = 100; % en Hz
varbw1_o = 25;
meanf1_o = 300; % en Hz
meanbw1_o = 70;

% Second formant f1
varf2_o = 700;
varbw2_o = 25;
meanf2_o = 2800;
meanbw2_o = 70;


for mm = 1:numsamples
    f_bw_1_o(1,mm) = randn(1,1) * varf1_o + meanf1_o;
    f_bw_1_o(2,mm) = randn(1,1) * varbw1_o + meanbw1_o;
    f_bw_2_o(1,mm) = randn(1,1) * varf2_o + meanf2_o;
    f_bw_2_o(2,mm) = randn(1,1) * varbw2_o + meanbw2_o;
end

%% SNF para f0 y f1, /o/

signal_1_ovoice = zeros(numsamples, length(signal_1));

for nn = 1:numsamples
    [signal_1_o b1 a1] = SNF_Inv(signal_1, f_bw_1_o(1,nn), Fs, f_bw_1_o(2,nn), 1); % SNF_Inv(signal_2, Fc, Fs, B, 1)
    [signal_1_o b1 a1] = SNF_Inv(signal_1_o, f_bw_2_o(1,nn), Fs, f_bw_2_o(2,nn), 1);
    for aa = 1:length(signal_1)
        signal_1_ovoice(nn,aa) = signal_1_o(1,aa);
    end
end

%sound(signal_1_ovoice(randi(numsamples),1:end), Fs);
figure
freqz(b1,a1)


%% Vowel /u/

%Parametros para SNF respecto a formant chart

% First formant f0
varf1_u = 100; % en Hz
varbw1_u = 25;
meanf1_u = 300; % en Hz
meanbw1_u = 70;

% Second formant f1
varf2_u = 700;
varbw2_u = 25;
meanf2_u = 2800;
meanbw2_u = 70;


for mm = 1:numsamples
    f_bw_1_u(1,mm) = randn(1,1) * varf1_u + meanf1_u;
    f_bw_1_u(2,mm) = randn(1,1) * varbw1_u + meanbw1_u;
    f_bw_2_u(1,mm) = randn(1,1) * varf2_u + meanf2_u;
    f_bw_2_u(2,mm) = randn(1,1) * varbw2_u + meanbw2_u;
end

%% SNF para f0 y f1, /u/

signal_1_uvoice = zeros(numsamples, length(signal_1));

for nn = 1:numsamples
    [signal_1_u b1 a1] = SNF_Inv(signal_1, f_bw_1_u(1,nn), Fs, f_bw_1_u(2,nn), 1); % Primer Filtro SNF para f0
    [signal_1_u b1 a1] = SNF_Inv(signal_1_u, f_bw_2_u(1,nn), Fs, f_bw_2_u(2,nn), 1); %Segundo Filtro SNF para f1
    for aa = 1:length(signal_1)
        signal_1_uvoice(nn,aa) = signal_1_u(1,aa);
    end
end

%sound(signal_1_uvoice(randi(numsamples),1:end), Fs);
figure
freqz(b1,a1)



%% Graficos

Y_1 = fft(signal_1_avoice(randi(numsamples), 1:end));
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

% Vector

for vv = 1:length(signal_1_avoice)
    v_tiempo(1,vv) = vv * (1/Fs);
end

figure
plot(v_tiempo, (signal_1_avoice(randi(numsamples), 1:end)))
        