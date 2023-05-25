clear all;
close all;
clc;

%% Modelo de Rosenberg

Fs = 8000; %Frecuencia de muestreo modelo

Ts = 1/Fs; %Tiempo de cada muestra

a = 1; %Amplitud


% Tiempo apertura

T_milis = 100; % tiempo en milisegundos
T_muestras = T_milis*Fs/1000 + 1; % tiempo en samples
Tp = 0.40 * (T_muestras-1); % porcentaje del tiempo de apertura con pendiente positiva(0.40 sugerido)
Tn = 0.16 * (T_muestras-1); % porcentaje del tiempo de apertura con pendiente negativa(0.16 sugerido)


% Vector tiempo

for j = 1:1:T_muestras
    y(j) = (j-1) * Ts;
end


% Vectores de ceros

f_A = zeros(1 , T_muestras);
f_B = zeros(1 , T_muestras);
f_C = zeros(1 , T_muestras);
f_D = zeros(1 , T_muestras);
f_E = zeros(1 , T_muestras);
f_F = zeros(1 , T_muestras);


% Modelo matematico para 0 hasta Tp

for t = 1:1:Tp + 1
    f_A(1,t) = a*((t-1)/Tp); 
    f_B(1,t) = a*(3*((t-1)/Tp)^2 - 2*((t-1)/Tp)^3);
    f_C(1,t) = (a/2)*(1-cos(pi*(t-1)/Tp));
    f_D(1,t) = (a/2)*(1-cos(pi*((t-1)/Tp)));
    f_E(1,t) = a*(sin((pi/2)*(t-1)/Tp));
    f_F(1,t) = (3*a/2)*(t/Tp);
    if f_F(1,t) > 1
        f_F(1,t) = 1;
    else
        f_F(1,t) = f_F(1,t);
    end % falta que sea <=
end


% Modelo matematico para Tp+1 hasta T = Tp+Tn

for t = Tp + 2:1:Tp + Tn + 1
    f_A(1,t) = a*(1-(((t-1)-Tp)/Tn));
    f_B(1,t) = a*(1-(((t-1)-Tp)/Tn)^2);
    f_C(1,t) = a*(cos((pi/2)*((t-1)-Tp)/Tn));
    f_D(1,t) = (a/2)*(1+cos(pi*((t-1)-Tp)/Tn)); % Rosenberg propone f_D(1,t) = (a/2)*(1+cos((pi/2)*((t-1)-Tp)/Tn)) , sin embargo esta mal
                                                % En (V. Espinoza, 2012)
                                                % esta corregido
    f_E(1,t) = a*(cos((pi/2)*((t-1)-Tp)/Tn));
    f_F(1,t) = (3*a/2)*(1-(((t-1)-Tp)/Tn));
    if f_F(1,t) > 1
        f_F(1,t) = 1;
    else
        f_F(1,t) = f_F(1,t);
    end
end


%% Graficos

figure
stem(y,f_D)
box on
grid on

% Grafico todos los pulsos del modelo de Rosenberg

figure('Name','All Glottal Pulse Shapes','NumberTitle','off');
plot(y,f_A)
title('All Glottal Pulse Shapes')
hold on
plot(y,f_B,'r')
plot(y,f_C,'c')
plot(y,f_D,'m-.')
plot(y,f_E,'g')
plot(y,f_F,'k')
legend('f_A','f_B','f_C','f_D','f_E','f_F')
hold off

% Grafico solo del pulso f_D, que es el pulso a estudiar

figure('Name','Rosenberg f_D Pulse Shape','NumberTitle','off');
plot(y,f_D,'m')
title('Glottal Pulse "f_D"')

% Grafico de la primera derivada de f_D
% representa la tasa de cambio de flujo de aire o bien, aceleración de
% flujo de aire

figure
plot(y(1:end-1),diff(f_D))
title ('Aceleración de Flujo de Aire')

%% Referencias

% 1.- Effect of Glottal Pulse Shape on the Quality of Natural Vowels.
%     Rosenberg. 1971.
% 2.- Stationary and dynamic aerodynamic assessment of vocal hyperfunction
%     using enhanced supraglottal and subglottal inverse ﬁltering methods.
%     V. Espinoza. 2018.
