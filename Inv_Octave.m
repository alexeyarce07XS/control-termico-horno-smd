clc;
clear;
close all;

% CONTROL TERMICO DE UN HORNO SMD
% GNU OCTAVE - MODELO DE PRIMER ORDEN Y CONTROL PI

% 1. PARAMETROS
R = 0.80;          % Resistencia termica (°C/W)
C = 12000;         % Capacidad termica (J/°C)
P = 2000;          % Potencia maxima (W)
Ta = 25;           % Temperatura ambiente (°C)
Tref = 180;        % Temperatura deseada (°C)

Kp = 0.02;         % Ganancia proporcional
Ki = 0.00005;      % Ganancia integral
Qd = 0;            % Perturbacion termica constante (W)

% 2. FUNCION DE TRANSFERENCIA Y PARAMETROS
K = P * R;
tau = R * C;

fprintf('\n=== MODELO DEL HORNO SMD ===\n');
fprintf('Ganancia K = %.2f °C/u\n', K);
fprintf('Constante de tiempo = %.2f s (%.2f min)\n', tau, tau/60);
fprintf('G(s) = %.2f / (%.2f*s + 1)\n', K, tau);

% 3. FUNCION DE TRANSFERENCIA EN LAZO CERRADO CON PI
% T(s) = K(Kp*s + Ki) / (tau*s^2 + (1+K*Kp)*s + K*Ki)
numCL = [K*Kp, K*Ki];
denCL = [tau, 1 + K*Kp, K*Ki];

polos = roots(denCL);
ceros = roots(numCL);

fprintf('\n=== LAZO CERRADO ===\n');
fprintf('Numerador: ');
disp(numCL);
fprintf('Denominador: ');
disp(denCL);

fprintf('Polos:\n');
disp(polos);
fprintf('Ceros:\n');
disp(ceros);

if all(real(polos) < 0)
    fprintf('Estabilidad: ESTABLE (todos los polos tienen parte real negativa).\n');
else
    fprintf('Estabilidad: INESTABLE o marginal.\n');
end

% 4. SIMULACION NUMERICA (EULER)
dt = 2;                    % Paso de simulacion (s)
tiempo_final = 4*3600;     % Cuatro horas
t = 0:dt:tiempo_final;
n = length(t);

T_abierto = zeros(1,n);
T_cerrado = zeros(1,n);
u_abierto = zeros(1,n);
u_cerrado = zeros(1,n);

T_abierto(1) = Ta;
T_cerrado(1) = Ta;

% Perturbacion constante; usar Qd=0 para el caso nominal
Q = Qd*ones(1,n);

% Mando fijo de lazo abierto calculado para la referencia nominal
u_fijo = (Tref - Ta)/K;
u_fijo = max(0, min(1, u_fijo));

% Estado integral del controlador PI
I = 0;

for i = 1:n-1

    % LAZO ABIERTO: potencia fija, sin corregir el error
    u_abierto(i) = u_fijo;
    dT_abierto = (P*u_abierto(i) ...
        - (T_abierto(i)-Ta)/R + Q(i))/C;
    T_abierto(i+1) = T_abierto(i) + dt*dT_abierto;

    % LAZO CERRADO: controlador PI con anti-windup
    error = Tref - T_cerrado(i);
    u_previo = Kp*error + I;
    u_cerrado(i) = max(0, min(1, u_previo));

    % Integracion condicional para evitar windup
    if (u_previo >= 0 && u_previo <= 1) ...
        || (u_previo > 1 && error < 0) ...
        || (u_previo < 0 && error > 0)
        I = I + Ki*error*dt;
    end

    u_cerrado(i) = max(0, min(1, Kp*error + I));

    dT_cerrado = (P*u_cerrado(i) ...
        - (T_cerrado(i)-Ta)/R + Q(i))/C;
    T_cerrado(i+1) = T_cerrado(i) + dt*dT_cerrado;
end

u_abierto(n) = u_abierto(n-1);
u_cerrado(n) = u_cerrado(n-1);

% 5. GRAFICA DE TEMPERATURA
figure(1);
plot(t/60, T_abierto, 'r--', 'linewidth', 2);
hold on;
plot(t/60, T_cerrado, 'b', 'linewidth', 2);
plot(t/60, Tref*ones(1,n), 'k:', 'linewidth', 1.5);
grid on;
xlabel('Tiempo (min)');
ylabel('Temperatura (°C)');
title('Control termico del horno SMD');
legend('Lazo abierto', 'Lazo cerrado PI', 'Referencia');

% 6. GRAFICA DE LA SENAL DE CONTROL
figure(2);
plot(t/60, 100*u_abierto, 'r--', 'linewidth', 2);
hold on;
plot(t/60, 100*u_cerrado, 'b', 'linewidth', 2);
grid on;
xlabel('Tiempo (min)');
ylabel('Potencia del calentador (%)');
title('Senal de control');
legend('Lazo abierto', 'Lazo cerrado PI');

% 7. RESULTADOS FINALES
fprintf('\n=== RESULTADOS DE LA SIMULACION ===\n');
fprintf('Mando fijo de lazo abierto = %.2f %%\n', 100*u_fijo);
fprintf('Temperatura final en lazo abierto = %.2f °C\n', T_abierto(n));
fprintf('Temperatura final en lazo cerrado = %.2f °C\n', T_cerrado(n));
fprintf('Error final = %.2f °C\n', Tref-T_cerrado(n));
fprintf('\nSimulacion finalizada.\n');
