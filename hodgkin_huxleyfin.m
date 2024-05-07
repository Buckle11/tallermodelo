% Parámetros del modelo de Hodgkin-Huxley
Cm = 1;         % Capacitancia de la membrana (uF/cm^2)
gNa = 120;      % Conductancia máxima de sodio (mS/cm^2)
gK = 36;        % Conductancia máxima de potasio (mS/cm^2)
gL = 0.3;       % Conductancia máxima de fugas (mS/cm^2)
ENa = 115;       % Potencial de equilibrio del sodio (mV)
EK = -12;       % Potencial de equilibrio del potasio (mV)
EL = -10.6;     % Potencial de equilibrio de las fugas (mV)

% Paso de integración
dt = 0.01;      % Paso de tiempo (ms)
tmax = 50;      % Tiempo máximo de simulación (ms)
t = 0:dt:tmax;  % Vector de tiempos

% Inicialización de las variables de estado
V = zeros(size(t));  % Potencial de membrana (mV)
m = zeros(size(t));  % Variable de activación de sodio
h = zeros(size(t));  % Variable de inactivación de sodio
n = zeros(size(t));  % Variable de activación de potasio

% Condiciones iniciales
V(1) = -65;  % Potencial de membrana inicial (mV)
m(1) = 0.5;    % Inicialización de m
h(1) = 0.06;    % Inicialización de h
n(1) = 0.5;    % Inicialización de n

% Estímulo aplicado
segundo_inicio = input('Segundo en el que desea aplicar la corriente: ');
segundo_fin = input('Segundo en el que desea detener la corriente: ');
amperaje = input('¿Qué amperaje desea inducir (uA)?: ');

% Convertir los segundos a índices de tiempo
indice_inicio = find(t >= segundo_inicio, 1);
indice_fin = find(t <= segundo_fin, 1, 'last');

% Inicialización de la corriente
I = zeros(size(t));  % Inicialmente sin corriente

% Corregir la asignación de corriente estimulada
I(indice_inicio:indice_fin) = amperaje;  % Estímulo de amperaje entre segundo_inicio y segundo_fin

% Simulación del modelo de Hodgkin-Huxley
for i = 1:length(t)-1
    % Ecuaciones del modelo de Hodgkin-Huxley
    alpha_m = (0.1*(V(i)+40)) / (1 - exp(-(V(i)+40)/10));
    beta_m = 4 * exp(-(V(i)+65)/18);
    alpha_h = 0.07 * exp(-(V(i)+65)/20);
    beta_h = 1 / (1 + exp(-(V(i)+35)/10));
    alpha_n = (0.01*(V(i)+55)) / (1 - exp(-(V(i)+55)/10));
    beta_n = 0.125 * exp(-(V(i)+65)/80);
    
    % Actualización de las variables de estado
    m(i+1) = m(i) + dt * (alpha_m * (1 - m(i)) - beta_m * m(i));
    h(i+1) = h(i) + dt * (alpha_h * (1 - h(i)) - beta_h * h(i));
    n(i+1) = n(i) + dt * (alpha_n * (1 - n(i)) - beta_n * n(i));
    
    % Corriente iónica
    INa = gNa * m(i)^3 * h(i) * (V(i) - ENa);
    IK = gK * n(i)^4 * (V(i) - EK);
    IL = gL * (V(i) - EL);
    
    % Ecuación del modelo de la membrana
    V(i+1) = V(i) + dt * (1/Cm) * (I(i+1) - INa - IK - IL);
end

% Gráficos
figure;
subplot(2,1,1);
plot(t, V, 'b');
title('Potencial de membrana (mV)');
xlabel('Tiempo (ms)');
ylabel('Potencial (mV)');

subplot(2,1,2);
hold on;
plot(t, m, 'r', 'LineWidth', 1.5);
plot(t, h, 'g', 'LineWidth', 1.5);
plot(t, n, 'b', 'LineWidth', 1.5);
title('Variables de compuerta');
xlabel('Tiempo (ms)');
ylabel('Probabilidad');
legend('m', 'h', 'n');
hold off;

