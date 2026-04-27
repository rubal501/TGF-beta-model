clear; clc; close all;
 
% ── Parametros ────────────────────────────────────────────────────────────────
%       Simbolo   Descripcion
p(1)  = 2;    % N_1      : exponente Hill (TGFbeta_RA)
p(2)  = 1.0;  % TGFbeta  : concentracion de ligando (normalizada)
p(3)  = 2;  % K_2      : constante de Hill (rSMAD_I)
p(4)  = 1;  % K_1      : constante de Hill (TGFbeta_RA)
p(5)  = 0.10; % delta_1  : degradacion de TGFbeta_RA por SMURF/SMAD7
p(6)  = 0.30; % gamma_1  : inactivacion de rSMAD_A
p(7)  = 0.05; % delta_2  : degradacion de rSMAD_A por SMURF/SMAD7
p(8)  = 0.50; % alpha_1  : produccion de SMAD7
p(9)  = 0.10; % Omega_1  : degradacion de SMAD7
p(10) = 0.20; % Omega_2  : disociacion del complejo SMURF/SMAD7
p(11) = 0.30; % alpha_2  : formacion del complejo SMURF/SMAD7
p(12) = 0.40; % alpha_3  : produccion de SMURF
p(13) = 0.15; % Omega_3  : degradacion de SMURF
 
% ── Condiciones iniciales ─────────────────────────────────────────────────────
%  El receptor empieza inactivo; todo el rSMAD esta en forma inactiva;
%  hay un nivel basal pequeno de SMURF libre.
%         TGFb_RA  rSMAD_A  rSMAD_I  SMAD7  SMURF_SMAD7  SMURF
y0 = [    0.0,     0.0,     1.0,     0.0,   0.0,          0.1  ];
 
% ── Integracion ───────────────────────────────────────────────────────────────
tspan = [0, 300];   % 0 a 300 minutos (5 horas)
opts  = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'NonNegative', 1:6);
[t, Y] = ode45(@(t,y) tgfbeta_model(t, y, p), tspan, y0, opts);
 
% ── Graficas ──────────────────────────────────────────────────────────────────
nombres = {'TGF\beta_{RA}', 'rSMAD_A', 'rSMAD_I', ...
           'SMAD7', 'SMURF/SMAD7', 'SMURF'};
colores = {'#2471a3', '#27ae60', '#85c1e9', ...
           '#e74c3c', '#8e44ad', '#e67e22'};
 
figure('Position', [100 100 1000 600]);
 
for i = 1:6
    subplot(2, 3, i);
    plot(t, Y(:,i), 'Color', colores{i}, 'LineWidth', 2);
    xlabel('Tiempo (min)', 'FontSize', 10);
    ylabel('Concentracion (u.a.)', 'FontSize', 10);
    title(nombres{i}, 'Interpreter', 'tex', 'FontSize', 11);
    grid on;
    xlim([0, tspan(2)]);
    ylim([0, max(Y(:,i))*1.15 + 1e-6]);
end
 
sgtitle('Dinamica del modelo TGF\beta - SMAD7 - SMURF', ...
        'Interpreter', 'tex', 'FontSize', 13, 'FontWeight', 'bold');
 