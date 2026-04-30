%% tgfbeta_equilibrio.m
%  Calculo de puntos de equilibrio del modelo TGFbeta/rSMAD/SMAD7/SMURF.
%  Llama directamente a tgfbeta_model.m con variables simbolicas para que
%  cualquier cambio en el modelo se refleje automaticamente aqui.

clear; clc;

% ── Parametros (mismos que tgfbeta_sim.m) ────────────────────────────────────
p(1)  = 2;    % N_1
p(2)  = 1.0;  % TGFbeta
p(3)  = 0.5;  % K_2
p(4)  = 0.5;  % K_1
p(5)  = 0.10; % delta_1
p(6)  = 0.30; % gamma_1
p(7)  = 0.05; % delta_2
p(8)  = 0.50; % alpha_1
p(9)  = 0.10; % Omega_1
p(10) = 0.20; % Omega_2
p(11) = 0.30; % alpha_2
p(12) = 0.40; % alpha_3
p(13) = 0.15; % Omega_3

% ── Variables simbolicas ──────────────────────────────────────────────────────
syms RA rA rI S7 SS SF real
y_sym = [RA; rA; rI; S7; SS; SF];

% ── Obtener las ecuaciones simbolicas desde tgfbeta_model ────────────────────
%  Al pasar y_sym (simbolico) y p (numerico), MATLAB evalua la funcion
%  propagando las simbolicas. El resultado es un vector de expresiones
%  simbolicas que representan dy/dt.
fprintf('Construyendo sistema simbolico desde tgfbeta_model.m...\n');
f_sym = tgfbeta_model(0, y_sym, p);   % t=0 no afecta (sistema autonomo)

% ── Punto inicial (semilla) ───────────────────────────────────────────────────
%  Usar el valor final de la simulacion como semilla para guiar vpasolve
%  hacia el atractor biologicamente relevante.
y_semilla = [0.9, 0.5, 0.5, 1.5, 0.4, 0.8];

fprintf('Buscando equilibrio con semilla: [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f]\n\n', ...
        y_semilla);

% ── vpasolve ──────────────────────────────────────────────────────────────────
sol = vpasolve(f_sym == 0, y_sym);

% ── Mostrar resultados ────────────────────────────────────────────────────────
nombres  = {'TGFbeta_RA', 'rSMAD_A', 'rSMAD_I', 'SMAD7', 'SMURF_SMAD7', 'SMURF'};
campos   = fieldnames(sol);
valores  = zeros(1, 6);

fprintf('%-20s %15s %15s\n', 'Variable', 'Valor', 'Biologicamente valido');
fprintf('%s\n', repmat('-', 1, 52));

todos_validos = true;
for i = 1:6
    v = double(sol.(campos{i}));
    valores(i) = v;
    valido = isreal(v) && v >= 0;
    if ~valido; todos_validos = false; end
    marca = 'Si';
    if ~valido; marca = '*** NO ***'; end
    fprintf('%-20s %15.6f %15s\n', nombres{i}, v, marca);
end
fprintf('\n');

if todos_validos
    fprintf('Equilibrio valido (todas las concentraciones >= 0 y reales).\n\n');
else
    fprintf('ADVERTENCIA: equilibrio con valores negativos o complejos.\n');
    fprintf('Prueba con otra semilla en y_semilla.\n\n');
end

% ── Verificacion: norma del residuo ──────────────────────────────────────────
residuo = double(subs(f_sym, y_sym, valores));
fprintf('Norma del residuo ||f(x*)||: %.2e\n\n', norm(residuo));

% ── Estabilidad: Jacobiano en x* ─────────────────────────────────────────────
fprintf('Calculando Jacobiano simbolico...\n');
J_sym = jacobian(f_sym, y_sym);
J_num = double(subs(J_sym, y_sym, valores));
vaps  = eig(J_num);

fprintf('\nValores propios del Jacobiano en el equilibrio:\n');
for i = 1:length(vaps)
    fprintf('  lambda_%d = %8.4f + %8.4fi\n', i, real(vaps(i)), imag(vaps(i)));
end

max_real = max(real(vaps));
fprintf('\nParte real maxima: %.6f --> ', max_real);
if max_real < 0
    fprintf('equilibrio ESTABLE.\n');
elseif max_real > 0
    fprintf('equilibrio INESTABLE.\n');
else
    fprintf('caso no concluyente.\n');
end