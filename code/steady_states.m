% Solve for steady states by setting all time derivatives to zero.
p = zeros(1,13);

p(1)  = 2;    % N_1      : exponente Hill (TGFbeta_RA)
p(2)  = 1.0;  % TGFbeta  : concentracion de ligando (normalizada)
p(3)  = 0.5;  % K_2      : constante de Hill (rSMAD_I)
p(4)  = 0.5;  % K_1      : constante de Hill (TGFbeta_RA)
p(5)  = 0.10; % delta_1  : degradacion de TGFbeta_RA por SMURF/SMAD7
p(6)  = 0.30; % gamma_1  : inactivacion de rSMAD_A
p(7)  = 0.05; % delta_2  : degradacion de rSMAD_A por SMURF/SMAD7
p(8)  = 0.50; % alpha_1  : produccion de SMAD7
p(9)  = 0.10; % Omega_1  : degradacion de SMAD7
p(10) = 0.20; % Omega_2  : disociacion del complejo SMURF/SMAD7
p(11) = 0.30; % alpha_2  : formacion del complejo SMURF/SMAD7
p(12) = 0.40; % alpha_3  : produccion de SMURF
p(13) = 0.15; % Omega_3  : degradacion de SMURF


N_A             = p(1);   % N_1      : exponente Hill (TGFbeta_RA)
TGFbeta         = p(2);   % TGFbeta  : concentracion de ligando (constante)
k2              = p(3);   % K_2      : constante de Hill (rSMAD_I)
K_A             = p(4);   % K_1      : constante de Hill (TGFbeta_RA)
delta_ss7       = p(5);   % delta_1  : degradacion de TGFbeta_RA mediada por SMURF/SMAD7
gamma_rssa      = p(6);   % gamma_1  : inactivacion de rSMAD_A
delta_rsma      = p(7);   % delta_2  : degradacion de rSMAD_A mediada por SMURF/SMAD7
alpha_s7        = p(8);   % alpha_1  : produccion de SMAD7
omega_s7        = p(9);   % Omega_1  : degradacion de SMAD7
omega_smurfsmad = p(10);  % Omega_2  : tasa de disociacion del complejo SMURF/SMAD7
alpha_smurfsmad = p(11);  % alpha_2  : tasa de formacion del complejo SMURF/SMAD7
alpha_sf        = p(12);  % alpha_3  : produccion de SMURF
omega_SF        = p(13);  % Omega_3  : degradacion de SMURF

% Declare the symbolic state variables in the same order as the ODE model.
syms TGFbeta_RA rSMAD_A rSMAD_I SMAD7 SMURF_SMAD7 SMURF 



hill = @(x,K,N) (x.^N)./(K.^N + x.^N);



% Algebraic form of the steady-state equations.


% dy1: activacion del receptor por TGFbeta (Hill sobre TGFbeta_RA),
%      degradacion mediada por el complejo SMURF/SMAD7.
dy1 = hill(TGFbeta_RA, K_A, N_A) * TGFbeta ...
	- delta_ss7 * SMURF_SMAD7 * TGFbeta_RA;

% dy2: activacion de rSMAD inhibida por SMAD7,
%      inactivacion espontanea y degradacion por SMURF/SMAD7.
dy2 = (1/(1 + SMAD7)) * (TGFbeta_RA * hill(rSMAD_I, k2, 1)) ...
	- gamma_rssa * rSMAD_A ...
	- delta_rsma * rSMAD_A * SMURF_SMAD7;

% dy3: regeneracion de rSMAD_I por inactivacion de rSMAD_A,
%      consumo directo por activacion (sin factor inhibitorio).
dy3 = gamma_rssa * rSMAD_A ...
	- TGFbeta_RA * hill(rSMAD_I, k2, 1);

% dy4: produccion de SMAD7 por rSMAD_A, degradacion espontanea,
%      consumo al formar el complejo, liberacion al disociarse.
dy4 = alpha_s7 * rSMAD_A ...
	- omega_s7 * SMAD7 ...
	- alpha_smurfsmad * SMAD7 * SMURF ...
	+ omega_smurfsmad * SMURF_SMAD7;

% dy5: formacion y disociacion del complejo SMURF/SMAD7.
dy5 = alpha_smurfsmad * SMAD7 * SMURF ...
	- omega_smurfsmad * SMURF_SMAD7;

% dy6: produccion de SMURF por rSMAD_A, degradacion espontanea,
%      consumo al formar el complejo, liberacion al disociarse.
%      Se usa omega_smurfsmad (Omega_2) consistentemente con dy4 y dy5.
dy6 = alpha_sf * rSMAD_A ...
	- omega_SF * SMURF ...
	- alpha_smurfsmad * SMURF * SMAD7 ...
	+ omega_smurfsmad * SMURF_SMAD7;


% Solve the nonlinear system for the symbolic variables.
eqs = [dy1 == 0, dy2 == 0, dy3 == 0, dy4 == 0, dy5 == 0, dy6 == 0];
variables = [TGFbeta_RA rSMAD_A rSMAD_I SMAD7 SMURF_SMAD7 SMURF ];

vpasolve(eqs,variables)



