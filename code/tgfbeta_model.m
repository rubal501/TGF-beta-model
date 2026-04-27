function dy = tgfbeta_model(t, y, p)

% Unpack model parameters from the input vector.
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

% State vector:
% y(1) = TGFbeta_RA   : receptor de TGFbeta activo
% y(2) = rSMAD_A      : rSMAD activo
% y(3) = rSMAD_I      : rSMAD inactivo
% y(4) = SMAD7        : SMAD7 libre
% y(5) = SMURF_SMAD7  : complejo SMURF/SMAD7
% y(6) = SMURF        : SMURF libre
TGFbeta_RA  = y(1);
rSMAD_A     = y(2);
rSMAD_I     = y(3);
SMAD7       = y(4);
SMURF_SMAD7 = y(5);
SMURF       = y(6);

% ── ODE system ────────────────────────────────────────────────────────────

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

dy = [dy1; dy2; dy3; dy4; dy5; dy6];

end


function H = hill(x, K, N)
% Funcion de activacion de Hill estandar.
H = (x^N) / (K^N + x^N);
end
