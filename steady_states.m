p = [2,    ... % N_A
         1,    ... % TGFbeta
         0.5,  ... % k2
         1,    ... % K_A
         0.1,  ... % delta_ss7
         0.2,  ... % gamma_rssa
         0.05, ... % delta_rsma
         0.3,  ... % alpha_s7
         0.1,  ... % omega_s7
         0.05, ... % omega_smurfsmad
         0.02, ... % alpha_smurfsmad
         0.1,  ... % omega_s7sf
         0.05, ... % alpha_sf
         0.1,  ... % omega_SF
         0.01, ... % oemga_sf
         0.05];    % omega_S7SF


N_A               = p(1);
TGFbeta           = p(2);
k2                = p(3);
K_A               = p(4);
delta_ss7         = p(5);
gamma_rssa        = p(6);
delta_rsma        = p(7);
alpha_s7          = p(8);
omega_s7          = p(9);
omega_smurfsmad   = p(10);
alpha_smurfsmad   = p(11);
omega_s7sf        = p(12);
alpha_sf          = p(13);
omega_SF          = p(14);
oemga_sf          = p(15);
omega_S7SF        = p(16);




syms TGFbeta_RA rSMAD_A rSMAD_I SMAD7 SMURF_SMAD7 SMURF 

hill = @(x,K,N) (x.^N)./(K.^N + x.^N);
RegSMADA = (1/(1+SMAD7))* (TGFbeta_RA * hill(rSMAD_I, k2, 1) );


% d TGFB-Ra
dy1 = hill(TGFbeta,K_A,N_A)*TFGB_RA - delta_ss7 * SMURF_SMAD7 * TGFbeta_RA;
% d RSMAD_A
dy2 = RegSMADA - gamma_rssa * rSMAD_A - delta_rsma * rSMAD_A * SMURF;
% d RSMAD_I
dy3 = gamma_rssa * rSMAD_A - RegSMADA ;
% d SMAD7
dy4 = alpha_s7 * rSMAD_A - omega_s7 * SMAD7  ...
- alpha_smurfsmad * SMAD7 * SMURF + omega_smurfsmad * SMURF_SMAD7;
% d SMURF_SMAD7
dy5 = alpha_smurfsmad * SMAD7 * SMURF - omega_smurfsmad * SMURF_SMAD7;
% d SMURF 
dy6 = alpha_sf * rSMAD_A - omega_SF * SMURF  ...
- alpha_smurfsmad * SMURF * SMAD7 + omega_S7SF * SMURF_SMAD7;



% Define the system of equations
eqs = [dy1 == 0, dy2 == 0, dy3 == 0, dy4 == 0, dy5 == 0, dy6 == 0];
variables = [TGFbeta_RA rSMAD_A rSMAD_I SMAD7 SMURF_SMAD7 SMURF ];

vpasolve(eqs,variables)




