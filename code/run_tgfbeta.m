% Simulate the TGF-beta signaling model over a fixed time horizon.
tspan = [0 100];

% Parameter vector shared with tgfbeta_model.m.
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

% Initial conditions for:
% [TGFbeta_RA, rSMAD_A, rSMAD_I, SMAD7, SMURF_SMAD7, SMURF].
y0 = [0.01, 0, 1, 0.01, 0.0, 0];

% Pass the parameter vector through an anonymous function so ode45 receives
% the signature f(t, y) that it expects.
[t, y] = ode45(@(t, y) tgfbeta_model(t, y, p), tspan, y0);

% Plot the temporal trajectories of all state variables.
figure;
hold on;
plot(t, y(:, 1))
plot(t, y(:, 2))
plot(t, y(:, 3))
plot(t, y(:, 4))
plot(t, y(:, 5))
plot(t, y(:, 6))
legend({"TGFbeta_RA", "RSMAD_A", "RSMAD_I", "SMAD7", "SMURF/SMAD7", "SMURF"})
hold off;
