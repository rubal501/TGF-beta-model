function dy = tgfbeta_model(t,y,p)

    %#ok<INUSD>
    % Unpack model parameters from the input vector.
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

    % State vector:
    % y(1) = TGFbeta_RA
    % y(2) = rSMAD_A
    % y(3) = rSMAD_I
    % y(4) = SMAD7
    % y(5) = SMURF_SMAD7
    % y(6) = SMURF
    TGFbeta_RA = y(1);
    rSMAD_A = y(2);
    rSMAD_I = y(3);
    SMAD7 = y(4);
    SMURF_SMAD7 = y(5);
    SMURF = y(6);

    % Effective rSMAD activation is reduced by the inhibitory action of
    % SMAD7 and modulated by the inactive rSMAD pool.
    RegSMADA = (1/(1+SMAD7)) * (TGFbeta_RA * hill(rSMAD_I, k2, 1));


    % ODE system for receptor activation, SMAD cycling, and SMURF feedback.
    dy1 = hill(TGFbeta, K_A, N_A) - delta_ss7 * SMURF_SMAD7 * TGFbeta_RA;
    dy2 = RegSMADA - gamma_rssa * rSMAD_A - delta_rsma * rSMAD_A * SMURF;
    dy3 = gamma_rssa * rSMAD_A - RegSMADA ;
    dy4 = alpha_s7 * rSMAD_A - omega_s7 * SMAD7  ...
    - alpha_smurfsmad * SMAD7 * SMURF + omega_smurfsmad * SMURF_SMAD7;
    dy5 = alpha_smurfsmad * SMAD7 * SMURF - omega_smurfsmad * SMURF_SMAD7;
    dy6 = alpha_sf * rSMAD_A - omega_SF * SMURF  ...
    - alpha_smurfsmad * SMURF * SMAD7 + omega_S7SF * SMURF_SMAD7;


     
    dy = [ dy1 ; dy2 ; dy3 ; dy4 ; dy5; dy6 ];

end


function H = hill(x,K,N)
    % Standard Hill activation function.
    H = (x^N)/(K^N + x^N);
end
