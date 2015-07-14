function P = S_Exp_Para(mode)

switch mode
  case 'paper-2AFC-corr'
    P.G.number_orientations=2;
    P.G.prior_task=[0 1]; % [cardinal, oblique]

    P.G.number_locations=1;
    P.G.dimension_X=10;
    P.G.dimension_G=4;

    P.G.kappa_O=[1 0]; % attended and unattended
    P.G.kappa_G=3;
    P.G.delta=15; % strength of X-G coupling for corr & CPs
    P.I.stimulus_regime='static';
    P.I.stimulus_contrast=zeros(1,P.G.number_orientations);

    P.S.number_repetitions=5; % number of trials

    P.S.number_repetitions=16; % number of trials
    P.I.n_zero_signal=20; % number of frames before onset of stimulus

    P.S.number_burn_in=0; % number of burn-in samPles
    P.S.number_samples_to_use=100; %Number  of non-burn samples to be used for evidence
    P.S.number_samples_per_evidence=5; % for dynamic-switching-signal-blocked
    P.S.phi_O=[[0:2:2*(P.G.number_orientations-1)];[1:2:2*P.G.number_orientations-1]]*pi/2/P.G.number_orientations;
    P.S.pO=ones(1,P.G.number_orientations)/P.G.number_orientations;
    P.S.pL=ones(1,P.G.number_locations)/P.G.number_locations;
    P.S.tauStyle=0.1; % initial value for style
    P.S.sigmaStyle=0.1;
    P.S.odds_inc=1/20; % with 80 samples/sec, every 250ms independent signal
    P.S.nT = 2; % Number of possible tasks

  case 'paper-2AFC-PK'
    P.G.number_orientations=2;
    P.G.prior_task=[0 1]; % [cardinal, oblique]
    P.G.number_locations=1;
    P.G.dimension_X=1024;
    P.G.dimension_G=256;
    P.G.kappa_O=[1 0]; % attended and unattended
    P.G.kappa_G=3;
    P.G.delta=80; % strength of X-G coupling
    P.I.stimulus_regime='dynamic-switching-signal-blocked';
    P.I.stimulus_contrast=zeros(1,P.G.number_orientations);
    P.S.number_repetitions=16; % number of trials; 
    P.I.n_zero_signal=50; % number of frames before onset of stimulus
    P.S.number_burn_in=0; % number of burn-in samPles
    P.S.number_samples_to_use=100; %Number  of non-burn samples to be used for evidence
    P.S.number_samples_per_evidence=2; % for dynamic-switching-signal-blocked
    P.S.phi_O=[[0:2:2*(P.G.number_orientations-1)];[1:2:2*P.G.number_orientations-1]]*pi/2/P.G.number_orientations;
    P.S.pO=ones(1,P.G.number_orientations)/P.G.number_orientations;
    P.S.pL=ones(1,P.G.number_locations)/P.G.number_locations;
    P.S.tauStyle=1;
    P.S.sigmaStyle=1e-10;
    P.S.odds_inc=1/20; % with 80 samples/sec, every 250ms independent signal
    P.S.nT = 2; % Number of possible tasks
    
  case 'test-2AFC-corr'
    P.G.number_orientations=2;
    P.G.prior_task=[0 1]; % [cardinal, oblique]
    P.G.number_locations=1;
    P.dimension_X=128;
    P.dimension_G=16;
    P.G.kappa_O=[1 0]; % attended and unattended
    P.G.kappa_G=3;
    P.G.delta=0.01; % strength of X-G coupling for corr & CPs
    P.I.stimulus_regime='static';
    P.I.stimulus_contrast=zeros(1,P.G.number_orientations);
    P.S.number_repetitions=128; % number of trials
    P.I.n_zero_signal=20; % number of frames before onset of stimulus
    P.S.number_burn_in=0; % number of burn-in samPles
    P.S.number_samples_to_use=100; %Number  of non-burn samples to be used for evidence
    P.S.number_samples_per_evidence=5; % for dynamic-switching-signal-blocked
    P.S.phi_O=[[0:2:2*(P.G.number_orientations-1)];[1:2:2*P.G.number_orientations-1]]*pi/2/P.G.number_orientations;
    P.S.pO=ones(1,P.G.number_orientations)/P.G.number_orientations;
    P.S.pL=ones(1,P.G.number_locations)/P.G.number_locations;
    P.S.tauStyle=0.1; % initial value for style
    P.S.sigmaStyle=0.1;
    P.S.odds_inc=1/20; % with 80 samples/sec, every 250ms independent signal
    P.S.nT = 2; % Number of possible tasks

  case 'test-contrast'
    P.G.number_orientations=2;
    P.G.prior_task=[1 0]; % [cardinal, oblique]
    P.G.number_locations=1;
    P.dimension_X=256;
    P.dimension_G=64;
    P.G.kappa_O=[1 0]; % attended and unattended
    P.G.kappa_G=3;
    P.G.delta=80/1024; % strength of X-G coupling in paper
    P.I.stimulus_regime='static';
    P.I.stimulus_contrast=zeros(1,P.G.number_orientations);
    P.S.number_repetitions=1; % number of trials
    P.I.n_zero_signal=20; % number of frames before onset of stimulus
    P.S.number_burn_in=0; % number of burn-in samPles
    P.S.number_samples_to_use=100; %Number  of non-burn samples to be used for evidence
    P.S.number_samples_per_evidence=5; % for dynamic-switching-signal-blocked
    P.S.phi_O=[[0:2:2*(P.G.number_orientations-1)];[1:2:2*P.G.number_orientations-1]]*pi/2/P.G.number_orientations;
    P.S.pO=ones(1,P.G.number_orientations)/P.G.number_orientations;
    P.S.pL=ones(1,P.G.number_locations)/P.G.number_locations;
    P.S.tauStyle=0.1;
    P.S.sigmaStyle=0.1;
    P.S.odds_inc=1/20; % with 80 samples/sec, every 250ms independent signal
    P.S.nT = 2; % Number of possible tasks
    
  otherwise
    warning('invalid option');
end