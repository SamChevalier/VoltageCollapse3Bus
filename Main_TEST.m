%% Numerical Simulation For 4-bus System with Random Load Fluctuations
%  Gen --- XFM ---|
%           |    SVC ---|
%           |--------- Load
clear all; clc;
global struc stop_TD_loop

%% 0. Load the Fast and Slow Noise Data
load('Noise_Data')
n_sims = size(Slw_Ns_Mtrx,1); % Use the amount of load data to determine

%% 1. Which test do you wish to run? How many n_sims?
%*% RBC = 1
%*% MBC = 2
%*% VBC = 3
struc.intstep  = 0.01; % Integration Step Size
struc.var_crit = 0;    % Define Below
tbegin         = 0;
tfinal         = 800;

%% 2. Enter main loop

% Loop over 3 Tests: RBC, MBC, VBC
for hh = 1:3
    struc.Test = hh;
    
    % Params
    stop_TD_loop   = 0;                 % Reset
    struc.load_fac = 2.05;              % Initial Loading: 1.123 doesn't converge, 1.1229 does!
    struc.V_crit   = 0.985;             % Critical Voltage
    dyn_4;
    
    %*% Noise Stats
    struc.tcorr   = 1;                  % Correlation time of noise
    struc.D       = 0.0005;             % Diffusion Coefficient
    
    struc.VP_nos = Fst_Ns_Mtrx(1,:);
    struc.VQ_nos = Fst_Ns_Mtrx(1,:);
    struc.rnd_wv = Slw_Ns_Mtrx(1,:);
    
    %*% Simulation Parameters
    struc.ii       = 1;    % Index for Recording Load/Voltage Data
    struc.jj       = 1;    % Index for Recording Variances
    struc.bus_ind  = 1;
    struc.Kr       = 0.1;   % Reference
    struc.Km       = 1;     % Mean
    struc.Kv       = 1000;  % Variance
    
    struc.P        = zeros(length(tbegin:struc.intstep:tfinal),1);
    struc.Q        = zeros(length(tbegin:struc.intstep:tfinal),1);
    struc.Tw       = 3;      % Data Collection Window (s)
    struc.V3B      = zeros(length(0:struc.intstep:struc.Tw)+1,1);  % B3 VMag vector
    struc.V4B      = zeros(length(0:struc.intstep:struc.Tw)+1,1);  % B4 VMag vector
    
    struc.P0       = PQ.con(3,4);       % Initial value of load active power
    struc.Q0       = PQ.con(3,5);       % Initial value of load reactive power
    struc.P0_noise = struc.P0;          % This is constantly updated to new "Steady State" power (P)
    struc.Q0_noise = struc.Q0;          % This is constantly updated to new "Steady State" power (Q)
    
    initpsat;                           % Initialize PSAT global variables
    datafile = 'dyn_4';                 % Test case data file
    runpsat(datafile,'data');           % Initialize datafile
    runpsat('perturb','pert');          % "perturb": Perturbation file
    Settings.freq   = 60;               % Change System Freq from default to 60
    clpsat.readfile = 1;                % Read data file before running power flow
    
    % Power Flow
    runpsat('pf');
    
    if struc.var_crit == 0
        % Use HELM/FPP/Statistical Solver to get Critical Var on Bus 4
        s_max = 0.1229;    % HELM isn't actually used here since the calc is trivial :-)
        PC    = 0.0001;    % Probability of collapse (1 = 100%)
        D     = struc.D;   % Diffusion coeff (s = 1 means 100% load pocket increase)
        t_min = 1;         % Time in minutes to protect against collapse
        [Crit_Var,s_margin] = Calc_CritVars(t_min,D,s_max,PC);
        struc.var_crit = Crit_Var;
        
        % Now re-run power flow
        runpsat(datafile,'data');
        struc.load_fac = 2.05;
        runpsat('pf');
    end
    
    % PF Results
    voltages        = DAE.y(Bus.v);
    struc.Vr_SVC    = voltages(3);          % Use power flow results at B3 as reference signal
    struc.B3_Vind   = Bus.v(3);             % Bus 3 index
    struc.B4_Vind   = Bus.v(4);             % Bus 4 index
    % runpsat('sssa');                  % Find bad states: runpsat('eigrep')
    % mx_reval = max(real(SSSA.eigs));  % Eigenvalues
    
    % SETTINGS FOR TIME DOMAIN SIMULATION
    Settings.coi   = 1;                % Use center of inertia for synchronous machines
    Settings.t0    = tbegin;           % Initial simulation time
    Settings.tf    = tfinal;           % Final simulation time
    Settings.pq2z  = 0;                % Do not convert PQ loads to constant impedance (Z) loads after power flow
    Settings.fixt  = 1;                % Enable fixed time-step solver
    Settings.tstep = struc.intstep;    % Fixed time step value
    
    % Time Domain Setup
    Varname.idx  = [(DAE.n + Bus.n + 1) : (DAE.n + 2*Bus.n) ...
                    (DAE.n + DAE.m + 1) : (DAE.n + DAE.m + Bus.n)];
    
    % Define Output Variables
    struc.ix_Vm = DAE.n+Bus.n+1:DAE.n+2*Bus.n;
    struc.ix_Pi = DAE.n + DAE.m + 1:DAE.n + DAE.m + Bus.n;
    
    % ******** Time Domain Simulation Process ******** %
    runpsat('td'); % Run TD (dynamic) simulation
    close all
    
    % Output variable values
    Vm   = Varout.vars(:,1:4);   % Vm: Bus Voltage Magnitudes
    P4   = Varout.vars(:,8);     % Active Power Injections
    time = Varout.t;             % Simulation Times
    
    % Save time, voltage, and Power Data for each run
    if struc.Test == 1
        T_T1 = time;
        P_T1 = P4;
        V_T1 = Vm;
    elseif struc.Test == 2
        T_T2 = time;
        P_T2 = P4;
        V_T2 = Vm;
    else
        T_T3 = time;
        P_T3 = P4;
        V_T3 = Vm;
    end
end
