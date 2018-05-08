function [Crit_Var,s_margin] = Calc_CritVars(t_min,D,s_max,PC)
global struc Fl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will compute expected critical variance (bus 4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define a vector of loading parameters between 0 and s_max: the loading
% parameter on the edge of the stability will be in this vector.
s_vec     = 0:0.0001:s_max;

% Using this vector, define a new vector of collapse probabilites. We are
% pretty much sliding the mean until we find the correct value.
collapse_prob_vec = 1-erf(s_vec/sqrt(4*D*t_min));

% Find the value of "s" which corresponds to PC: the desired probability.
% In some situations, the probability PC cannot be achieved. If this is so,
% then s_margin
Diff_vec = collapse_prob_vec - PC;
[~,ind]  = (min(abs(Diff_vec)));
s_margin = s_vec(ind);

% We must subtract s_margin from the collapse distance in order to solve
% for the actual margin.
s_margin = s_max - s_margin;

% We now scale Load Factor by this value and compute the critical variances

%  Gen --- XFM ---|
%           |    SVC ---|
%           |--------- Load
    struc.load_fac = 2.05*(1 + s_margin);
    dyn_4;
    
    %*% Noise Stats
    gamma         = 1/struc.tcorr;
    nos_std       = 0.0025*sqrt(2*gamma); % Should be 0.01
    npq           = size(PQ.con,1);
    
    P0             = PQ.con(:,4);
    Q0             = PQ.con(:,5);
    
    initpsat;                           % Initialize PSAT global variables
    datafile = 'dyn_4';                 % Test case data file
    runpsat(datafile,'data');           % Initialize datafile
    runpsat('perturb','pert');          % "perturb": Perturbation file
    Settings.freq   = 60;               % Change System Freq from default to 60
    clpsat.readfile = 1;                % Read data file before running power flow
    
    % Power Flow
    runpsat('pf');
    
    % ******** Develop Alg and State Variable Variance Matrices ******** %
    n_fl = size(Fl.con,1); % How large is u?
    
    % Define the indices of the generator deltas in the state variable vec
    del_ix = [Syn.delta];
    
    % Assemble A - Can call fm_abcd; also
    As = full(DAE.Fx) - full(DAE.Fy)*inv(full(DAE.Gy))*full(DAE.Gx);
    
    % In the state matrix, set the deltas equal to the negative of delta_1
    % This will effectively set delta1 = 0
    As(del_ix(2:end),2) = -As(1,2);
    
    % Remove the first row and column
    Asr = As(2:end,2:end);
    
    % How many states (x) and algebraic (y) variables?
    ny = size(DAE.Gy,1);
    nx = size(Asr,1);
    
    % Build dG_du
    dG_du = zeros(ny,n_fl); 
    for jj = 1:3
        dG_du(Fl.con(jj,1),jj)   = P0(jj);   % Theta is realted to P
        dG_du(Fl.con(jj,1)+4,jj) = Q0(jj);   % V is realted to Q
        % The other differential variables are not important to us
    end
    
    % Build A2
    A2 = - (DAE.Fy(2:end,:)/DAE.Gy) * dG_du;
    
    % Build A3
    A3 = zeros(n_fl,nx);
    
    % Build E (A4)
    gamma = 1/struc.tcorr;
    E     = gamma*eye(n_fl);
    
    % Assemble the final A matrix
    A = [Asr A2;
        A3  -E];
    
    % Assemble B
    nos_std       = 0.0025*sqrt(2*gamma); % Should be 0.01
    nos_vec = nos_std*ones(3,1);
    B_Sub   = diag(nos_vec);
    B       = [zeros(nx,npq); B_Sub];
    
    % Create covariance matrix for the differential variables (sigma_x)
    Asp   = -sparse(A);
    Bsp   = sparse(B);
    sigma = lyap(Asp,-Bsp*Bsp.');
    
    % Convert to algebraic Covaraince Matrix
    K           = [-DAE.Gy\DAE.Gx(:,2:end)  -DAE.Gy\dG_du];
    sigma_y     = K * sigma * K.';
    
    % Grab the voltage variances (diagonal elements) from the Cov Mat
    sigma_VA    = sigma_y(1:8,1:8);
    V_var_sigma = diag(sigma_VA(5:8,5:8));
    Crit_Var    = V_var_sigma(4);
end

