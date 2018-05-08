function perturb(t) 
% Adding correlated noise (Ornstein-Uhlenbeck process) to load P & Q
global struc Fl Shunt DAE stop_TD_loop

%*% Don't let the load go less than 25%
if struc.rnd_wv(struc.ii) < -0.75
    struc.rnd_wv(struc.ii) = -0.75;
end

%*% Fast Noise: Random Load Fluctuations %*%
        Fl.con(3,2) = struc.P0 + struc.rnd_wv(struc.ii) * struc.P0 + struc.P0_noise.*(struc.VP_nos(struc.ii));
        Fl.con(3,5) = struc.Q0 + struc.rnd_wv(struc.ii) * struc.Q0 + struc.Q0_noise.*(struc.VQ_nos(struc.ii));
        
%*% After taking a random walk, record the new loads %*%
        struc.P0_noise = struc.P0 + struc.rnd_wv(struc.ii) * struc.P0;
        struc.Q0_noise = struc.Q0 + struc.rnd_wv(struc.ii) * struc.Q0;

% Record the loads
struc.P(struc.ii) = Fl.con(3,2);
struc.Q(struc.ii) = Fl.con(3,5);
struc.ii          = struc.ii + 1;

% Save Voltage Magnitude at Buses 3 (SVC) and 4 (Load) for struc.T secs
struc.V3B(struc.bus_ind) = DAE.y(struc.B3_Vind);
struc.V4B(struc.bus_ind) = DAE.y(struc.B4_Vind);
struc.bus_ind            = struc.bus_ind + 1;

% Bifurcation Test: if any voltage is above 1.35 or below 0.75, then
% bifurcation has probably occured

    if (DAE.y(struc.B3_Vind) > 1.35) || (DAE.y(struc.B3_Vind) < 0.75) || (DAE.y(struc.B4_Vind) > 1.35) || (DAE.y(struc.B4_Vind) < 0.75)
        stop_TD_loop = 1; % Kill the simulation!
    end
 
if (rem(t,struc.Tw) < (struc.intstep)) && (t > 0.1) && (struc.bus_ind > 200)
    % If this condition is met, process data and take control action!
    struc.V3B(struc.V3B==0) = []; % Remove 0's if necessary
    struc.V4B(struc.V4B==0) = []; % Remove 0's if necessary
 
    % Downsample to ~30 Hz
    struc.V3B = struc.V3B(1:3:end);
    struc.V4B = struc.V4B(1:3:end);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%*% Controller Tests %*%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%*% RBC %*%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if struc.Test == 1     % RBC
        disp(t)
        % Find the mean of the buffered data at bus 3
        B3_mean    = mean(struc.V3B);
        Delta_B3   = struc.Vr_SVC - B3_mean; % Subtract mean from reference
        
        % This is the control
        Delta_Bsvc = Delta_B3*struc.Kr;
        
        % Apply control to Shunt
        Shunt.con(6) = Shunt.con(6) + Delta_Bsvc;
   

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%*% MBC %*%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    elseif struc.Test == 2 % MBC
        disp(t)
        % Find the mean of the buffered data at bus 3
        B3_mean    = mean(struc.V3B);
        % Find the mean of the buffered data at bus 4
        B4_mean    = mean(struc.V4B);
        
        Delta_B3   = struc.Vr_SVC - B3_mean; % Subtract mean from reference
        Delta_B4   = struc.V_crit - B4_mean; % Subtract mean from reference
        if Delta_B4 < 0
            Delta_B4 = 0;
        end
        
        % This is the control
        Delta_Bsvc = Delta_B3*struc.Kr+Delta_B4*struc.Km;
        
        % Apply control to Shunt
        Shunt.con(6) = Shunt.con(6) + Delta_Bsvc;
        
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%*% VBC %*%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    else
        disp(t)
        % Find the mean of the buffered data at bus 3
        B3_mean    = mean(struc.V3B);
        % Find the mean of the buffered data at bus 4
        B4_mean    = mean(struc.V4B);
        
        Delta_B3   = struc.Vr_SVC - B3_mean; % Subtract mean from reference
        Delta_B4   = struc.V_crit - B4_mean; % Subtract mean from reference
        if Delta_B4 < 0
            Delta_B4 = 0;
        end
        
        % Now filter and take the variance of BV4
        SGF_4      = sgolayfilt(struc.V4B,2,35);
        var_v4     = var(SGF_4 - struc.V4B);
        Delta_var  = var_v4 - struc.var_crit;
        if Delta_var < 0
            Delta_var = 0;
        end
        
        % This is the control
        Delta_Bsvc = Delta_var*struc.Kv + Delta_B3*struc.Kr + Delta_B4*struc.Km;
        
        % Apply control to Shunt
        Shunt.con(6) = Shunt.con(6) + Delta_Bsvc;
        
        % Record the Variance
        struc.vars(struc.jj) = var_v4;
        struc.jj             = struc.jj + 1;
    end

    % Reset
    struc.bus_ind = 1;
    struc.V3B     = zeros(length(0:struc.intstep:struc.Tw)+1,1);  % B3 VMag vector
    struc.V4B     = zeros(length(0:struc.intstep:struc.Tw)+1,1);  % B4 VMag vector
end
end