% 3 Bus System: Gen, Load, SVC :)
global struc Fl Shunt

% Bus Num, V Base, V Guess, Phase Guess, Area, Region
Bus.con = [ ...
   1  18;    % Bus 2
   2  230;   % Bus 7
   3  230;   % Bus 8
   4  230];  % Bus 5 (load bus)

% Bus Num, P Rating, V Rating, V Mag, Ref Angle
SW.con = [ ...
   1  100  18  1.04  0];

% Bus Num, Pow Rating, V Rating, P, Q, Max V, Min V, Conv. to Imp.
PQ.con = [ ...
   2  100.0  230  0                    0;
   3  100.0  230  0                    0;
   4  100.0  230  struc.load_fac*0.9   struc.load_fac*0.3];

% Bus Num, P%, P-V coeff, P-F coeff, Q%, Q-V coeff, Q-F coeff, tau, conn.
Fl.con = [ ...
  2  100  0  0  100  0  0  0.0001  1;
  3  100  0  0  100  0  0  0.0001  1;
  4  100  0  0  100  0  0  0.0001  1];

% From Bus, To Bus, P Rating, V Rating, F Rating, Length, NA, R, X, B, NA,
% NA, I Max, P Max, S Max, conn status
Line.con = [ ... %                 R       X       B
1  2  100  18   60  0  0.0783      0       0.125   0      0  0  0  0  0  1;
2  3  100  230  60  0  0           0.0085  0.072   0.3    0  0  0  0  0  1;
2  4  100  230  60  0  0           0.0461  0.561   0.3    0  0  0  0  0  1;
3  4  100  230  60  0  0           0.0461  0.0561  0.3    0  0  0  0  0  1];

% *Controllable* Shunt Device
Shunt.con = [ ...
    3 100 230 60 0 .45 1];

% Gen
M       = 10;
D       = 1;
Syn.con = [ ...              8        9                 13       14
  1  100  18    60  4  0  0  0.08958  0.01198  0  6  0  0.08645  0.01969  0  0.535  0  M  D  0  0  1  1  0.002  0  0];

% Exc
Exc.con = [ ... 
  1  2  5  -5  40  0.2  0.063  0.35  1  0.314  0.001  0.0039  1.555];

% TG
Tg.con = [ ...
   1   1   1   0.02   20    0.1     0.1   0.45   0.00   12.0   50.0];

% Bus names
Varname.bus = { ...
  'BUS-1 '; 'BUS-2 '; 'BUS-3 '; 'BUS-4 '};
