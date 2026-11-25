%% =====================================================================
%  AXIAL COMPRESSOR STAGE-BY-STAGE MEAN-LINE DESIGN
%  Method: Cohen, Rogers & Saravanamuttoo ("Gas Turbine Theory")
%  Features: Exact RPM kinematics, dT0 Bisection Solver, First-Power Law
% =====================================================================
clear; clc; close all;

%% ---------------------------------------------------------------------
%  1. USER INPUT BLOCK
% ------------------------------------------------------------------------
gam   = 1.4;                 % ratio of specific heats
R     = 287;                 % J/kg.K, specific gas constant (air)
Cp    = gam*R/(gam-1);       % J/kg.K

mdot     = 131.2447;         % kg/s, mass flow rate at inlet
N_rpm    = 5819;             % rpm, LP spool speed (EXACT INPUT)
T01      = 389.86;           % K, inlet stagnation temperature
P01      = 73.66e3;          % Pa, inlet stagnation pressure
eta_p    = 0.8782;           % polytropic efficiency
nu_ht    = 0.27;             % hub-to-tip radius ratio at inlet
M1       = 0.549;            % inlet absolute Mach no.
PR_total = 4.1;              % overall design pressure ratio

N_stages = 7;                

% Prescribed per-stage design parameters
dT0_1 = 20.0;                % K, Stage 1 prescribed temp rise
dT0_2 = 25.0;                % K, Stage 2 prescribed temp rise
lambda = [0.98, 0.93, 0.88, 0.83, 0.83, 0.83, 0.83];   % work-done factors
Lambda = [ NaN, 0.70, 0.50, 0.50, 0.50, 0.50, 0.50 ];  % degree of reaction

VORTEX_LAW = 'first_power';  % 'first_power' -> Cw = a*r (solid body)
DEHALLER_LIMIT = 0.72;       % stall guideline for plotting

fprintf('=== INPUT SUMMARY ===\n');
fprintf('Stages: %d | mdot=%.4f kg/s | N = %d rpm | PR_total=%.3f\n\n', ...
    N_stages, mdot, N_rpm, PR_total);

%% ---------------------------------------------------------------------
%  2. STEP 1 : INLET KINEMATICS & ANNULUS (DRIVEN BY RPM)
% ------------------------------------------------------------------------
T1   = T01 / (1 + (gam-1)/2 * M1^2);      % static temperature, K
P1   = P01 / (T01/T1)^(gam/(gam-1));      % static pressure, Pa
rho1 = P1 / (R*T1);                       % static density, kg/m3
Ca   = M1 * sqrt(gam*R*T1);               % mean axial velocity (constant for mean-line)

% Annulus sizing based on mass flow
A1  = mdot / (rho1 * Ca);                 % Inlet area, m^2
rt1 = sqrt( A1 / (pi*(1 - nu_ht^2)) );    % tip radius at inlet
rr1 = nu_ht * rt1;                        % root radius at inlet
rm1 = (rr1 + rt1)/2;                      % mean radius at inlet

omega = N_rpm * (2*pi) / 60;              % rad/s
U_tip1 = omega * rt1;                     
U = omega * rm1;                          % U_mean, held constant for every stage
phi = Ca / U;                             % Actual computed flow coefficient

fprintf('=== STEP 1: INLET KINEMATICS ===\n');
fprintf('T1 = %.2f K | P1 = %.2f kPa | Ca (=C1) = %.2f m/s\n', T1, P1/1000, Ca);
fprintf('rt1 = %.4f m | rr1 = %.4f m | rm1 = %.4f m\n', rt1, rr1, rm1);
fprintf('U_tip1 = %.2f m/s | U_mean = %.2f m/s | phi = %.3f\n\n', U_tip1, U, phi);

%% ---------------------------------------------------------------------
%  3. BISECTION SOLVER FOR STAGE 3-7 TEMPERATURE DISTRIBUTION
% ------------------------------------------------------------------------
dT_low = 10; dT_high = 60; tol = 1e-4;
dT0 = zeros(1, N_stages);
dT0(1) = dT0_1; 
dT0(2) = dT0_2;
dT_guess = 0;

while (dT_high - dT_low) > tol
    dT_guess = (dT_low + dT_high) / 2;
    for i = 3:N_stages
        dT0(i) = dT_guess;
    end
    
    PR_calc = 1.0;
    T0_curr = T01;
    for i = 1:N_stages
        PR_calc = PR_calc * (1 + eta_p * dT0(i) / T0_curr)^(gam/(gam-1));
        T0_curr = T0_curr + dT0(i);
    end
    
    if PR_calc < PR_total
        dT_low = dT_guess;
    else
        dT_high = dT_guess;
    end
end
fprintf('=== BISECTION SOLVER ===\n');
fprintf('Balanced Stages 3-7 at dT0 = %.3f K to hit PR = %.3f\n\n', dT_guess, PR_total);

%% ---------------------------------------------------------------------
%  4. STEP 2 : STAGE-BY-STAGE VELOCITY TRIANGLES
% ------------------------------------------------------------------------
beta1 = zeros(1,N_stages); beta2 = zeros(1,N_stages);   
alpha1 = zeros(1,N_stages); alpha2 = zeros(1,N_stages);   
Cw1 = zeros(1,N_stages); Cw2 = zeros(1,N_stages);   
PR_st = zeros(1,N_stages);                                
P01_st = zeros(1,N_stages); T01_st = zeros(1,N_stages);   
P03_st = zeros(1,N_stages); T03_st = zeros(1,N_stages);   
Lambda_out = zeros(1,N_stages);                           
DeHaller = zeros(1,N_stages);
Deflection_deg = zeros(1,N_stages);

fprintf('=== STEP 2: STAGE-BY-STAGE DESIGN ===\n');
for i = 1:N_stages
    if i == 1
        T01_st(i) = T01;  P01_st(i) = P01;
    else
        T01_st(i) = T03_st(i-1);  P01_st(i) = P03_st(i-1);
    end

    dT0i = dT0(i);
    lam_i = lambda(i);

    if isnan(Lambda(i))
        % Stage 1: Zero inlet swirl
        alpha1(i) = 0;  Cw1(i) = 0;
        dCw = Cp*dT0i / (lam_i*U);
        Cw2(i)  = dCw;
        beta1(i)  = atan(U/Ca);
        beta2(i)  = atan( (U-Cw2(i))/Ca );
        alpha2(i) = atan( Cw2(i)/Ca );
        Lambda_out(i) = 1 - (Cw2(i)+Cw1(i))/(2*U);
    else
        % Stages 2-7: Prescribed Reaction
        Lam_i = Lambda(i);
        rhs_diff = dT0i*Cp/(lam_i*U*Ca);
        rhs_sum  = Lam_i*2*U/Ca;
        tanb1 = (rhs_sum + rhs_diff)/2;
        tanb2 = (rhs_sum - rhs_diff)/2;
        beta1(i) = atan(tanb1);  beta2(i) = atan(tanb2);

        tana1 = U/Ca - tanb1;                      
        tana2 = tana1 + (tanb2 - tanb1);            
        alpha1(i) = atan(tana1);  alpha2(i) = atan(tana2);

        Cw1(i) = Ca*tana1;  Cw2(i) = Ca*tana2;
        Lambda_out(i) = Lam_i;
    end

    PR_st(i) = (1 + eta_p*dT0i/T01_st(i))^(gam/(gam-1));
    P03_st(i) = PR_st(i)*P01_st(i);
    T03_st(i) = T01_st(i) + dT0i;

    DeHaller(i) = cos(beta1(i))/cos(beta2(i));
    Deflection_deg(i) = rad2deg(beta1(i) - beta2(i));

    fprintf('  Stage %d: b1=%6.2f  b2=%6.2f  a1=%6.2f  a2=%6.2f deg | PR=%.3f | DeH=%.3f\n', ...
        i, rad2deg(beta1(i)), rad2deg(beta2(i)), rad2deg(alpha1(i)), rad2deg(alpha2(i)), PR_st(i), DeHaller(i));
end
fprintf('\n');

%% ---------------------------------------------------------------------
%  5. STEP 3a : EXIT ANNULUS SIZING
% ------------------------------------------------------------------------
r_mean = zeros(1,N_stages); r_root = zeros(1,N_stages); r_tip = zeros(1,N_stages);
h_blade = zeros(1,N_stages); A3_st = zeros(1,N_stages);

fprintf('=== STEP 3a: EXIT ANNULUS SIZING ===\n');
for i = 1:N_stages
    C3sq = Ca^2 + Cw2(i)^2;
    T3_st  = T03_st(i) - C3sq/(2*Cp);
    P3_st  = P03_st(i) * (T3_st/T03_st(i))^(gam/(gam-1));
    rho3   = P3_st/(R*T3_st);
    A3_st(i) = mdot/(rho3*Ca);

    r_mean(i) = rm1;  % constant mean-line design
    h_blade(i)= A3_st(i)/(2*pi*r_mean(i));
    r_root(i) = r_mean(i) - h_blade(i)/2;
    r_tip(i)  = r_mean(i) + h_blade(i)/2;

    fprintf('  Stage %d: r_root=%.4f | r_tip=%.4f | h=%.4f m\n', ...
        i, r_root(i), r_tip(i), h_blade(i));
end
fprintf('\n');

%% ---------------------------------------------------------------------
%  6. STEP 3b : AIR ANGLE VARIATION WITH RADIUS (FIRST-POWER LAW)
% ------------------------------------------------------------------------
station_names = {'Root','Mean','Tip'};
beta1_rmt = zeros(N_stages,3); beta2_rmt = zeros(N_stages,3);
alpha1_rmt = zeros(N_stages,3); alpha2_rmt = zeros(N_stages,3);

fprintf('=== STEP 3b: AIR ANGLES (%s) ===\n', VORTEX_LAW);
for i = 1:N_stages
    radii = [r_root(i), r_mean(i), r_tip(i)];
    for k = 1:3
        r = radii(k);
        % First-Power Law: Swirl velocity proportional to radius
        Cw1_r = Cw1(i) * (r / r_mean(i));
        Cw2_r = Cw2(i) * (r / r_mean(i));
        
        % SREE Integration Check: Trapping imaginary Ca logic
        radicand1 = Ca^2 + 2*(Cw1(i)^2 - Cw1_r^2);
        radicand2 = Ca^2 + 2*(Cw2(i)^2 - Cw2_r^2);
        
        if radicand1 < 0
            Ca1_r = NaN; % Flow recirculation / stagnation
            fprintf('*** WARNING: Stage %d %s Ca1 stagnation (imaginary flow)\n', i, station_names{k});
        else
            Ca1_r = sqrt(radicand1);
        end
        
        if radicand2 < 0
            Ca2_r = NaN;
            fprintf('*** WARNING: Stage %d %s Ca2 stagnation (imaginary flow)\n', i, station_names{k});
        else
            Ca2_r = sqrt(radicand2);
        end
        
        U_r = omega * r;
        
        beta1_rmt(i,k)  = atan( (U_r - Cw1_r)/Ca1_r );
        beta2_rmt(i,k)  = atan( (U_r - Cw2_r)/Ca2_r );
        alpha1_rmt(i,k) = atan( Cw1_r/Ca1_r );
        alpha2_rmt(i,k) = atan( Cw2_r/Ca2_r );
    end
end
fprintf('\n');

%% =======================================================================
%  7. PLOTS
% =========================================================================
stages_x = 1:N_stages;

% Plot 1: Annulus
figure('Name','Annulus Sizing');
plot(stages_x, r_tip,  '-o','LineWidth',1.8); hold on;
plot(stages_x, r_mean, '-s','LineWidth',1.8);
plot(stages_x, r_root, '-^','LineWidth',1.8);
xlabel('Stage'); ylabel('Radius (m)'); title('Annulus Radius'); legend('Tip','Mean','Root'); grid on;

% Plot 2: De Haller
figure('Name','De Haller Number');
bar(stages_x, DeHaller, 'FaceColor',[0.30 0.55 0.80]); hold on;
plot([0.5, N_stages+0.5], [DEHALLER_LIMIT, DEHALLER_LIMIT], 'r--','LineWidth',1.6);
xlabel('Stage'); ylabel('V_2/V_1'); title('De Haller Number (red = 0.72 limit)'); grid on;

% Plot 3: Air Angle Spread (Handles NaNs gracefully)
figure('Name','Root vs Tip Angles');
subplot(2,1,1);
plot(stages_x, rad2deg(beta1_rmt(:,1)),'-^', stages_x, rad2deg(beta1_rmt(:,3)),'-o');
title('\beta_1 Variation (\beta_1 Root vs Tip)'); legend('Root','Tip'); ylabel('Deg'); grid on;
subplot(2,1,2);
plot(stages_x, rad2deg(beta2_rmt(:,1)),'-^', stages_x, rad2deg(beta2_rmt(:,3)),'-o');
title('\beta_2 Variation (\beta_2 Root vs Tip)'); legend('Root','Tip'); xlabel('Stage'); ylabel('Deg'); grid on;