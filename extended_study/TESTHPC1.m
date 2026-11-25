%% =====================================================================
%  AXIAL COMPRESSOR STAGE-BY-STAGE MEAN-LINE DESIGN
%  Method: Cohen, Rogers & Saravanamuttoo, Gas Turbine Theory
%  Reverse-engineered from RR Olympus Mk.531 LPC data (7 stages)
%
%  Edit the INPUT block only, then press Run.
%  Ten figure windows open: 5 plots and 5 tables.
%
%  One blade speed: U = omega * r_mean, from the given RPM and the
%  inlet annulus. The flow coefficient phi = Ca/U is an OUTPUT.
%  Stage pressure ratios use the true polytropic relation.
%  Annulus areas use stator-exit velocity (station 3), not rotor-exit.
%  first_power is Dixon's constant-reaction law, not solid-body swirl.
% =====================================================================
clear; clc; close all;

%% ---------------------------------------------------------------------
%  1. INPUT
% ------------------------------------------------------------------------
gam   = 1.4;
R     = 287;                          % J/kg.K
Cp    = gam*R/(gam-1);

mdot     = 131.2447;                  % kg/s
N_rpm    = 8382;                      % rpm, LP spool (authoritative)
T01      = 610.294;                    % K
P01      = 289.988e3;                   % Pa
eta_p    = 0.8402;                    % polytropic efficiency
nu_ht    = 0.58;                      % hub-to-tip at inlet
M1       = 0.499;                     % inlet Mach, no swirl so C1 = Ca
PR_total = 2.9;
N_stages = 7;

% Stage 1 kept low (tip is already slightly transonic).
% Stage 2 ramps. Stages 3..N are equalised so the last stage is not a dump.
dT0_stage1 = 20.0;                    % K
dT0_stage2 = 25.0;                    % K
dT0    = [dT0_stage1, dT0_stage2, nan, nan, nan, nan, nan];
lambda = [0.98, 0.93, 0.88, 0.83, 0.83, 0.83, 0.83];   % work-done factor
% Stage 1 has no IGV: Lambda is an OUTPUT there.
Lambda = [ nan, 0.70, 0.50, 0.50, 0.50, 0.50, 0.50 ];

VORTEX_LAW = 'exponential';
%   'free'        : Cw * r = const, Ca = const.
%   'first_power' : Cw1 = a*r - b/r , Cw2 = a*r + b/r  (constant reaction)
%   'exponential' : Cw1 = a - b/r   , Cw2 = a + b/r    (Dixon zero-power)
%                   Constant work. Milder twist than free vortex.
%                   Stage 1 follows the assignment sheet (U = r*omega).

DEHALLER_LIMIT = 0.72;

assert(numel(dT0)==N_stages && numel(lambda)==N_stages && numel(Lambda)==N_stages, ...
    'dT0, lambda and Lambda must each have N_stages entries');

fprintf('=== INPUT SUMMARY ===\n');
fprintf('Stages: %d | mdot=%.4f kg/s | N=%.0f rpm | PR_total=%.3f\n', ...
    N_stages, mdot, N_rpm, PR_total);
fprintf('eta_p=%.4f (polytropic) | nu_ht=%.3f | M1=%.3f | vortex=%s\n\n', ...
    eta_p, nu_ht, M1, VORTEX_LAW);

%% ---------------------------------------------------------------------
%  2. STEP 1 : INLET
% ------------------------------------------------------------------------
% No IGV => C1 = Ca, taken from the given inlet Mach number.
% Annulus from continuity and nu_ht. Blade speed from the shaft:
%     U = omega * r_mean
% phi = Ca/U is then known. It is not an independent input.

P1 = P01 * (1 + (gam - 1)/2 * M1^2)^(-gam / (gam - 1));
T1   = T01 / (1 + (gam-1)/2 * M1^2);
rho1 = P1 / (R*T1);
a1   = sqrt(gam*R*T1);
Ca   = M1 * a1;
A1   = mdot / (rho1*Ca);

rt1 = sqrt( A1 / (pi*(1-nu_ht^2)) );
rr1 = nu_ht * rt1;
rm1 = 0.5*(rr1 + rt1);

omega  = N_rpm*2*pi/60;
U      = omega * rm1;                 % the only mean blade speed
phi    = Ca / U;
U_tip1 = omega * rt1;
Mrel_tip1 = sqrt(Ca^2 + U_tip1^2) / a1;

fprintf('=== STEP 1: INLET ===\n');
fprintf('T1 = %.2f K | rho1 = %.4f kg/m3 | a1 = %.2f m/s\n', T1, rho1, a1);
fprintf('Ca (=C1, no IGV) = %.2f m/s\n', Ca);
fprintf('A1 = %.4f m2 | rt1 = %.4f m | rr1 = %.4f m | rm1 = %.4f m\n', A1, rt1, rr1, rm1);
fprintf('omega = %.2f rad/s\n', omega);
fprintf('U_mean = omega*rm1 = %.2f m/s   (used in every later step)\n', U);
fprintf('U_tip  = omega*rt1 = %.2f m/s\n', U_tip1);
fprintf('phi    = Ca/U      = %.3f     (output, not an input)\n', phi);
fprintf('Tip relative Mach (stage-1 LE, no swirl) = %.3f\n\n', Mrel_tip1);

%% ---------------------------------------------------------------------
%  3. STEP 1b : EQUALISE STAGES 3..N
% ------------------------------------------------------------------------
% Unknown X = common dT0 for stages 3..N-1.
% Choose X so the last stage, reverse-engineered from the leftover PR,
% comes out equal to X as well.

X_lo = 1;  X_hi = 150;
for iter = 1:60
    X_try = 0.5*(X_lo + X_hi);
    T_try = T01;  PR_prod = 1;
    for k = 1:N_stages-1
        if k==1,      dTk = dT0_stage1;
        elseif k==2,  dTk = dT0_stage2;
        else,         dTk = X_try;
        end
        PR_prod = PR_prod * pr_poly(dTk, T_try, eta_p, gam);
        T_try   = T_try + dTk;
    end
    dT_last = dT_from_pr(PR_total/PR_prod, T_try, eta_p, gam);
    if dT_last - X_try > 0, X_lo = X_try; else, X_hi = X_try; end
end
dT0_balanced = 0.5*(X_lo + X_hi);
dT0 = [dT0_stage1, dT0_stage2, repmat(dT0_balanced,1,N_stages-3), NaN];

fprintf('=== STEP 1b: BALANCED STAGE LOADING ===\n');
fprintf('Common dT0 for stages 3-%d = %.3f K\n', N_stages-1, dT0_balanced);
fprintf('(stage %d reverse-engineers to the same value; overall PR is hit exactly)\n\n', N_stages);

%% ---------------------------------------------------------------------
%  4. STEP 2 : MEAN-LINE VELOCITY TRIANGLES
% ------------------------------------------------------------------------
beta1 = zeros(1,N_stages);  beta2 = zeros(1,N_stages);
alpha1 = zeros(1,N_stages); alpha2 = zeros(1,N_stages);
Cw1 = zeros(1,N_stages);    Cw2 = zeros(1,N_stages);
PR_st = zeros(1,N_stages);
P01_st = zeros(1,N_stages); T01_st = zeros(1,N_stages);
P03_st = zeros(1,N_stages); T03_st = zeros(1,N_stages);
Lambda_out = zeros(1,N_stages);
DeHaller = zeros(1,N_stages);
Defl_deg = zeros(1,N_stages);
dT0_used = dT0;
psi = zeros(1,N_stages);

fprintf('=== STEP 2: MEAN-LINE STAGE DESIGN ===\n');
for i = 1:N_stages

    if i == 1
        T01_st(i) = T01;  P01_st(i) = P01;
    else
        T01_st(i) = T03_st(i-1);  P01_st(i) = P03_st(i-1);
    end

    if i == N_stages && isnan(dT0(i))
        PR_left = PR_total / prod(PR_st(1:i-1));
        dT0_used(i) = dT_from_pr(PR_left, T01_st(i), eta_p, gam);
        fprintf('  Stage %d [REVERSE-ENGINEERED]: PR_left = %.4f  -> dT0 = %.3f K\n', ...
            i, PR_left, dT0_used(i));
    end

    dT0i = dT0_used(i);
    lam  = lambda(i);
    psi(i) = Cp*dT0i / (lam*U*U);

    if isnan(Lambda(i))
        % no IGV
        alpha1(i) = 0;  Cw1(i) = 0;
        Cw2(i)    = Cp*dT0i / (lam*U);
        beta1(i)  = atan(U/Ca);
        beta2(i)  = atan( (U-Cw2(i))/Ca );
        alpha2(i) = atan( Cw2(i)/Ca );
        Lambda_out(i) = 1 - (Cw2(i)+Cw1(i))/(2*U);
    else
        % prescribed reaction:  tan b1 - tan b2 = dT0*Cp/(lam*U*Ca)
        %                       tan b1 + tan b2 = Lambda*2U/Ca
        Lam = Lambda(i);
        rhs_diff = dT0i*Cp / (lam*U*Ca);
        rhs_sum  = Lam*2*U / Ca;
        tanb1 = 0.5*(rhs_sum + rhs_diff);
        tanb2 = 0.5*(rhs_sum - rhs_diff);
        beta1(i) = atan(tanb1);  beta2(i) = atan(tanb2);
        tana1 = U/Ca - tanb1;
        tana2 = tana1 + (tanb1 - tanb2);
        alpha1(i) = atan(tana1);  alpha2(i) = atan(tana2);
        Cw1(i) = Ca*tana1;  Cw2(i) = Ca*tana2;
        Lambda_out(i) = Lam;
    end

    PR_st(i)  = pr_poly(dT0i, T01_st(i), eta_p, gam);
    P03_st(i) = PR_st(i)*P01_st(i);
    T03_st(i) = T01_st(i) + dT0i;

    W1 = hypot(Ca, U-Cw1(i));
    W2 = hypot(Ca, U-Cw2(i));
    DeHaller(i) = W2/W1;
    Defl_deg(i) = rad2deg(beta1(i) - beta2(i));

    fprintf(['  Stage %d: b1=%6.2f  b2=%6.2f  a1=%6.2f  a2=%6.2f deg | ', ...
             'Lambda=%.3f | PR=%.4f | DeHaller=%.3f | Defl=%.2f deg | ', ...
             'dT0=%.2f K | psi=%.3f\n'], ...
        i, rad2deg(beta1(i)), rad2deg(beta2(i)), rad2deg(alpha1(i)), rad2deg(alpha2(i)), ...
        Lambda_out(i), PR_st(i), DeHaller(i), Defl_deg(i), dT0i, psi(i));
end

PR_achieved = prod(PR_st);
fprintf('\nOverall PR achieved = %.4f  (target = %.4f)\n', PR_achieved, PR_total);
fprintf('Overall T03 = %.2f K | Overall P03 = %.2f kPa\n', T03_st(end), P03_st(end)/1000);
fprintf('Sum dT0 = %.2f K | phi = %.3f (constant Ca, constant rm)\n\n', sum(dT0_used), phi);

if any(DeHaller < DEHALLER_LIMIT)
    bad = find(DeHaller < DEHALLER_LIMIT);
    fprintf('*** WARNING: De Haller < %.2f at stage(s) %s\n', DEHALLER_LIMIT, mat2str(bad));
    fprintf('*** Rear stages are highly loaded: U is the real blade speed (%.0f m/s),\n', U);
    fprintf('*** so psi is high. Raise N_rpm, add a stage, or lower PR_total to ease this.\n\n');
end

%% ---------------------------------------------------------------------
%  5. STEP 3a : EXIT ANNULUS  (stator exit = station 3)
% ------------------------------------------------------------------------
% C3 is the absolute speed AFTER the stator.
%   stages 1..N-1 : stator sets up the next rotor, so Cw3(i) = Cw1(i+1)
%   last stage    : axial discharge, Cw3 = 0
% Constant mean radius (Cohen). Hub-to-tip rises because density rises.
% That is expected.

r_mean  = zeros(1,N_stages);
r_root  = zeros(1,N_stages);
r_tip   = zeros(1,N_stages);
h_blade = zeros(1,N_stages);
T3_st   = zeros(1,N_stages);
P3_st   = zeros(1,N_stages);
rho3_st = zeros(1,N_stages);
A3_st   = zeros(1,N_stages);
nu_ht_st = zeros(1,N_stages);

fprintf('=== STEP 3a: EXIT ANNULUS (stator exit) ===\n');
fprintf('Constant mean radius. nu_ht rises with density — expected, not an error.\n');
for i = 1:N_stages
    if i < N_stages,  Cw3 = Cw1(i+1);  else,  Cw3 = 0;  end
    C3sq = Ca^2 + Cw3^2;
    T3_st(i)   = T03_st(i) - C3sq/(2*Cp);
    assert(T3_st(i) > 0, 'T3 <= 0 at stage %d', i);
    P3_st(i)   = P03_st(i) * (T3_st(i)/T03_st(i))^(gam/(gam-1));
    rho3_st(i) = P3_st(i) / (R*T3_st(i));
    A3_st(i)   = mdot / (rho3_st(i)*Ca);

    r_mean(i)  = rm1;
    h_blade(i) = A3_st(i) / (2*pi*r_mean(i));
    r_root(i)  = r_mean(i) - h_blade(i)/2;
    r_tip(i)   = r_mean(i) + h_blade(i)/2;
    assert(r_root(i) > 0, 'hub radius <= 0 at stage %d — annulus collapsed', i);
    nu_ht_st(i) = r_root(i)/r_tip(i);

    fprintf('  Stage %d: T3=%.2f K | P3=%.2f kPa | A=%.4f m2 | h=%.4f m | r_root=%.4f | r_tip=%.4f | nu_ht=%.3f\n', ...
        i, T3_st(i), P3_st(i)/1000, A3_st(i), h_blade(i), r_root(i), r_tip(i), nu_ht_st(i));
end
fprintf('  Inlet :                          A=%.4f m2 | h=%.4f m | r_root=%.4f | r_tip=%.4f | nu_ht=%.3f\n\n', ...
    A1, rt1-rr1, rr1, rt1, nu_ht);

%% ---------------------------------------------------------------------
%  6. STEP 3b : ROOT / MEAN / TIP
% ------------------------------------------------------------------------
% U(r) = omega * r. No other blade-speed scale exists.
% Station 1 is evaluated at the stage-inlet radius, station 2 at the
% stage-exit radius (the hub and casing both move in a tapered annulus).

station_names = {'Root','Mean','Tip'};
beta1_rmt  = nan(N_stages,3); beta2_rmt  = nan(N_stages,3);
alpha1_rmt = nan(N_stages,3); alpha2_rmt = nan(N_stages,3);
Ca1_rmt = nan(N_stages,3);    Ca2_rmt = nan(N_stages,3);
Lam_rmt = nan(N_stages,3);

r1_st = zeros(N_stages,3);    r2_st = zeros(N_stages,3);
r1_st(1,:) = [rr1, rm1, rt1];
r2_st(1,:) = [r_root(1), r_mean(1), r_tip(1)];
for i = 2:N_stages
    r1_st(i,:) = [r_root(i-1), r_mean(i-1), r_tip(i-1)];
    r2_st(i,:) = [r_root(i),   r_mean(i),   r_tip(i)];
end

DH_rmt = nan(N_stages,3);  Defl_rmt = nan(N_stages,3);

fprintf('=== STEP 3b: AIR ANGLES (%s) ===\n', VORTEX_LAW);
fprintf('Stage 1 (sheet): alpha1=0, Ca=const, Cw2=mean-line, U=r*omega.\n');
fprintf('Later stages: %s. Station 1 at inlet radius, station 2 at exit radius.\n', VORTEX_LAW);
n_nan = 0;

for i = 1:N_stages
    if i==1, law_i = 'sheet'; else, law_i = VORTEX_LAW; end
    fprintf('  Stage %d  [%s]:\n', i, law_i);

    for k = 1:3
        r1 = r1_st(i,k);  r2 = r2_st(i,k);
        U1 = omega*r1;    U2 = omega*r2;

        if i == 1
            U1 = omega*r1;  U2 = U1;
            Cw1_r = 0;  Cw2_r = Cw2(1);
            Ca1_r = Ca; Ca2_r = Ca;
        elseif strcmp(VORTEX_LAW,'first_power') && ~isnan(Lambda(i))
            [a_fp, b_fp] = first_power_ab(Cw1(i), Cw2(i), rm1);
            Cw1_r = a_fp*r1 - b_fp/r1;
            Cw2_r = a_fp*r2 + b_fp/r2;
            Ca1_r = first_power_ca(a_fp, b_fp, rm1, r1, Ca, 1);
            Ca2_r = first_power_ca(a_fp, b_fp, rm1, r2, Ca, 2);
        elseif strcmp(VORTEX_LAW,'exponential') && ~isnan(Lambda(i))
            [a_e, b_e] = exp_ab(Cw1(i), Cw2(i), rm1);
            Cw1_r = a_e - b_e/r1;
            Cw2_r = a_e + b_e/r2;
            Ca1_r = exp_ca(a_e, b_e, rm1, r1, Ca, 1);
            Ca2_r = exp_ca(a_e, b_e, rm1, r2, Ca, 2);
        else
            Cw1_r = Cw1(i)*(rm1/r1);
            Cw2_r = Cw2(i)*(rm1/r2);
            Ca1_r = Ca;  Ca2_r = Ca;
        end
        if isnan(Ca1_r) || isnan(Ca2_r)
            n_nan = n_nan + 1;
            fprintf('    *** WARNING: %s: Ca^2<=0 under %s. NaN.\n', station_names{k}, VORTEX_LAW);
        end

        Ca1_rmt(i,k) = Ca1_r;  Ca2_rmt(i,k) = Ca2_r;
        if ~isnan(Ca1_r) && Ca1_r > 0
            beta1_rmt(i,k)  = atan( (U1-Cw1_r)/Ca1_r );
            alpha1_rmt(i,k) = atan( Cw1_r/Ca1_r );
        end
        if ~isnan(Ca2_r) && Ca2_r > 0
            beta2_rmt(i,k)  = atan( (U2-Cw2_r)/Ca2_r );
            alpha2_rmt(i,k) = atan( Cw2_r/Ca2_r );
        end
        if ~isnan(Ca1_r) && ~isnan(Ca2_r)
            Lam_rmt(i,k) = 1 - (Cw1_r+Cw2_r)/(U1+U2);
            DH_rmt(i,k)  = hypot(Ca2_r, U2-Cw2_r) / hypot(Ca1_r, U1-Cw1_r);
            Defl_rmt(i,k)= rad2deg(beta1_rmt(i,k) - beta2_rmt(i,k));
        end

        fprintf('    %-5s: b1=%6.2f b2=%6.2f a1=%6.2f a2=%6.2f deg | U=%6.1f | DH=%.3f | defl=%5.2f deg\n', ...
            station_names{k}, ...
            rad2deg(beta1_rmt(i,k)), rad2deg(beta2_rmt(i,k)), ...
            rad2deg(alpha1_rmt(i,k)), rad2deg(alpha2_rmt(i,k)), ...
            U1, DH_rmt(i,k), Defl_rmt(i,k));
    end
end

fprintf('\n=== STAGE 1 TABLE (assignment sheet) ===\n');
fprintf('  %-5s %8s %8s %8s %8s %8s %8s %8s\n', ...
    '', 'beta1', 'beta2', 'alpha1', 'alpha2', 'U', 'DH', 'defl');
for k = 1:3
    fprintf('  %-5s %8.2f %8.2f %8.2f %8.2f %8.1f %8.3f %8.2f\n', ...
        station_names{k}, rad2deg(beta1_rmt(1,k)), rad2deg(beta2_rmt(1,k)), ...
        rad2deg(alpha1_rmt(1,k)), rad2deg(alpha2_rmt(1,k)), ...
        omega*r1_st(1,k), DH_rmt(1,k), Defl_rmt(1,k));
end
fprintf('\n=== DE HALLER & DEFLECTION (mean line) ===\n');
for i = 1:N_stages
    fprintf('  Stage %d:  De Haller = %.3f   deflection = %.2f deg\n', i, DeHaller(i), Defl_deg(i));
end
fprintf('\n');
if n_nan == 0
    fprintf('No unrealisable stations. Vortex law is closed at root, mean and tip of every stage.\n\n');
end

fprintf('=== SELF-CHECK: Step 3b MEAN row vs Step 2 ===\n');
ok = true;
for i = 1:N_stages
    db1 = abs(rad2deg(beta1_rmt(i,2)) - rad2deg(beta1(i)));
    db2 = abs(rad2deg(beta2_rmt(i,2)) - rad2deg(beta2(i)));
    if db1 > 0.05 || db2 > 0.05
        ok = false;
        fprintf('  Stage %d: |db1|=%.3f  |db2|=%.3f deg  *** MISMATCH\n', i, db1, db2);
    else
        fprintf('  Stage %d: |db1|=%.3f  |db2|=%.3f deg\n', i, db1, db2);
    end
end
if ok
    fprintf('  Mean row reproduces Step 2. One U, one design.\n\n');
end

%% =====================================================================
%  7. PLOTS AND TABLES  (every figure opens in MATLAB)
% =====================================================================
stages_x = 1:N_stages;
x_m = 0:N_stages;
rt_m = [rt1, r_tip];
rm_m = [rm1, r_mean];
rr_m = [rr1, r_root];

% ---- 1. Annulus -------------------------------------------------------
figure('Name','1. Annulus (meridional)');
plot(x_m, rt_m, '-o', 'LineWidth', 1.8, 'DisplayName', 'Tip'); hold on;
plot(x_m, rm_m, '-s', 'LineWidth', 1.8, 'DisplayName', 'Mean');
plot(x_m, rr_m, '-^', 'LineWidth', 1.8, 'DisplayName', 'Root');
xlabel('Station (0 = inlet, n = stage-n exit)'); ylabel('Radius (m)');
title('Annulus (constant mean radius)');
legend('Location','eastoutside'); grid on; xlim([0 N_stages]);

% ---- 2+3. De Haller and deflection, numbers on the bars ---------------
figure('Name','2. De Haller and deflection','Position',[80 80 1100 420]);
subplot(1,2,1);
bh = bar(stages_x, DeHaller, 'FaceColor', [0.30 0.55 0.80]); hold on;
plot([0.5 N_stages+0.5], [DEHALLER_LIMIT DEHALLER_LIMIT], 'r--', 'LineWidth', 1.6);
for i = 1:N_stages
    text(i, DeHaller(i)+0.02, sprintf('%.3f', DeHaller(i)), ...
        'HorizontalAlignment','center', 'FontSize', 9);
end
xlabel('Stage'); ylabel('De Haller W_2/W_1');
title('De Haller number (mean line)');
ylim([0 1.15]); grid on; legend('De Haller', sprintf('%.2f limit', DEHALLER_LIMIT));

subplot(1,2,2);
bar(stages_x, Defl_deg, 'FaceColor', [0.80 0.45 0.25]); hold on;
for i = 1:N_stages
    text(i, Defl_deg(i)+0.8, sprintf('%.1f^\\circ', Defl_deg(i)), ...
        'HorizontalAlignment','center', 'FontSize', 9);
end
xlabel('Stage'); ylabel('Deflection \beta_1 - \beta_2 (deg)');
title('Rotor deflection (mean line)');
ylim([0 max(Defl_deg)*1.25]); grid on;

% ---- 4. Velocity triangles --------------------------------------------
figure('Name','3. Velocity triangles','Position',[50 50 1400 700]);
ncols = ceil(sqrt(N_stages)); nrows = ceil(N_stages/ncols);
for i = 1:N_stages
    subplot(nrows, ncols, i); hold on;
    xoff = 1.3*max(U, Ca);
    quiver(0,    0, U,          0,  0, 'Color', [0.6 0.6 1], 'LineWidth', 1.2);
    quiver(0,    0, Cw1(i),     Ca, 0, 'b',                 'LineWidth', 1.6);
    quiver(U,    0, Cw1(i)-U,   Ca, 0, 'b--',               'LineWidth', 1.3);
    quiver(xoff, 0, U,          0,  0, 'Color', [1 0.6 0.6],'LineWidth', 1.2);
    quiver(xoff, 0, Cw2(i),     Ca, 0, 'r',                 'LineWidth', 1.6);
    quiver(xoff+U,0,Cw2(i)-U,   Ca, 0, 'r--',               'LineWidth', 1.3);
    axis equal; grid on;
    title(sprintf('Stage %d', i));
    if i == 1
        legend({'U','C_1','W_1','U','C_2','W_2'}, ...
            'Location','southoutside', 'Orientation','horizontal');
    end
end
put_sgtitle('Velocity triangles (blue = rotor inlet, red = rotor exit, solid = abs, dashed = rel)');

% ---- 5. Air angles vs radius ------------------------------------------
figure('Name','4. Air angles vs radius','Position',[50 50 1200 800]);
angle_sets = {beta1_rmt, '\beta_1 (rotor inlet, relative)'; ...
              beta2_rmt, '\beta_2 (rotor exit, relative)'; ...
              alpha1_rmt,'\alpha_1 (rotor inlet, absolute)'; ...
              alpha2_rmt,'\alpha_2 (rotor exit, absolute)'};
for p = 1:4
    subplot(2,2,p);
    data_deg = rad2deg(angle_sets{p,1});
    plot(stages_x, data_deg(:,1), '-^', 'LineWidth', 1.6, 'DisplayName', 'Root'); hold on;
    plot(stages_x, data_deg(:,2), '-s', 'LineWidth', 1.6, 'DisplayName', 'Mean');
    plot(stages_x, data_deg(:,3), '-o', 'LineWidth', 1.6, 'DisplayName', 'Tip');
    xlabel('Stage'); ylabel('Angle (deg)');
    title(angle_sets{p,2});
    legend('Location','best'); grid on;
end
put_sgtitle(sprintf('Air angle variation with radius (%s)', strrep(VORTEX_LAW,'_',' ')));

% ---- 6. Stage coefficients --------------------------------------------
figure('Name','5. Stage coefficients');
plot(stages_x, lambda,     '-o', 'LineWidth', 1.6, 'DisplayName', '\lambda work-done'); hold on;
plot(stages_x, Lambda_out, '-s', 'LineWidth', 1.6, 'DisplayName', '\Lambda reaction');
plot(stages_x, psi,        '-^', 'LineWidth', 1.6, 'DisplayName', '\psi loading');
plot([0.5 N_stages+0.5], [phi phi], 'k:', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('\\phi = %.3f', phi));
xlabel('Stage'); title('Stage coefficients'); xlim([0.5 N_stages+0.5]);
legend('Location','best'); grid on;

% ---- 7. TABLE: mean-line results --------------------------------------
mean_cols = {'dT0_K','PR','Lambda','psi','beta1','beta2','alpha1','alpha2','DeHaller','defl_deg'};
mean_data = [dT0_used(:), PR_st(:), Lambda_out(:), psi(:), ...
             rad2deg(beta1(:)), rad2deg(beta2(:)), ...
             rad2deg(alpha1(:)), rad2deg(alpha2(:)), ...
             DeHaller(:), Defl_deg(:)];
show_table('6. TABLE  Mean-line stages', ...
    mean_cols, arrayfun(@(k)sprintf('St %d',k), 1:N_stages, 'uni', 0), mean_data);

% ---- 8. TABLE: stage 1 root / mean / tip (assignment sheet) -----------
s1_cols = {'r_m','U_mps','beta1','beta2','alpha1','alpha2','DeHaller','defl_deg'};
s1_data = [r1_st(1,:)', (omega*r1_st(1,:))', ...
           rad2deg(beta1_rmt(1,:))', rad2deg(beta2_rmt(1,:))', ...
           rad2deg(alpha1_rmt(1,:))', rad2deg(alpha2_rmt(1,:))', ...
           DH_rmt(1,:)', Defl_rmt(1,:)'];
show_table('7. TABLE  Stage 1 root / mean / tip', ...
    s1_cols, station_names, s1_data);

% ---- 9. TABLE: exit annulus -------------------------------------------
ann_cols = {'T3_K','P3_kPa','A_m2','h_m','r_root','r_tip','nu_ht'};
ann_data = [T3_st(:), P3_st(:)/1000, A3_st(:), h_blade(:), r_root(:), r_tip(:), nu_ht_st(:)];
show_table('8. TABLE  Exit annulus', ...
    ann_cols, arrayfun(@(k)sprintf('St %d',k), 1:N_stages, 'uni', 0), ann_data);

% ---- 10. TABLE: every stage, root / mean / tip ------------------------
rmt_cols = {'beta1','beta2','alpha1','alpha2','U_mps','Lambda','DeHaller','defl_deg'};
rmt_rows = {};
rmt_data = zeros(N_stages*3, 8);
n = 0;
for i = 1:N_stages
    for k = 1:3
        n = n + 1;
        rmt_rows{n} = sprintf('S%d %s', i, station_names{k}); %#ok<SAGROW>
        rmt_data(n,:) = [rad2deg(beta1_rmt(i,k)), rad2deg(beta2_rmt(i,k)), ...
                         rad2deg(alpha1_rmt(i,k)), rad2deg(alpha2_rmt(i,k)), ...
                         omega*r1_st(i,k), Lam_rmt(i,k), DH_rmt(i,k), Defl_rmt(i,k)];
    end
end
show_table('9. TABLE  Root / mean / tip all stages', rmt_cols, rmt_rows, rmt_data);

% ---- 11. TABLE: inlet summary -----------------------------------------
in_cols = {'value'};
in_rows = {'T1_K','rho1','Ca_mps','U_mean','U_tip','phi','rt1_m','rr1_m','rm1_m','Mrel_tip'};
in_data = [T1; rho1; Ca; U; U_tip1; phi; rt1; rr1; rm1; Mrel_tip1];
show_table('10. TABLE  Inlet', in_cols, in_rows, in_data);

fprintf('=== DONE. 10 figure windows (plots + tables). ===\n');
fprintf('Use the Windows menu in MATLAB to flip between them.\n');

% ---- 12. Exponential vs first-power ----------------------------------
fp_b1 = nan(N_stages,3); fp_b2 = nan(N_stages,3);
fp_DH = nan(N_stages,3); fp_Ca2 = nan(N_stages,3);
for i = 2:N_stages
    [af, bf] = first_power_ab(Cw1(i), Cw2(i), rm1);
    for k = 1:3
        r1 = r1_st(i,k); r2 = r2_st(i,k);
        U1 = omega*r1;   U2 = omega*r2;
        Cw1f = af*r1 - bf/r1;  Cw2f = af*r2 + bf/r2;
        Ca1f = first_power_ca(af, bf, rm1, r1, Ca, 1);
        Ca2f = first_power_ca(af, bf, rm1, r2, Ca, 2);
        fp_Ca2(i,k) = Ca2f;
        if ~isnan(Ca1f) && Ca1f > 0, fp_b1(i,k) = atan((U1-Cw1f)/Ca1f); end
        if ~isnan(Ca2f) && Ca2f > 0, fp_b2(i,k) = atan((U2-Cw2f)/Ca2f); end
        if ~isnan(Ca1f) && ~isnan(Ca2f)
            fp_DH(i,k) = hypot(Ca2f, U2-Cw2f) / hypot(Ca1f, U1-Cw1f);
        end
    end
end

figure('Name','11. Exponential vs first-power angles','Position',[60 60 1100 720]);
subplot(2,1,1);
plot(stages_x, rad2deg(beta1_rmt(:,1)), '-^', 'LineWidth', 1.6, 'DisplayName', 'exp root'); hold on;
plot(stages_x, rad2deg(beta1_rmt(:,3)), '-o', 'LineWidth', 1.6, 'DisplayName', 'exp tip');
plot(stages_x, rad2deg(fp_b1(:,1)), '--^', 'LineWidth', 1.4, 'DisplayName', '1st-p root');
plot(stages_x, rad2deg(fp_b1(:,3)), '--o', 'LineWidth', 1.4, 'DisplayName', '1st-p tip');
plot(stages_x, rad2deg(beta1_rmt(:,2)), ':s', 'Color', [0.4 0.4 0.4], 'DisplayName', 'mean (same)');
ylabel('\beta_1 (deg)'); grid on; legend('Location','best');
title('Relative inlet angle');
subplot(2,1,2);
plot(stages_x, rad2deg(beta2_rmt(:,1)), '-^', 'LineWidth', 1.6, 'DisplayName', 'exp root'); hold on;
plot(stages_x, rad2deg(beta2_rmt(:,3)), '-o', 'LineWidth', 1.6, 'DisplayName', 'exp tip');
plot(stages_x, rad2deg(fp_b2(:,1)), '--^', 'LineWidth', 1.4, 'DisplayName', '1st-p root');
plot(stages_x, rad2deg(fp_b2(:,3)), '--o', 'LineWidth', 1.4, 'DisplayName', '1st-p tip');
plot(stages_x, rad2deg(beta2_rmt(:,2)), ':s', 'Color', [0.4 0.4 0.4], 'DisplayName', 'mean (same)');
ylabel('\beta_2 (deg)'); xlabel('Stage'); grid on; legend('Location','best');
title('Relative exit angle');
put_sgtitle('Exponential (solid) vs first-power (dashed). Mean line is identical.');

figure('Name','12. De Haller root/tip, both laws');
plot(2:N_stages, DH_rmt(2:end,1), '-^', 'LineWidth', 1.6, 'DisplayName', 'exp root'); hold on;
plot(2:N_stages, DH_rmt(2:end,3), '-o', 'LineWidth', 1.6, 'DisplayName', 'exp tip');
plot(2:N_stages, fp_DH(2:end,1), '--^', 'LineWidth', 1.4, 'DisplayName', '1st-p root');
plot(2:N_stages, fp_DH(2:end,3), '--o', 'LineWidth', 1.4, 'DisplayName', '1st-p tip');
plot([1.5 N_stages+0.5], [0.72 0.72], 'k:', 'DisplayName', '0.72');
xlabel('Stage'); ylabel('De Haller'); title('De Haller at root and tip');
legend('Location','best'); grid on; ylim([0 1.15]);

cmp_cols = {'b1_exp','b1_fp','b2_exp','b2_fp','DH_exp','DH_fp','Ca2_exp','Ca2_fp'};
cmp_rows = {};
cmp_data = zeros((N_stages-1)*3, 8);
n = 0;
for i = 2:N_stages
    for k = 1:3
        n = n + 1;
        cmp_rows{n} = sprintf('S%d %s', i, station_names{k}); %#ok<SAGROW>
        cmp_data(n,:) = [rad2deg(beta1_rmt(i,k)), rad2deg(fp_b1(i,k)), ...
                         rad2deg(beta2_rmt(i,k)), rad2deg(fp_b2(i,k)), ...
                         DH_rmt(i,k), fp_DH(i,k), Ca2_rmt(i,k), fp_Ca2(i,k)];
    end
end
show_table('13. TABLE  Exponential vs first-power', cmp_cols, cmp_rows, cmp_data);


%% =====================================================================
%  Local functions
% =====================================================================
function PR = pr_poly(dT, Tin, eta_p, gam)
    PR = (1 + dT/Tin)^(gam*eta_p/(gam-1));
end

function dT = dT_from_pr(PR, Tin, eta_p, gam)
    dT = Tin * (PR^((gam-1)/(gam*eta_p)) - 1);
end

function [a, b] = first_power_ab(Cw1m, Cw2m, rm)
    a = (Cw1m + Cw2m) / (2*rm);
    b = (Cw2m - Cw1m) * rm / 2;
end

function Ca_r = first_power_ca(a, b, rm, r, Ca_m, station)
    lnrr = log(r/rm);
    dr2  = r^2 - rm^2;
    if station == 1
        Casq = Ca_m^2 - 2*a^2*dr2 + 4*a*b*lnrr;
    else
        Casq = Ca_m^2 - 2*a^2*dr2 - 4*a*b*lnrr;
    end
    if Casq > 0, Ca_r = sqrt(Casq); else, Ca_r = NaN; end
end

function [a, b] = exp_ab(Cw1m, Cw2m, rm)
    a = (Cw1m + Cw2m) / 2;
    b = (Cw2m - Cw1m) * rm / 2;
end

function Ca_r = exp_ca(a, b, rm, r, Ca_m, station)
    lnrr = log(r/rm);
    inv  = 1/r - 1/rm;
    if station == 1, s = -1; else, s = 1; end
    Casq = Ca_m^2 - 2*a^2*lnrr + 2*a*s*b*inv;
    if Casq > 0, Ca_r = sqrt(Casq); else, Ca_r = NaN; end
end

function put_sgtitle(str)
    if exist('sgtitle','file') || exist('sgtitle','builtin')
        sgtitle(str);
    else
        annotation('textbox', [0 0.95 1 0.05], 'String', str, ...
            'HorizontalAlignment','center', 'EdgeColor','none', ...
            'FontWeight','bold', 'FontSize', 11);
    end
end

function show_table(fig_name, col_names, row_names, data)
    nR = size(data,1);
    nC = size(data,2);
    fig = figure('Name', fig_name, 'Position', [120 120 max(720, 90*nC+160) max(220, 28*nR+80)]);
    uit = uitable(fig, ...
        'Data', data, ...
        'ColumnName', col_names, ...
        'RowName', row_names, ...
        'Units', 'normalized', ...
        'Position', [0.02 0.04 0.96 0.92], ...
        'FontSize', 10, ...
        'ColumnWidth', 'auto');
    try
        uit.ColumnFormat = repmat({'bank'}, 1, nC);
    catch
    end
end

