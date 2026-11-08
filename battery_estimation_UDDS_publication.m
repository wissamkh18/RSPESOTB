%% Robust State-Parameter Estimation in Second-Order Thvenin Batteries
%  Second-order Thevenin model:
%     x1' = -theta*u,  x2' = -a*x2 + g1*u,  x3' = -b*x3 + g2*u,
%     y   = p(x1) - x2 - x3 - Rs*u,
%  theta = 1/Q, a = 1/(R1*C1), b = 1/(R2*C2), g1 = 1/C1, g2 = 1/C2.
%
%  Case 1: identification of R1, R2, C1, C2 over rests. The simulation is
%          stopped shortly after the second usable rest.
%  Case 2: online estimator over one real EPA UDDS cycle with R1, R2, C1, C2
%          known, and initial estimates different from the true state:
%          - polarization observers (open loop),
%          - Rs updated at each current jump (hybrid least squares),
%          - x1 and theta = 1/Q propagated by Coulomb counting and updated
%            every Tu by a sliding-window regression (hybrid estimator).
%  Both nominal (onoff=0) and disturbed (onoff=1) online cases are run.
%  Identification in Case 1 remains noise-free.
%  Publication outputs: one 2x2 identification figure, one 2x2 nominal
%  figure, one nominal SOC figure, one 1x2 disturbed figure, and one UDDS
%  current-profile figure. All figures are exported as EPS.
clear; clc; close all;
TuChoice = 60;      % update period [s]; tune this in practice
TChoice  = 300;     % sliding regression window [s]

% Consistent publication formatting for all figures.
lineWidth = 1.35;
fontSize = 9.5;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex', ...
           'defaultTextInterpreter',          'latex', ...
           'defaultLegendInterpreter',        'latex', ...
           'defaultLegendBox',                'off', ...
           'defaultAxesFontSize',             fontSize, ...
           'defaultAxesBox',                  'on', ...
           'defaultAxesXGrid',                'on', ...
           'defaultAxesYGrid',                'on', ...
           'defaultLineLineWidth',            lineWidth);

%% ======================= Common settings =======================
% True parameters
P.Rs = 0.0013;  P.R1 = 0.0042;  P.R2 = 0.0020;      % ohm
P.C1 = 17111;   P.C2 = 440.5675;                    % F
P.Q  = 72000;                                       % C (20 Ah)
P.a  = 1/(P.R1*P.C1);  P.b = 1/(P.R2*P.C2);
P.g1 = 1/P.C1;         P.g2 = 1/P.C2;

% OCV map p (descending powers of s), from the paper
P.ocv  = [8.4073 -19.892 11.497 4.161 -4.5533 0.34365 0.64685 3.5016];
P.docv = polyder(P.ocv);
fprintf('min p''(s) on [0,1] = %.4f (must be > 0)\n', ...
        min(polyval(P.docv, linspace(0,1,2001))));

dt   = 0.1;          % measurement sampling period [s]
% No artificial Imax scaling: current comes from the vehicle/pack model.
x0   = [0.90; 0.02; 0.01];              % true initial state (not relaxed)

% Required input: epa_udds_speed.csv, columns [time_s, speed_mph].
% Genuine EPA UDDS trace at 1 Hz, converted to representative cell current.
% Use exactly one original EPA UDDS cycle for both identification and estimation.
[uOne, ~, tOne] = udds_current_epa();
assert(abs(dt*round(1/dt)-1)<1e-10, 'dt must divide 1 s exactly');
% Each CSV speed sample represents a 1 s ZOH interval.
u1Hz = uOne(:);
nUp = round(1/dt);
u   = repelem(u1Hz,nUp);
N   = numel(u);
t   = (0:N-1)'*dt;
[x, y, yMinus] = simulate_plant(P, u, x0, dt, zeros(N,3), zeros(N,1));
if any(x(1,:) < 0 | x(1,:) > 1)
    error('SOC leaves [0,1]: adjust the operational scenario; do not clip SOC silently.');
end
fprintf('Using one original EPA UDDS cycle: duration %.0f s, current [%.2f, %.2f] A\n', ...
    N*dt,min(u),max(u));

%% ======================= Case 1: identification over rests =======================
minRest = 2;                 % minimum rest length used [s]
% Values shown before identification (only for display, e.g., catalogue data)
R1d = 2.0*P.R1;  R2d = 0.5*P.R2;  C1d = 0.5*P.C1;  C2d = 2.0*P.C2;

idLog = nan(N,4);            % [R1 R2 C1 C2]
cur   = [R1d R2d C1d C2d];
rests = struct('s',{},'e',{},'a',{},'b',{},'x2s',{},'x3s',{},'x2e',{},'x3e',{});
inRest = false;  kStop = N;  nUsable = 0;
for k = 1:N
    if u(k) == 0 && ~inRest, inRest = true; kRest = k; end
    if inRest && (u(k) ~= 0 || k == N)
        inRest = false;

        kEnd = k - (u(k) ~= 0);
        idx  = kRest:kEnd;
        if numel(idx)*dt >= minRest && numel(idx) >= 5
            r = identify_rest(y(idx), dt);
            if r.ok
                nR = numel(idx);
                rI = struct('s',kRest,'e',kEnd+1,'a',r.a,'b',r.b, ...
                    'x2s',r.x2s,'x3s',r.x3s, ...
                    'x2e',r.x2s*exp(-r.a*nR*dt),'x3e',r.x3s*exp(-r.b*nR*dt));
                if ~isempty(rests)
                    pr = rests(end);
                    Fa = load_integral(u, dt, pr.e, rI.s, r.a);
                    Fb = load_integral(u, dt, pr.e, rI.s, r.b);
                    L  = (rI.s - pr.e)*dt;
                    g1 = (rI.x2s - exp(-r.a*L)*pr.x2e)/Fa;
                    g2 = (rI.x3s - exp(-r.b*L)*pr.x3e)/Fb;
                    cur = [g1/r.a, g2/r.b, 1/g1, 1/g2];
                end
                rests(end+1) = rI;            %#ok<SAGROW>
                nUsable = nUsable + 1;
                if nUsable == 2
                    tIdentification = t(kEnd);  % End of the second usable rest
                    kStop = min(N, k + round(30/dt));
                end
            end
        end
    end
    idLog(k,:) = cur;
    if k >= kStop, break; end
end
t1 = t(1:kStop);  idLog = idLog(1:kStop,:);
fprintf('\n===== IDENTIFICATION TIME =====\n');
fprintf('Exact identification time = %.3f s\n', tIdentification);
fprintf('Simulation stopping time  = %.3f s\n', t1(end));
fprintf('\nCase 1: identified after two rests\n');
fprintf('  R1 = %.6g (true %.6g)   R2 = %.6g (true %.6g)\n', cur(1), P.R1, cur(2), P.R2);
fprintf('  C1 = %.6g (true %.6g)   C2 = %.6g (true %.6g)\n', cur(3), P.C1, cur(4), P.C2);

%% Identification accuracy
trueParams = [P.R1 P.R2 P.C1 P.C2];

relErrors = 100*abs(cur-trueParams)./trueParams;

fprintf('\nIdentification results:\n');
fprintf('Stopping time: %.1f s\n',t1(end));
fprintf('Maximum relative error: %.20f %%\n',max(relErrors));

fprintf('R1 error: %.20f %%\n',relErrors(1));
fprintf('R2 error: %.20f %%\n',relErrors(2));
fprintf('C1 error: %.20f %%\n',relErrors(3));
fprintf('C2 error: %.20f %%\n',relErrors(4));

% One 2x2 identification figure: R1, R2, C1, and C2.
blue = [0.00 0.36 0.67];
red  = [0.82 0.24 0.18];
figId = plot_identification_2x2(t1,idLog,[P.R1 P.R2 P.C1 P.C2],blue,red);

%% ======================= Case 2: both online estimation scenarios =======================
T  = TChoice;        % 300-s regression window (s)
Tu = TuChoice;       % 60-s update period (s); adjustable in practice
mu = 0;               % No forgetting of past jumps

fprintf('\nRunning the nominal case (onoff=0)...\n');
RESnom = run_online_case(0,P,u,x0,dt,T,Tu,mu);
fprintf('Running the disturbed/noisy case (onoff=1)...\n');
RESunc = run_online_case(1,P,u,x0,dt,T,Tu,mu);

%% Steady-state errors: disturbed case
% Evaluate errors over the final part of the UDDS cycle
idx = t >= 900;

% Series resistance: absolute error in mOhm
errRs = 1000*abs(RESunc.lg.Rs(idx) - P.Rs);

% SOC: absolute error in percentage points
x1true = RESunc.x(1,:)';
errSOC = 100*abs(RESunc.lg.x1(idx) - x1true(idx));

% Capacity: relative error in percent
errQ = 100*abs(RESunc.lg.Q(idx) - P.Q)/P.Q;

fprintf('\n===== DISTURBED CASE ERRORS (t >= 900 s) =====\n');
fprintf('Maximum Rs error  = %.10f mOhm\n', max(errRs));
fprintf('Maximum SOC error = %.10f percentage points\n', max(errSOC));
fprintf('Maximum Q error   = %.10f %%\n', max(errQ));
fprintf('Nominal: %d SOC/capacity updates; first at %.1f s\n', ...
    numel(RESnom.tUpd),first_or_nan(RESnom.tUpd));
fprintf('Disturbed: %d SOC/capacity updates; first at %.1f s\n', ...
    numel(RESunc.tUpd),first_or_nan(RESunc.tUpd));

%% ======================= Publication figures for online estimation =======================
% Nominal case: one 2x2 figure for x2, x3, Rs, Q, plus a separate x1 figure.
figNom = plot_nominal_summary(t,RESnom,P,blue,red);
figNomX1 = plot_soc_figure(t,RESnom,blue,red);

% Disturbed/noisy case: one 1x2 figure for Rs and Q, as requested.
figUnc = plot_uncertainty_summary(t,RESunc,P,blue,red);

% Publication UDDS current figure: ONE original EPA UDDS cycle only.
% This figure belongs to Section II-B (system analysis).
% The online simulation above also uses the same single original cycle.
% TOP = original UDDS, BOTTOM LEFT = rest, BOTTOM RIGHT = jump.
% No (a), (b), (c) labels inside the figure.
ucOne = uOne(:);
tcOne = tOne(:);

% Representative rest interval in the original UDDS cycle
restFlag = abs(ucOne(1:end-1)) < 1e-10;
changes = diff([false; restFlag; false]);
restBeg = find(changes == 1);
restEnd = find(changes == -1);
durations = restEnd - restBeg;
valid = find(durations >= 10 & restBeg > 1);
assert(~isempty(valid), 'The EPA current profile has no qualifying rest.');
[~,iLongest] = max(durations(valid));
ridx = valid(iLongest);
restZoom = [max(tcOne(1),tcOne(restBeg(ridx))-5), ...
            min(tcOne(end),tcOne(min(restEnd(ridx),numel(tcOne)))+5)];

% Representative current jump close to 194 s in the original cycle
jumpCandidates = find(tcOne(2:end)>=193 & tcOne(2:end)<=195);
assert(~isempty(jumpCandidates), 'No UDDS samples near 194 s');
[~,imax] = max(abs(ucOne(jumpCandidates+1)-ucOne(jumpCandidates)));
jumpTime = tcOne(jumpCandidates(imax)+1);
jumpZoom = [jumpTime-1, jumpTime+1];

% Single-column publication figure (original cycle only)
figUDDS = figure('Color','w','Units','centimeters','Position',[3 3 17.5 9.5]);
tlUDDS = tiledlayout(figUDDS,2,2,'TileSpacing','compact','Padding','compact');
axTop = nexttile(tlUDDS,[1 2]);
stairs(axTop,tcOne,ucOne,'Color',blue,'LineWidth',1.2);
hold(axTop,'on');
yline(axTop,0,'k:','LineWidth',0.7);
xlim(axTop,[tcOne(1),tcOne(end)]);
yRange = [floor(min(ucOne)/5)*5, ceil(max(ucOne)/5)*5];
if diff(yRange)==0, yRange=yRange+[-1 1]; end
ylim(axTop,yRange);
rectangle(axTop,'Position',[restZoom(1),yRange(1),diff(restZoom),diff(yRange)], ...
    'EdgeColor',[0.45 0.45 0.45],'LineStyle','--','LineWidth',0.9);
rectangle(axTop,'Position',[jumpZoom(1),yRange(1),diff(jumpZoom),diff(yRange)], ...
    'EdgeColor',[0.45 0.45 0.45],'LineStyle','-.','LineWidth',0.9);
xlabel(axTop,'$t$ (s)'); ylabel(axTop,'$u$ (A)');
style_pub_axes(axTop);

axRest = nexttile(tlUDDS,3);
stairs(axRest,tcOne,ucOne,'Color',blue,'LineWidth',1.2);
yline(axRest,0,'k:','LineWidth',0.6);
xlim(axRest,restZoom);
rr = ucOne(tcOne>=restZoom(1) & tcOne<=restZoom(2));
ylim(axRest,[min(rr)-2,max(rr)+2]);
xlabel(axRest,'$t$ (s)'); ylabel(axRest,'$u$ (A)');
style_pub_axes(axRest);

axJump = nexttile(tlUDDS,4);
stairs(axJump,tcOne,ucOne,'Color',blue,'LineWidth',1.2);
xlim(axJump,jumpZoom);
jj = ucOne(tcOne>=jumpZoom(1) & tcOne<=jumpZoom(2));
ylim(axJump,[min(jj)-2,max(jj)+2]);
xlabel(axJump,'$t$ (s)'); ylabel(axJump,'$u$ (A)');
style_pub_axes(axJump);

% Statistics for the original UDDS cycle
uStat = u1Hz(:);

% Rest intervals lasting at least 10 seconds
isRest = abs(uStat) <= 1e-10;
edges = diff([false; isRest; false]);

restStart = find(edges == 1);
restEnd = find(edges == -1);

restDuration = restEnd - restStart; % seconds at 1 Hz
nRest = sum(restDuration >= 10);

% Current jumps exceeding 0.5 A
duStat = diff(uStat);
nJumps = sum(abs(duStat) > 0.5);

fprintf('Number of rest intervals >= 10 s: %d\n', nRest);
fprintf('Number of current jumps > 0.5 A: %d\n', nJumps);

fprintf('Minimum current: %.2f A\n', min(uStat));
fprintf('Maximum current: %.2f A\n', max(uStat));

% Single-cycle statistics for Section II-B (same current)
uOriginal = uOne(:);
restOriginal = abs(uOriginal) <= 1e-10;
edgeOriginal = diff([false; restOriginal; false]);
begOriginal = find(edgeOriginal==1);
endOriginal = find(edgeOriginal==-1);
fprintf('Original UDDS: %d rests >= 10 s, %d jumps > 0.5 A\n', ...
    sum(endOriginal-begOriginal >= 10),sum(abs(diff(uOriginal))>0.5));

%% ======================= Publication-ready EPS export =======================
outdir = 'figures_eps';
if ~isfolder(outdir)
    mkdir(outdir);
end

% Identification: one 2x2 figure for R1, R2, C1, and C2.
export_pub_eps(figId,fullfile(outdir,'identification_R1_R2_C1_C2.eps'));

% Nominal case: one 2x2 figure for x2, x3, Rs, Q, and one separate x1 figure.
export_pub_eps(figNom,fullfile(outdir,'nominal_x23_Rs_Q.eps'));
export_pub_eps(figNomX1,fullfile(outdir,'nominal_x1.eps'));

% Disturbed/noisy case: one 1x2 figure for Rs and Q.
export_pub_eps(figUnc,fullfile(outdir,'uncertainties_Rs_Q.eps'));

% Original EPA UDDS current profile (one cycle, not seven).
export_pub_eps(figUDDS,fullfile(outdir,'UDDS_original_with_zooms.eps'));

fprintf('\nSaved 5 publication-ready EPS figures in %s\n',fullfile(pwd,outdir));
fprintf('  identification_R1_R2_C1_C2.eps\n');
fprintf('  nominal_x23_Rs_Q.eps\n');
fprintf('  nominal_x1.eps\n');
fprintf('  uncertainties_Rs_Q.eps\n');
fprintf('  UDDS_original_with_zooms.eps\n');

%% ======================= Local functions =======================
function [u, v, ts] = udds_current_epa()
% Reads the EPA speed cycle from epa_udds_speed.csv [time_s, speed_mph].
% EPA provides vehicle SPEED, not battery current; the latter is modeled.
filename = 'epa_udds_speed.csv';
% Try MATLAB's current folder, then the folder holding this .m script.
if ~isfile(filename)
    scriptFolder = fileparts(mfilename('fullpath'));
    candidate = fullfile(scriptFolder, filename);
    if isfile(candidate)
        filename = candidate;
    else
        % Current Folder can differ from the folder shown in File Explorer.
        [selectedFile, selectedFolder] = uigetfile( ...
            {'*.csv','CSV files (*.csv)'}, ...
            'Select your epa_udds_speed.csv');
        if isequal(selectedFile,0)
            error('UDDS CSV not selected. Select epa_udds_speed.csv to run.');
        end
        filename = fullfile(selectedFolder, selectedFile);
    end
end
fprintf('Reading UDDS file: %s\n', filename);
D = readmatrix(filename);
assert(size(D,2)>=2, 'CSV must have at least time_s and speed_mph.');
D = D(:,1:2);
D = D(all(isfinite(D),2),:);  % remove text-header / malformed rows
assert(size(D,1)>=1000, 'EPA UDDS should contain at least 1000 data points.');
ts = D(:,1);
v_mph = D(:,2);
assert(all(diff(ts)>0), 'UDDS time must be strictly increasing.');
assert(abs(ts(1))<1e-9 && all(abs(diff(ts)-1)<1e-9), ...
    'This script expects a 1 Hz EPA UDDS CSV starting at t=0 s.');

v = 0.44704*v_mph; % mph -> m/s
m = 1500;
g = 9.81;
Cr = 0.010;
rho = 1.225;
Cd = 0.29;
Af = 2.2;
grade = 0;
eta_drive = 0.90;
eta_regen = 0.65;
P_aux = 0;
Vpack_nom = 350;
Ns = 100; %#ok<NASGU> % pack configuration 100s3p
Np = 3;

acc = [0; diff(v)./diff(ts)];
F_total = m*acc + m*g*Cr*cos(grade) ...
          + 0.5*rho*Cd*Af.*v.^2 + m*g*sin(grade);
Pwheel = F_total.*v;
Pbatt = zeros(size(Pwheel));
traction = Pwheel>=0;
Pbatt(traction) = Pwheel(traction)/eta_drive + P_aux;
Pbatt(~traction) = eta_regen*Pwheel(~traction) + P_aux;
u = (Pbatt/Vpack_nom)/Np;
% Positive current: cell discharge; negative current: regenerative charge.
% u(k) is held constant for [ts(k),ts(k+1)); exact zeros remain rests.
assert(all(isfinite(u)), 'Non-finite UDDS-derived battery current.');
end

function [x, y, yMinus] = simulate_plant(P, u, x0, dt, D, w)
% Exact discretization (u and D piecewise constant on each sample).
% D = [d1 d2 d3] (N x 3): process disturbances; w (N x 1): voltage noise,
% piecewise constant, so that w(t_k^-) = w(k-1) and w(t_k^+) = w(k).
N  = numel(u);
ea = exp(-P.a*dt);  eb = exp(-P.b*dt);
x  = zeros(3,N);  x(:,1) = x0;
for k = 1:N-1
    x(1,k+1) = x(1,k) + (-u(k)/P.Q + D(k,1))*dt;
    x(2,k+1) = ea*x(2,k) + (1-ea)*(P.g1*u(k) + D(k,2))/P.a;
    x(3,k+1) = eb*x(3,k) + (1-eb)*(P.g2*u(k) + D(k,3))/P.b;
end
pv     = polyval(P.ocv, x(1,:))';
uPrev  = [u(1); u(1:end-1)];
wPrev  = [w(1); w(1:end-1)];
y      = pv - x(2,:)' - x(3,:)' - P.Rs*u     + w;       % y(t_k^+)
yMinus = pv - x(2,:)' - x(3,:)' - P.Rs*uPrev + wPrev;   % y(t_k^-)
if any(x(1,:) < 0 | x(1,:) > 1), warning('SOC left [0,1]'); end
end

function r = identify_rest(yr, h)
% Prony-type least squares on a rest interval (u = 0)
r.ok = false;
d = diff(yr(:));
if numel(d) < 4 || max(abs(d)) < 1e-12, return; end
Phi = [d(2:end-1), -d(1:end-2)];
Y   = d(3:end);
if rcond(Phi'*Phi) < 1e-14, return; end
sig  = (Phi'*Phi)\(Phi'*Y);
disc = sig(1)^2 - 4*sig(2);
if disc <= 0, return; end
rr = sort((sig(1) + [1; -1]*sqrt(disc))/2, 'descend');
if any(rr <= 0 | rr >= 1), return; end
r.a = -log(rr(1))/h;  r.b = -log(rr(2))/h;
AB  = [1 1; rr(1) rr(2)] \ d(1:2);           % amplitudes of the two modes
r.x2s = AB(1)/(1-rr(1));
r.x3s = AB(2)/(1-rr(2));
r.ok  = true;
end

function F = load_integral(u, dt, k0, k1, a)
% int_{t_{k0}}^{t_{k1}} exp(-a (t_{k1}-tau)) u(tau) dtau, u piecewise constant
kk = (k0:k1-1)';
T1 = (k1-1)*dt;  tau0 = (kk-1)*dt;  tau1 = kk*dt;
F  = sum(u(kk).*(exp(-a*(T1-tau1)) - exp(-a*(T1-tau0)))/a);
end

function s = ocv_inverse(v, ocv, docv)
% p^{-1}(Pi_p(v)) by safeguarded Newton iterations
p0 = polyval(ocv,0);  p1 = polyval(ocv,1);
v  = min(max(v, p0), p1);
s = 0.5; lo = 0; hi = 1;
for it = 1:50
    f = polyval(ocv,s) - v;
    if f > 0, hi = s; else, lo = s; end
    sn = s - f/polyval(docv,s);
    if sn <= lo || sn >= hi, sn = (lo+hi)/2; end
    if abs(sn - s) < 1e-14, s = sn; break; end
    s = sn;
end
end
function RES = run_online_case(onoff,P,u,x0,dt,T,Tu,mu)
% Exact same online estimator for both scenarios; only noise changes.
N = numel(u);
d1bar = 5e-7;       % SOC process disturbance (1/s)
d23bar = 1e-4;      % RC state disturbances (V/s)
wbar = 1e-3;        % measured voltage noise (V)
duMin = 2*onoff;    % jump threshold: 0 A nominal, 2 A noisy
rng(2);            % reproducibility of uncertain case
Dist = onoff*[d1bar*ones(N,1),d23bar*(2*rand(N,2)-1)];
w = onoff*wbar*(2*rand(N,1)-1);
[x,y,yMinus] = simulate_plant(P,u,x0,dt,Dist,w);

ea = exp(-P.a*dt); eb = exp(-P.b*dt);
NT = round(T/dt);
q = [0; cumsum(u(1:end-1))*dt]/3600; % Ah for numerical conditioning
xh = [0;0]; Rsh = 2*P.Rs; x1h = 0.5; thh = 1/(1.3*P.Q);
S = 0; firstJump = NaN; tNext = Inf; tUpd = [];
s0 = zeros(N,1); sq = s0; sq2 = s0; sx = s0; sqx = s0;
lg.x2 = nan(N,1); lg.x3 = nan(N,1); lg.Rs = nan(N,1);
lg.x1 = nan(N,1); lg.Q = nan(N,1); lg.x1m = nan(N,1);

for k = 1:N
    if k>1
        S = exp(-mu*dt)*S;
        du = u(k)-u(k-1);
        if du ~= 0 && abs(du)>=duMin
            % One-sided limits: state continuity isolates the ohmic jump.
            dy = y(k)-yMinus(k);
            S = S + du^2;
            Rsh = Rsh - du/S*(dy+Rsh*du);
            if isnan(firstJump)
                firstJump = (k-1)*dt;
                tNext = firstJump+T+dt;
            end
        end
    end
    ybar = y(k)+xh(1)+xh(2)+Rsh*u(k);
    x1m = ocv_inverse(ybar,P.ocv,P.docv);
    if k==1
        s0(k)=1; sq(k)=q(k); sq2(k)=q(k)^2; sx(k)=x1m; sqx(k)=q(k)*x1m;
    else
        s0(k)=s0(k-1)+1;
        sq(k)=sq(k-1)+q(k);
        sq2(k)=sq2(k-1)+q(k)^2;
        sx(k)=sx(k-1)+x1m;
        sqx(k)=sqx(k-1)+q(k)*x1m;
    end
    tk = (k-1)*dt;
    if tk>=tNext-1e-9
        j=k-NT;
        n=s0(k)-s0(j); Sq=sq(k)-sq(j); Sq2=sq2(k)-sq2(j);
        Sx=sx(k)-sx(j); Sqx=sqx(k)-sqx(j); qt=q(k);
        M=[n,n*qt-Sq; n*qt-Sq,n*qt^2-2*qt*Sq+Sq2]/n;
        be=[Sx;qt*Sx-Sqx]/n;
        if rcond(M)>1e-12
            vv=M\be;
            x1h=vv(1); thh=vv(2)/3600;  % return to 1/C
            tUpd(end+1)=tk; %#ok<AGROW>
        else
            fprintf('Skipped ill-conditioned regression at t = %.1f s\n',tk);
        end
        tNext=tNext+Tu;
    end
    lg.x2(k)=xh(1); lg.x3(k)=xh(2); lg.Rs(k)=Rsh;
    lg.x1(k)=x1h; lg.Q(k)=1/thh; lg.x1m(k)=x1m;
    xh(1)=ea*xh(1)+(1-ea)*P.g1/P.a*u(k);
    xh(2)=eb*xh(2)+(1-eb)*P.g2/P.b*u(k);
    x1h=x1h-thh*u(k)*dt;
end
RES.x=x; RES.y=y; RES.yMinus=yMinus;
RES.lg=lg; RES.tUpd=tUpd; RES.firstJump=firstJump;
end

function fig = plot_identification_2x2(tid,idEst,trueParams,blue,red)
% One 2x2 figure for R1, R2, C1, and C2.
fig = figure('Color','w','Units','centimeters','Position',[2 2 18 11.5]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

yLabels = {'$R_1$ ($\Omega$)','$R_2$ ($\Omega$)','$C_1$ (F)','$C_2$ (F)'};
trueLabels = {'$R_1$','$R_2$','$C_1$','$C_2$'};
estLabels = {'$\hat{R}_1$','$\hat{R}_2$','$\hat{C}_1$','$\hat{C}_2$'};

for i = 1:4
    ax = nexttile(tl,i);
    pTrue = plot(ax,tid,trueParams(i)*ones(size(tid)), ...
        '-','Color',blue,'LineWidth',1.45);
    hold(ax,'on');
    pEst = plot(ax,tid,idEst(:,i),'--','Color',red,'LineWidth',1.45);
    xlabel(ax,'$t$ (s)');
    ylabel(ax,yLabels{i});
    legend(ax,[pTrue pEst],{trueLabels{i},estLabels{i}}, ...
        'Location','best','FontSize',9.5);
    style_pub_axes(ax);
end
end

function fig = plot_nominal_summary(t,RES,P,blue,red)
% Nominal case: 2x2 panels for x2, x3, Rs, and Q.
fig = figure('Color','w','Units','centimeters','Position',[2 2 18 11.5]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax = gobjects(1,4);

ax(1) = nexttile(tl,1);
plot_pair_panel(ax(1),t,RES.x(2,:)',RES.lg.x2,'$x_2$ (V)', ...
    '$x_2$','$\hat{x}_2$',blue,red,'northwest');

ax(2) = nexttile(tl,2);
plot_pair_panel(ax(2),t,RES.x(3,:)',RES.lg.x3,'$x_3$ (V)', ...
    '$x_3$','$\hat{x}_3$',blue,red,'northwest');

ax(3) = nexttile(tl,3);
plot_pair_panel(ax(3),t,P.Rs*ones(size(t)),RES.lg.Rs,'$R_s$ ($\Omega$)', ...
    '$R_s$','$\hat{R}_s$',blue,red,'northwest');

ax(4) = nexttile(tl,4);
plot_pair_panel(ax(4),t,P.Q*ones(size(t)),RES.lg.Q,'$Q$ (C)', ...
    '$Q$','$\hat{Q}$',blue,red,'northwest');

drawnow;
add_pub_inset(ax(1),t,[RES.x(2,:)' RES.lg.x2],[0 250],[blue;red],{'-','--'});
add_pub_inset(ax(2),t,[RES.x(3,:)' RES.lg.x3],[0 12],[blue;red],{'-','--'});
add_pub_inset(ax(3),t,[P.Rs*ones(size(t)) RES.lg.Rs],[0 100],[blue;red],{'-','--'});
add_pub_inset(ax(4),t,[P.Q*ones(size(t)) RES.lg.Q],[300 1000],[blue;red],{'-','--'});
end

function fig = plot_soc_figure(t,RES,blue,red)
% Separate SOC figure for the nominal case.
fig = figure('Color','w','Units','centimeters','Position',[2 2 18 7.5]);
ax = axes('Parent',fig);
plot_pair_panel(ax,t,RES.x(1,:)',RES.lg.x1,'$x_1$', ...
    '$x_1$','$\hat{x}_1$',blue,red,'northwest');
drawnow;
add_pub_inset(ax,t,[RES.x(1,:)' RES.lg.x1],[300 700],[blue;red],{'-','--'});
end

function fig = plot_uncertainty_summary(t,RES,P,blue,red)
% Disturbed/noisy case: side-by-side Rs and Q plots.
fig = figure('Color','w','Units','centimeters','Position',[2 2 18 7.5]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

axRs = nexttile(tl,1);
plot_pair_panel(axRs,t,P.Rs*ones(size(t)),RES.lg.Rs,'$R_s$ ($\Omega$)', ...
    '$R_s$','$\hat{R}_s$',blue,red,'northwest');

axQ = nexttile(tl,2);
plot_pair_panel(axQ,t,P.Q*ones(size(t)),RES.lg.Q,'$Q$ (C)', ...
    '$Q$','$\hat{Q}$',blue,red,'northwest');

drawnow;
zoomWindow = [max(t(1),900) min(t(end),1200)];
add_pub_inset(axRs,t,[P.Rs*ones(size(t)) RES.lg.Rs], ...
    zoomWindow,[blue;red],{'-','--'});
add_pub_inset(axQ,t,[P.Q*ones(size(t)) RES.lg.Q], ...
    zoomWindow,[blue;red],{'-','--'});
end

function plot_pair_panel(ax,t,yTrue,yEstimate,yLabel,trueLabel,estLabel,blue,red,legendLocation)
pTrue = plot(ax,t,yTrue,'-','Color',blue,'LineWidth',1.45);
hold(ax,'on');
pEst = plot(ax,t,yEstimate,'--','Color',red,'LineWidth',1.45);
xlim(ax,[t(1) t(end)]);
xlabel(ax,'$t$ (s)');
ylabel(ax,yLabel);
legend(ax,[pTrue pEst],{trueLabel,estLabel}, ...
    'Location',legendLocation,'FontSize',9.5);
style_pub_axes(ax);
end

function style_pub_axes(ax)
set(ax,'FontSize',9.5,'LineWidth',0.9,'Box','on', ...
    'TickLabelInterpreter','latex','XGrid','on','YGrid','on', ...
    'GridAlpha',0.18);
ax.XLabel.FontSize = 10.5;
ax.YLabel.FontSize = 10.5;
end

function add_pub_inset(axMain,t,Y,tz,colors,styles)
% Draw a dashed data-space rectangle and a vector inset within this panel.
mask=t>=tz(1) & t<=tz(2);
assert(any(mask),'Zoom window has no data.');
yseg=Y(mask,:);
lo=min(yseg(:)); hi=max(yseg(:));
ypad=max(0.07*(hi-lo),max(abs([lo hi]))*0.003);
if ypad==0, ypad=1e-5; end
yz=[lo-ypad hi+ypad];
rectangle(axMain,'Position',[tz(1) yz(1) diff(tz) diff(yz)], ...
    'LineStyle','--','EdgeColor',[0.25 0.25 0.25],'LineWidth',0.8, ...
    'HandleVisibility','off');
% Use a figure-normalized inset positioned inside the upper-right of the panel.
set(axMain,'Units','normalized');
p=get(axMain,'Position');
pInset=[p(1)+0.56*p(3),p(2)+0.49*p(4),0.40*p(3),0.43*p(4)];
axIn=axes('Parent',ancestor(axMain,'figure'),'Units','normalized', ...
    'Position',pInset,'Box','on','Color','w','FontSize',8.5, ...
    'TickLabelInterpreter','latex');
hold(axIn,'on');
for i=1:size(Y,2)
    plot(axIn,t(mask),Y(mask,i),'Color',colors(i,:), ...
        'LineStyle',styles{i},'LineWidth',1.05);
end
xlim(axIn,tz); ylim(axIn,yz);
set(axIn,'XGrid','on','YGrid','on','LineWidth',0.8);
end

function export_pub_eps(fig,filename)
set(fig,'Color','w','Renderer','painters','PaperPositionMode','auto');
print(fig,filename,'-depsc2','-painters');
end

function v = first_or_nan(a)
if isempty(a), v = NaN; else, v = a(1); end
end
