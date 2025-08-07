%% Vibration Analysis of a Spring-Mass-Damper System
% Author:    Syed Arsalan Gul
% Date:      06/08/2025
% About:     This script simulates a 1-DOF spring-mass-damper system subjected to
%            harmonic excitation. It solves the motion numerically, and performs
%            FFT-based frequency analysis

%% ------------------ INPUT SECTION ------------------
% System Parameters
V.m = 1.0;          % Mass [kg]
V.c = 2.0;          % Damping Coefficient [Ns/m]
V.k = 20.0;         % Stiffness [N/m]

% Excitation Force
V.F = @(t) 5 * sin(3*t);  % External harmonic force [N]

% Initial Conditions
V.x0 = [0; 0];       % [Initial displacement; Initial velocity]

% Time Domain
V.tspan = [0 20];    % Simulation duration [s]

%% ------------ COMPUTATIONAL ANALYSIS --------------
% Define ODE system
V.odeFunc = @(t, x)[
    x(2); 
    (V.F(t) - V.c * x(2) - V.k * x(1)) / V.m
];

% Solve ODE using ode45
[t, x] = ode45(V.odeFunc, V.tspan, V.x0);

% Extract Displacement and Velocity
V.x = x(:,1);     % Displacement [m]
V.v = x(:,2);     % Velocity [m/s]

% System characteristics
V.omega_n = sqrt(V.k/V.m);                % Natural frequency [rad/s]
V.zeta = V.c / (2 * sqrt(V.k*V.m));       % Damping ratio

% Frequency-Domain Analysis
V.Fs = 1 / mean(diff(t));         % Sampling frequency [Hz]
L = length(V.x);                  % Signal length
Y = fft(V.x);                     % FFT of displacement
P2 = abs(Y / L);                  % Two-sided spectrum
P1 = P2(1:L/2+1);                 % Single-sided spectrum
P1(2:end-1) = 2 * P1(2:end-1);
f = V.Fs * (0:(L/2)) / L;         % Frequency axis [Hz]

% Dominant frequency (peak in spectrum)
[maxMag, idxMax] = max(P1);
dominantFreq = f(idxMax);

%% ---------------- OUTPUT SECTION ------------------
fprintf('>> Vibration Analysis Results:\n');
fprintf('-----------------------------------------------------\n');
fprintf('Mass (m)                         : %.2f kg\n', V.m);
fprintf('Damping Coefficient (c)          : %.2f Ns/m\n', V.c);
fprintf('Stiffness (k)                    : %.2f N/m\n', V.k);
fprintf('Excitation Frequency             : %.2f rad/s\n', 3);
fprintf('Natural Frequency (ω_n)          : %.2f rad/s\n', V.omega_n);
fprintf('Damping Ratio (ζ)                : %.3f\n', V.zeta);
fprintf('Dominant Response Frequency      : %.2f Hz\n', dominantFreq);
fprintf('-----------------------------------------------------\n');

%% ----------- INTERPRETATION OF RESULTS -------------
fprintf('>> INTERPRETATION:\n');

% Interpretation: Damping Type
if V.zeta < 1
    dampingType = 'Underdamped';
elseif V.zeta == 1
    dampingType = 'Critically Damped';
else
    dampingType = 'Overdamped';
end
fprintf('System is %s (ζ = %.3f).\n', dampingType, V.zeta);

% Interpretation: Resonance Possibility
exc_freq = 3; % Excitation frequency in rad/s
freq_ratio = exc_freq / V.omega_n;
if abs(freq_ratio - 1) < 0.1
    fprintf('Excitation frequency is close to the natural frequency → Possible resonance!\n');
else
    fprintf('Excitation frequency is NOT close to the natural frequency → Resonance unlikely.\n');
end

% Interpretation: Steady-State Achievement
steady_state_threshold = 0.01;
if max(abs(V.x(end-100:end))) < steady_state_threshold
    fprintf('System response has settled → Overdamped or fast decay.\n');
else
    fprintf('System exhibits sustained oscillations → Due to continuous forcing.\n');
end
fprintf('-----------------------------------------------------\n\n');

%% ---------------- PLOTTING SECTION ------------------

% === 1. Displacement vs Time Plot ===
figure;
plot(t, V.x, 'b-', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Displacement [m]');
title('Time-Domain Response: Displacement');
grid on;

% Add centered top annotation
yl = ylim;
xl = xlim;
text(mean(xl), yl(2)*0.95, ...
    ['→ Steady-state oscillation due to continuous harmonic forcing.' newline ...
     '→ Amplitude remains bounded due to damping.'], ...
    'FontSize', 10, 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontWeight', 'bold', 'BackgroundColor', 'white');

% === 2. Velocity vs Time Plot ===
figure;
plot(t, V.v, 'r--', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Velocity [m/s]');
title('Time-Domain Response: Velocity');
grid on;

% Add centered top annotation
yl = ylim;
xl = xlim;
text(mean(xl), yl(2)*0.95, ...
    ['→ Velocity oscillates at same frequency as displacement.' newline ...
     '→ Peaks occur when displacement crosses zero (out of phase).'], ...
    'FontSize', 10, 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontWeight', 'bold', 'BackgroundColor', 'white');

% === 3. Frequency Spectrum Plot ===
figure;
plot(f, P1, 'k-', 'LineWidth', 2);
xlabel('Frequency [Hz]');
ylabel('|X(f)|');
title('Frequency Spectrum of Displacement');
grid on;

% Mark dominant frequency
hold on;
plot(dominantFreq, maxMag, 'ro', 'MarkerSize', 8, 'LineWidth', 2);

% Add centered top annotation
yl = ylim;
xl = xlim;
text(mean(xl), yl(2)*0.95, ...
    sprintf(['→ Dominant Frequency: %.2f Hz — matches excitation.' newline ...
             '→ Indicates steady-state driven response.'], dominantFreq), ...
    'FontSize', 10, 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontWeight', 'bold', 'BackgroundColor', 'white');

%% ---------------- SAVE ANALYSIS ------------------

save('SpringMassDamper_Analysis', 'V', 't', 'x', 'f', 'P1');
