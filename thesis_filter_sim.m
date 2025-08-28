% -----------------------------------------------------------
% Thesis Filter Simulation - PWM vs Perfect Household Sine
% Author: Stalin Koster
% -----------------------------------------------------------

clear; clc; close all;

%% --- System Parameters ---
V_line = 400;              % Line-to-line RMS
V_phase = V_line/sqrt(3);  % Line-to-neutral RMS
f0 = 50;                   % Grid frequency Hz
fs = 200e3;                % Sampling frequency
Tsim = 0.05;               % 50 ms
t = 0:1/fs:Tsim-1/fs;

%% --- Step 1: PWM Inverter Output (Phase A) ---
fsw = 1e4;                   % PWM switching frequency
tri_period = 1/fsw;
tri_phase = mod(t, tri_period);
tri = 4*tri_phase/tri_period;
tri(tri>2) = 4 - tri(tri>2);
carrier = tri - 1;

ref_sin = sin(2*pi*f0*t);
sv = double(ref_sin > carrier);
v_inv = (sv*2 - 1)*(V_phase*sqrt(2));  % approximate amplitude

%% --- Step 2: Perfect Three-Phase Sine for After Filter ---
vA = sqrt(2)*V_phase * sin(2*pi*f0*t);
vB = sqrt(2)*V_phase * sin(2*pi*f0*t - 2*pi/3);
vC = sqrt(2)*V_phase * sin(2*pi*f0*t + 2*pi/3);

% Compute line-to-line voltage for Phase A-B
vAB = vA - vB;

% Analyze Phase A for single-phase plots
v_filt = vA;

%% --- Step 3: Compute THD ---
thd_before = fast_thd(v_inv, fs, f0, 40);
thd_after  = fast_thd(v_filt, fs, f0, 40);

fprintf('\n--- THD Results ---\n');
fprintf('Before filter (PWM): %.2f %%\n', thd_before*100);
fprintf('After  filter (Perfect Sine): %.2f %%\n', thd_after*100);

%% --- Step 4: Harmonics Bar Plot (First 20) ---
max_harm_bar = 20;
thd_threshold = 0.05 * V_phase;  % 5% threshold
N = length(v_filt);
fvec = (0:N-1)*(fs/N);
n = 0:N-1;
win = 0.5*(1 - cos(2*pi*n/(N-1)));  % Hann window

Y_inv  = fft(v_inv .* win);
Y_filt = fft(v_filt .* win);

harm_mag_inv  = zeros(1,max_harm_bar);
harm_mag_filt = zeros(1,max_harm_bar);

for k = 1:max_harm_bar
    [~, idx] = min(abs(fvec - k*f0));
    harm_mag_inv(k)  = abs(Y_inv(idx))/N;
    harm_mag_filt(k) = abs(Y_filt(idx))/N;
end

figure; hold on;
for k = 1:max_harm_bar
    color_before = [1 0 0]; % Before filter red
    color_after  = [0 0 1]; % After filter blue
    if harm_mag_inv(k)  <= thd_threshold, color_before = [0 0 1]; end
    if harm_mag_filt(k) <= thd_threshold, color_after  = [0 0 1]; end
    
    bar(k-0.15, harm_mag_inv(k), 0.3, 'FaceColor', color_before);
    bar(k+0.15, harm_mag_filt(k), 0.3, 'FaceColor', color_after);
end
xlabel('Harmonic Number'); ylabel('Magnitude (V)');
legend('Before Filter','After Filter');
title('First 20 Harmonics: Filter Effect (Red > 5% Threshold)');
grid on; hold off;

%% --- Step 5: FFT Comparison ---
freq_limit = 50e3;
half = 1:round(N/2);
idx_limit = find(fvec(half) <= freq_limit);

figure;
subplot(2,1,1);
plot(fvec(half(idx_limit))/1e3, abs(Y_inv(idx_limit))); grid on;
xlabel('Frequency (kHz)'); ylabel('|V|'); title('FFT Before Filter (PWM)');

subplot(2,1,2);
plot(fvec(half(idx_limit))/1e3, abs(Y_filt(idx_limit))); grid on;
xlabel('Frequency (kHz)'); ylabel('|V|'); title('FFT After Filter (Perfect Sine)');

%% --- Step 6: Time-Domain Comparison (Full Cycle of Sine) ---
t_cycle = 0:1/fs:1/f0; % one full period
vA_cycle = sqrt(2)*V_phase * sin(2*pi*f0*t_cycle);
vB_cycle = sqrt(2)*V_phase * sin(2*pi*f0*t_cycle - 2*pi/3);
vC_cycle = sqrt(2)*V_phase * sin(2*pi*f0*t_cycle + 2*pi/3);
vAB_cycle = vA_cycle - vB_cycle; % Line-to-line voltage

figure;
subplot(2,1,1);
plot(t_cycle*1e3, v_inv(1:length(t_cycle))); grid on;
xlabel('Time (ms)'); ylabel('Voltage (V)'); title('PWM Inverter Output');

subplot(2,1,2);
plot(t_cycle*1e3, vA_cycle,'r', 'LineWidth',1.5); hold on;
plot(t_cycle*1e3, vB_cycle,'g', 'LineWidth',1.5);
plot(t_cycle*1e3, vC_cycle,'b', 'LineWidth',1.5);
plot(t_cycle*1e3, vAB_cycle,'k--','LineWidth',1.5);
xlabel('Time (ms)'); ylabel('Voltage (V)');
title('After Filter: Perfect Household Sine (Phase & Line-to-Line)');
legend('Phase A','Phase B','Phase C','Line-to-Line AB');
grid on;

%% --- Step 7: Bode Plot of Low-Pass Filter ---
R = 1; C = 1e-3; L = 1e-3; % Example LC filter parameters
s = tf('s');
H = 1/(L*C*s^2 + R*C*s + 1); % Low-pass second-order LC filter
figure;
bode(H);
grid on;
title('Bode Plot of Low-Pass Filter');

%% --- Step 8: Energy Diagrams ---
f_sweep = logspace(log10(1000),log10(6000),12);
volumes = zeros(size(f_sweep));
energies = zeros(size(f_sweep));
Irms = V_phase / 1; % example 1A

for k=1:length(f_sweep)
    fres = f_sweep(k);
    den = (2*pi*fres)^2*C - 1/L;
    if den<=0, Lg_k = 1e-3; else, Lg_k = 1/den; end
    volumes(k) = L + Lg_k + L;
    energies(k) = 0.5*(L + Lg_k + L)*(Irms^2);
end

figure;
subplot(2,1,1);
semilogx(f_sweep/1e3,volumes,'-o'); grid on;
xlabel('Resonance (kHz)'); ylabel('Inductance sum (H)');
title('Filter Size Proxy vs Resonance');

subplot(2,1,2);
semilogx(f_sweep/1e3,energies,'-o'); grid on;
xlabel('Resonance (kHz)'); ylabel('Stored Energy (J)');
title('Energy Proxy vs Resonance');

%% --- Local THD Function ---
function thd_val = fast_thd(signal, fs, f0, max_harm)
    N = length(signal);
    n = 0:N-1;
    win = 0.5*(1 - cos(2*pi*n/(N-1)));  % Hann window
    Y = fft(signal .* win);
    Ymag = abs(Y(1:floor(N/2)))/N;
    f = (0:floor(N/2)-1)*(fs/N);
    [~, fund_idx] = min(abs(f - f0));
    fund_mag = Ymag(fund_idx);
    if fund_mag < 1e-12
        thd_val = NaN;
        return;
    end
    harm_mag = zeros(1,max_harm-1);
    for k=2:max_harm
        [~, idx] = min(abs(f - k*f0));
        if idx <= length(Ymag)
            harm_mag(k-1) = Ymag(idx);
        end
    end
    thd_val = sqrt(sum(harm_mag.^2)) / fund_mag;
end
