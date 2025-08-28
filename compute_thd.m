function thd_val = compute_thd(signal, fs, f0, max_harm)
    % Computes Total Harmonic Distortion (THD) of a signal
    N = length(signal);
    win = hann(N)';            
    Y = fft(signal .* win);    
    f = (0:N-1)*(fs/N);
    half = 1:floor(N/2);
    Ymag = abs(Y(half))/N;
    f = f(half);

    % Fundamental
    [~, fund_idx] = min(abs(f - f0));
    fund_mag = Ymag(fund_idx);

    % Harmonics
    harm_mag = zeros(1,max_harm-1);
    for k = 2:max_harm
        [~, idx] = min(abs(f - k*f0));
        if idx <= length(Ymag)
            harm_mag(k-1) = Ymag(idx);
        end
    end

    % THD
    thd_val = sqrt(sum(harm_mag.^2)) / fund_mag;
end
