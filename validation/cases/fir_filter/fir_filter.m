% fir_filter.m -- FIR filter response via Octave filter().
% 5-tap moving average applied to a step + sinusoid signal.
%
% Usage: octave --no-gui fir_filter.m > fir_filter_octave.csv

% 5-tap moving average
taps = [0.2, 0.2, 0.2, 0.2, 0.2];
n_samples = 100;

% Input: step at sample 10 + sinusoid
x = zeros(1, n_samples);
for i = 1:n_samples
    if i >= 10
        x(i) = 1.0;
    end
    x(i) = x(i) + 0.3 * sin(2 * pi * 0.05 * (i-1));
end

y = filter(taps, 1, x);

fprintf('sample,input,output\n');
for i = 1:n_samples
    fprintf('%.15e,%.15e,%.15e\n', i-1, x(i), y(i));
end
