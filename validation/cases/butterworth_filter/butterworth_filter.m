% butterworth_filter.m -- 4th-order Butterworth lowpass filter step response.
% Cascaded biquad sections matching ctrlpp::make_butterworth<4>().
%
% Usage: octave --no-gui butterworth_filter.m > butterworth_filter_octave.csv

pkg load signal;

% Design 4th-order Butterworth at 50 Hz, sampled at 1000 Hz
[b, a] = butter(4, 50 / 500);  % normalized cutoff = fc / (fs/2)

% Step response
n_samples = 200;
x = ones(1, n_samples);
y = filter(b, a, x);

fprintf('sample,response\n');
for i = 1:n_samples
    fprintf('%.15e,%.15e\n', i-1, y(i));
end
