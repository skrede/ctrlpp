% generate_summary.m -- Aggregate per-case reports into top-level summaries.
% Reads all cases/*/analysis/report.md and produces:
%   analysis/summary.md         (consumer-facing)
%   analysis/detailed_report.md (reviewer-facing)
%
% Usage: octave --no-gui generate_summary.m <validation_dir>

function generate_summary()
    warning('off', 'Octave:graphics-toolkit');
    graphics_toolkit('gnuplot');

    args = argv();
    if numel(args) < 1
        val_dir = '.';
    else
        val_dir = args{1};
    end

    analysis_dir = fullfile(val_dir, 'analysis');
    mkdir(analysis_dir);

    cases_dir = fullfile(val_dir, 'cases');
    case_dirs = dir(cases_dir);

    % Collect all case results
    results = {};
    for i = 1:numel(case_dirs)
        d = case_dirs(i);
        if ~d.isdir || strcmp(d.name, '.') || strcmp(d.name, '..')
            continue;
        end
        report_path = fullfile(cases_dir, d.name, 'analysis', 'report.md');
        if ~exist(report_path, 'file')
            continue;
        end

        r.name = d.name;
        r.report_path = report_path;

        % Parse report for verdict and worst digits
        fid = fopen(report_path, 'r');
        r.verdict = 'UNKNOWN';
        r.worst_digits = Inf;
        r.best_digits = 0;
        r.worst_signal = '';
        r.best_signal = '';
        r.worst_max_abs = 0;
        r.n_signals = 0;
        r.n_rows = 0;

        while ~feof(fid)
            ln = fgetl(fid);
            if ~ischar(ln) || isempty(ln); continue; end

            % Parse verdict
            if ~isempty(strfind(ln, '**PASS**'))
                r.verdict = 'PASS';
            elseif ~isempty(strfind(ln, '**FAIL**'))
                r.verdict = 'FAIL';
            end

            % Parse rows
            tok = regexp(ln, '\*\*Rows\*\*:\s*(\d+)', 'tokens');
            if ~isempty(tok)
                r.n_rows = str2double(tok{1}{1});
            end

            % Parse table rows: | signal | max_abs | ... | digits | ... | result |
            if ln(1) == '|' && ~isempty(regexp(ln, '\|\s*\d'))
                parts = strtrim(strsplit(ln, '|'));
                parts = parts(~cellfun(@isempty, parts));
                if numel(parts) >= 8
                    sig_name = parts{1};
                    max_abs = str2double(parts{2});
                    digits_str = parts{7};

                    if strcmp(sig_name, 'time') || strcmp(sig_name, 'step') || strcmp(sig_name, 'sample')
                        continue;
                    end

                    r.n_signals = r.n_signals + 1;

                    if strcmp(digits_str, 'exact')
                        digits_val = 16;
                    else
                        digits_val = str2double(digits_str);
                    end

                    if digits_val < r.worst_digits
                        r.worst_digits = digits_val;
                        r.worst_signal = sig_name;
                        r.worst_max_abs = max_abs;
                    end
                    if digits_val > r.best_digits
                        r.best_digits = digits_val;
                        r.best_signal = sig_name;
                    end
                end
            end
        end
        fclose(fid);

        results{end+1} = r;
    end

    n_cases = numel(results);
    n_pass = sum(cellfun(@(r) strcmp(r.verdict, 'PASS'), results));
    n_fail = n_cases - n_pass;

    % Sort by worst digits (worst case first for reviewer)
    worst_digits_arr = cellfun(@(r) r.worst_digits, results);
    [~, sort_idx] = sort(worst_digits_arr);

    % Find best and worst cases for plots
    best_idx = sort_idx(end);
    worst_idx = sort_idx(1);

    % ===== Consumer summary =====
    fid = fopen(fullfile(analysis_dir, 'summary.md'), 'w');

    fprintf(fid, '# ctrlpp Validation Summary\n\n');
    fprintf(fid, '- **Date**: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '- **Cases**: %d tested, %d passed, %d failed\n', n_cases, n_pass, n_fail);
    if n_fail == 0
        fprintf(fid, '- **Result**: All cases pass within numerical tolerance\n\n');
    else
        fprintf(fid, '- **Result**: %d case(s) failed\n\n', n_fail);
    end

    fprintf(fid, 'All computations validated against GNU Octave %s with the Control and Signal packages.\n\n', version());

    fprintf(fid, '## Validated components\n\n');
    fprintf(fid, '| Case | Module | Digits | Verdict |\n');
    fprintf(fid, '|------|--------|--------|---------|\n');

    for i = 1:n_cases
        r = results{i};
        if isinf(r.worst_digits)
            digits_str = 'exact';
        else
            digits_str = sprintf('%.0f', r.worst_digits);
        end
        module = classify_module(r.name);
        fprintf(fid, '| %s | %s | %s | %s |\n', r.name, module, digits_str, r.verdict);
    end

    fprintf(fid, '\n## Not validated\n\n');
    fprintf(fid, 'The following components have no Octave toolbox equivalent and are not included:\n\n');
    fprintf(fid, '- `pid` policies: saturation, anti-windup, rate limiting, gain scheduling (nonlinear; `pid()` is linear-only)\n');
    fprintf(fid, '- `mrac` (no MRAC in any Octave package)\n');
    fprintf(fid, '- `rls`, `recursive_arx` (no RLS in any Octave package)\n');
    fprintf(fid, '- `biquad::low_pass`, `biquad::notch`, `biquad::dirty_derivative` (RBJ cookbook coefficients)\n');
    fprintf(fid, '- `ekf`, `ukf`, `particle_filter`, `manifold_ukf` (nonlinear estimators)\n');
    fprintf(fid, '- `mekf`, `complementary_filter` (SO(3) attitude estimators)\n');
    fprintf(fid, '- `mpc`, `nmpc`, `mhe`, `nmhe` (optimization-based controllers/estimators)\n');
    fprintf(fid, '- `trapezoidal_trajectory`, `double_s_trajectory`, `modified_sin_trajectory`, `modified_trap_trajectory` (velocity profiles)\n');
    fprintf(fid, '- `online_planner_2nd`, `online_planner_3rd` (real-time planners)\n');
    fprintf(fid, '- `bspline_trajectory`, `smoothing_spline` (advanced spline types)\n');
    fprintf(fid, '- `time_scaling`, `synchronize` (trajectory utilities)\n\n');
    fprintf(fid, 'These require custom Octave implementations or external toolboxes not available.\n');

    fclose(fid);

    % ===== Detailed reviewer report =====
    fid = fopen(fullfile(analysis_dir, 'detailed_report.md'), 'w');

    fprintf(fid, '# ctrlpp Validation: Detailed Report\n\n');
    fprintf(fid, '- **Date**: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '- **Reference**: GNU Octave %s, Control %s, Signal %s\n', version(), pkg_version('control'), pkg_version('signal'));
    fprintf(fid, '- **Cases**: %d tested, %d passed, %d failed\n', n_cases, n_pass, n_fail);
    fprintf(fid, '- **Default tolerances**: atol=1e-10, rtol=1e-8\n\n');

    fprintf(fid, '## Overview\n\n');
    fprintf(fid, '| Case | Module | Type | Signals | Rows | Worst digits | Worst signal | Max abs err | Verdict |\n');
    fprintf(fid, '|------|--------|------|---------|------|-------------|-------------|------------|--------|\n');

    for i = 1:n_cases
        idx = sort_idx(i);
        r = results{idx};
        module = classify_module(r.name);
        if r.n_rows > 1
            case_type = 'time-series';
        else
            case_type = 'algebraic';
        end
        if isinf(r.worst_digits)
            digits_str = 'exact';
        else
            digits_str = sprintf('%.1f', r.worst_digits);
        end
        fprintf(fid, '| %s | %s | %s | %d | %d | %s | %s | %.2e | %s |\n', ...
                r.name, module, case_type, r.n_signals, r.n_rows, ...
                digits_str, r.worst_signal, r.worst_max_abs, r.verdict);
    end

    % Best and worst case highlights
    r_best = results{best_idx};
    r_worst = results{worst_idx};

    fprintf(fid, '\n## Highlights\n\n');
    fprintf(fid, '### Best case: %s\n\n', r_best.name);
    fprintf(fid, 'Minimum %.0f digits of agreement on signal `%s`.\n\n', r_best.worst_digits, r_best.worst_signal);

    % Copy best case plots
    best_plots = dir(fullfile(cases_dir, r_best.name, 'analysis', '*.png'));
    if ~isempty(best_plots)
        [~] = copyfile(fullfile(cases_dir, r_best.name, 'analysis', best_plots(1).name), analysis_dir);
        fprintf(fid, '![%s](%s)\n\n', r_best.name, best_plots(1).name);
    end

    fprintf(fid, '### Worst (acceptable) case: %s\n\n', r_worst.name);
    fprintf(fid, 'Minimum %.1f digits of agreement on signal `%s` (max abs err = %.2e).\n\n', ...
            r_worst.worst_digits, r_worst.worst_signal, r_worst.worst_max_abs);

    worst_plots = dir(fullfile(cases_dir, r_worst.name, 'analysis', '*.png'));
    if ~isempty(worst_plots)
        [~] = copyfile(fullfile(cases_dir, r_worst.name, 'analysis', worst_plots(1).name), analysis_dir);
        fprintf(fid, '![%s](%s)\n\n', r_worst.name, worst_plots(1).name);
    end

    % Omissions section
    fprintf(fid, '## Omitted components\n\n');
    fprintf(fid, 'Components without Octave toolbox equivalents (would require reimplementing the algorithm under test):\n\n');
    fprintf(fid, '| Component | Reason |\n');
    fprintf(fid, '|-----------|--------|\n');
    fprintf(fid, '| `pid` policies (saturation, anti-windup, etc.) | Nonlinear; `pid()` is linear-only |\n');
    fprintf(fid, '| `mrac` | No MRAC in any Octave package |\n');
    fprintf(fid, '| `rls`, `recursive_arx` | No RLS in any Octave package |\n');
    fprintf(fid, '| `biquad::low_pass/notch/dirty_derivative` | RBJ cookbook design; no toolbox function |\n');
    fprintf(fid, '| `ekf` | Nonlinear; no Octave class |\n');
    fprintf(fid, '| `ukf` | Nonlinear; no Octave class |\n');
    fprintf(fid, '| `particle_filter` | Stochastic; no Octave class |\n');
    fprintf(fid, '| `manifold_ukf` | SO(3)-specific; no toolbox |\n');
    fprintf(fid, '| `mekf` | Multiplicative EKF on SO(3); no toolbox |\n');
    fprintf(fid, '| `complementary_filter` | Mahony filter; no toolbox |\n');
    fprintf(fid, '| `mpc`, `nmpc` | Optimization-based; requires solver integration |\n');
    fprintf(fid, '| `mhe`, `nmhe` | Optimization-based estimation |\n');
    fprintf(fid, '| `trapezoidal_trajectory` | No robotics trajectory package |\n');
    fprintf(fid, '| `double_s_trajectory` | 7-segment S-curve; no toolbox |\n');
    fprintf(fid, '| `modified_sin/trap_trajectory` | Specialized profiles; no toolbox |\n');
    fprintf(fid, '| `online_planner_*` | Real-time planners; no toolbox |\n');
    fprintf(fid, '| `bspline_trajectory` | B-spline evaluation; no toolbox |\n');
    fprintf(fid, '| `smoothing_spline` | Regularized spline; no toolbox |\n');
    fprintf(fid, '| `time_scaling`, `synchronize` | Trajectory utilities |\n');

    fprintf(fid, '\n## Per-case details\n\n');
    fprintf(fid, 'Full per-signal statistics and plots are in each `cases/<name>/analysis/report.md`.\n');

    fclose(fid);

    fprintf('Summary written to %s\n', fullfile(analysis_dir, 'summary.md'));
    fprintf('Detailed report written to %s\n', fullfile(analysis_dir, 'detailed_report.md'));
end

function module = classify_module(name)
    if ~isempty(regexp(name, 'dare|lqr|lqi|pid|mrac|pole'))
        module = 'control';
    elseif ~isempty(regexp(name, 'c2d|analysis|tf_ss|propagate'))
        module = 'model';
    elseif ~isempty(regexp(name, 'kalman|luenberger|ekf|ukf'))
        module = 'estimation';
    elseif ~isempty(regexp(name, 'biquad|butterworth|fir|filter'))
        module = 'dsp';
    elseif ~isempty(regexp(name, 'spline|trajectory|trap'))
        module = 'trajectory';
    elseif ~isempty(regexp(name, 'so3'))
        module = 'lie';
    elseif ~isempty(regexp(name, 'rls|arx|n4sid'))
        module = 'sysid';
    else
        module = 'other';
    end
end

function v = pkg_version(name)
    pkgs = pkg('list', name);
    if ~isempty(pkgs)
        v = pkgs{1}.version;
    else
        v = 'n/a';
    end
end

generate_summary();
