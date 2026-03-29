% validate_compare.m -- Universal CSV signal comparator for ctrlpp validation.
%
% Usage: octave --no-gui validate_compare.m <ref.csv> <cand.csv> [atol] [rtol] [output_dir]
%
% Compares two CSV files column-by-column. Both must have identical column
% headers (first row). Comparison uses combined absolute + relative tolerance:
%   pass iff |ref - cand| <= atol + rtol * |ref|   for every element.
%
% When output_dir is given, writes:
%   - Per-signal overlay + error plots (PNG)
%   - report.md with statistics and embedded images
%
% Exits with code 0 on PASS, 1 on FAIL, 2 on error.

function validate_compare()
    warning('off', 'Octave:graphics-toolkit');
    graphics_toolkit('gnuplot');

    args = argv();

    if numel(args) < 2
        fprintf(2, 'usage: validate_compare.m <ref.csv> <cand.csv> [atol] [rtol] [output_dir]\n');
        exit(2);
    end

    ref_file  = args{1};
    cand_file = args{2};

    atol = 1e-10;
    rtol = 1e-8;
    output_dir = '';
    if numel(args) >= 3, atol = str2double(args{3}); end
    if numel(args) >= 4, rtol = str2double(args{4}); end
    if numel(args) >= 5, output_dir = args{5}; end

    [ref_hdr,  ref_data]  = load_csv(ref_file);
    [cand_hdr, cand_data] = load_csv(cand_file);

    if numel(ref_hdr) ~= numel(cand_hdr)
        fprintf(2, 'FAIL: column count mismatch (%d vs %d)\n', numel(ref_hdr), numel(cand_hdr));
        exit(1);
    end

    for i = 1:numel(ref_hdr)
        if ~strcmp(ref_hdr{i}, cand_hdr{i})
            fprintf(2, 'FAIL: column %d header mismatch: "%s" vs "%s"\n', i, ref_hdr{i}, cand_hdr{i});
            exit(1);
        end
    end

    if size(ref_data, 1) ~= size(cand_data, 1)
        fprintf(2, 'FAIL: row count mismatch (%d vs %d)\n', size(ref_data, 1), size(cand_data, 1));
        exit(1);
    end

    n_rows = size(ref_data, 1);
    n_cols = numel(ref_hdr);
    is_timeseries = n_rows > 1;

    % Detect time column
    time_col = 0;
    if is_timeseries
        for i = 1:n_cols
            if strcmp(ref_hdr{i}, 'time') || strcmp(ref_hdr{i}, 't')
                time_col = i;
                break;
            end
        end
    end

    % Compute per-column statistics
    all_pass = true;
    col_stats = struct();

    for col = 1:n_cols
        r = ref_data(:, col);
        c = cand_data(:, col);
        err = abs(r - c);
        tol = atol + rtol * abs(r);
        violations = err > tol;

        s.name = ref_hdr{col};
        s.max_abs_err = max(err);
        s.mean_abs_err = mean(err);
        s.std_err = std(err);

        % Use a robust threshold to avoid division by near-zero reference values
        ref_scale = max(abs(r));
        nz_thresh = max(eps, ref_scale * 1e-12);
        nz = abs(r) > nz_thresh;
        rel_err = zeros(size(r));
        if any(nz)
            rel_err(nz) = err(nz) ./ abs(r(nz));
            s.max_rel_err = max(rel_err(nz));
            s.mean_rel_err = mean(rel_err(nz));
        else
            s.max_rel_err = 0;
            s.mean_rel_err = 0;
        end

        % Digits of agreement: -log10(max relative error)
        if s.max_rel_err > 0
            s.digits = -log10(s.max_rel_err);
        else
            s.digits = Inf;
        end

        % Relative L2 norm: ||e||_2 / ||ref||_2
        ref_norm = norm(r);
        if ref_norm > eps
            s.rel_l2 = norm(err) / ref_norm;
        else
            s.rel_l2 = 0;
        end

        % Error growth rate (time-series only)
        s.growth_rate = NaN;
        s.growth_class = 'n/a';
        if is_timeseries && n_rows > 10 && col ~= time_col
            nz_err = err > eps;
            if sum(nz_err) > 10
                if time_col > 0
                    t_nz = ref_data(nz_err, time_col);
                else
                    t_nz = find(nz_err);
                end
                log_err = log10(err(nz_err));
                % Linear fit: log10(err) = a*t + b
                p = polyfit(t_nz, log_err, 1);
                s.growth_rate = p(1);
                if abs(p(1)) < 0.01
                    s.growth_class = 'stable';
                elseif p(1) < 0.1
                    s.growth_class = 'linear';
                else
                    s.growth_class = 'exponential';
                end
            else
                s.growth_class = 'negligible';
            end
        end

        s.pass = ~any(violations);

        if any(violations)
            worst_idx = find(err == max(err(violations)), 1);
            fprintf('FAIL  %-20s  max_err=%.3e  at row %d  (ref=%.6e  cand=%.6e  tol=%.3e)\n', ...
                    ref_hdr{col}, max(err(violations)), worst_idx, r(worst_idx), c(worst_idx), tol(worst_idx));
            all_pass = false;
        else
            fprintf('PASS  %-20s  max_err=%.3e  digits=%.1f\n', ref_hdr{col}, max(err), s.digits);
        end

        col_stats(col).s = s;
    end

    % Generate analysis outputs if output_dir given
    if ~isempty(output_dir)
        if is_timeseries
            generate_timeseries_plots(ref_hdr, ref_data, cand_data, time_col, output_dir);
        else
            generate_algebraic_plot(ref_hdr, ref_data, cand_data, output_dir);
        end

        write_report(ref_hdr, col_stats, all_pass, atol, rtol, n_rows, is_timeseries, time_col, output_dir);
    end

    if all_pass
        exit(0);
    else
        exit(1);
    end
end

function generate_timeseries_plots(headers, ref_data, cand_data, time_col, output_dir)
    if time_col > 0
        x_axis = ref_data(:, time_col);
        x_label = headers{time_col};
    else
        x_axis = (1:size(ref_data, 1))';
        x_label = 'sample';
    end

    for col = 1:numel(headers)
        if col == time_col
            continue;
        end

        r = ref_data(:, col);
        c = cand_data(:, col);
        err = abs(r - c);

        f = figure('visible', 'off');

        subplot(2, 1, 1);
        plot(x_axis, r, 'b-', 'linewidth', 1.5);
        hold on;
        plot(x_axis, c, 'r--', 'linewidth', 1.5);
        hold off;
        ylabel(headers{col});
        legend('octave', 'ctrlpp', 'location', 'best');
        title(sprintf('%s: overlay', headers{col}));
        grid on;

        subplot(2, 1, 2);
        semilogy(x_axis, max(err, eps), 'k-', 'linewidth', 1.0);
        xlabel(x_label);
        ylabel('|error|');
        title(sprintf('%s: absolute error', headers{col}));
        grid on;

        png_path = fullfile(output_dir, [headers{col}, '.png']);
        print(f, png_path, '-dpng', '-r150');
        close(f);
    end
end

function generate_algebraic_plot(headers, ref_data, cand_data, output_dir)
    n = numel(headers);
    ref_vals = ref_data(1, :);
    cand_vals = cand_data(1, :);

    f = figure('visible', 'off');

    x = 1:n;
    bar_width = 0.35;
    bar(x - bar_width/2, ref_vals, bar_width, 'facecolor', [0.2 0.4 0.8]);
    hold on;
    bar(x + bar_width/2, cand_vals, bar_width, 'facecolor', [0.8 0.2 0.2]);
    hold off;

    set(gca, 'xtick', x, 'xticklabel', headers);
    legend('octave', 'ctrlpp', 'location', 'best');
    title('Algebraic comparison');
    ylabel('value');
    grid on;

    % Add error annotations
    for i = 1:n
        err = abs(ref_vals(i) - cand_vals(i));
        y_pos = max(abs(ref_vals(i)), abs(cand_vals(i)));
        text(i, y_pos * 1.05, sprintf('err=%.1e', err), 'horizontalalignment', 'center', 'fontsize', 7);
    end

    png_path = fullfile(output_dir, 'comparison.png');
    print(f, png_path, '-dpng', '-r150');
    close(f);
end

function write_report(headers, col_stats, all_pass, atol, rtol, n_rows, is_timeseries, time_col, output_dir)
    report_path = fullfile(output_dir, 'report.md');
    fid = fopen(report_path, 'w');

    % Extract case name from output_dir path
    [parent_dir, ~] = fileparts(output_dir);
    [~, case_name] = fileparts(parent_dir);

    if all_pass
        verdict = 'PASS';
    else
        verdict = 'FAIL';
    end

    fprintf(fid, '# Validation: %s\n\n', case_name);
    fprintf(fid, '- **Date**: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '- **Verdict**: **%s**\n', verdict);
    fprintf(fid, '- **Rows**: %d\n', n_rows);
    fprintf(fid, '- **Tolerances**: atol=%.0e, rtol=%.0e\n\n', atol, rtol);

    fprintf(fid, '## Per-signal statistics\n\n');
    fprintf(fid, '| Signal | Max abs | Mean abs | Std err | Max rel | Mean rel | Digits | Rel L2 |');
    if is_timeseries
        fprintf(fid, ' Growth | Class |');
    end
    fprintf(fid, ' Result |\n');

    fprintf(fid, '|--------|---------|----------|---------|---------|----------|--------|--------|');
    if is_timeseries
        fprintf(fid, '--------|-------|');
    end
    fprintf(fid, '--------|\n');

    for col = 1:numel(headers)
        s = col_stats(col).s;
        if s.pass
            result_str = 'PASS';
        else
            result_str = 'FAIL';
        end

        if isinf(s.digits)
            digits_str = 'exact';
        else
            digits_str = sprintf('%.1f', s.digits);
        end

        fprintf(fid, '| %s | %.2e | %.2e | %.2e | %.2e | %.2e | %s | %.2e |', ...
                s.name, s.max_abs_err, s.mean_abs_err, s.std_err, ...
                s.max_rel_err, s.mean_rel_err, digits_str, s.rel_l2);

        if is_timeseries
            if isnan(s.growth_rate)
                fprintf(fid, ' n/a | %s |', s.growth_class);
            else
                fprintf(fid, ' %.2e | %s |', s.growth_rate, s.growth_class);
            end
        end

        fprintf(fid, ' %s |\n', result_str);
    end

    fprintf(fid, '\n## Plots\n\n');

    if is_timeseries
        for col = 1:numel(headers)
            if col == time_col
                continue;
            end
            fprintf(fid, '### %s\n\n', headers{col});
            fprintf(fid, '![%s](%s.png)\n\n', headers{col}, headers{col});
        end
    else
        fprintf(fid, '![comparison](comparison.png)\n\n');
    end

    fclose(fid);
end

function [headers, data] = load_csv(filename)
    fid = fopen(filename, 'r');
    if fid == -1
        fprintf(2, 'ERROR: cannot open %s\n', filename);
        exit(2);
    end

    hdr_line = fgetl(fid);
    hdr_line = strtrim(hdr_line);
    if hdr_line(1) == '#'
        hdr_line = hdr_line(2:end);
    end
    headers = strtrim(strsplit(hdr_line, ','));

    data = [];
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line) || isempty(strtrim(line))
            continue;
        end
        vals = str2double(strsplit(strtrim(line), ','));
        data = [data; vals];
    end

    fclose(fid);
end

validate_compare();
