#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{

struct benchmark_record
{
    std::string system;
    std::string solver;
    std::string algorithm;
    std::string warm_start;
    int nx{0};
    int horizon{0};
    double objective{0.0};
    double max_violation{0.0};
    double gradient_norm{0.0};
    bool success{false};
    int iterations{0};
    double solve_time_ms{0.0};
};

auto split_csv_line(const std::string& line) -> std::vector<std::string>
{
    std::vector<std::string> fields;
    std::istringstream stream(line);
    std::string field;
    while(std::getline(stream, field, ','))
        fields.push_back(field);
    return fields;
}

auto is_numeric_string(const std::string& s) -> bool
{
    if(s.empty())
        return false;
    for(char c : s)
    {
        if(c < '0' || c > '9')
            return false;
    }
    return true;
}

auto parse_csv(const std::string& path) -> std::vector<benchmark_record>
{
    std::vector<benchmark_record> records;
    std::ifstream file(path);
    if(!file.is_open())
    {
        std::cerr << "Warning: cannot open " << path << '\n';
        return records;
    }

    std::string header_line;
    if(!std::getline(file, header_line))
        return records;

    auto headers = split_csv_line(header_line);

    // Find column indices by name
    int col_system = -1;
    int col_solver = -1;
    int col_algorithm = -1;
    int col_warm_start = -1;
    int col_nx = -1;
    int col_horizon = -1;
    int col_objective = -1;
    int col_max_violation = -1;
    int col_gradient_norm = -1;
    int col_success = -1;
    int col_iterations = -1;
    int col_solve_time_ms = -1;

    for(int i = 0; i < static_cast<int>(headers.size()); ++i)
    {
        const auto& h = headers[static_cast<std::size_t>(i)];
        if(h == "system") col_system = i;
        else if(h == "solver") col_solver = i;
        else if(h == "algorithm") col_algorithm = i;
        else if(h == "warm_start") col_warm_start = i;
        else if(h == "nx") col_nx = i;
        else if(h == "horizon") col_horizon = i;
        else if(h == "objective") col_objective = i;
        else if(h == "max_violation") col_max_violation = i;
        else if(h == "gradient_norm") col_gradient_norm = i;
        else if(h == "success") col_success = i;
        else if(h == "iterations") col_iterations = i;
        else if(h == "solve_time_ms") col_solve_time_ms = i;
    }

    // Check that all required columns were found
    if(col_system < 0 || col_solver < 0 || col_algorithm < 0 ||
       col_nx < 0 || col_horizon < 0 || col_objective < 0 ||
       col_max_violation < 0 || col_success < 0 || col_solve_time_ms < 0)
    {
        std::cerr << "Warning: " << path << " missing required columns\n";
        return records;
    }

    // The step-budget CSV uses "budget" instead of "warm_start";
    // if warm_start column is missing, check for "budget" column
    int col_budget = -1;
    if(col_warm_start < 0)
    {
        for(int i = 0; i < static_cast<int>(headers.size()); ++i)
        {
            if(headers[static_cast<std::size_t>(i)] == "budget")
            {
                col_budget = i;
                break;
            }
        }
    }

    std::string line;
    while(std::getline(file, line))
    {
        if(line.empty())
            continue;

        auto fields = split_csv_line(line);
        auto n = static_cast<int>(fields.size());

        auto get = [&](int col) -> const std::string&
        {
            static const std::string empty;
            return (col >= 0 && col < n) ? fields[static_cast<std::size_t>(col)] : empty;
        };

        benchmark_record rec;
        rec.system = get(col_system);
        rec.solver = get(col_solver);
        rec.algorithm = get(col_algorithm);

        if(col_warm_start >= 0)
            rec.warm_start = get(col_warm_start);
        else if(col_budget >= 0)
            rec.warm_start = get(col_budget);
        else
            rec.warm_start = "unknown";

        try
        {
            rec.nx = std::stoi(get(col_nx));
            rec.horizon = std::stoi(get(col_horizon));
            rec.objective = std::stod(get(col_objective));
            rec.max_violation = std::stod(get(col_max_violation));
            rec.solve_time_ms = std::stod(get(col_solve_time_ms));

            if(col_gradient_norm >= 0)
                rec.gradient_norm = std::stod(get(col_gradient_norm));
            if(col_success >= 0)
                rec.success = (get(col_success) == "1");
            if(col_iterations >= 0)
                rec.iterations = std::stoi(get(col_iterations));
        }
        catch(...)
        {
            continue;
        }

        records.push_back(std::move(rec));
    }

    return records;
}

// Problem identity: (system, nx, horizon, warm_start)
struct problem_key
{
    std::string system;
    int nx;
    int horizon;
    std::string warm_start;

    auto operator<=>(const problem_key&) const = default;
};

// Solver identity: (solver, algorithm)
struct solver_key
{
    std::string solver;
    std::string algorithm;

    auto operator<=>(const solver_key&) const = default;

    auto label() const -> std::string
    {
        return solver + "_" + algorithm;
    }
};

// Metric extraction functions
using metric_fn = double (*)(const benchmark_record&);

auto metric_solve_time(const benchmark_record& r) -> double
{
    return r.solve_time_ms;
}

// Grouped data: problem -> solver -> record
using grouped_data = std::map<problem_key, std::map<solver_key, benchmark_record>>;

auto group_records(const std::vector<benchmark_record>& records) -> grouped_data
{
    grouped_data grouped;
    for(const auto& rec : records)
    {
        problem_key pk{rec.system, rec.nx, rec.horizon, rec.warm_start};
        solver_key sk{rec.solver, rec.algorithm};
        grouped[pk][sk] = rec;
    }
    return grouped;
}

auto collect_solvers(const grouped_data& data) -> std::vector<solver_key>
{
    std::set<solver_key> solvers;
    for(const auto& [pk, solver_map] : data)
    {
        for(const auto& [sk, rec] : solver_map)
            solvers.insert(sk);
    }
    return {solvers.begin(), solvers.end()};
}

auto log_space(double lo, double hi, int n) -> std::vector<double>
{
    std::vector<double> result;
    result.reserve(static_cast<std::size_t>(n));
    double log_lo = std::log10(lo);
    double log_hi = std::log10(hi);
    for(int i = 0; i < n; ++i)
    {
        double t = (n > 1) ? static_cast<double>(i) / (n - 1) : 0.0;
        result.push_back(std::pow(10.0, log_lo + t * (log_hi - log_lo)));
    }
    return result;
}

constexpr double infinity_ratio = 1e20;
constexpr double epsilon = 1e-12;

enum class profile_type
{
    time,
    objective,
    violation
};

auto compute_ratios(const grouped_data& data,
                    const std::vector<solver_key>& solvers,
                    profile_type type)
    -> std::map<solver_key, std::vector<double>>
{
    std::map<solver_key, std::vector<double>> ratios;
    for(const auto& sk : solvers)
        ratios[sk] = {};

    for(const auto& [pk, solver_map] : data)
    {
        // Find the best value for this problem
        double best = std::numeric_limits<double>::max();

        for(const auto& [sk, rec] : solver_map)
        {
            if(!rec.success)
                continue;

            double val = 0.0;
            switch(type)
            {
            case profile_type::time:
                val = rec.solve_time_ms + epsilon;
                break;
            case profile_type::objective:
                val = std::abs(rec.objective) + epsilon;
                break;
            case profile_type::violation:
                val = rec.max_violation + epsilon;
                break;
            }
            best = std::min(best, val);
        }

        // For objective accuracy: compute cost as |obj - best_obj| + epsilon
        // We need the best objective separately
        double best_obj = std::numeric_limits<double>::max();
        if(type == profile_type::objective)
        {
            for(const auto& [sk, rec] : solver_map)
            {
                if(rec.success)
                    best_obj = std::min(best_obj, rec.objective);
            }
        }

        for(const auto& sk : solvers)
        {
            auto it = solver_map.find(sk);
            if(it == solver_map.end() || !it->second.success)
            {
                ratios[sk].push_back(infinity_ratio);
                continue;
            }

            const auto& rec = it->second;
            double val = 0.0;
            switch(type)
            {
            case profile_type::time:
                val = rec.solve_time_ms + epsilon;
                break;
            case profile_type::objective:
                val = std::abs(rec.objective - best_obj) + epsilon;
                best = epsilon;
                break;
            case profile_type::violation:
                val = rec.max_violation + epsilon;
                break;
            }

            double ratio = val / best;
            ratios[sk].push_back(ratio);
        }
    }

    return ratios;
}

void write_profile_csv(const std::string& filename,
                       const std::vector<solver_key>& solvers,
                       const std::map<solver_key, std::vector<double>>& ratios)
{
    auto tau_values = log_space(1.0, 100.0, 50);
    auto n_problems = ratios.begin()->second.size();

    std::ofstream out(filename);
    if(!out.is_open())
    {
        std::cerr << "Error: cannot write " << filename << '\n';
        return;
    }

    // Header
    out << "tau";
    for(const auto& sk : solvers)
        out << ',' << sk.label();
    out << '\n';

    // Rows
    for(double tau : tau_values)
    {
        out << tau;
        for(const auto& sk : solvers)
        {
            const auto& solver_ratios = ratios.at(sk);
            int count = 0;
            for(double r : solver_ratios)
            {
                if(r <= tau)
                    ++count;
            }
            double rho = static_cast<double>(count) / static_cast<double>(n_problems);
            out << ',' << rho;
        }
        out << '\n';
    }

    std::cout << "  Wrote " << filename << " (" << n_problems << " problems, "
              << solvers.size() << " solvers)\n";
}

auto find_quality_csvs() -> std::vector<std::string>
{
    std::vector<std::string> paths;
    for(const auto& entry : std::filesystem::directory_iterator("."))
    {
        auto name = entry.path().filename().string();
        if(name.starts_with("bench_") && name.ends_with("_quality.csv"))
            paths.push_back(entry.path().string());
    }
    std::ranges::sort(paths);
    return paths;
}

}

int main(int argc, char* argv[])
{
    std::vector<std::string> csv_paths;

    if(argc > 1)
    {
        for(int i = 1; i < argc; ++i)
            csv_paths.emplace_back(argv[i]);
    }
    else
    {
        csv_paths = find_quality_csvs();
    }

    if(csv_paths.empty())
    {
        std::cerr << "No quality CSV files found. Pass paths as arguments or run "
                     "from a directory containing bench_*_quality.csv files.\n";
        return 1;
    }

    // Parse all CSVs
    std::vector<benchmark_record> all_records;
    for(const auto& path : csv_paths)
    {
        auto records = parse_csv(path);
        std::cout << "Read " << records.size() << " records from " << path << '\n';
        all_records.insert(all_records.end(), records.begin(), records.end());
    }

    // Filter: skip rows where warm_start is a pure numeric string (budget points)
    std::erase_if(all_records, [](const benchmark_record& r)
    {
        return is_numeric_string(r.warm_start);
    });

    std::cout << "Total records after filtering: " << all_records.size() << '\n';

    if(all_records.empty())
    {
        std::cerr << "No valid records to process.\n";
        return 1;
    }

    // Group and compute profiles
    auto grouped = group_records(all_records);
    auto solvers = collect_solvers(grouped);

    std::cout << "Problems: " << grouped.size() << ", Solvers: " << solvers.size() << '\n';

    // Compute and write three performance profiles
    auto time_ratios = compute_ratios(grouped, solvers, profile_type::time);
    write_profile_csv("profile_time.csv", solvers, time_ratios);

    auto obj_ratios = compute_ratios(grouped, solvers, profile_type::objective);
    write_profile_csv("profile_objective.csv", solvers, obj_ratios);

    auto viol_ratios = compute_ratios(grouped, solvers, profile_type::violation);
    write_profile_csv("profile_violation.csv", solvers, viol_ratios);

    std::cout << "Done. Generated 3 performance profile CSVs.\n";
    return 0;
}
