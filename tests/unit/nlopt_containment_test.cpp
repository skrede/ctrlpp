#include "ctrlpp/mpc/nlopt_solver.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>

#include <new>
#include <span>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <utility>
#include <optional>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

#if defined(__GLIBC__) && defined(__GLIBCXX__)
#include <array>
#include <cstdio>
#include <string>
#include <csignal>

#include <unistd.h>
#include <pthread.h>
#include <sys/wait.h>
#endif

namespace
{

std::atomic<bool> fail_next_allocation{false};

auto injected_failure() -> bool
{
    return fail_next_allocation.exchange(false, std::memory_order_relaxed);
}

auto try_allocate(std::size_t size) noexcept -> void*
{
    if(injected_failure())
        return nullptr;
    return std::malloc(size > 0 ? size : 1);
}

auto try_allocate(std::size_t size, std::align_val_t alignment) noexcept -> void*
{
    if(injected_failure())
        return nullptr;
    const auto align = static_cast<std::size_t>(alignment);
    const std::size_t padded = align * ((size + align - 1) / align);
    const std::size_t request = padded > 0 ? padded : align;
#if defined(_MSC_VER)
    return _aligned_malloc(request, align);
#else
    return std::aligned_alloc(align, request);
#endif
}

auto allocate(std::size_t size) -> void*
{
    if(void* pointer = try_allocate(size))
        return pointer;
    throw std::bad_alloc{};
}

auto allocate(std::size_t size, std::align_val_t alignment) -> void*
{
    if(void* pointer = try_allocate(size, alignment))
        return pointer;
    throw std::bad_alloc{};
}

void aligned_free(void* pointer) noexcept
{
#if defined(_MSC_VER)
    _aligned_free(pointer);
#else
    std::free(pointer);
#endif
}

}

// Every form is replaced because a sanitizer runtime supplies its own definition of any form left
// out, and its memory would then reach these std::free-based deletes as an allocator mismatch.
void* operator new(std::size_t size)
{
    return allocate(size);
}

void* operator new[](std::size_t size)
{
    return allocate(size);
}

void* operator new(std::size_t size, std::align_val_t alignment)
{
    return allocate(size, alignment);
}

void* operator new[](std::size_t size, std::align_val_t alignment)
{
    return allocate(size, alignment);
}

void* operator new(std::size_t size, const std::nothrow_t&) noexcept
{
    return try_allocate(size);
}

void* operator new[](std::size_t size, const std::nothrow_t&) noexcept
{
    return try_allocate(size);
}

void* operator new(std::size_t size, std::align_val_t alignment, const std::nothrow_t&) noexcept
{
    return try_allocate(size, alignment);
}

void* operator new[](std::size_t size, std::align_val_t alignment, const std::nothrow_t&) noexcept
{
    return try_allocate(size, alignment);
}

void operator delete(void* pointer) noexcept
{
    std::free(pointer);
}

void operator delete[](void* pointer) noexcept
{
    std::free(pointer);
}

void operator delete(void* pointer, std::size_t) noexcept
{
    std::free(pointer);
}

void operator delete[](void* pointer, std::size_t) noexcept
{
    std::free(pointer);
}

void operator delete(void* pointer, std::align_val_t) noexcept
{
    aligned_free(pointer);
}

void operator delete[](void* pointer, std::align_val_t) noexcept
{
    aligned_free(pointer);
}

void operator delete(void* pointer, std::size_t, std::align_val_t) noexcept
{
    aligned_free(pointer);
}

void operator delete[](void* pointer, std::size_t, std::align_val_t) noexcept
{
    aligned_free(pointer);
}

void operator delete(void* pointer, const std::nothrow_t&) noexcept
{
    std::free(pointer);
}

void operator delete[](void* pointer, const std::nothrow_t&) noexcept
{
    std::free(pointer);
}

void operator delete(void* pointer, std::align_val_t, const std::nothrow_t&) noexcept
{
    aligned_free(pointer);
}

void operator delete[](void* pointer, std::align_val_t, const std::nothrow_t&) noexcept
{
    aligned_free(pointer);
}

namespace
{

using NloptSolver = ctrlpp::nlopt_solver<double>;

struct call_counts
{
    std::int32_t constrained{0};
    std::int32_t plain{0};
};

auto bowl(std::span<const double> z) -> double
{
    return z[0] * z[0] + z[1] * z[1];
}

void bowl_gradient(std::span<const double> z, std::span<double> g)
{
    g[0] = 2.0 * z[0];
    g[1] = 2.0 * z[1];
}

void three_rows(std::span<const double> z, std::span<double> c)
{
    c[0] = z[0];
    c[1] = z[1];
    c[2] = z[0] + z[1];
}

auto three_row_problem() -> ctrlpp::nlp_problem<double>
{
    return {
        .n_vars = 2,
        .n_constraints = 3,
        .cost = bowl,
        .gradient = bowl_gradient,
        .constraints = three_rows,
        .constraint_jacobian = [](std::span<const double>, std::span<double>) {},
        .x_lower = {},
        .x_upper = {},
        .c_lower = Eigen::VectorXd::Constant(3, -1.0),
        .c_upper = Eigen::VectorXd::Constant(3, 1.0),
    };
}

auto three_row_problem(std::int32_t& calls) -> ctrlpp::nlp_problem<double>
{
    auto problem = three_row_problem();
    problem.constraints = [&calls](std::span<const double> z, std::span<double> c)
    {
        ++calls;
        three_rows(z, c);
    };
    return problem;
}

auto rowless_problem(std::int32_t& calls) -> ctrlpp::nlp_problem<double>
{
    return {
        .n_vars = 2,
        .n_constraints = 0,
        .cost = bowl,
        .gradient = bowl_gradient,
        .constraints = [&calls](std::span<const double>, std::span<double>) { ++calls; },
        .constraint_jacobian = [](std::span<const double>, std::span<double>) {},
        .x_lower = {},
        .x_upper = {},
        .c_lower = {},
        .c_upper = {},
    };
}

struct resetup_fixture
{
    resetup_fixture() : counts{}, constrained{three_row_problem(counts.constrained)}, plain{rowless_problem(counts.plain)} {}

    resetup_fixture(const resetup_fixture&) = delete;
    resetup_fixture& operator=(const resetup_fixture&) = delete;

    call_counts counts;
    ctrlpp::nlp_problem<double> constrained;
    ctrlpp::nlp_problem<double> plain;
};

auto interior_start() -> ctrlpp::nlp_update<double>
{
    return {.x0 = Eigen::Vector2d{0.5, 0.5}};
}

auto accepted(ctrlpp::solve_status status) -> bool
{
    return status == ctrlpp::solve_status::optimal || status == ctrlpp::solve_status::solved_inaccurate;
}

void require_rowless_solve(NloptSolver& solver, resetup_fixture& fixture)
{
    fixture.counts = {};
    const auto start = interior_start();
    ctrlpp::nlp_result<double> result{};
    REQUIRE_NOTHROW(result = solver.solve(start));
    CHECK(accepted(result.status));
    CHECK(result.objective < bowl({start.x0.data(), 2}));
    CHECK(fixture.counts.plain == 0);
    CHECK(fixture.counts.constrained == 0);
}

struct throw_at
{
    void tick()
    {
        if(++calls == at)
            throw at;
    }

    std::int32_t calls{0};
    std::int32_t at{0};
};

template <typename Trigger>
auto counted_cost(Trigger& trigger) -> ctrlpp::nlp_problem<double>
{
    auto problem = three_row_problem();
    problem.cost = [&trigger](std::span<const double> z) -> double
    {
        trigger.tick();
        return bowl(z);
    };
    return problem;
}

template <typename Trigger>
auto counted_rows(Trigger& trigger) -> ctrlpp::nlp_problem<double>
{
    auto problem = three_row_problem();
    problem.constraints = [&trigger](std::span<const double> z, std::span<double> c)
    {
        trigger.tick();
        three_rows(z, c);
    };
    return problem;
}

// SLSQP stops round-off limited on a rank-deficient equality Jacobian, here two identical rows.
template <typename Trigger>
auto duplicate_equality(Trigger& trigger) -> ctrlpp::nlp_problem<double>
{
    return {
        .n_vars = 2,
        .n_constraints = 2,
        .cost =
            [&trigger](std::span<const double> z) -> double
        {
            trigger.tick();
            return (z[0] - 1.0) * (z[0] - 1.0) + (z[1] - 1.0) * (z[1] - 1.0);
        },
        .gradient =
            [](std::span<const double> z, std::span<double> g)
        {
            g[0] = 2.0 * (z[0] - 1.0);
            g[1] = 2.0 * (z[1] - 1.0);
        },
        .constraints =
            [](std::span<const double> z, std::span<double> c)
        {
            c[0] = z[0] + z[1];
            c[1] = z[0] + z[1];
        },
        .constraint_jacobian = [](std::span<const double>, std::span<double>) {},
        .x_lower = {},
        .x_upper = {},
        .c_lower = Eigen::Vector2d{1.0, 1.0},
        .c_upper = Eigen::Vector2d{1.0, 1.0},
    };
}

auto solve_fresh(const ctrlpp::nlp_problem<double>& problem, const ctrlpp::nlp_update<double>& start) -> ctrlpp::nlp_result<double>
{
    NloptSolver solver;
    REQUIRE(solver.setup(problem).has_value());
    ctrlpp::nlp_result<double> result{};
    REQUIRE_NOTHROW(result = solver.solve(start));
    return result;
}

void require_replayed_throw_is_contained(throw_at& trigger, const ctrlpp::nlp_problem<double>& problem, const ctrlpp::nlp_update<double>& start)
{
    const std::int32_t recorded = trigger.calls;
    REQUIRE(recorded > 0);
    trigger = {.calls = 0, .at = recorded};
    const auto result = solve_fresh(problem, start);
    CHECK(trigger.calls == recorded);
    CHECK(result.status == ctrlpp::solve_status::error);
}

}

TEST_CASE("nlopt_solver re-set up on a problem without constraints evaluates no constraint", "[nlopt][resetup]")
{
    resetup_fixture fixture;
    NloptSolver solver;
    REQUIRE(solver.setup(fixture.constrained).has_value());
    REQUIRE(solver.setup(fixture.plain).has_value());
    require_rowless_solve(solver, fixture);
}

TEST_CASE("nlopt_solver copy-assigned from a solver without constraints evaluates no constraint", "[nlopt][resetup]")
{
    resetup_fixture fixture;
    NloptSolver target;
    NloptSolver source;
    REQUIRE(target.setup(fixture.constrained).has_value());
    REQUIRE(source.setup(fixture.plain).has_value());
    REQUIRE_NOTHROW(target = source);
    require_rowless_solve(target, fixture);
}

TEST_CASE("nlopt_solver move-assigned from a solver without constraints evaluates no constraint", "[nlopt][resetup]")
{
    resetup_fixture fixture;
    NloptSolver target;
    NloptSolver source;
    REQUIRE(target.setup(fixture.constrained).has_value());
    REQUIRE(source.setup(fixture.plain).has_value());
    REQUIRE_NOTHROW(target = std::move(source));
    require_rowless_solve(target, fixture);
}

TEST_CASE("nlopt_solver copy-constructed from a solver re-set up without constraints evaluates no constraint", "[nlopt][resetup]")
{
    resetup_fixture fixture;
    NloptSolver source;
    REQUIRE(source.setup(fixture.constrained).has_value());
    REQUIRE(source.setup(fixture.plain).has_value());
    std::optional<NloptSolver> copy;
    REQUIRE_NOTHROW(copy.emplace(source));
    require_rowless_solve(*copy, fixture);
}

TEST_CASE("nlopt_solver move-constructed from a solver re-set up without constraints evaluates no constraint", "[nlopt][resetup]")
{
    resetup_fixture fixture;
    NloptSolver source;
    REQUIRE(source.setup(fixture.constrained).has_value());
    REQUIRE(source.setup(fixture.plain).has_value());
    std::optional<NloptSolver> moved;
    REQUIRE_NOTHROW(moved.emplace(std::move(source)));
    require_rowless_solve(*moved, fixture);
}

TEST_CASE("nlopt_solver contains a non-standard exception from the violation evaluation after optimize", "[nlopt][containment]")
{
    throw_at trigger;
    const auto problem = counted_rows(trigger);
    const auto dry = solve_fresh(problem, interior_start());
    REQUIRE(dry.status == ctrlpp::solve_status::optimal);
    require_replayed_throw_is_contained(trigger, problem, interior_start());
}

TEST_CASE("nlopt_solver contains a non-standard exception from the round-off cost re-evaluation", "[nlopt][containment]")
{
    throw_at trigger;
    const auto problem = duplicate_equality(trigger);
    const ctrlpp::nlp_update<double> start{.x0 = Eigen::Vector2d::Zero()};
    const auto dry = solve_fresh(problem, start);
    REQUIRE(dry.status == ctrlpp::solve_status::solved_inaccurate);
    require_replayed_throw_is_contained(trigger, problem, start);
}

TEST_CASE("nlopt_solver contains a non-standard exception from the cost inside optimize", "[nlopt][containment]")
{
    throw_at trigger;
    const auto problem = counted_cost(trigger);
    const auto dry = solve_fresh(problem, interior_start());
    REQUIRE(dry.status == ctrlpp::solve_status::optimal);
    require_replayed_throw_is_contained(trigger, problem, interior_start());
}

TEST_CASE("nlopt_solver setup reports an allocation failure as allocation_failed and recovers", "[nlopt][containment]")
{
    resetup_fixture fixture;
    NloptSolver solver;

    std::optional<ctrlpp::expected<void, ctrlpp::nlopt_setup_error>> outcome;
    fail_next_allocation.store(true, std::memory_order_relaxed);
    REQUIRE_NOTHROW(outcome.emplace(solver.setup(fixture.constrained)));
    const bool injection_consumed = !fail_next_allocation.exchange(false, std::memory_order_relaxed);
    REQUIRE(injection_consumed);
    REQUIRE_FALSE(outcome->has_value());
    CHECK(outcome->error() == ctrlpp::nlopt_setup_error::allocation_failed);

    REQUIRE(solver.setup(fixture.constrained).has_value());
    ctrlpp::nlp_result<double> result{};
    REQUIRE_NOTHROW(result = solver.solve(interior_start()));
    CHECK(accepted(result.status));
}

#if defined(__GLIBC__) && defined(__GLIBCXX__)

namespace
{

struct cancel_at
{
    void tick()
    {
        if(++calls == at)
        {
            pthread_cancel(pthread_self());
            pthread_testcancel();
        }
    }

    std::int32_t calls{0};
    std::int32_t at{0};
};

struct cancellation_run
{
    const ctrlpp::nlp_problem<double>& problem;
    const ctrlpp::nlp_update<double>& start;
};

struct child_outcome
{
    int status;
    std::string diagnostics;
};

void* solve_on_thread(void* argument)
{
    auto const& run = *static_cast<const cancellation_run*>(argument);
    NloptSolver solver;
    if(solver.setup(run.problem).has_value())
        static_cast<void>(solver.solve(run.start));
    return nullptr;
}

// The test framework's abort handler is inherited across fork and is reset so the abort reaches the exit status.
// LIBC_FATAL_STDERR_ makes glibc write its fatal message to the redirected stderr rather than to the terminal.
void run_child(int diagnostics_pipe, cancellation_run run)
{
    std::signal(SIGABRT, SIG_DFL);
    dup2(diagnostics_pipe, STDERR_FILENO);
    setenv("LIBC_FATAL_STDERR_", "1", 1);
    pthread_t thread{};
    if(pthread_create(&thread, nullptr, solve_on_thread, &run) != 0)
        _exit(2);
    void* exit_value = nullptr;
    if(pthread_join(thread, &exit_value) != 0)
        _exit(3);
    _exit(exit_value == PTHREAD_CANCELED ? 0 : 1);
}

// The solve runs in a forked child that reports through its exit status, so a process abort is observed rather than suffered.
auto cancel_in_child(const ctrlpp::nlp_problem<double>& problem, const ctrlpp::nlp_update<double>& start) -> child_outcome
{
    std::array<int, 2> pipe_ends{};
    REQUIRE(pipe(pipe_ends.data()) == 0);
    const pid_t child = fork();
    REQUIRE(child >= 0);
    if(child == 0)
    {
        close(pipe_ends[0]);
        run_child(pipe_ends[1], {.problem = problem, .start = start});
    }
    close(pipe_ends[1]);
    child_outcome outcome{.status = 0, .diagnostics = {}};
    std::array<char, 256> buffer{};
    for(ssize_t n = 0; (n = read(pipe_ends[0], buffer.data(), buffer.size())) > 0;)
        outcome.diagnostics.append(buffer.data(), static_cast<std::size_t>(n));
    close(pipe_ends[0]);
    REQUIRE(waitpid(child, &outcome.status, 0) == child);
    std::fputs(outcome.diagnostics.c_str(), stderr);
    return outcome;
}

auto unwound(const child_outcome& outcome) -> bool
{
    return WIFEXITED(outcome.status) && WEXITSTATUS(outcome.status) == 0;
}

auto aborted_unrethrown(const child_outcome& outcome) -> bool
{
    return WIFSIGNALED(outcome.status) && WTERMSIG(outcome.status) == SIGABRT && outcome.diagnostics.find("FATAL: exception not rethrown") != std::string::npos;
}

void replay_last_call(cancel_at& trigger)
{
    const std::int32_t recorded = trigger.calls;
    REQUIRE(recorded > 0);
    trigger = {.calls = 0, .at = recorded};
}

}

TEST_CASE("nlopt_solver lets a thread cancelled in the round-off cost re-evaluation unwind", "[nlopt][cancellation]")
{
    cancel_at trigger;
    const auto problem = duplicate_equality(trigger);
    const ctrlpp::nlp_update<double> start{.x0 = Eigen::Vector2d::Zero()};
    REQUIRE(solve_fresh(problem, start).status == ctrlpp::solve_status::solved_inaccurate);
    replay_last_call(trigger);
    CHECK(unwound(cancel_in_child(problem, start)));
}

TEST_CASE("nlopt_solver lets a thread cancelled in the violation evaluation after optimize unwind", "[nlopt][cancellation]")
{
    cancel_at trigger;
    const auto problem = counted_rows(trigger);
    REQUIRE(solve_fresh(problem, interior_start()).status == ctrlpp::solve_status::optimal);
    replay_last_call(trigger);
    CHECK(unwound(cancel_in_child(problem, interior_start())));
}

// NLopt's own callback wrapper catches the unwind inside optimize and does not re-raise it, so glibc aborts.
TEST_CASE("A thread cancelled in a callback NLopt invokes during optimize still aborts the process", "[nlopt][cancellation]")
{
    cancel_at trigger{.calls = 0, .at = 1};
    const auto problem = counted_cost(trigger);
    CHECK(aborted_unrethrown(cancel_in_child(problem, interior_start())));
}

#endif
