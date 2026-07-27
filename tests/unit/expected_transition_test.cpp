// Cross-state assignment of the owned expected has to reinitialize its union:
// one member is destroyed and the other constructed in the same storage. If the
// construction throws, the active-member flag must still describe the member the
// union really holds, otherwise the destructor destroys an object whose lifetime
// already ended. This translation unit drives that transition matrix
// exhaustively with instrumented members that throw on demand.
//
// It lives in the exceptions carve-out tree: throwing fixtures need exceptions,
// and the default -fno-exceptions tree does not add this target at all, so the
// fixtures need no preprocessor gates. The carve-out tree is the one built under
// the address and undefined-behavior sanitizers, which is what turns a surviving
// double destruction into a hard failure on top of the live-instance ledger
// asserted here.
#include "ctrlpp/expected.h"

#include <catch2/catch_test_macros.hpp>

#include <utility>
#include <type_traits>

namespace
{

// Live-instance ledger shared by every probe instantiation. A member destroyed
// twice drives this below its entry value, so the absence of a double
// destruction is asserted positively here rather than left to the sanitizer.
int live_probes = 0;

struct probe_failure
{
};

// Which constructor of a copy made FROM an instance throws. Armed per instance,
// so the destination of an assignment can stay quiet while the source trips.
enum class trip
{
    never,
    on_copy,
    on_move,
};

// Keeps the value-side and error-side probes distinct types, so an expected
// specialization can arm one side without arming the other.
enum class probe_role
{
    value,
    error,
};

// Instrumented member type. NothrowCopy and NothrowMove decide which
// reinitialization branch the enclosing expected selects for this type, and the
// per-instance trip decides whether that branch actually throws.
//
// The payload lives on the heap so that the two ways a defective
// reinitialization goes wrong are sanitizer-visible: reading a member whose
// lifetime already ended touches freed memory, and destroying it a second time
// frees the same block twice. The expected object itself is a stack local, and
// the sanitizer sees nothing at all in stack storage that is merely destroyed.
//
// The cell is an owning raw pointer rather than a unique_ptr on purpose. A
// unique_ptr destructor nulls its stored pointer after deleting, so a
// double-destroyed member would delete a null pointer and a read of it would
// return the moved-from sentinel: the fixture would hide the very memory error
// it exists to expose. The raw cell leaves the freed pointer in place, which is
// what a real double destruction looks like.
template <bool NothrowCopy, bool NothrowMove, probe_role Role>
class probe
{
public:
    explicit probe(int payload)
        : m_cell(new int(payload))
    {
        ++live_probes;
    }

    // Both constructors throw from the member-initializer list, before this
    // probe's own cell is initialized, so a throwing construction over a
    // destroyed member writes nothing into the union and the freed pointer stays
    // observable. Throwing from the constructor body would first store this
    // probe's own cell over it and mask the error.
    //
    // The allocation is the one remaining throw source on the nothrow
    // instantiations; exhausting the heap in a unit test terminates, which is an
    // acceptable outcome for a fixture.
    probe(const probe& other) noexcept(NothrowCopy)
        : m_trip(other.trip_after_copy_check())
        , m_cell(new int(other.payload()))
    {
        ++live_probes;
    }

    probe(probe&& other) noexcept(NothrowMove)
        : m_trip(other.trip_after_move_check())
        , m_cell(std::exchange(other.m_cell, nullptr))
    {
        ++live_probes;
    }

    // Assignment never throws: the same-state control cells assign rather than
    // reinitialize, and they must stay a clean baseline against which the
    // cross-state cells are read.
    probe& operator=(const probe& other)
    {
        if (this != &other)
        {
            int* replacement = new int(other.payload());
            delete m_cell;
            m_trip = other.m_trip;
            m_cell = replacement;
        }
        return *this;
    }

    probe& operator=(probe&& other) noexcept
    {
        if (this != &other)
        {
            delete m_cell;
            m_trip = other.m_trip;
            m_cell = std::exchange(other.m_cell, nullptr);
        }
        return *this;
    }

    ~probe()
    {
        delete m_cell;
        --live_probes;
    }

    int payload() const
    {
        return m_cell ? *m_cell : moved_from_payload;
    }

    void arm(trip which) noexcept
    {
        m_trip = which;
    }

private:
    trip trip_after_copy_check() const
    {
        if constexpr (!NothrowCopy)
        {
            if (m_trip == trip::on_copy)
                throw probe_failure{};
        }
        return m_trip;
    }

    trip trip_after_move_check() const
    {
        if constexpr (!NothrowMove)
        {
            if (m_trip == trip::on_move)
                throw probe_failure{};
        }
        return m_trip;
    }

    // Reported by a probe whose cell was moved out. No assertion below reads the
    // payload of a moved-from probe; the sentinel keeps the accessor total rather
    // than dereferencing a null cell.
    static constexpr int moved_from_payload = -1;

    trip m_trip = trip::never;
    int* m_cell  = nullptr;
};

using nothrow_value_probe = probe<true, true, probe_role::value>;
using nothrow_error_probe = probe<true, true, probe_role::error>;

using throwing_copy_value_probe = probe<false, true, probe_role::value>;
using throwing_copy_error_probe = probe<false, true, probe_role::error>;

using throwing_both_value_probe = probe<false, false, probe_role::value>;
using throwing_both_error_probe = probe<false, false, probe_role::error>;

// Which reinitialization branch a pair selects is a property of these traits, so
// pin them down here: the fixture families below are named after the branch these
// assertions prove they reach.
static_assert(std::is_nothrow_copy_constructible_v<nothrow_value_probe>);
static_assert(std::is_nothrow_move_constructible_v<nothrow_value_probe>);
static_assert(std::is_nothrow_copy_constructible_v<nothrow_error_probe>);
static_assert(std::is_nothrow_move_constructible_v<nothrow_error_probe>);

static_assert(!std::is_nothrow_copy_constructible_v<throwing_copy_value_probe>);
static_assert(std::is_nothrow_move_constructible_v<throwing_copy_value_probe>);
static_assert(!std::is_nothrow_copy_constructible_v<throwing_copy_error_probe>);
static_assert(std::is_nothrow_move_constructible_v<throwing_copy_error_probe>);

static_assert(!std::is_nothrow_copy_constructible_v<throwing_both_value_probe>);
static_assert(!std::is_nothrow_move_constructible_v<throwing_both_value_probe>);
static_assert(!std::is_nothrow_copy_constructible_v<throwing_both_error_probe>);
static_assert(!std::is_nothrow_move_constructible_v<throwing_both_error_probe>);

enum class held
{
    value,
    error,
};

enum class assign_kind
{
    copy,
    move,
};

constexpr int destination_payload = 11;
constexpr int source_payload      = 22;

const char* describe(held which)
{
    return which == held::value ? "value" : "error";
}

const char* describe(assign_kind which)
{
    return which == assign_kind::copy ? "copy-assign" : "move-assign";
}

// Every builder below returns a prvalue, so the expected is initialized in place
// and no probe is moved while building a fixture. An armed probe would otherwise
// trip during setup instead of during the assignment under test.
template <typename T, typename E>
ctrlpp::expected<T, E> make_in_state(held which, int payload)
{
    if (which == held::value)
        return ctrlpp::expected<T, E>{T{payload}};
    return ctrlpp::expected<T, E>{ctrlpp::unexpected(E{payload})};
}

template <typename E>
ctrlpp::expected<void, E> make_void_in_state(held which, int payload)
{
    if (which == held::value)
        return ctrlpp::expected<void, E>{};
    return ctrlpp::expected<void, E>{ctrlpp::unexpected(E{payload})};
}

template <typename T, typename E>
void arm_active_member(ctrlpp::expected<T, E>& subject, trip which)
{
    if (subject.has_value())
        (*subject).arm(which);
    else
        subject.error().arm(which);
}

template <typename T, typename E>
int active_payload(const ctrlpp::expected<T, E>& subject)
{
    return subject.has_value() ? (*subject).payload() : subject.error().payload();
}

// A cell in which nothing is armed: the transition must complete and land on the
// source's state and payload.
template <typename T, typename E>
void check_transition_completes(const char* label, held from, held to, assign_kind kind)
{
    INFO(label << ": " << describe(from) << " -> " << describe(to) << " by " << describe(kind));

    const int baseline = live_probes;
    {
        auto destination = make_in_state<T, E>(from, destination_payload);
        auto source      = make_in_state<T, E>(to, source_payload);

        if (kind == assign_kind::copy)
            destination = source;
        else
            destination = std::move(source);

        CHECK(destination.has_value() == (to == held::value));
        CHECK(active_payload(destination) == source_payload);
    }
    CHECK(live_probes == baseline);
}

// A cell in which the constructed member is armed. Four assertions: the
// exception propagates, the destination still reports its ORIGINAL state, its
// original payload survives, and no instance was created or destroyed on net.
template <typename T, typename E>
void check_transition_rolls_back(const char* label, held from, held to, assign_kind kind)
{
    INFO(label << ": throwing " << describe(to) << " construction during " << describe(kind) << " from a " << describe(from) << " state");

    const int baseline = live_probes;
    {
        auto destination = make_in_state<T, E>(from, destination_payload);
        auto source      = make_in_state<T, E>(to, source_payload);
        arm_active_member(source, kind == assign_kind::copy ? trip::on_copy : trip::on_move);

        const int before  = live_probes;
        bool propagated   = false;
        try
        {
            if (kind == assign_kind::copy)
                destination = source;
            else
                destination = std::move(source);
        }
        catch (const probe_failure&)
        {
            propagated = true;
        }

        CHECK(propagated);
        CHECK(destination.has_value() == (from == held::value));
        CHECK(active_payload(destination) == destination_payload);
        CHECK(live_probes == before);
    }
    CHECK(live_probes == baseline);
}

// Starting state x target state x assignment flavor, then the armed cells. An
// armed cell exists exactly where the constructing operation is potentially
// throwing for this pair, which is also the condition that selects the
// reinitialization branch, so the enumeration below adapts to the fixture family
// instead of being hand-picked per family.
template <typename T, typename E>
void drive_transition_matrix(const char* label)
{
    check_transition_completes<T, E>(label, held::value, held::value, assign_kind::copy);
    check_transition_completes<T, E>(label, held::value, held::value, assign_kind::move);
    check_transition_completes<T, E>(label, held::error, held::error, assign_kind::copy);
    check_transition_completes<T, E>(label, held::error, held::error, assign_kind::move);

    check_transition_completes<T, E>(label, held::error, held::value, assign_kind::copy);
    check_transition_completes<T, E>(label, held::error, held::value, assign_kind::move);
    check_transition_completes<T, E>(label, held::value, held::error, assign_kind::copy);
    check_transition_completes<T, E>(label, held::value, held::error, assign_kind::move);

    if constexpr (!std::is_nothrow_copy_constructible_v<T>)
        check_transition_rolls_back<T, E>(label, held::error, held::value, assign_kind::copy);
    if constexpr (!std::is_nothrow_move_constructible_v<T>)
        check_transition_rolls_back<T, E>(label, held::error, held::value, assign_kind::move);
    if constexpr (!std::is_nothrow_copy_constructible_v<E>)
        check_transition_rolls_back<T, E>(label, held::value, held::error, assign_kind::copy);
    if constexpr (!std::is_nothrow_move_constructible_v<E>)
        check_transition_rolls_back<T, E>(label, held::value, held::error, assign_kind::move);
}

template <typename E>
void check_void_transition_completes(const char* label, held from, held to, assign_kind kind)
{
    INFO(label << ": " << describe(from) << " -> " << describe(to) << " by " << describe(kind));

    const int baseline = live_probes;
    {
        auto destination = make_void_in_state<E>(from, destination_payload);
        auto source      = make_void_in_state<E>(to, source_payload);

        if (kind == assign_kind::copy)
            destination = source;
        else
            destination = std::move(source);

        CHECK(destination.has_value() == (to == held::value));
        if (!destination.has_value())
            CHECK(destination.error().payload() == source_payload);
    }
    CHECK(live_probes == baseline);
}

template <typename E>
void check_void_transition_rolls_back(const char* label, assign_kind kind)
{
    INFO(label << ": throwing error construction during " << describe(kind) << " from a value state");

    const int baseline = live_probes;
    {
        auto destination = make_void_in_state<E>(held::value, destination_payload);
        auto source      = make_void_in_state<E>(held::error, source_payload);
        source.error().arm(kind == assign_kind::copy ? trip::on_copy : trip::on_move);

        const int before = live_probes;
        bool propagated  = false;
        try
        {
            if (kind == assign_kind::copy)
                destination = source;
            else
                destination = std::move(source);
        }
        catch (const probe_failure&)
        {
            propagated = true;
        }

        CHECK(propagated);
        CHECK(destination.has_value());
        CHECK(live_probes == before);
    }
    CHECK(live_probes == baseline);
}

}

// Both members construct without throwing from either flavor's argument, so
// every reinitialization takes the destroy-then-construct branch. That branch has
// no throwing cell by construction; it is here so the branch is exercised and its
// ledger is proven balanced.
TEST_CASE("expected assignment transitions with nothrow members", "[expected][transition]")
{
    drive_transition_matrix<nothrow_value_probe, nothrow_error_probe>("nothrow pair");
}

// Both members throw on copy construction but move without throwing, so a
// cross-state copy-assignment takes the temporary-first branch and a cross-state
// move-assignment falls back to the destroy-then-construct branch. Covers the
// value type's copy constructor and the error type's copy constructor as
// throwing members.
TEST_CASE("expected assignment transitions with a throwing copy and a nothrow move", "[expected][transition]")
{
    drive_transition_matrix<throwing_copy_value_probe, throwing_copy_error_probe>("throwing-copy pair");
}

// The value type throws from both constructors, so building it can never be
// staged and the error-to-value transitions take the staged-rollback branch, with
// the error staged out of the union first. Covers the value type's copy
// constructor and its move constructor as throwing members.
TEST_CASE("expected assignment transitions with a value type that always throws", "[expected][transition]")
{
    drive_transition_matrix<throwing_both_value_probe, throwing_copy_error_probe>("throwing value, nothrow-movable error");
}

// Mirror image: the error type throws from both constructors, so the
// value-to-error transitions take the staged-rollback branch with the value
// staged out first. Covers the error type's copy constructor and its move
// constructor as throwing members.
TEST_CASE("expected assignment transitions with an error type that always throws", "[expected][transition]")
{
    drive_transition_matrix<throwing_copy_value_probe, throwing_both_error_probe>("nothrow-movable value, throwing error");
}

// The void specialization constructs its error before touching the active-member
// flag and so was already correct. It is asserted here as a no-regression
// control: these cases must pass both before and after the fix.
TEST_CASE("expected<void, E> assignment transitions are a no-regression control", "[expected][transition]")
{
    check_void_transition_completes<throwing_both_error_probe>("void control", held::value, held::value, assign_kind::copy);
    check_void_transition_completes<throwing_both_error_probe>("void control", held::value, held::value, assign_kind::move);
    check_void_transition_completes<throwing_copy_error_probe>("void control", held::error, held::error, assign_kind::copy);
    check_void_transition_completes<throwing_copy_error_probe>("void control", held::error, held::error, assign_kind::move);
    check_void_transition_completes<throwing_copy_error_probe>("void control", held::error, held::value, assign_kind::copy);
    check_void_transition_completes<throwing_copy_error_probe>("void control", held::error, held::value, assign_kind::move);
    check_void_transition_completes<throwing_copy_error_probe>("void control", held::value, held::error, assign_kind::copy);
    check_void_transition_completes<throwing_copy_error_probe>("void control", held::value, held::error, assign_kind::move);

    check_void_transition_rolls_back<throwing_both_error_probe>("void control", assign_kind::copy);
    check_void_transition_rolls_back<throwing_both_error_probe>("void control", assign_kind::move);
}
