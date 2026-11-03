"""ctrlpp's data for the style lint: scopes, guard prefixes, include roots, carve-outs and word lists.

Every carve-out carries its reason as data, so a reader can check it against the code it exempts.
"""

CPP_SCOPE_PREFIXES = ("lib/", "examples/")
CPP_SUFFIXES = (".h", ".hpp", ".cpp", ".cc", ".inl", ".ipp")
HEADER_SUFFIXES = (".h", ".hpp")
SPELLING_NAMES = frozenset({"CMakeLists.txt"})
SPELLING_SUFFIXES = (".cmake", ".cmake.in")

ALLOWED_ATTRIBUTES = frozenset({"noreturn", "fallthrough", "maybe_unused"})

# (path prefix, attribute, reason)
ATTRIBUTE_CARVE_OUTS = (
    ("examples/embedded/", "gnu::noinline",
     "the on-target stack and allocation probes measure a call's own frame, which inlining folds into "
     "the caller's and so hides"),
)

# (path prefix, guard prefix): the guard is the prefix plus the path under the path prefix, uppercased,
# with every character other than a letter or digit written as '_'.
GUARD_PREFIXES = (
    ("lib/ctrlpp/include/", "HPP_GUARD_"),
    ("examples/", "HPP_GUARD_CTRLPP_EXAMPLES_"),
)

# Directories a quoted include resolves against after the includer's own directory, as the builds pass them.
INCLUDE_ROOTS = ("lib/ctrlpp/include", "examples/embedded/shared")

# (path glob, earlier include, later include glob, reason): the later include may sit wherever the earlier
# one's required position puts it, and no include-order finding is raised for the pair.
ORDER_DEPENDENCIES = (
    ("examples/embedded/*", "eigen_alloc_sentinel.h", "*",
     "the sentinel's runtime-no-malloc switch and eigen_assert override take effect only if defined before "
     "Eigen is parsed, and the header #errors when Eigen was parsed first"),
    ("examples/embedded/esp32/*", "freertos/FreeRTOS.h", "freertos/*",
     "ESP-IDF's FreeRTOS headers #error unless FreeRTOS.h is included before them"),
)

STANDARD_HEADERS = frozenset("""
    algorithm any array atomic barrier bit bitset charconv chrono codecvt compare complex concepts
    condition_variable coroutine deque exception execution filesystem format forward_list fstream functional
    future initializer_list iomanip ios iosfwd iostream istream iterator latch limits list locale map
    memory memory_resource mutex new numbers numeric optional ostream queue random ranges ratio regex
    scoped_allocator semaphore set shared_mutex source_location span sstream stack stdexcept stop_token
    streambuf string string_view strstream syncstream system_error thread tuple type_traits typeindex
    typeinfo unordered_map unordered_set utility valarray variant vector version
    expected flat_map flat_set generator mdspan print spanstream stacktrace stdfloat
    cassert cctype cerrno cfenv cfloat cinttypes climits clocale cmath csetjmp csignal cstdarg cstddef
    cstdint cstdio cstdlib cstring ctime cuchar cwchar cwctype
    assert.h ctype.h errno.h fenv.h float.h inttypes.h limits.h locale.h math.h setjmp.h signal.h stdarg.h
    stddef.h stdint.h stdio.h stdlib.h string.h time.h uchar.h wchar.h wctype.h
""".split())

ISE_STEMS = (
    "optim", "normal", "initial", "stabil", "standard", "canonical", "central", "organ", "discret", "minim",
    "maxim", "real", "general", "parametr", "parameter", "character", "linear", "regular", "factor", "vector",
    "symmetr", "local", "final", "util", "priorit", "categor", "custom", "summar", "recogn", "visual",
    "serial", "synchron", "special", "orthogonal", "orthonormal", "diagonal", "equal", "penal", "random",
    "sanit", "token", "amort", "author", "critic", "emphas", "quant", "digit", "global", "ideal", "nominal",
    "rational", "decimal", "capital", "modern", "neutral", "polar", "harmon", "memor", "apolog",
)

# (pattern, American replacement), matched case-insensitively against the word the pattern ends on.
BRITISH_PATTERNS = (
    (r"(" + "|".join(ISE_STEMS) + r")is(e|ed|es|ing|ation|ations|er|ers|able)\b", r"\1iz\2"),
    (r"\b(behavi|col|fav|hon|neighb|flav|lab|harb|hum|rum|vap)ou(r\w*)\b", r"\1o\2"),
    (r"\b(cent|met|lit|fib)re(s?)\b", r"\1er\2"),
    (r"\b(cent)r(ed|ing)\b", r"\1er\2"),
    (r"\b(anal|paral|catal)ys(e|ed|ing|er|ers)\b", r"\1yz\2"),
    (r"\b(model)l(er|ers)\b", r"\1\2"),
    (r"\b(model|label|cancel|signal|travel|level|channel|fuel|tunnel|dial|total|equal)l(ed|ing)\b", r"\1\2"),
    (r"\b(judg)ement(s?)\b", r"\1ment\2"),
    (r"\bwhilst\b", r"while"),
    (r"\bamongst\b", r"among"),
    (r"\bprogramme(s?)\b", r"program\1"),
    (r"\bartefact(s?)\b", r"artifact\1"),
    (r"\blicence(s|d)?\b", r"license\1"),
    (r"\b(def|off|pret)ence(s?)\b", r"\1ense\2"),
    (r"\bgrey(s|ed|ish|scale)?\b", r"gray\1"),
    (r"\baluminium\b", r"aluminum"),
    (r"\bscepti(c|cs|cal|cism)\b", r"skepti\1"),
    (r"\bfulfil(s|ment)?\b", r"fulfill\1"),
)

# (phrase, reason): exact text the spelling rule skips, such as the title of a British publication.
ALLOWED_PHRASES = ()

# Matched from the start of a data member's type, after any leading static, inline, constexpr, mutable or const; a
# function-pointer member's type is its whole declarator.
CALLBACK_TYPE_PATTERNS = (
    r"std::function\s*<",
    r"std::move_only_function\s*<",
    r"std::optional\s*<\s*std::(?:move_only_)?function\s*<",
    r"[^()]*\(\s*\*\s*\w+\s*\)\s*\(",
)

CARVE_OUT_SUFFIXES = ("_config", "_options", "_settings")
CARVE_OUT_NESTED = "config"

# (path, function, reason): a deduced or trailing return type that cannot be written in front.
UNSPELLABLE_RETURNS = ()

FALLIBLE_CALL_NAMES = (
    "create", "setup", "place_observer", "place", "lqr_gain_continuous", "lqr_gain", "partition_lqi_gain", "lqi_gain",
    "validate_biquad_design", "can_rescale_to", "rescale_to", "dare", "care", "set_params", "set_integral",
    "lqr_tv_gains", "lqr_cost", "validate_resolution_factor", "check_coefficient_resolution", "prewarped_frequency",
    "validated_prewarp", "resolved_coefficients", "chebyshev1_section", "moesp", "moesp_singular_values", "tf2ss",
    "ss2tf", "try_create", "make_particle_filter", "set_noise", "reset_covariance",
)

# Each name is fallible on some surfaces and infallible on others (the online planners' update, the constant LQR
# gain's compute, the trajectories' evaluate, the solver adapters' solve, Eigen's in-place normalize), and a call
# site's name cannot tell them apart, so a discarded result of one of them goes unflagged.
EXCLUDED_CALL_NAMES = frozenset({"update", "compute", "evaluate", "solve", "normalize"})
