#include "ctrlpp/expected.h"

auto make_result() -> ctrlpp::expected<int, int>
{
    return 0;
}

int main()
{
    make_result();
    return 0;
}
