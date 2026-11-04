#include "stack_watermark_board.h"

// Under -nostartfiles crti.o and crtn.o are not linked, but __libc_init_array, which the CMSIS
// Reset_Handler calls to run the static constructors, references _init and _fini. It calls _init before it
// walks the constructor list, which is why the stack is painted here.

extern "C" void _init(void)
{
    ctrlpp::paint_stack_at_startup();
}

extern "C" void _fini(void)
{
}

// A function-local static with a destructor registers that destructor against
// this symbol, normally defined by crtbegin.o. The image never returns from
// main, so the registration is never run; the definition only lets it link.
extern "C" {
void* __dso_handle = nullptr;
}
