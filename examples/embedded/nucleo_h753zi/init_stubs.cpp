#include "stack_watermark_board.h"

// __libc_init_array (called by the CMSIS Reset_Handler to run C++ static
// constructors) references _init/_fini, normally supplied by crti.o/crtn.o.
// Under -nostartfiles those CRT objects are not linked, so provide the stubs
// -- the .init_array/.fini_array sections carry the actual ctor/dtor lists and
// are walked by __libc_init_array itself. It calls _init before it walks the
// constructor list, which is why the stack is painted here.

extern "C" void _init(void)
{
    ctrlpp::paint_stack_at_startup();
}

extern "C" void _fini(void) {}

// A function-local static with a destructor registers that destructor against
// this symbol, normally defined by crtbegin.o. The image never returns from
// main, so the registration is never run; the definition only lets it link.
extern "C" {
void *__dso_handle = nullptr;
}
