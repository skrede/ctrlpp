// __libc_init_array (called by the CMSIS Reset_Handler to run C++ static
// constructors) references _init/_fini, normally supplied by crti.o/crtn.o.
// Under -nostartfiles those CRT objects are not linked, so provide empty stubs
// -- the .init_array/.fini_array sections carry the actual ctor/dtor lists and
// are walked by __libc_init_array itself.

extern "C" void _init(void) {}
extern "C" void _fini(void) {}
