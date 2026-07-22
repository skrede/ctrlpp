/* Hand-written linker script for STM32H753ZI (NUCLEO-H753ZI), single-core
 * Cortex-M7. CMSIS's cmsis-device-h7 pack ships ST linker templates only for
 * the dual-core H745/747/755/757 parts, so a hand-written script is the only
 * option for the single-core H753.
 *
 * .data/.bss/stack and a small heap live in DTCM (zero-wait-state, CPU-local);
 * code and read-only data live in flash.
 *
 * The symbol names below (_estack, _sdata/_edata/_sidata, _sbss/_ebss,
 * end/_end) are dictated verbatim by Source/Templates/gcc/startup_stm32h753xx.s
 * from the cmsis-device-h7 pack (end/_end are the heap base newlib's _sbrk
 * consumes). Renaming any of them silently breaks the startup copy/zero loops.
 */

MEMORY
{
  FLASH (rx)  : ORIGIN = 0x08000000, LENGTH = 2048K
  DTCM  (rwx) : ORIGIN = 0x20000000, LENGTH = 128K
}

ENTRY(Reset_Handler)

/* Initial stack pointer = top of DTCM; the stack grows down from here. */
_estack = ORIGIN(DTCM) + LENGTH(DTCM);

/* This leg is a 2x2 Riccati design plus a K*x product: no solver-sized frame.
 * newlib-nano's %f printf path locks a small heap, so keep a modest reserve;
 * the stack reserve is generous but far below argmin's solver budget. */
_Min_Heap_Size  = 0x2000;   /* 8 KiB */
_Min_Stack_Size = 0x8000;   /* 32 KiB */

SECTIONS
{
  /* Vector table must be first in flash (CPU reads _estack + Reset_Handler
   * from the first two words at 0x08000000 at reset). */
  .isr_vector :
  {
    . = ALIGN(4);
    KEEP(*(.isr_vector))
    . = ALIGN(4);
  } >FLASH

  .text :
  {
    . = ALIGN(4);
    *(.text)
    *(.text*)
    *(.glue_7)
    *(.glue_7t)
    *(.eh_frame)
    KEEP(*(.init))
    KEEP(*(.fini))
    . = ALIGN(4);
    _etext = .;
  } >FLASH

  .rodata :
  {
    . = ALIGN(4);
    *(.rodata)
    *(.rodata*)
    . = ALIGN(4);
  } >FLASH

  .ARM.extab : { *(.ARM.extab* .gnu.linkonce.armextab.*) } >FLASH
  .ARM :
  {
    __exidx_start = .;
    *(.ARM.exidx*)
    __exidx_end = .;
  } >FLASH

  .preinit_array :
  {
    PROVIDE_HIDDEN(__preinit_array_start = .);
    KEEP(*(.preinit_array*))
    PROVIDE_HIDDEN(__preinit_array_end = .);
  } >FLASH
  .init_array :
  {
    PROVIDE_HIDDEN(__init_array_start = .);
    KEEP(*(SORT(.init_array.*)))
    KEEP(*(.init_array*))
    PROVIDE_HIDDEN(__init_array_end = .);
  } >FLASH
  .fini_array :
  {
    PROVIDE_HIDDEN(__fini_array_start = .);
    KEEP(*(SORT(.fini_array.*)))
    KEEP(*(.fini_array*))
    PROVIDE_HIDDEN(__fini_array_end = .);
  } >FLASH

  /* .data: loaded from flash (LMA), resident in DTCM (VMA); startup copies it
   * from _sidata to [_sdata,_edata). */
  _sidata = LOADADDR(.data);
  .data :
  {
    . = ALIGN(4);
    _sdata = .;
    *(.data)
    *(.data*)
    . = ALIGN(4);
    _edata = .;
  } >DTCM AT>FLASH

  /* .bss: zeroed by startup over [_sbss,_ebss). */
  .bss :
  {
    . = ALIGN(4);
    _sbss = .;
    __bss_start__ = _sbss;
    *(.bss)
    *(.bss*)
    *(COMMON)
    . = ALIGN(4);
    _ebss = .;
    __bss_end__ = _ebss;
  } >DTCM

  /* Heap + stack reservation. end/_end mark the heap start (consumed by
   * newlib's _sbrk); the stack reserve below it makes the linker error out if
   * DTCM cannot hold data+bss+heap+stack. */
  ._user_heap_stack :
  {
    . = ALIGN(8);
    PROVIDE(end = .);
    PROVIDE(_end = .);
    . = . + _Min_Heap_Size;
    . = . + _Min_Stack_Size;
    . = ALIGN(8);
  } >DTCM

  /DISCARD/ : { libc.a(*) libm.a(*) libgcc.a(*) }

  .ARM.attributes 0 : { *(.ARM.attributes) }
}
