/* Hand-written linker script for STM32H753ZI (NUCLEO-H753ZI), single-core
 * Cortex-M7. CMSIS's cmsis-device-h7 pack ships ST linker templates only for
 * the dual-core H745/747/755/757 parts, so a hand-written script is the only
 * option for the single-core H753.
 *
 * .data/.bss and the stack live in DTCM (zero-wait-state, CPU-local); code and
 * read-only data live in flash. The heap lives in AXI SRAM: the predictive
 * family's setup needs more of it than DTCM can hold beside its static instance
 * and the stack. With the D-cache off, as this leg runs, every access to a
 * heap-resident object pays the AXI SRAM's wait states.
 *
 * The symbol names below (_estack, _sdata/_edata/_sidata, _sbss/_ebss,
 * end/_end) are dictated verbatim by Source/Templates/gcc/startup_stm32h753xx.s
 * from the cmsis-device-h7 pack (end/_end are the heap base _sbrk consumes).
 * Renaming any of them silently breaks the startup copy/zero loops.
 */

MEMORY
{
  FLASH (rx)  : ORIGIN = 0x08000000, LENGTH = 2048K
  DTCM  (rwx) : ORIGIN = 0x20000000, LENGTH = 128K
  AXISRAM (rwx) : ORIGIN = 0x24000000, LENGTH = 512K
}

ENTRY(Reset_Handler)

/* Initial stack pointer = top of DTCM; the stack grows down from here. */
_estack = ORIGIN(DTCM) + LENGTH(DTCM);

/* The heap reserve is the high-water the board reported for this image with
 * the predictive family linked in, measured under a 256 KiB reserve: 76,720 B
 * after the family's setup, and the same at the end of the run, since the armed
 * steady state allocates nothing. No margin is added: the figure is exact and
 * repeats for a given image, and a requirement that grows is refused by name at
 * setup instead of running into anything. */
_Min_Heap_Size  = 76720;
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

  /* The linker's --build-id note is placed explicitly rather than left to orphan
   * placement: objcopy -O binary pads or reorders an orphan section, and the
   * vector table above must stay at 0x08000000. KEEP holds it against
   * --gc-sections; the firmware reads the descriptor back through the two
   * symbols below and prints it, which is what makes a stale resident image
   * detectable. */
  .note.gnu.build-id :
  {
    . = ALIGN(4);
    PROVIDE(__start_note_gnu_build_id = .);
    KEEP(*(.note.gnu.build-id))
    PROVIDE(__stop_note_gnu_build_id = .);
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

  /* The stack reserve makes the linker error out if DTCM cannot hold
   * data+bss+stack. */
  ._user_stack :
  {
    . = ALIGN(8);
    . = . + _Min_Stack_Size;
    . = ALIGN(8);
  } >DTCM

  /* end/_end mark the heap start (consumed by the _sbrk in sbrk_ceiling.cpp)
   * and _eheap its ceiling, so a break request past the reserve is refused by
   * _sbrk rather than handed memory nothing reserved. */
  ._user_heap (NOLOAD) :
  {
    . = ALIGN(8);
    PROVIDE(end = .);
    PROVIDE(_end = .);
    . = . + _Min_Heap_Size;
    PROVIDE(_eheap = .);
  } >AXISRAM

  /DISCARD/ : { libc.a(*) libm.a(*) libgcc.a(*) }

  .ARM.attributes 0 : { *(.ARM.attributes) }
}
