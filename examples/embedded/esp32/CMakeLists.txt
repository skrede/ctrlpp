# ctrlpp ESP32 on-target embeddability check: compiles the core ctrlpp control
# surface under the xtensa-esp32-elf toolchain and, when flashed, runs a live
# closed-loop LQR on the chip. `idf.py build` proves language-level
# embeddability; on a connected board `idf.py flash monitor` runs app_main, which
# designs an infinite-horizon LQR gain for a discrete double-integrator with
# float scalars and streams a tagged report over UART2 (the telemetry adapter):
# the reset that started the run, the closed-loop trajectory as binary32 bit
# patterns, and the estimation, signal-processing, trajectory and control
# families each judged against its host-double reference.
#
# Two images from the same sources, the stack probe's aside. ctrlpp_esp32 is
# the allocation-evidence image: it carries the allocation sensor and proves it
# can see an allocation, then runs the stack probe, which ends in a deliberate
# overflow and a restart.
# ctrlpp_esp32_timing is built without the sensor and without assertions, so the
# code it times carries none of their checks; it stops the task watchdog,
# proves the core's cycle counter live, times the families, and runs no stack
# probe, so it reports once per boot. Each prints its kind and its own
# fingerprint.
#
# To run:
#   1. Install ESP-IDF v5.3+ (https://docs.espressif.com/projects/esp-idf/).
#   2. `. $IDF_PATH/export.sh`   (e.g. `. /opt/esp-idf/export.sh`)
#   3. `idf.py -C examples/embedded/esp32 set-target esp32 && idf.py -C examples/embedded/esp32 build`
#   4. To validate on hardware: `idf.py -p <flash-port> flash monitor`.
#      The report is streamed on UART2 (ESP TX GPIO5 -> telemetry USB RX).
#   5. For the timing image, after step 4 has written the bootloader and the
#      partition table: `python -m esptool --chip esp32 -p <flash-port>
#      write_flash 0x10000 build/ctrlpp_esp32_timing.bin` (the factory app
#      offset in partitions.csv).
#
# Exercises ctrlpp control (lqr_gain, DARE-backed) with Scalar=float, exceptions
# and RTTI off, size-optimized -- the realistic embedded posture.

cmake_minimum_required(VERSION 3.16)

set(CTRLPP_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../..")

# Expose ctrlpp's public include tree to every component in this project as a
# system include, so the strict warning surface does not trip on Eigen or the
# ctrlpp template depth.
list(APPEND CMAKE_INCLUDE_PATH "${CTRLPP_ROOT}/lib/ctrlpp/include")

include($ENV{IDF_PATH}/tools/cmake/project.cmake)
project(ctrlpp_esp32)

include("${CMAKE_CURRENT_LIST_DIR}/timing_image.cmake")
ctrlpp_esp32_timing_image()
