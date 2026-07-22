# ctrlpp ESP32 on-target embeddability check: compiles the core ctrlpp control
# surface under the xtensa-esp32-elf toolchain and, when flashed, runs a live
# closed-loop LQR on the chip. `idf.py build` proves language-level
# embeddability; on a connected board `idf.py flash monitor` runs app_main, which
# designs an infinite-horizon LQR gain for a discrete double-integrator with
# float scalars, streams the closed-loop trajectory as CSV over UART2 (the
# telemetry adapter), and reports gain deviation against the host-double
# reference.
#
# To run:
#   1. Install ESP-IDF v5.3+ (https://docs.espressif.com/projects/esp-idf/).
#   2. `. $IDF_PATH/export.sh`   (e.g. `. /opt/esp-idf/export.sh`)
#   3. `idf.py -C examples/embedded/esp32 set-target esp32 && idf.py -C examples/embedded/esp32 build`
#   4. To validate on hardware: `idf.py -p <flash-port> flash monitor`.
#      The CSV trajectory is streamed on UART2 (ESP TX GPIO17 -> telemetry USB RX).
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
