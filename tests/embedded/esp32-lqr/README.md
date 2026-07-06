# ctrlpp ESP32 LQR embeddability check

An on-target embeddability check that runs the core ctrlpp control surface on a
physical ESP32 (Xtensa LX6). It designs an infinite-horizon LQR gain for a
discrete double-integrator with `float` scalars on the chip, runs the closed
loop from a unit position error, and reports the result three ways:

- the on-device float gain and its deviation from an independent host-double
  reference (the "never trust the self-report" cross-check);
- per-step feedback-law timing in microseconds and CPU cycles;
- the full closed-loop trajectory streamed as CSV over UART2.

It exercises `ctrlpp::lqr_gain` (DARE-backed) with C++ exceptions and RTTI off
and size optimization on, i.e. a realistic embedded posture. The `idf.py build`
leg is a standalone language-level embeddability check; the flash/monitor legs
add on-hardware validation when a board is connected.

## Wiring

Two USB-serial links to the board:

| Host port      | Adapter            | Role                                          |
| -------------- | ------------------ | --------------------------------------------- |
| `/dev/ttyUSB0` | CP2102 (onboard)   | Flash + console log (`ESP_LOGI` diagnostics)  |
| `/dev/ttyUSB1` | FT232R (external)  | Telemetry: the CSV trajectory stream          |

The telemetry adapter is cross-wired to the ESP's UART2: **ESP TX GPIO17 -> USB
RX**, **ESP RX GPIO16 <- USB TX**, at 115200 8N1. The console log uses the
onboard CP2102 (UART0), also 115200.

## Build

```sh
. /opt/esp-idf/export.sh                 # ESP-IDF v5.1+ (tested on v6.0.1)
export CTRLPP_EIGEN_DIR=../../../build-tests/_deps/eigen3-src   # any Eigen >= 3.4
idf.py set-target esp32
idf.py build
```

Eigen is header-only and sits outside the ESP-IDF component graph; point
`CTRLPP_EIGEN_DIR` at any Eigen checkout (the ctrlpp host test build already
fetches one under `build-tests/_deps/eigen3-src`).

## Run on hardware

The demo runs once per boot: `app_main` prints the diagnostics and streams the
201-sample CSV, then returns. Press the board's EN/reset button to replay.

Flash and watch the console diagnostics (gain, deviation, timing, settling):

```sh
idf.py -p /dev/ttyUSB0 flash monitor
```

Capture the CSV trajectory from the telemetry adapter in a second terminal
(start it before resetting to catch the one-shot stream):

```sh
idf.py -p /dev/ttyUSB1 monitor -b 115200
#   header: k,t_s,x0,x1,u
```

## Reference

The host-double reference constants baked into `main.cpp` come from the
host-double build of the identical problem:

```
K = [9.467093672, 5.281739638]
final |x| = 2.02e-05   (from x0 = [1, 0], 200 steps at 50 Hz)
```

The on-device float gain should match to within float precision; the closed loop
should settle to a comparable final norm.

## Measured results (ESP32, ESP-IDF v6.0.1, -Os, exceptions/RTTI off)

First on-hardware run, captured over the console (UART0) and telemetry (UART2):

```
gain    K = [9.4692507, 5.2824945]   (host = [9.4670937, 5.2817396])
float-vs-host gain dev = [2.157e-03, 7.549e-04]
feedback law: 0.4415 us/eval, 71 cycles/eval (avg of 2000)
closed loop: 200 steps @ 50 Hz, peak |u| = -9.4693
settling: final |x| = 2.018e-05 (host = 2.021e-05, ratio = 0.999)
```

Full 201-sample CSV vs the host-double trajectory: max |dx0| = 6.0e-5,
max |dx1| = 1.75e-4, max |du| = 2.2e-3 -- all consistent with a float Riccati
solve, dominated by the gain's 3rd-4th significant-figure difference.

Footprint (`idf.py size`): app image ~390 KB (62% of the partition free), static
DRAM 12.8 KB (7%), IRAM 41 KB (31%), flash code 133 KB + data 193 KB. The float
LQR/DARE path and Eigen link and fit with generous headroom.
