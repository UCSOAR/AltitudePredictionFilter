# Port the Altitude Prediction Filter (APF) into this board repo

You are porting the Altitude Prediction Filter from the DAQBoard repo into THIS repo (a new flight
board, STM32H7 @ 480 MHz, FreeRTOS, SoarOS). On the DAQBoard it was fully working on the bench: tare
on the pad, launch detection, UKF in flight, apogee detection/prediction, output published on the
DataBroker and logged to flash. The goal is the same behaviour here.

## Ground rules
- The user builds in STM32CubeIDE. Do NOT create CMake, Makefiles or compiler/toolchain files.
- Do not redesign the filter. Port it, then adapt only what this board needs (pins, sensor drivers,
  task names). Ask before changing filter maths or tuning.
- Keep line endings as each file already has them (several files are CRLF).
- At the end, list every file you changed or added.

## What you have
The user copies the whole `AltitudePredictionFilter/` folder into the root of this repo. It contains:
- `Input/ Process/ Estimate/ Prediction/ Output/` - the flight code (all `.cpp` + `Inc/*.hpp`)
- `Legacy/` - the old HALO code, reference only. **Exclude from the build.**
- `Tools/` - PC test tools (replay, Taber test, viewer). `Tools/replay/apf_replay.cpp` has its own
  `main()`. **Exclude from the build.**
- `Cases/ testSuite/ covarianceCalc/` - data, not code
- `Porting/` - this prompt plus reference patches of every change made OUTSIDE the filter folder
  on the DAQBoard (read them, do not `git apply` them blindly, paths and context differ here):
  - `daqboard_components.patch` - DAQBoard `Components/` + `parser.py` + `.gitignore`
  - `daqboard_cproject.patch` - DAQBoard `.cproject` (include paths / exclusions)
  - `soaros.patch` - SoarOS submodule
  - `peripheral_drivers.patch` - PeripheralDriversSubmodule

**IMPORTANT:** on the DAQBoard, `AltitudePredictionFilter/` was in `.gitignore`. Do NOT add that
ignore line here. The folder must be tracked by git in this repo (except `Tools/*/build/`,
`Tools/*/results/`, `Tools/**/__pycache__/` which are generated).

## Pipeline (for orientation)
SensorFrame (IMU0, IMU1, BARO0, BARO1 from the DataBroker) -> PreProcess (tare, calibration, axis
alignment) -> Dynamite (Madgwick per IMU + baro/IMU trust-weighted fusion) -> LaunchDetection ->
UKF (scenario model, KD-tree) -> ApogeeDetection + ApogeePredictor -> FilterOutput (`FilterData`,
published + logged). It all runs in one task: `Input/SensorFrameTask.cpp` (`SensorFrameTask`).
Scenario set / covariances / process model are chosen in `Input/Inc/FilterProfile.hpp`
(`APF_PROFILE`, default `APF_PROFILE_STORMFLY`; leave the default).

## Step 1 - Survey this repo first (no edits yet)
Report back:
1. Does it have SoarOS / PeripheralDriversSubmodule / SoarCommunications as submodules, and which
   branch/commit? Is SoarOS here the same as DAQBoard's (`Jad/DAQ`)?
2. Where are the equivalents of DAQBoard's `Components/main_system.cpp`, `DebugTask.cpp`,
   `SystemDefines.hpp`, `LoggingTask/LoggingTask.cpp`, `LoggingModule/LoggingService.cpp`,
   `LoggingModule/Inc/Log.hpp`, and the sensor polling task (`PollingTask`)?
3. Which IMUs and barometers this board has, their SPI buses / CS pins, and what types they publish
   (`IMUData`, `BaroData` in `SoarOS/Components/SystemTypes/SensorDataTypes.hpp`), with their ids.
4. Whether the DataBroker already has `FILTER_DATA` (`DataBrokerMessageTypes.hpp`) and the logger
   already has `LoggingData::FILTER` (`Log.hpp`). On DAQBoard both already existed.
5. Free heap / RAM after the current build (the KD-tree allocates on the heap, see Step 6).
Then propose the plan and wait for the user's OK.

## Step 2 - SoarOS changes (`Porting/soaros.patch`)
SoarOS is its own repo. Make these changes in the SoarOS submodule used here, then tell the user
they must be committed and pushed in the SoarOS repo and the submodule pointer committed here.
1. `Components/SystemTypes/SensorDataTypes.hpp` - replace the old `FilterData` (uint32 alt/velo/...)
   with the new one: `timestampMs, altitude, velocity, acceleration, altitudeStd, velocityStd,
   predictedApogee, timeToEvent, groundAltitude, apfState, flags` (floats except the first and
   the two uint8). The filter does not compile without it.
2. `Core/Inc/Queue.hpp`, `Core/Queue.cpp` - add `bool Send(Command&, bool reportFull, uint32_t
   waitTicks)`; the old `Send(cmd, reportFull)` calls it with `DEFAULT_QUEUE_SEND_WAIT_TICKS`.
   The failure print shows sender task, queue, used/depth, cmd, size, time, dropped count, and is
   rate limited (first drop, then every 100th). Needs `#include "task.h"`.
3. `Components/DataBroker/Inc/Publisher.hpp` - publish with `Send(brokerData, true, 0)` (never
   block). **This is the fix for a real lockup**: with the blocking send, a slow LoggingTask made
   PollingTask block 15 ms per message, the SensorFrame queue filled, and the whole board stalled.
   Do not skip it.

## Step 3 - Sensor drivers (`Porting/peripheral_drivers.patch`)
Only apply the parts that match this board's sensors; check each with the user.
- LSM6DSO32: accel +/-32 g (CTRL1_XL 0b01010100, 0.976 mg/LSB), gyro +/-2000 dps
  (CTRL2_G 0b01011100) stored as **0.1 dps** (`scaleToInt16(raw, 7, 10)`).
- LSM6DSO: accel +/-16 g, gyro +/-2000 dps, also stored as 0.1 dps.
  The filter's clip limits assume these ranges: `Process/Inc/Dynamite.hpp`
  `DYNAMITE_ACCEL_RANGE_G = {16, 32}`, `DYNAMITE_GYRO_RANGE_DPS = 2000`, gyro data in 0.1 dps.
  If this board's IMUs or ranges differ, update those constants and tell the user.
- MS5607 / MS5611: read the PROM once (`loadCalibration()`), pressure at OSR 2048 (0x46, osDelay 6),
  temperature at OSR 1024 (0x54, osDelay 4) refreshed every 10 samples.
- PollingTask: ground poll 20 ms (was 300) so the tare and launch detection run at flight rate.
  (Magnetometer read-before-trigger swap is unrelated to the filter, ask the user.)

## Step 4 - Board integration (`Porting/daqboard_components.patch`)
Make the equivalent changes in this repo's files:
1. `SystemDefines.hpp`:
   `TASK_SENSORFRAME_PRIORITY = 3`, `TASK_SENSORFRAME_DEPTH_OBJS = 16`,
   `TASK_SENSORFRAME_DEPTH_WORDS = 1024`, and `TASK_LOGGING_QUEUE_DEPTH_OBJS = 32` (was 10).
   The filter task must not sit below logging in priority.
2. `main_system.cpp`: `#include "SensorFrameTask.hpp"` and `SensorFrameTask::Inst().InitTask();`
   after the polling task.
3. `DebugTask.cpp`: route commands starting `APF_` to `SensorFrameTask::HandleApfCommand(msg)`,
   print "unknown APF command" if it returns false. (`APF_Help` lists them.)
4. `LoggingTask.cpp`: subscribe to `FilterData`; on `FILTER_DATA` write a 20 byte slot:
   `buf[0] = LoggingData::FILTER`, bytes 1-4 `timestampMs`, bytes 5-19 `FilterLog::Pack(data, buf+5)`
   (`Output/Inc/FilterLog.hpp`), priority SECOND.
5. `LoggingService.cpp` flash dump: parse the FILTER slot and print
   `FILTER Timestamp=.. AltCm=.. VelDms=.. AccCms2=.. SdCm=.. ApogeeM=.. EventS=.. State=.. Flags=..`
   (byte layout in the patch / `FilterLog.hpp`).
6. `parser.py` (if this repo has the flash-dump parser): add the `filter_pattern` and `filter_rows`.
7. `SensorFrameTask` subscribes to the IMU and baro data and maps them by `id` to IMU0/IMU1 and
   BARO0/BARO1. Check `Input/SensorFrameTask.cpp` against this board's ids and units; the IMU
   axis mapping is in `Input/AxisAlignment.cpp` (body +Z must point up the rocket) - confirm it
   with the user for this board's mounting.

## Step 5 - CubeIDE project (`Porting/daqboard_cproject.patch`)
Edit `.cproject` for every build configuration (Debug and Release, C and C++):
- Include paths (workspace-relative, NOT the absolute `C:\Users\...` path the DAQBoard used):
  `AltitudePredictionFilter/Input/Inc`, `Process/Inc`, `Estimate/Inc`, `Prediction/Inc`, `Output/Inc`.
- Source exclusions: `AltitudePredictionFilter/Legacy` and `AltitudePredictionFilter/Tools`.
Tell the user to refresh the project in CubeIDE afterwards. Do not touch anything else in `.cproject`.

## Step 6 - Things to check / known issues
- **KD-tree heap:** `Estimate/KDTree.cpp` (legacy port) builds `std::shared_ptr` nodes, one per
  scenario row, on the first match: thousands of nodes for the 9 Stormfly cases. On the DAQBoard
  only ~14 KB heap was free, which is NOT enough. Check `configTOTAL_HEAP_SIZE` / free heap here
  and report. If it does not fit, stop and ask (options: bigger heap, or a static-array KD-tree).
- **CPU:** on DAQBoard the predictor took ~200 ms per run and polling dropped to ~4 Hz in flight.
  At 480 MHz it should be better; measure the SensorFrameTask loop time and report it.
- `Estimate/KDTree.cpp` has one deliberate fix vs legacy (`KDNode()` initialised to infinity) -
  keep it, it prevents a garbage-index crash.
- Leave `APF_PROFILE` at Stormfly. The Taber profiles are for PC tests.

## Step 7 - Verify
1. PC check (no board needed): `python AltitudePredictionFilter/Tools/replay/run_tests.py` needs
   clang++. Expected on DAQBoard: 17/20 PASS; CHECK on `oldconfig_offnominal`, `oldconfig_coning`,
   `taber` (expected). The script's include list has `../SoarOS/Components/SystemTypes` - fix
   the path if SoarOS lives elsewhere here.
2. The user builds in CubeIDE. Fix compile errors without changing filter behaviour.
3. Bench test with the user, in order:
   - Boot, `APF_Help`, check the state prints and the tare completes on the pad.
   - No "Could not send data to queue" storm while logging (one rate-limited line is fine).
   - Lift the board ~1 m quickly: launch detection should fire (floor `LAUNCH_ALT_MIN_M` in
     `Input/Inc/LaunchDetection.hpp`).
   - Flash dump shows `FILTER ...` lines; `parser.py` turns them into a CSV.
4. Final message: every changed/added file (this repo and each submodule separately), what must be
   committed/pushed in which repo, and anything left for the user.
