# STELAR-X repository identity

The semantic migration and the product rebrand are both complete. Current
implementation surfaces use one canonical spelling:

| Surface | Canonical name |
|---|---|
| Command | `./stelarx` |
| Java entry point | `stelarx.Main` |
| Java source root | `src/stelarx` |
| Java system properties | `stelarx.*` |
| Environment variables | `STELARX_*` |
| Portable JAR | `stelarx.jar` |
| Native libraries | `libstelarx_weight`, `libstelarx_dp`, `libstelarx_dist`, `libstelarx_sim` |
| JNI namespace | `Java_stelarx_*` |
| Monitor command | `./run-stelarx-with-monitor.sh` |
| Simulation command | `./test-stelarx-simulated.sh` |
| Experiment outputs | `stelarx_outputs`, `out-stelarx.tre`, `stat-stelarx.csv` |
| Crash reports | `crash_logs/stelarx-crash-*`, `crash_logs/stelarx-hotspot-crash-*` |

`ASTRAL-X` still appears where it is historically meaningful: migration notes
describe the source implementation, ASTRAL-MP comparison material keeps the
name of that separate method, and `raw-prev` contains archived code. Those are
provenance labels, not active STELAR-X interfaces.

The checkout directory can be renamed or moved freely. Active scripts resolve
their resources relative to their own location, and explicit root/dataset
options accept absolute paths when an external location is needed.
