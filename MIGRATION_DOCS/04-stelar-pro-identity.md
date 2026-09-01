# STELAR-Pro repository identity

The semantic migration and the product rebrand are both complete. Current
implementation surfaces use one canonical spelling:

| Surface | Canonical name |
|---|---|
| Command | `./stelar-pro` |
| Java entry point | `stelarx.Main` |
| Java source root | `src/stelarx` |
| Java system properties | `stelarpro.*` |
| Environment variables | `STELAR_PRO_*` |
| Portable JAR | `stelar-pro.jar` |
| Native libraries | `libstelar_pro_weight`, `libstelar_pro_dp`, `libstelar_pro_dist`, `libstelar_pro_sim` |
| JNI namespace | `Java_stelarx_*` |
| Monitor command | `./run-stelar-pro-with-monitor.sh` |
| Simulation command | `./test-stelar-pro-simulated.sh` |
| Experiment outputs | `stelar-pro-outputs`, `out-stelar-pro.tre`, `stat-stelar-pro.csv` |
| Crash reports | `crash_logs/stelar-pro-crash-*`, `crash_logs/stelar-pro-hotspot-crash-*` |

The `stelarx` Java package and `Java_stelarx_*` JNI symbols remain internal ABI
identifiers so existing native entry points do not require an unrelated package
migration. Public commands, artifacts, messages, properties, and library base
names use STELAR-Pro.

`ASTRAL-X` still appears where it is historically meaningful: migration notes
describe the source implementation, ASTRAL-MP comparison material keeps the
name of that separate method, and `raw-prev` contains archived code. Those are
provenance labels, not active STELAR-Pro interfaces.

The checkout directory can be renamed or moved freely. Active scripts resolve
their resources relative to their own location, and explicit root/dataset
options accept absolute paths when an external location is needed.
