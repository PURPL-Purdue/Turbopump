#include "failure_monitor.h"

/*
 * Keep sensor limits and fault policy in this module. The control task should
 * consume the result instead of knowing individual safety thresholds. The valid
 * redlines can differ by ECU state, so the current state is a required input.
 */
bool FailureMonitor_Evaluate(const sensor_data_t *sensors,
                             ecu_state_t current_state)
{
  (void)sensors;
  (void)current_state;
  /* TODO: redline checks per state.
   * TODO: noise / bad-reading handling - use sensors->valid[] and decide how
   * many bad or out-of-range readings (consecutive or within a window) count
   * as a FAULT, rather than tripping on a single noisy sample. */
  return false;
}
