#include "failure_monitor.h"

/*
 * Keep sensor limits and fault policy in this module. The control task should
 * consume the result instead of knowing individual safety thresholds. The valid
 * redlines can differ by ECU state, so the current state is a required input.
 */
fault_flags_t FailureMonitor_Evaluate(const sensor_data_t *sensors,
                                      ecu_state_t current_state)
{
  (void)sensors;
  (void)current_state;
  return FAULT_NONE;
}
