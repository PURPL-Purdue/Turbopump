#ifndef FAILURE_MONITOR_H
#define FAILURE_MONITOR_H

/* Safety evaluation interface; no limits are defined in this skeleton. */
#include "ecu_types.h"

/* Evaluate a sample and return the sensor faults active for this ECU state. */
fault_flags_t FailureMonitor_Evaluate(const sensor_data_t *sensors,
                                      ecu_state_t current_state);

#endif
