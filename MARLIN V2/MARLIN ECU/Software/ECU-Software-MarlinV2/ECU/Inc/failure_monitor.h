#ifndef FAILURE_MONITOR_H
#define FAILURE_MONITOR_H

/* Safety evaluation interface; no limits are defined in this skeleton. */
#include "ecu_types.h"

/* True if the sample breaks a limit for this ECU state. */
bool FailureMonitor_Evaluate(const sensor_data_t *sensors,
                             ecu_state_t current_state);

#endif
