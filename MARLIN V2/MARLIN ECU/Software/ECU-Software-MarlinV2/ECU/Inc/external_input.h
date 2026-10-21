#ifndef EXTERNAL_INPUT_H
#define EXTERNAL_INPUT_H

/*
 * Mission Control receive path.  Whatever receives frames (the future
 * Ethernet code) decodes them into ecu_inbound_frame_t and hands them to
 * AppTasks_SubmitInbound(), which queues them for the control task.
 */

#include "ecu_types.h"

/* Prepare the receive transport. */
void ExternalInput_Initialize(void);

/* Check the receive buffers for new frames and submit each one.  Polled from
 * the telemetry task; no receive interrupt. */
void ExternalInput_Poll(void);

#endif
