#ifndef TELEMETRY_H
#define TELEMETRY_H

/* External telemetry interface; Ethernet serialization belongs in telemetry.c. */
#include "ecu_types.h"

/* Prepare the future telemetry transport. */
void Telemetry_Initialize(void);

/* Send one outbound frame produced by the control task. */
void Telemetry_Process(const ecu_outbound_frame_t *frame);

#endif
