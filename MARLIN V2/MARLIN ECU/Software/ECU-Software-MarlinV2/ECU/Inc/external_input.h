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

#endif
