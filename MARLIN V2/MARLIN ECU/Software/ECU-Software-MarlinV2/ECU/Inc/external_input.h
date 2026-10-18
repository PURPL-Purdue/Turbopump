#ifndef EXTERNAL_INPUT_H
#define EXTERNAL_INPUT_H

/*
 * Mission Control command path.  Whatever receives commands (the future
 * Ethernet code, or a UART/button test stub) decodes them into mc_command_t
 * and calls ExternalInput_Submit().  The control task drains them each cycle
 * with ExternalInput_Receive().  A queue in between keeps the two sides
 * thread-safe and means two commands in the same 1 ms are never merged.
 */

#include <stdbool.h>

#include "ecu_types.h"

/* Most commands the control task will take in one cycle (and queue depth). */
#define MC_COMMAND_QUEUE_LENGTH 16U

/* Create the command queue.  Call before the scheduler starts. */
void ExternalInput_Initialize(void);

/* Queue one decoded command.  Never blocks; false if the queue is full.
 * Callable from a task or an ISR. */
bool ExternalInput_Submit(const mc_command_t *command);

/* Take the next waiting command.  Never blocks; false if none is waiting. */
bool ExternalInput_Receive(mc_command_t *command);

#endif
