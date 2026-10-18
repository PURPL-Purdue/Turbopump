#ifndef ACTUATORS_H
#define ACTUATORS_H

/* Drives the physical outputs (valves, igniters) from an ecu_command_t. */
#include "ecu_types.h"

/* Set up the output pins and turn every output off. */
void Actuators_Initialize(void);

/* Drive every output: 0 is off, anything else is on. */
void Actuators_Apply(const ecu_command_t *command);

/* Turn every output off.  Also called from the fault handlers in Core/. */
void Actuators_ForceSafe(void);

#endif
