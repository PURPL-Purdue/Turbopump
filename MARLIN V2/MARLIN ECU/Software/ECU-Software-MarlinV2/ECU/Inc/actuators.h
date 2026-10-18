#ifndef ACTUATORS_H
#define ACTUATORS_H

/* Drives the physical outputs (valves, igniters) from an ecu_command_t. */
#include "ecu_types.h"

typedef enum
{
  /* PLACEHOLDER types. */
  ACTUATOR_TYPE_SOLENOID = 0,
  ACTUATOR_TYPE_SERVO,
  ACTUATOR_TYPE_IGNITER
} actuator_type_t;

/* One actuator's properties plus its measured position. */
typedef struct
{
  actuator_type_t type;
  bool normally_open;   /* Safe state is de-energized.  */
  uint16_t safe_value;  /* Value ForceSafe drives.       */
} actuator_t;

/* Set up the output pins and turn every output off. */
void Actuators_Initialize(void);

/* Drive every output from the command. */
void Actuators_Apply(const ecu_command_t *command);

/* Turn every output off.  Also called from the fault handlers in Core/. */
void Actuators_ForceSafe(void);

/* Read-only access to one actuator's properties and position. */
const actuator_t *Actuators_Get(actuator_id_t id);

#endif
