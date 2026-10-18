#include "actuators.h"

/* TODO: private hardware binding (GPIO port/pin, timer channel) per actuator. */

/* PLACEHOLDER properties until the valve/actuator list is final. */
static actuator_t actuators[ACT_COUNT] = {
  [ACT_VALVE_0]   = { .type = ACTUATOR_TYPE_SOLENOID },
  [ACT_VALVE_1]   = { .type = ACTUATOR_TYPE_SOLENOID },
  [ACT_VALVE_2]   = { .type = ACTUATOR_TYPE_SOLENOID },
  [ACT_VALVE_3]   = { .type = ACTUATOR_TYPE_SOLENOID },
  [ACT_VALVE_4]   = { .type = ACTUATOR_TYPE_SOLENOID },
  [ACT_VALVE_5]   = { .type = ACTUATOR_TYPE_SOLENOID },
  [ACT_IGNITER_0] = { .type = ACTUATOR_TYPE_IGNITER },
  [ACT_IGNITER_1] = { .type = ACTUATOR_TYPE_IGNITER },
};

void Actuators_Initialize(void)
{
  /* TODO: configure the output pins. */
}

void Actuators_Apply(const ecu_command_t *command)
{
  (void)command;
  /* TODO: drive the outputs. */
}

void Actuators_ForceSafe(void)
{
  /* TODO: drive every output to actuators[i].safe_value. */
}

const actuator_t *Actuators_Get(actuator_id_t id)
{
  return &actuators[id];
}
