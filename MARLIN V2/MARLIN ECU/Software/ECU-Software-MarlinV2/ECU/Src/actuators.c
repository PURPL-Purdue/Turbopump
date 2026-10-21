#include "actuators.h"

/*
 * PLACEHOLDER properties until the valve/actuator list is final.  Hardware
 * binding fields (port, pin, timer, timer_channel) are left NULL/0 for now.
 */
static actuator_t actuators[ACT_COUNT] = {
  [ACT_VALVE_0] = { .type = ACTUATOR_TYPE_SOLENOID },
  [ACT_VALVE_1] = { .type = ACTUATOR_TYPE_SOLENOID },
  [ACT_VALVE_2] = { .type = ACTUATOR_TYPE_SOLENOID },
  [ACT_VALVE_3] = { .type = ACTUATOR_TYPE_SOLENOID },
  [ACT_VALVE_4] = { .type = ACTUATOR_TYPE_SOLENOID },
  [ACT_VALVE_5] = { .type = ACTUATOR_TYPE_SOLENOID },
  [ACT_VALVE_6] = { .type = ACTUATOR_TYPE_SOLENOID },
  [ACT_VALVE_7] = { .type = ACTUATOR_TYPE_SOLENOID },
};

void Actuators_Initialize(void)
{
  /* TODO: configure each actuators[i] port/pin/timer, then drive safe. */
}

void Actuators_Apply(const ecu_command_t *command)
{
  (void)command;
  /* TODO: drive the outputs. */
}

const actuator_t *Actuators_Get(actuator_id_t id)
{
  return &actuators[id];
}
