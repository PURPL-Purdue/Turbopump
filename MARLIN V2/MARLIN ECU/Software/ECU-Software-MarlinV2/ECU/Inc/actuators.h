#ifndef ACTUATORS_H
#define ACTUATORS_H

/* Drives the physical outputs (valves) from an ecu_command_t. */
#include "ecu_types.h"
#include "stm32h7xx_hal.h"

typedef enum
{
  /* PLACEHOLDER types. */
  ACTUATOR_TYPE_SOLENOID = 0,
  ACTUATOR_TYPE_SERVO
} actuator_type_t;

/* One actuator's properties and hardware binding.  No position feedback:
 * an actuator's state is only known from what was last commanded. */
typedef struct
{
  actuator_type_t type;
  bool normally_open;          /* Safe state is de-energized.            */

  /* Hardware binding, unassigned (NULL/0) until the pinout is final. */
  GPIO_TypeDef *port;          /* Drive pin (solenoid).                  */
  uint16_t pin;
  TIM_HandleTypeDef *timer;    /* PWM output (servo); NULL if unused.    */
  uint32_t timer_channel;
} actuator_t;

/* Configure each actuator's hardware from its table entry and turn every
 * output off.  The actuator objects themselves are static (actuators.c). */
void Actuators_Initialize(void);

/* Drive every output from the command. */
void Actuators_Apply(const ecu_command_t *command);

/* Read-only access to one actuator's properties. */
const actuator_t *Actuators_Get(actuator_id_t id);

#endif
