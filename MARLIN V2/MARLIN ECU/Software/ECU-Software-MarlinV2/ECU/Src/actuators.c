#include "actuators.h"

#include "stm32h7xx_hal.h"
#include "stm32h7xx_nucleo.h"

typedef struct
{
  GPIO_TypeDef *port; /* NULL: not wired to a pin yet, so it is skipped. */
  uint16_t pin;
} actuator_pin_t;

/*
 * PLACEHOLDER pin map: the first three valves drive the Nucleo LEDs so the
 * state machine can be watched on the bench.  Replace with the real ECU pins.
 */
static const actuator_pin_t actuator_pins[ACT_COUNT] = {
  [ACT_VALVE_0] = { LED1_GPIO_PORT, LED1_PIN },
  [ACT_VALVE_1] = { LED2_GPIO_PORT, LED2_PIN },
  [ACT_VALVE_2] = { LED3_GPIO_PORT, LED3_PIN },
};

void Actuators_Initialize(void)
{
  GPIO_InitTypeDef init = {0};

  LED1_GPIO_CLK_ENABLE();
  LED2_GPIO_CLK_ENABLE();
  LED3_GPIO_CLK_ENABLE();

  /* Pins start low (off) when they become outputs. */
  Actuators_ForceSafe();

  init.Mode = GPIO_MODE_OUTPUT_PP;
  init.Pull = GPIO_NOPULL;
  init.Speed = GPIO_SPEED_FREQ_LOW;
  for (uint32_t id = 0U; id < ACT_COUNT; id++)
  {
    if (actuator_pins[id].port != NULL)
    {
      init.Pin = actuator_pins[id].pin;
      HAL_GPIO_Init(actuator_pins[id].port, &init);
    }
  }
}

void Actuators_Apply(const ecu_command_t *command)
{
  for (uint32_t id = 0U; id < ACT_COUNT; id++)
  {
    if (actuator_pins[id].port != NULL)
    {
      HAL_GPIO_WritePin(actuator_pins[id].port, actuator_pins[id].pin,
                        (command->out[id] != 0U) ? GPIO_PIN_SET : GPIO_PIN_RESET);
    }
  }
}

void Actuators_ForceSafe(void)
{
  const ecu_command_t all_off = {0};

  Actuators_Apply(&all_off);
}
