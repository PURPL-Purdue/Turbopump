#include "sensors.h"

/*
 * ADC, SPI, DFSDM, DMA, calibration, and conversion code can be added here
 * without exposing those hardware details to the control task.
 */

/* PLACEHOLDER properties until the instrumentation list is final. */
static sensor_t sensors[SNS_COUNT] = {
  [SNS_PT_0] = { .type = SENSOR_TYPE_PRESSURE,    .channel = 0U, .scale = 1.0f },
  [SNS_PT_1] = { .type = SENSOR_TYPE_PRESSURE,    .channel = 1U, .scale = 1.0f },
  [SNS_PT_2] = { .type = SENSOR_TYPE_PRESSURE,    .channel = 2U, .scale = 1.0f },
  [SNS_PT_3] = { .type = SENSOR_TYPE_PRESSURE,    .channel = 3U, .scale = 1.0f },
  [SNS_TC_0] = { .type = SENSOR_TYPE_TEMPERATURE, .channel = 4U, .scale = 1.0f },
  [SNS_TC_1] = { .type = SENSOR_TYPE_TEMPERATURE, .channel = 5U, .scale = 1.0f },
};

void Sensors_Initialize(void)
{
  /* TODO: configure sensor peripherals. */
}

void Sensors_StartDma(void)
{
  /* TODO: start the DMA transfers. */
}

bool Sensors_Read(sensor_data_t *data)
{
  (void)data;
  /* TODO: copy sensors[i].value / .valid into the snapshot. */
  return false;
}

const sensor_t *Sensors_Get(sensor_id_t id)
{
  return &sensors[id];
}
