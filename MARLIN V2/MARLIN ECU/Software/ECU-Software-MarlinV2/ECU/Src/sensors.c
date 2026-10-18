#include "sensors.h"

/*
 * ADC, SPI, DFSDM, DMA, calibration, and conversion code can be added here
 * without exposing those hardware details to the control task.
 */
void Sensors_Initialize(void)
{
  /* Configure sensor peripherals here. */
}

void Sensors_StartDma(void)
{
  /* Start the selected DMA transfers here. */
}

bool Sensors_Read(sensor_data_t *data)
{
  /* No real sensor source exists yet, so no sample is reported. */
  (void)data;
  return false;
}