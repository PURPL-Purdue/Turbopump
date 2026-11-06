#include "autosequencer.h"

#include <stddef.h>

const autosequencer_step_t *Autosequencer_GetStep(const autosequencer_table_t *table,
                                                  uint32_t index)
{
  if (index >= table->length)
  {
    return NULL;
  }

  return &table->steps[index];
}

bool Autosequencer_StepComplete(const autosequencer_step_t *step,
                                const sensor_data_t *sensors)
{
  (void)step;
  (void)sensors;
  /* TODO: evaluate the step's exit conditions. */
  return false;
}
