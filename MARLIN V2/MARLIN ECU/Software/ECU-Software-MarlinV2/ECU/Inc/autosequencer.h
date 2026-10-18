#ifndef AUTOSEQUENCER_H
#define AUTOSEQUENCER_H

/*
 * The autosequencer is a list of steps.  Each step holds its outputs until its
 * exit conditions are met (values that must be reached), then the sequence
 * advances to the next step.  Steps are not tied to a fixed time.
 */

#include <stdbool.h>
#include <stdint.h>

#include "ecu_types.h"

typedef struct
{
  ecu_command_t command; /* Outputs held while this step is active. */
  /* TODO: exit conditions (values that must be reached to advance). */
} autosequencer_step_t;

typedef struct
{
  const autosequencer_step_t *steps;
  uint32_t length;
} autosequencer_table_t;

/* The tables themselves (autosequencer_table.c). */
extern const autosequencer_table_t autosequencer_main_sequence;
extern const autosequencer_table_t autosequencer_fault_sequence;

/* Step `index` of the table, or NULL once past the end. */
const autosequencer_step_t *Autosequencer_GetStep(const autosequencer_table_t *table,
                                                  uint32_t index);

/* True once the step's exit conditions are met. */
bool Autosequencer_StepComplete(const autosequencer_step_t *step,
                                const sensor_data_t *sensors);

#endif
