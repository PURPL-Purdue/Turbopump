#include <stddef.h>

#include "autosequencer.h"

/*
 * Sequence data, kept apart from the control code so the sequence can change
 * on its own.  These will come from the autosequencer spreadsheet.
 */

/* Main sequence: EMPTY until the real sequence exists.  START then just goes
 * straight back to STANDBY. */
const autosequencer_table_t autosequencer_main_sequence = {
  .rows = NULL,
  .length = 0U,
  .running_start_row = 0U,
  .shutdown_start_row = 0U,
};

/* Fault sequence: PLACEHOLDER single row with every output off. */
static const ecu_command_t fault_rows[] = {
  { .out = { 0U } },
};

const autosequencer_table_t autosequencer_fault_sequence = {
  .rows = fault_rows,
  .length = sizeof(fault_rows) / sizeof(fault_rows[0]),
};
