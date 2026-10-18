#ifndef AUTOSEQUENCER_H
#define AUTOSEQUENCER_H

/*
 * The autosequencer table: row N is the complete output command for control
 * tick N (N * CONTROL_PERIOD_MS after the sequence starts).  The tables are
 * const, so they live in flash and are read directly; nothing is copied.
 */

#include <stdbool.h>
#include <stdint.h>

#include "ecu_types.h"

typedef struct
{
  const ecu_command_t *rows;
  uint32_t length;

  /* Main sequence only: rows [0, running_start_row) are STARTUP,
   * [running_start_row, shutdown_start_row) are RUNNING, the rest SHUTDOWN. */
  uint32_t running_start_row;
  uint32_t shutdown_start_row;
} autosequencer_table_t;

/* The tables themselves (autosequencer_table.c). */
extern const autosequencer_table_t autosequencer_main_sequence;
extern const autosequencer_table_t autosequencer_fault_sequence;

/* Copy row `tick` into `command`.  Returns false once past the end. */
bool Autosequencer_GetRow(const autosequencer_table_t *table,
                          uint32_t tick,
                          ecu_command_t *command);

#endif
