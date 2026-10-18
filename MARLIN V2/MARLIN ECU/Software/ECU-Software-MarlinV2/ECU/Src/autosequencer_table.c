#include <stddef.h>

#include "autosequencer.h"

/*
 * Sequence data, kept apart from the control code so the sequence can change
 * on its own.  EMPTY until the sequences are defined.
 */

const autosequencer_table_t autosequencer_main_sequence = {
  .steps = NULL,
  .length = 0U,
};

const autosequencer_table_t autosequencer_fault_sequence = {
  .steps = NULL,
  .length = 0U,
};
