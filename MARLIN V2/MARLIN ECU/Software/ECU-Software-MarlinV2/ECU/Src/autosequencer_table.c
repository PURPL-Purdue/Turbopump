#include <stddef.h>

#include "autosequencer.h"

/*
 * Sequence data, kept apart from the control code so the sequence can change
 * on its own.  EMPTY until the sequence is defined.
 */

const autosequencer_table_t autosequencer_sequence = {
  .steps = NULL,
  .length = 0U,
};

/* TODO: section start indices (main, fault, ...) once the sequence exists. */
