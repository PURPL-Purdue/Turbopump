#include "autosequencer.h"

bool Autosequencer_GetRow(const autosequencer_table_t *table,
                          uint32_t tick,
                          ecu_command_t *command)
{
  if (tick >= table->length)
  {
    return false;
  }

  *command = table->rows[tick];
  return true;
}
