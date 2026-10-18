#include "telemetry.h"

/*
 * The telemetry task calls this module after receiving a complete snapshot.
 * Ethernet buffers, encoding, and the wire protocol belong here later.
 */
void Telemetry_Initialize(void)
{
  /* Configure the future telemetry transport here. */
}

void Telemetry_Process(const telemetry_frame_t *frame)
{
  (void)frame;
  /* Serialize and transmit frames through the Ethernet driver here. */
}