#include "external_input.h"

#include "app_tasks.h"

void ExternalInput_Initialize(void)
{
  /* TODO: configure the receive transport. */
}

void ExternalInput_Poll(void)
{
  /* TODO: read each received frame from the Ethernet RX buffers, decode it
   * into an ecu_inbound_frame_t and AppTasks_SubmitInbound() it. */
}
