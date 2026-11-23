#ifndef APP_TASKS_H
#define APP_TASKS_H

/* Public entry points for the application layer. */
#include <stdbool.h>

#include "ecu_types.h"

/* Called before the RTOS scheduler starts. */
void AppTasks_Initialize(void);

/* Called by the two CubeMX-created RTOS task wrappers. */
void AppControlTask_Run(void *argument);
void AppTelemetryTask_Run(void *argument);

/* Queue one inbound frame for the control task.  Never blocks; false if the
 * queue is full.  Callable from a task or an ISR. */
bool AppTasks_SubmitInbound(const ecu_inbound_frame_t *frame);

#endif
