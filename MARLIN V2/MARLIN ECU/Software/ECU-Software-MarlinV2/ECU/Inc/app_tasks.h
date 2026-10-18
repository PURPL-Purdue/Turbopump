#ifndef APP_TASKS_H
#define APP_TASKS_H

/* Public entry points for the application layer. */
#include "ecu_types.h"

/* Called before the RTOS scheduler starts. */
void AppTasks_Initialize(void);

/* Called by the two CubeMX-created RTOS task wrappers. */
void AppControlTask_Run(void *argument);
void AppTelemetryTask_Run(void *argument);

#endif