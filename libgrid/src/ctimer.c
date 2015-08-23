/*
 * CPU usage time routines.
 *
 */

#include "grid.h"

EXPORT void grid_timer_start(grid_timer *timer) {

  gettimeofday(&timer->zero_time, 0);
  timer->zero_clock = clock();
}

EXPORT double grid_timer_wall_clock_time(grid_timer *timer) {

  struct timeval now;

  gettimeofday(&now, 0);

  return 1.0 * (now.tv_sec - timer->zero_time.tv_sec)
    + 1e-6 * (now.tv_usec - timer->zero_time.tv_usec);
}

EXPORT double grid_timer_cpu_time(grid_timer *timer) {

  clock_t now = clock();

  return (now - timer->zero_clock) / ((double) CLOCKS_PER_SEC);
}

