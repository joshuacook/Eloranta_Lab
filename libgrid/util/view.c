/*
 * Show grid axis cut animation using xmgrace.
 *
 * Usage: view file1 file2...
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <grace_np.h>

void my_error_function(const char *msg) {
  fprintf(stderr, "library message: \"%s\"\n", msg);
}

int main(int argc, char **argv) {

  int i, j, npts;
  FILE *fp;
  double x, x_prev, y, begin, end, step = -1.0, ymax = -1.0, ymin = 0.0;

  if(argc < 2) {
    fprintf(stderr, "Usage: view file1 file2...\n");
    exit(1);
  }

  if(!(fp = fopen(argv[1], "r"))) {
    fprintf(stderr, "Can't open file %s.\n", argv[1]);
    exit(1);
  }
  /* First pass to figure out grid step length, begin and end. */
  for (i = 0; ; i++) {
    if(fscanf(fp, " %le %le", &x, &y) != 2) break;
    if(y > ymax) ymax = 1.5*y;
    if(!i) begin = x;
    if(i && step == -1.0) step = x - x_prev;
    if(feof(fp)) break;
    x_prev = x;
  }
  end = x_prev;
  step = fabs(step);
  fclose(fp);

  npts = 1 + (int) (0.5 + ((end - begin) / step));

  printf("begin = %le, end = %le, step = %le, ymax = %le, ymin = %le.\n", begin, end, step, ymax, ymin);

  GraceRegisterErrorFunction(my_error_function);
  /* Start Grace with a buffer size of 2048 and open the pipe */
  if (GraceOpenVA("/usr/bin/xmgrace", 2048, "-free", "-nosigcatch","-geometry", "930x730", NULL) == -1) {
    fprintf(stderr, "Can't start Grace. \n");
    exit(1);
  }

  /* Send some initialization commands to Grace */
  GracePrintf("world xmin %le", begin);
  GracePrintf("world xmax %le", end);
  GracePrintf("world ymin %le", ymin);
  GracePrintf("world ymax %le", ymax);

  GracePrintf("xaxis tick major %le", (end - begin) / 10.0);
  GracePrintf("xaxis tick minor %le", (end - begin) / 20.0);
  GracePrintf("yaxis tick major %le", ymax / 10.0);
  GracePrintf("yaxis tick minor %le", ymax / 20.0);

  for (i = 1; i < argc; i++) {

    if(!(fp = fopen(argv[i], "r"))) {
      fprintf(stderr, "Can't open file %s.\n", argv[i]);
      exit(1);
    }

    if(i > 1)
      GracePrintf("kill g0.s0");
    for (j = 0; j < npts; j++) {
      char buf[128];
      if(fscanf(fp, " %le %le", &x, &y) != 2) {
	fprintf(stderr, "File format error.\n");
	exit(1);
      }
      sprintf(buf, "File = %d", i);
      GracePrintf("title \"%s\"", buf);
      GracePrintf("g0.s0 point %le, %le", x, y);
    }
    GracePrintf("redraw");
    if(i == 1) sleep(5);
    else usleep(300000);
    fclose(fp);
  }
}
