/*
///////////////////////////////////////////////////////////////////////////////
//        Copyright (c) 2012-2020 Oscar Riveros. all rights reserved.        //
//                        oscar.riveros@peqnp.science                        //
//                                                                           //
//   without any restriction, Oscar Riveros reserved rights, patents and     //
//  commercialization of this knowledge or derived directly from this work.  //
///////////////////////////////////////////////////////////////////////////////

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <memory.h>
#include <complex.h>

struct data {
    int size;
    double complex *xy;
    int *seq;
    int *opt;
    double master;
};

double oracle (struct data *data)
{
  int i = 0;
  double local = 0.0;
  for (i = 0; i < data->size; i++) {
    local += round(cabs (data->xy[data->seq[(i + 1) % data->size]] - data->xy[data->seq[i]]));
    if (local > data->master) {
      return local;
    }
  }
  return local;
}

void invert (int i, int j, struct data *data)
{
  int aux = 0;
  while (i < j) {
    aux = data->seq[j];
    data->seq[j] = data->seq[i];
    data->seq[i] = aux;
    i++;
    j--;
  }
}

void hess (struct data *data)
{
  int i = 0, j = 0, anchor = 0;
  double local = 0, global = 0;
  data->master = DBL_MAX;
  for (;;) {
    global = DBL_MAX;
    anchor = 1;
    for (i = 0; i < data->size; i++) {
      for (j = 0; j < data->size; j++) {
        oo:
        invert (i, j, data);
        local = oracle (data);
        if (local < global) {
          global = local;
          if (global < data->master) {
            data->master = global;
            printf ("c %lf\n", data->master);
            memcpy(data->opt, data->seq, data->size * sizeof (int));
            anchor = 0;
          }
          goto oo;
        } else if (local > data->master) {
          goto oo;
        }
      }
    }
    if (anchor) {
      break;
    }
  }
}

int main (int argc, char **argv)
{
  int i = 0, aux = 0;
  double x = 0.0, y = 0.0;
  struct data data;

  FILE *instance = fopen (argv[1], "r");
  fscanf (instance, "%i", &data.size);

  data.xy = (double complex *) malloc (data.size * sizeof (double complex));
  data.seq = (int *) malloc (data.size * sizeof (int));
  data.opt = (int *) malloc (data.size * sizeof (int));

  for (i = 0; i < data.size; i++) {
    fscanf (instance, "%i %lf %lf", &aux, &x, &y);
    data.xy[i] = x + I * y;
    data.seq[i] = i;
  }

  hess (&data);

  for (i = 0; i < data.size; i++) {
    printf ("%i ", data.opt[i]);
  }
  printf ("\n");

  return 0;
}
