/*
newtonsMethodForQuadraticSystem.h
*/

void differential_for_quadratic(
    double t,
     double *x,
      double *xp);

void trace_one_period(
    double *initial_xz,
     double *end_xz,
      double **Jacobian,
       double period);

void find_poincare_fixed_point(
    double *initial_xz, //also the output fixed points
     double period);