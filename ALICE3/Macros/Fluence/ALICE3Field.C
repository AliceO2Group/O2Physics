std::function<void(const double *, double *)> field() {
    return [](const double *x, double *b) {
      double Rc;
      double R1;
      double R2;
      double B1;
      double B2;
      double beamStart = 500.; //[cm]
      double tokGauss = 1. / 0.1; // conversion from Tesla to kGauss
  
      bool isMagAbs = true;
  
      // ***********************
      // LAYOUT 1
      // ***********************
  
      // RADIUS
      Rc = 185.; //[cm]
      R1 = 220.; //[cm]
      R2 = 290.; //[cm]
  
      // To set the B2
      B1 = 2.;                                    //[T]
      B2 = -Rc * Rc * B1 / ((R2 * R2 - R1 * R1)); //[T]
  
      if ((abs(x[2]) <= beamStart) && (sqrt(x[0] * x[0] + x[1] * x[1]) < Rc)) {
        b[0] = 0.;
        b[1] = 0.;
        b[2] = B1 * tokGauss;
      } else if ((abs(x[2]) <= beamStart) &&
                 (sqrt(x[0] * x[0] + x[1] * x[1]) >= Rc &&
                  sqrt(x[0] * x[0] + x[1] * x[1]) < R1)) {
        b[0] = 0.;
        b[1] = 0.;
        b[2] = 0.;
      } else if ((abs(x[2]) <= beamStart) &&
                 (sqrt(x[0] * x[0] + x[1] * x[1]) >= R1 &&
                  sqrt(x[0] * x[0] + x[1] * x[1]) < R2)) {
        b[0] = 0.;
        b[1] = 0.;
        if (isMagAbs) {
          b[2] = B2 * tokGauss;
        } else {
          b[2] = 0.;
        }
      } else {
        b[0] = 0.;
        b[1] = 0.;
        b[2] = 0.;
      }
    };
  }  