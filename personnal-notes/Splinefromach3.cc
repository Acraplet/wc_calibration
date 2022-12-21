 
class Monotone_Spline: public TSpline3_red {
// ************************
// closely follows TSpline3_red class to fit in easily with existing machinery
// Monotone spline is similar to regular cubic spline but enforce the condition that the interpolated value at any point
// must be between its two nearest knots, DOES NOT make the entire spline monotonic, only the segments

  public:
    // Empty constructor
    Monotone_Spline()
    :TSpline3_red()
    {
    }

    // The constructor that takes a TSpline3 pointer and copies in to memory
    Monotone_Spline(TSpline3* &spline, int Param = -1)
    { // need to override this so that daughter class SetFunc Gets called instead of TSpline3_red's version
      Par = NULL;
      XPos = NULL;
      YResp = NULL;
      SetFunc(spline, Param);
    }

    // Empty destructor
    ~Monotone_Spline() {
      delete[] Secants;
      delete[] Tangents;
      // this should also call base class destructor automatically
    }

    // Set a function
    inline void SetFunc(TSpline3* &spline, int Param = -1) {
      nPoints = spline->GetNp();
      ParamNo = Param;
      if (Par != NULL) {
        for (int i = 0; i < nPoints; ++i) {
          delete[] Par[i];
          Par[i] = NULL;
        }
        delete[] Par;
        Par = NULL;
      }
      if (XPos != NULL) delete[] XPos;
      if (YResp != NULL) delete[] YResp;
      // Save the parameters for each knot
      Par = new __float__*[nPoints];
      // Save the positions of the knots
      XPos = new __float__[nPoints];
      // Save the y response at each knot
      YResp = new __float__[nPoints];
      // values of the secants at each point (for calculating monotone spline)
      Secants = new __float__[nPoints -1];
      // values of the tangens at each point (for calculating monotone spline)
      Tangents = new __float__[nPoints];

      // get the knot values for the spline
      for (int i = 0; i < nPoints; ++i) {
        // 3 is the size of the TSpline3 coefficients
        Par[i] = new __float__[3];

        double x = -999.99, y = -999.99;
        spline->GetKnot(i, x, y);

        XPos[i]   = x;
        YResp[i]  = y;

        Tangents[i] = 0.0;
      }

      // deal with the case of two points (just do linear interpolation between them)
      if (nPoints ==2){
          Par[0][0] = (YResp[1] - YResp[0]) / ((XPos[1] - XPos[0]) * (XPos[1] - XPos[0]));
          Par[0][1] = 0.0;
          Par[0][2] = 0.0;
          // extra "virtual" segment at end to make Par array shape fit with knot arrays shapes
          Par[1][1] = 0.0;
          Par[1][2] = 0.0;

          return;
      } // if nPoints !=2 do full monotonic spline treatment:

      // first pass over knots to calculate the secants
      for (int i = 0; i < nPoints-1; ++i) {
        Secants[i] = (YResp[i+1] - YResp[i]) / (XPos[i+1] - XPos[i]);
        //std::cout<<"secant "<<i<<": "<<Secants[i]<<std::endl;
      }

      Tangents[0] = Secants[0];
      Tangents[nPoints-1] = Secants[nPoints -2];

      __float__ alpha;
      __float__ beta;

      // second pass over knots to calculate tangents
      for (int i = 1; i < nPoints-1; ++i) {
        if ((Secants[i-1] >= 0.0 && Secants[i] >= 0.0) | (Secants[i-1] < 0.0 && Secants[i] < 0.0)){ //check for same sign
          Tangents[i] = (Secants[i-1] + Secants[i]) /2.0;
        }
      }

      // third pass over knots to rescale tangents
      for (int i = 0; i < nPoints-1; ++i) {
        if (Secants[i] == 0.0){
          Tangents[i] = 0.0;
          Tangents[i+1] = 0.0;
        }

        else{
          alpha = Tangents[i]  / Secants[i];
          beta = Tangents[i+1] / Secants[i];

          if (alpha <0.0){
            Tangents[i] = 0.0;
          }
          if (beta < 0.0){
            Tangents[i+1] = 0.0;
          }

          if (alpha * alpha + beta * beta >9.0){
            __float__ tau = 3.0 / sqrt(alpha * alpha + beta * beta);
            Tangents[i]   = tau * alpha * Secants[i];
            Tangents[i+1] = tau * beta  * Secants[i];
          }
        }
        //std::cout<<"alpha, beta : "<<alpha<<", "<<beta<<std::endl;
        //std::cout<<"tangent "<<i<<": "<<Tangents[i]<<std::endl;

      } // finished rescaling tangents

      //std::cout<<"tangent "<<nPoints-1<<": "<<Tangents[nPoints-1]<<std::endl;

      // fourth pass over knots to calculate the coefficients for the spline
      __float__ dx;
      for(int i = 0; i <nPoints-1; i++){
        double b, c, d = -999.999;
        dx = XPos[i+1] - XPos[i];

        b = Tangents[i] * dx;
        c = 3.0* (YResp[i+1] - YResp[i]) -2.0 *dx * Tangents[i] - dx * Tangents[i +1];
        d = 2.0* (YResp[i] - YResp[i+1]) + dx * (Tangents[i] + Tangents[i+1]);

        Par[i][0] = b /  dx;
        Par[i][1] = c / (dx * dx);
        Par[i][2] = d / (dx * dx * dx);

        if((Par[i][0] == -999) | (Par[i][1] == -999) | (Par[i][2] ==-999) | (Par[i][0] == -999.999) | (Par[i][1] == -999.999) | (Par[i][2] ==-999.999)){
            std::cout<<"bad spline parameters for segment "<<i<<", will cause problems with GPU: (b, c, d) = "<<Par[i][0]<<", "<<Par[i][1]<<", "<<Par[i][2]<<std::endl;
        }
        //std::cout<<"b : "<<b<<std::endl;
        //std::cout<<"dx: "<<dx<<", x_0: "<<XPos[i]<<", x_1: "<<XPos[i+1]<<std::endl;
        //std::cout<<"    "<<" , y_0: "<<YResp[i]<<", y_1: "<<YResp[i+1]<<std::endl;
      }

      // include params for final "segment" outside of the spline so that par array fits with x and y arrays,
      // should never actually get used but if not set then the GPU code gets very angry
      Par[nPoints-1][0] = 0.0;
      Par[nPoints-1][1] = 0.0;
      Par[nPoints-1][2] = 0.0;

      // check the input spline for linear segments, if there are any then overwrite the calculated coefficients
      // this will pretty much only ever be the case if they are set to be linear in samplePDFND i.e. the user wants it to be linear
      for(int i = 0; i <nPoints-1; i++){
        double x = -999.99, y = -999.99, b = -999.99, c = -999.99, d = -999.99;
        spline->GetCoeff(i, x, y, b, c, d);

        if((c == 0.0 && d == 0.0)){
          Par[i][0] = b;
          Par[i][1] = 0.0;
          Par[i][2] = 0.0;
        }
      }

      delete spline;
      spline = NULL;
    }
    // finished calculating coeffs

    protected: //changed to protected from private so can be accessed by derived classes
    // values of the secants at each point (for calculating monotone spline)
    __float__ *Secants;
    // values of the tangents at each point (for calculating monotone spline)
    __float__ *Tangents;
};
