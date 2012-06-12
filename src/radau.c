  void Radau15::init() {
    
    type = RA15;
    
    h[0] = 0.0;
    h[1] = 0.05626256053692215;
    h[2] = 0.18024069173689236;
    h[3] = 0.35262471711316964;
    h[4] = 0.54715362633055538;
    h[5] = 0.73421017721541053;
    h[6] = 0.88532094683909577;
    h[7] = 0.97752061356128750;
    
    xc[0] = 0.5;
    xc[1] = 0.16666666666666667;
    xc[2] = 0.08333333333333333;
    xc[3] = 0.05;
    xc[4] = 0.03333333333333333;
    xc[5] = 0.02380952380952381;
    xc[6] = 0.01785714285714286;
    xc[7] = 0.01388888888888889;
    
    vc[0] = 0.5;
    vc[1] = 0.3333333333333333;
    vc[2] = 0.25;
    vc[3] = 0.2;
    vc[4] = 0.1666666666666667;
    vc[5] = 0.1428571428571429;
    vc[6] = 0.125;
    
    // r.resize(28);
    //
    int j,k,l;
    l=0;
    for (j=1;j<8;++j) {
      for(k=0;k<j;++k) {
      r[l] = 1.0 / (h[j] - h[k]);
      ++l;
      }
    }
    
    /* 
       for(k=0;k<28;++k) {
       printf("r[%02i]= %e\n",k,r[k]);
       }
    */
    
    // c.resize(21);
    // d.resize(21);
    //
    c[0] = -h[1];
    d[0] =  h[1];
    l=0;
    for (j=2;j<7;++j) {
      ++l;
      c[l] = -h[j] * c[l-j+1];
      d[l] =  h[1] * d[l-j+1];
      for(k=2;k<j;++k) {
      ++l;
      c[l] = c[l-j] - h[j] * c[l-j+1];
      d[l] = d[l-j] + h[k] * d[l-j+1];
      }
      ++l;
      c[l] = c[l-j] - h[j];
      d[l] = d[l-j] + h[j]; 
    }
    
    /*
      for(k=0;k<21;++k) {
      printf("c[%02i]= %e    d[%02i]= %e\n",k,c[k],k,d[k]);
      }
    */
    
    nv    = 0;
    niter = 6;
    
    // s.resize(9);
    
    size = 0;
  }
  
  Radau15::~Radau15() {
    
  }
  
  void Radau15::Bodies_Mass_or_N_Bodies_Changed(const Frame &frame) {
    
    // cerr << " *** Bodies_Mass_or_N_Bodies_Changed() called!! ***" << endl;
    
    nv = 3*frame.size();
    //
    if (nv > x.size()) {
      g.resize(7);
      b.resize(7);
      e.resize(7);
      //
      unsigned int l;
      for (l=0;l<7;++l) {
      g[l].resize(nv);
      b[l].resize(nv);
      e[l].resize(nv);    
      }
      //
      x.resize(nv);
      v.resize(nv);
      a.resize(nv);
      //
      x1.resize(nv);
      v1.resize(nv);
      a1.resize(nv);
      //
      acc.resize(frame.size());
      mass.resize(frame.size());
    }
    // reset (long... bad style... may use memset...)
    /* unsigned int j,k;
       for (j=0;j<7;++j) {
       for(k=0;k<nv;++k) {
       b[j][k] = 0.0;
       e[j][k] = 0.0;
       }    
       }
    */
    // better
    memset(&b[0][0],0,7*nv);
    memset(&e[0][0],0,7*nv);
    
    /* 
       {
       // test
       for (unsigned int j=0;j<7;++j) {
       for(unsigned int k=0;k<nv;++k) {
       std::cerr << "zero? j:" << j << " k:" << k << " b:" << b[j][k] << " e:" << e[j][k] << std::endl;
       }    
       }
       }
    */
    
    for(unsigned int k=0;k<frame.size();++k) {
      mass[k] = frame[k].mass();
    }
    
    size = frame.size();
  }
  
  void Radau15::Step(const Frame & frame_in, Frame & frame_out, Interaction * interaction) {
    
    // cerr << "-> inside  Radau15::Step()..." << endl;
    
    // cerr << "RA15: initial timestep: " << timestep << endl;
    
    // N.B. Input/output must be in coordinates with respect to the central body.
    
    // frames...
    frame_out = frame_in;
    
    niter = 2;
    // if (frame_out.size() != mass.size()) {
    if (frame_out.size() != size) {
      Bodies_Mass_or_N_Bodies_Changed(frame_out);
      niter = 6;
    } else {
      unsigned int l=0;
      while (l != frame_out.size()) {
      if (frame_out[l].mass() != mass[l]) {
        Bodies_Mass_or_N_Bodies_Changed(frame_out);
        niter = 6;
        break;
      }
      ++l;
      }
    }
    
    // cerr << "niter: " << niter << endl;
    
    interaction->Acceleration(frame_out,acc);
    
    unsigned int j,k;
    for(k=0;k<frame_in.size();++k) {
      x1[3*k]   = frame_in[k].position().x;
      x1[3*k+1] = frame_in[k].position().y;
      x1[3*k+2] = frame_in[k].position().z;
      //
      v1[3*k]   = frame_in[k].velocity().x;
      v1[3*k+1] = frame_in[k].velocity().y;
      v1[3*k+2] = frame_in[k].velocity().z;
      //
      a1[3*k]   = acc[k].x;
      a1[3*k+1] = acc[k].y;  
      a1[3*k+2] = acc[k].z;
    }
    
    for(k=0;k<nv;++k) {
      g[0][k] = b[6][k]*d[15] + b[5][k]*d[10] + b[4][k]*d[6] + b[3][k]*d[3]  + b[2][k]*d[1]  + b[1][k]*d[0]  + b[0][k];
      g[1][k] = b[6][k]*d[16] + b[5][k]*d[11] + b[4][k]*d[7] + b[3][k]*d[4]  + b[2][k]*d[2]  + b[1][k];
      g[2][k] = b[6][k]*d[17] + b[5][k]*d[12] + b[4][k]*d[8] + b[3][k]*d[5]  + b[2][k];
      g[3][k] = b[6][k]*d[18] + b[5][k]*d[13] + b[4][k]*d[9] + b[3][k];
      g[4][k] = b[6][k]*d[19] + b[5][k]*d[14] + b[4][k];
      g[5][k] = b[6][k]*d[20] + b[5][k];
      g[6][k] = b[6][k];
    }
    
    double tmp,gk;
    double q1,q2,q3,q4,q5,q6,q7;
    
    unsigned int main_loop_counter;
    for (main_loop_counter=0;main_loop_counter<niter;++main_loop_counter) {
      for(j=1;j<8;++j) {
      
      // s[0] = timestep * h[j];
      s[0] = timestep.GetDouble() * h[j];
      s[1] = s[0] * s[0] * 0.5;
      s[2] = s[1] * h[j] * 0.3333333333333333;
      s[3] = s[2] * h[j] * 0.5;
      s[4] = s[3] * h[j] * 0.6;
      s[5] = s[4] * h[j] * 0.6666666666666667;
      s[6] = s[5] * h[j] * 0.7142857142857143;
      s[7] = s[6] * h[j] * 0.75;
      s[8] = s[7] * h[j] * 0.7777777777777778;
      
      for(k=0;k<nv;++k) {
        x[k] = ( s[8]*b[6][k] +
               s[7]*b[5][k] + 
               s[6]*b[4][k] + 
               s[5]*b[3][k] + 
               s[4]*b[2][k] + 
               s[3]*b[1][k] + 
               s[2]*b[0][k] ) +
          s[1]*a1[k] + 
          s[0]*v1[k] + 
          x1[k];
      }
      
      // needed only if using a velocity-dependent interaction...
      if (interaction->depends_on_velocity()) { 
        // s[0] = timestep * h[j];
        s[0] = timestep.GetDouble() * h[j];
        s[1] = s[0] * h[j] * 0.5;
        s[2] = s[1] * h[j] * 0.6666666666666667;
        s[3] = s[2] * h[j] * 0.75;
        s[4] = s[3] * h[j] * 0.8;
        s[5] = s[4] * h[j] * 0.8333333333333333;
        s[6] = s[5] * h[j] * 0.8571428571428571;
        s[7] = s[6] * h[j] * 0.875;
        
        for(k=0;k<nv;++k) {
          v[k] = ( s[7]*b[6][k] + 
                 s[6]*b[5][k] + 
                 s[5]*b[4][k] + 
                 s[4]*b[3][k] + 
                 s[3]*b[2][k] + 
                 s[2]*b[1][k] +
                 s[1]*b[0][k] ) +
            s[0]*a1[k] + 
            v1[k];
        }
      }
      
      {
        Vector rr,vv,drr,dvv;
        for(k=0;k<frame_out.size();++k) {
          
          frame_out[k] = frame_in[k];
          
          rr.x = x[3*k];
          rr.y = x[3*k+1];
          rr.z = x[3*k+2];
          
          drr = rr - frame_in[k].position();
          frame_out[k].AddToPosition(drr);
          
          vv.x = v[3*k];
          vv.y = v[3*k+1];
          vv.z = v[3*k+2];
          
          dvv = vv - frame_in[k].velocity();
          frame_out[k].AddToVelocity(dvv);
        }
      }
      
      if (interaction->IsSkippingJPLPlanets()) {
        frame_out.SetTime(frame_in+timestep*h[j]);
        frame_out.ForceJPLEphemerisData();
      }
      //
      interaction->Acceleration(frame_out,acc);
      
      for(k=0;k<frame_out.size();++k) {
        a[3*k]   = acc[k].x;
        a[3*k+1] = acc[k].y;  
        a[3*k+2] = acc[k].z;
      }
      
      switch (j) {
      case 1: 
        for(k=0;k<nv;++k) {
          tmp = g[0][k];
          g[0][k]  = (a[k] - a1[k]) * r[0];
          b[0][k] += g[0][k] - tmp;
        } break;
      case 2: 
        for(k=0;k<nv;++k) {
          tmp = g[1][k];
          gk = a[k] - a1[k];
          g[1][k] = (gk*r[1] - g[0][k])*r[2];
          tmp = g[1][k] - tmp;
          b[0][k] += tmp * c[0];
          b[1][k] += tmp;
        } break;
      case 3: 
        for(k=0;k<nv;++k) {
          tmp = g[2][k];
          gk = a[k] - a1[k];
          g[2][k] = ((gk*r[3] - g[0][k])*r[4] - g[1][k])*r[5];
          tmp = g[2][k] - tmp;
          b[0][k] += tmp * c[1];
          b[1][k] += tmp * c[2];
          b[2][k] += tmp;
        } break;
      case 4:
        for(k=0;k<nv;++k) {
          tmp = g[3][k];
          gk = a[k] - a1[k];
          g[3][k] = (((gk*r[6] - g[0][k])*r[7] - g[1][k])*r[8] - g[2][k])*r[9];
          tmp = g[3][k] - tmp;
          b[0][k] += tmp * c[3];
          b[1][k] += tmp * c[4];
          b[2][k] += tmp * c[5];
          b[3][k] += tmp;
        } break;
      case 5:
        for(k=0;k<nv;++k) {
          tmp = g[4][k];
          gk = a[k] - a1[k];
          g[4][k] = ((((gk*r[10] - g[0][k])*r[11] - g[1][k])*r[12] - g[2][k])*r[13] - g[3][k])*r[14];
          tmp = g[4][k] - tmp;
          b[0][k] += tmp * c[6];
          b[1][k] += tmp * c[7];
          b[2][k] += tmp * c[8];
          b[3][k] += tmp * c[9];
          b[4][k] += tmp;
        } break;
      case 6:
        for(k=0;k<nv;++k) {
          tmp = g[5][k];
          gk = a[k] - a1[k];
          g[5][k] = (((((gk*r[15] - g[0][k])*r[16] - g[1][k])*r[17] - g[2][k])*r[18] - g[3][k])*r[19] - g[4][k])*r[20];
          tmp = g[5][k] - tmp;
          b[0][k] += tmp * c[10];
          b[1][k] += tmp * c[11];
          b[2][k] += tmp * c[12];
          b[3][k] += tmp * c[13];
          b[4][k] += tmp * c[14];
          b[5][k] += tmp;
        } break;
      case 7:
        for(k=0;k<nv;++k) {
          tmp = g[6][k];
          gk = a[k] - a1[k];
          g[6][k] = ((((((gk*r[21] - g[0][k])*r[22] - g[1][k])*r[23] - g[2][k])*r[24] - g[3][k])*r[25] - g[4][k])*r[26] - g[5][k])*r[27];
          tmp = g[6][k] - tmp;
          b[0][k] += tmp * c[15];
          b[1][k] += tmp * c[16];
          b[2][k] += tmp * c[17];
          b[3][k] += tmp * c[18];
          b[4][k] += tmp * c[19];
          b[5][k] += tmp * c[20];
          b[6][k] += tmp;
        } break;
      default:
        ORSA_LOGIC_ERROR("aieeee!!!");
      }
      
      }
    }
    
    timestep_done = timestep;
    
    // Estimate suitable sequence size for the next call
    tmp = 0.0;
    for(k=0;k<nv;++k) {
      tmp = MAX(tmp,fabs(b[6][k]));
    }
    // if (tmp!=0.0) tmp /= (72.0 * secure_pow(fabs(timestep),7));
    // if (tmp!=0.0) tmp /= (72.0 * secure_pow(fabs(timestep.GetDouble()),7));
    if (tmp!=0.0) tmp /= (72.0 * pow(fabs(timestep.GetDouble()),7));
    
    if (tmp < 1.0e-50) { // is equal to zero?
      // timestep = timestep_done * 1.4;
      timestep = timestep_done * 1.4;
    } else {
      
      // old rule...
      // timestep = copysign(secure_pow(accuracy/tmp,0.1111111111111111),timestep_done); // 1/9=0.111...
      // timestep = copysign(secure_pow(accuracy/tmp,0.1111111111111111),timestep_done.GetDouble()); // 1/9=0.111...
      timestep = copysign(pow(accuracy/tmp,0.1111111111111111),timestep_done.GetDouble()); // 1/9=0.111...
      
    }
    //
    // if (fabs(timestep/timestep_done) < 1.0) {
    if (fabs(timestep.GetDouble()/timestep_done.GetDouble()) < 1.0) {
      timestep = timestep_done * 0.8;
      // std::cerr << "Radau: step rejected! New proposed timestep: " << timestep.GetDouble() << std::endl;
      frame_out = frame_in;
      niter = 6;
      return;
    }
    //
    if (fabs(timestep.GetDouble()/timestep_done.GetDouble()) > 1.4) timestep = timestep_done * 1.4;
    
    // std::cerr << "RA15: new timestep: " << timestep.GetDouble() << std::endl;
    
    // Find new position and velocity values at end of the sequence
    tmp = timestep_done.GetDouble() * timestep_done.GetDouble();
    for(k=0;k<nv;++k) {
      x1[k] = ( xc[7]*b[6][k] +
            xc[6]*b[5][k] + 
            xc[5]*b[4][k] + 
            xc[4]*b[3][k] + 
            xc[3]*b[2][k] + 
            xc[2]*b[1][k] + 
            xc[1]*b[0][k] + 
            xc[0]*a1[k]   ) * tmp + v1[k]*timestep_done.GetDouble() + x1[k];
      
      v1[k] = ( vc[6]*b[6][k] + 
            vc[5]*b[5][k] + 
            vc[4]*b[4][k] +
            vc[3]*b[3][k] + 
            vc[2]*b[2][k] + 
            vc[1]*b[1][k] +
            vc[0]*b[0][k] + 
            a1[k])        * timestep_done.GetDouble() + v1[k];
    }
    
    {
      Vector rr,vv,drr,dvv;
      for(k=0;k<frame_out.size();++k) {
      
      frame_out[k] = frame_in[k];
      
      rr.x = x1[3*k];
      rr.y = x1[3*k+1];
      rr.z = x1[3*k+2];
      
      drr = rr - frame_in[k].position();  
      frame_out[k].AddToPosition(drr);
      
      vv.x = v1[3*k];
      vv.y = v1[3*k+1];
      vv.z = v1[3*k+2];
      
      dvv = vv - frame_in[k].velocity();
      frame_out[k].AddToVelocity(dvv);
      }
    }
    
    // frame_out += timestep_done;
    // frame_out.SetTime(frame_in + timestep_done);
    frame_out.SetTime(frame_in);
    frame_out += timestep_done;
    
    // Predict new B values to use at the start of the next sequence. The predicted
    // values from the last call are saved as E. The correction, BD, between the
    // actual and predicted values of B is applied in advance as a correction.
    //
    q1 = timestep.GetDouble() / timestep_done.GetDouble();
    q2 = q1 * q1;
    q3 = q1 * q2;
    q4 = q2 * q2;
    q5 = q2 * q3;
    q6 = q3 * q3;
    q7 = q3 * q4;
    
    for(k=0;k<nv;++k) {
       
      s[0] = b[0][k] - e[0][k];
      s[1] = b[1][k] - e[1][k];
      s[2] = b[2][k] - e[2][k];
      s[3] = b[3][k] - e[3][k];
      s[4] = b[4][k] - e[4][k];
      s[5] = b[5][k] - e[5][k];
      s[6] = b[6][k] - e[6][k];
      
      // Estimate B values for the next sequence
      
      e[0][k] = q1*(b[6][k]* 7.0 + b[5][k]* 6.0 + b[4][k]* 5.0 + b[3][k]* 4.0 + b[2][k]* 3.0 + b[1][k]*2.0 + b[0][k]);
      e[1][k] = q2*(b[6][k]*21.0 + b[5][k]*15.0 + b[4][k]*10.0 + b[3][k]* 6.0 + b[2][k]* 3.0 + b[1][k]);
      e[2][k] = q3*(b[6][k]*35.0 + b[5][k]*20.0 + b[4][k]*10.0 + b[3][k]* 4.0 + b[2][k]);
      e[3][k] = q4*(b[6][k]*35.0 + b[5][k]*15.0 + b[4][k]* 5.0 + b[3][k]);
      e[4][k] = q5*(b[6][k]*21.0 + b[5][k]* 6.0 + b[4][k]);
      e[5][k] = q6*(b[6][k]* 7.0 + b[5][k]);
      e[6][k] = q7* b[6][k];
      
      b[0][k] = e[0][k] + s[0];
      b[1][k] = e[1][k] + s[1];
      b[2][k] = e[2][k] + s[2];
      b[3][k] = e[3][k] + s[3];
      b[4][k] = e[4][k] + s[4];
      b[5][k] = e[5][k] + s[5];
      b[6][k] = e[6][k] + s[6];
      
    }
    
    // cerr << "-> out of Radau15::Step()..." << endl;
    
  }
  
} // namespace orsa

