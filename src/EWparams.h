//======== Electroweak quantities =================

// W and Z masses, widths and leptonic branching ratios:
const double mW = 80.398;   const double mZ = 91.1876;
const double Gamma_W = 2.118;  const double Gamma_Z = 2.4952;
const double Br_l_W = 0.1082;  const double Br_l_Z = 0.03363;  // Br_l_W was 0.1057 !!!
// sin^2 theta_W(M_Z)_{MS-bar}, "effective" (Z pole asymmetry) value
const double sstW = 0.23143 ;    // modern LEP/SLC value of Z "leptonic" angle 
//  const double sstW = 0.227 ;  // value used by van Neerven et al.
const double cstW = 1. - sstW ;
const double GF = 1.16639e-5 ; // Fermi constant
const double alpha_QED_0 = 1./137.036 ;  // fine structure constant
// For denominator factors, compute with sin^2theta_W defined by:
const double sstWden = PI * alpha_QED_0 /sqrt(2)/mW/mW/GF;
const double cstWden = 1. - sstWden ;
// fermion vector/axial couplings to gamma/Z, in van Neerven's normalization:
const double v_l_Z = - 1. + 4.*sstW ;   const double a_l_Z = 1. ;
const double Q_u = 2./3.;    const double Q_d = - 1./3.;
const double v_u_Z = 1. - 8./3.*sstW ;  const double v_d_Z = - 1. + 4./3.*sstW ;
const double a_u_Z = - 1. ;  const double a_d_Z = 1. ;
// quark vector/axial couplings to W, van Neerven's normalization 
// (same for u as d, so we drop u,d subscript):
const double v_W = 1./sqrt(2);  const double a_W = -1./sqrt(2);
// CKM data:
const double V_ud = 0.975;  const double V_us = 0.222;
const double V_cd = 0.222; const double V_cs = 0.974;
const double V_ud_sq = V_ud*V_ud;     const double V_us_sq = V_us*V_us;
const double V_cd_sq = V_cd*V_cd;     const double V_cs_sq = V_cs*V_cs;
const double V_ub_sq = 1. - V_ud_sq - V_us_sq;
const double V_cb_sq = 1. - V_cd_sq - V_cs_sq;
const double V_td_sq = 1. - V_ud_sq - V_cd_sq;
const double V_ts_sq = 1. - V_us_sq - V_cs_sq;
const double V_tb_sq = 1. - V_ub_sq - V_cb_sq;
