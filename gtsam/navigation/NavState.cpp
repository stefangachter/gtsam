/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file    NavState.cpp
 * @brief   Navigation state composing of attitude, position, and velocity
 * @author  Frank Dellaert
 * @date    July 2015
 **/

#include <gtsam/navigation/NavState.h>

#include <string>

namespace gtsam {

//------------------------------------------------------------------------------
NavState NavState::Create(const Rot3& R, const Point3& t, const Velocity3& v,
    OptionalJacobian<9, 3> H1, OptionalJacobian<9, 3> H2,
    OptionalJacobian<9, 3> H3) {
  if (H1)
    *H1 << I_3x3, Z_3x3, Z_3x3;
  if (H2)
    *H2 << Z_3x3, R.transpose(), Z_3x3;
  if (H3)
    *H3 << Z_3x3, Z_3x3, R.transpose();
  return NavState(R, t, v);
}
//------------------------------------------------------------------------------
NavState NavState::FromPoseVelocity(const Pose3& pose, const Vector3& vel,
    OptionalJacobian<9, 6> H1, OptionalJacobian<9, 3> H2) {
  if (H1)
    *H1 << I_3x3, Z_3x3, Z_3x3, I_3x3, Z_3x3, Z_3x3;
  if (H2)
    *H2 << Z_3x3, Z_3x3, pose.rotation().transpose();
  return NavState(pose, vel);
}

//------------------------------------------------------------------------------
const Rot3& NavState::attitude(OptionalJacobian<3, 9> H) const {
  if (H)
    *H << I_3x3, Z_3x3, Z_3x3;
  return R_;
}

//------------------------------------------------------------------------------
const Point3& NavState::position(OptionalJacobian<3, 9> H) const {
  if (H)
    *H << Z_3x3, R(), Z_3x3;
  return t_;
}

//------------------------------------------------------------------------------
const Vector3& NavState::velocity(OptionalJacobian<3, 9> H) const {
  if (H)
    *H << Z_3x3, Z_3x3, R();
  return v_;
}

//------------------------------------------------------------------------------
Vector3 NavState::bodyVelocity(OptionalJacobian<3, 9> H) const {
  const Rot3& nRb = R_;
  const Vector3& n_v = v_;
  Matrix3 D_bv_nRb;
  Vector3 b_v = nRb.unrotate(n_v, H ? &D_bv_nRb : 0);
  if (H)
    *H << D_bv_nRb, Z_3x3, I_3x3;
  return b_v;
}

//------------------------------------------------------------------------------
Matrix5 NavState::matrix() const {
  Matrix3 R = this->R();

  Matrix5 T = Matrix5::Identity();
  T.block<3, 3>(0, 0) = R;
  T.block<3, 1>(0, 3) = t_;
  T.block<3, 1>(0, 4) = v_;
  return T;
}

//------------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const NavState& state) {
  os << "R: " << state.attitude() << "\n";
  os << "p: " << state.position().transpose() << "\n";
  os << "v: " << state.velocity().transpose();
  return os;
}

//------------------------------------------------------------------------------
void NavState::print(const std::string& s) const {
  std::cout << (s.empty() ? s : s + " ") << *this << std::endl;
}

//------------------------------------------------------------------------------
bool NavState::equals(const NavState& other, double tol) const {
  return R_.equals(other.R_, tol) && traits<Point3>::Equals(t_, other.t_, tol)
      && equal_with_abs_tol(v_, other.v_, tol);
}

//------------------------------------------------------------------------------
NavState NavState::inverse() const {
  Rot3 Rt = R_.inverse();
  return NavState(Rt, Rt * (-t_), Rt * -(v_));
}

//------------------------------------------------------------------------------
NavState NavState::Expmap(const Vector9& xi, OptionalJacobian<9, 9> Hxi) {
  // Get angular velocity w and components rho (for t) and nu (for v) from xi
  Vector3 w = xi.head<3>(), rho = xi.segment<3>(3), nu = xi.tail<3>();

  // Compute rotation using Expmap
  Matrix3 Jr;
  Rot3 R = Rot3::Expmap(w, Hxi ? &Jr : nullptr);

  // Compute translations and optionally their Jacobians Q in w
  // The Jacobians with respect to rho and nu are equal to Jr
  Matrix3 Qt, Qv;
  Vector3 t = Pose3::ExpmapTranslation(w, rho, Hxi ? &Qt : nullptr);
  Vector3 v = Pose3::ExpmapTranslation(w, nu, Hxi ? &Qv : nullptr);

  if (Hxi) {
    *Hxi << Jr, Z_3x3, Z_3x3,  //
        Qt, Jr, Z_3x3,         //
        Qv, Z_3x3, Jr;
  }

  return NavState(R, t, v);
}

//------------------------------------------------------------------------------
Vector9 NavState::Logmap(const NavState& state, OptionalJacobian<9, 9> Hstate) {
  if (Hstate) *Hstate = LogmapDerivative(state);

  const Vector3 phi = Rot3::Logmap(state.rotation());
  const Vector3& p = state.position();
  const Vector3& v = state.velocity();
  const double t = phi.norm();
  if (t < 1e-8) {
    Vector9 log;
    log << phi, p, v;
    return log;

  } else {
    const Matrix3 W = skewSymmetric(phi / t);

    const double Tan = tan(0.5 * t);
    const Vector3 Wp = W * p;
    const Vector3 Wv = W * v;
    const Vector3 rho = p - (0.5 * t) * Wp + (1 - t / (2. * Tan)) * (W * Wp);
    const Vector3 nu = v - (0.5 * t) * Wv + (1 - t / (2. * Tan)) * (W * Wv);
    Vector9 log;
    // Order is ω, p, v
    log << phi, rho, nu;
    return log;
  }
}

//------------------------------------------------------------------------------
Matrix9 NavState::AdjointMap() const {
  const Matrix3 R = R_.matrix();
  Matrix3 A = skewSymmetric(t_.x(), t_.y(), t_.z()) * R;
  Matrix3 B = skewSymmetric(v_.x(), v_.y(), v_.z()) * R;
  // Eqn 2 in Barrau20icra
  Matrix9 adj;
  adj << R, Z_3x3, Z_3x3, A, R, Z_3x3, B, Z_3x3, R;
  return adj;
}

//------------------------------------------------------------------------------
Vector9 NavState::Adjoint(const Vector9& xi_b, OptionalJacobian<9, 9> H_state,
                          OptionalJacobian<9, 9> H_xib) const {
  const Matrix9 Ad = AdjointMap();

  // Jacobians
  if (H_state) *H_state = -Ad * adjointMap(xi_b);
  if (H_xib) *H_xib = Ad;

  return Ad * xi_b;
}

//------------------------------------------------------------------------------
Matrix9 NavState::adjointMap(const Vector9& xi) {
  Matrix3 w_hat = skewSymmetric(xi(0), xi(1), xi(2));
  Matrix3 v_hat = skewSymmetric(xi(3), xi(4), xi(5));
  Matrix3 a_hat = skewSymmetric(xi(6), xi(7), xi(8));
  Matrix9 adj;
  adj << w_hat, Z_3x3, Z_3x3, v_hat, w_hat, Z_3x3, a_hat, Z_3x3, w_hat;
  return adj;
}

//------------------------------------------------------------------------------
Vector9 NavState::adjoint(const Vector9& xi, const Vector9& y,
                          OptionalJacobian<9, 9> Hxi,
                          OptionalJacobian<9, 9> H_y) {
  if (Hxi) {
    Hxi->setZero();
    for (int i = 0; i < 9; ++i) {
      Vector9 dxi;
      dxi.setZero();
      dxi(i) = 1.0;
      Matrix9 Gi = adjointMap(dxi);
      Hxi->col(i) = Gi * y;
    }
  }

  const Matrix9& ad_xi = adjointMap(xi);
  if (H_y) *H_y = ad_xi;

  return ad_xi * y;
}

//------------------------------------------------------------------------------
Matrix9 NavState::ExpmapDerivative(const Vector9& xi) {
  Matrix9 J;
  Expmap(xi, J);
  return J;
}

//------------------------------------------------------------------------------
Matrix9 NavState::LogmapDerivative(const NavState& state) {
  const Vector9 xi = Logmap(state);
  const Vector3 w = xi.head<3>();
  Vector3 rho = xi.segment<3>(3);
  Vector3 nu = xi.tail<3>();
  
  Matrix3 Qt, Qv;
  const Rot3 R = Rot3::Expmap(w);
  Pose3::ExpmapTranslation(w, rho, Qt);
  Pose3::ExpmapTranslation(w,  nu, Qv);
  const Matrix3 Jw = Rot3::LogmapDerivative(w);
  const Matrix3 Qt2 = -Jw * Qt * Jw;
  const Matrix3 Qv2 = -Jw * Qv * Jw;

  Matrix9 J;
  J <<  Jw, Z_3x3, Z_3x3, 
       Qt2,    Jw, Z_3x3,
       Qv2, Z_3x3,    Jw;
  return J;
}


//------------------------------------------------------------------------------
NavState NavState::ChartAtOrigin::Retract(const Vector9& xi,
                                          ChartJacobian Hxi) {
  return Expmap(xi, Hxi);
}

//------------------------------------------------------------------------------
Vector9 NavState::ChartAtOrigin::Local(const NavState& state,
                                       ChartJacobian Hstate) {
  return Logmap(state, Hstate);
}

//------------------------------------------------------------------------------
NavState NavState::retract(const Vector9& xi, //
    OptionalJacobian<9, 9> H1, OptionalJacobian<9, 9> H2) const {
  Rot3 nRb = R_;
  Point3 n_t = t_, n_v = v_;
  Matrix3 D_bRc_xi, D_R_nRb, D_t_nRb, D_v_nRb;
  const Rot3 bRc = Rot3::Expmap(dR(xi), H2 ? &D_bRc_xi : 0);
  const Rot3 nRc = nRb.compose(bRc, H1 ? &D_R_nRb : 0);
  const Point3 t = n_t + nRb.rotate(dP(xi), H1 ? &D_t_nRb : 0);
  const Point3 v = n_v + nRb.rotate(dV(xi), H1 ? &D_v_nRb : 0);
  if (H1) {
    *H1 << D_R_nRb, Z_3x3, Z_3x3, //
    // Note(frank): the derivative of n_t with respect to xi is nRb
    // We pre-multiply with nRc' to account for NavState::Create
    // Then we make use of the identity nRc' * nRb = bRc'
    nRc.transpose() * D_t_nRb, bRc.transpose(), Z_3x3,
    // Similar reasoning for v:
    nRc.transpose() * D_v_nRb, Z_3x3, bRc.transpose();
  }
  if (H2) {
    *H2 << D_bRc_xi, Z_3x3, Z_3x3, //
    Z_3x3, bRc.transpose(), Z_3x3, //
    Z_3x3, Z_3x3, bRc.transpose();
  }
  return NavState(nRc, t, v);
}

//------------------------------------------------------------------------------
Vector9 NavState::localCoordinates(const NavState& g, //
    OptionalJacobian<9, 9> H1, OptionalJacobian<9, 9> H2) const {
  Matrix3 D_dR_R, D_dt_R, D_dv_R;
  const Rot3 dR = R_.between(g.R_, H1 ? &D_dR_R : 0);
  const Point3 dP = R_.unrotate(g.t_ - t_, H1 ? &D_dt_R : 0);
  const Vector dV = R_.unrotate(g.v_ - v_, H1 ? &D_dv_R : 0);

  Vector9 xi;
  Matrix3 D_xi_R;
  xi << Rot3::Logmap(dR, (H1 || H2) ? &D_xi_R : 0), dP, dV;
  if (H1) {
    *H1 << D_xi_R * D_dR_R, Z_3x3, Z_3x3,  //
        D_dt_R, -I_3x3, Z_3x3,             //
        D_dv_R, Z_3x3, -I_3x3;
  }
  if (H2) {
    *H2 << D_xi_R, Z_3x3, Z_3x3,    //
        Z_3x3, dR.matrix(), Z_3x3,  //
        Z_3x3, Z_3x3, dR.matrix();
  }

  return xi;
}

//------------------------------------------------------------------------------
// sugar for derivative blocks
#define D_R_R(H) (H)->block<3,3>(0,0)
#define D_R_t(H) (H)->block<3,3>(0,3)
#define D_R_v(H) (H)->block<3,3>(0,6)
#define D_t_R(H) (H)->block<3,3>(3,0)
#define D_t_t(H) (H)->block<3,3>(3,3)
#define D_t_v(H) (H)->block<3,3>(3,6)
#define D_v_R(H) (H)->block<3,3>(6,0)
#define D_v_t(H) (H)->block<3,3>(6,3)
#define D_v_v(H) (H)->block<3,3>(6,6)

//------------------------------------------------------------------------------
NavState NavState::update(const Vector3& b_acceleration, const Vector3& b_omega,
    const double dt, OptionalJacobian<9, 9> F, OptionalJacobian<9, 3> G1,
    OptionalJacobian<9, 3> G2) const {

  Vector9 xi;
  Matrix39 D_xiP_state;
  Vector3 b_v = bodyVelocity(F ? &D_xiP_state : 0);
  double dt22 = 0.5 * dt * dt;

  // Integrate on tangent space. TODO(frank): coriolis?
  dR(xi) << dt * b_omega;
  dP(xi) << dt * b_v + dt22 * b_acceleration;
  dV(xi) << dt * b_acceleration;

  // Bring back to manifold
  Matrix9 D_newState_xi;
  NavState newState = retract(xi, F, G1 || G2 ? &D_newState_xi : 0);

  // Derivative wrt state is computed by retract directly
  // However, as dP(xi) also depends on state, we need to add that contribution
  if (F) {
    F->middleRows<3>(3) += dt * D_t_t(F) * D_xiP_state;
  }
  // derivative wrt acceleration
  if (G1) {
    // D_newState_dPxi = D_newState_xi.middleCols<3>(3)
    // D_dPxi_acc = dt22 * I_3x3
    // D_newState_dVxi = D_newState_xi.rightCols<3>()
    // D_dVxi_acc = dt * I_3x3
    // *G2 = D_newState_acc = D_newState_dPxi * D_dPxi_acc + D_newState_dVxi * D_dVxi_acc
    *G1 = D_newState_xi.middleCols<3>(3) * dt22
        + D_newState_xi.rightCols<3>() * dt;
  }
  // derivative wrt omega
  if (G2) {
    // D_newState_dRxi = D_newState_xi.leftCols<3>()
    // D_dRxi_omega = dt * I_3x3
    // *G1 = D_newState_omega = D_newState_dRxi * D_dRxi_omega
    *G2 = D_newState_xi.leftCols<3>() * dt;
  }
  return newState;
}

//------------------------------------------------------------------------------
Vector9 NavState::coriolis(double dt, const Vector3& omega, bool secondOrder,
    OptionalJacobian<9, 9> H) const {
  const Rot3& nRb = R_;
  const Velocity3& n_v = v_;
  const Point3& n_t = t_;

  const double dt2 = dt * dt;
  const Vector3 omega_cross_vel = omega.cross(n_v);

  Vector9 xi;
  Matrix3 D_dP_R, D_dP_dR, D_dV_dR;
  dR(xi) << nRb.unrotate((-dt) * omega, H ? &D_dP_R : 0);
  dP(xi) << nRb.unrotate((-dt2) * omega_cross_vel, H ? &D_dP_dR : 0);
  dV(xi) << nRb.unrotate((-2.0 * dt) * omega_cross_vel, H ? &D_dV_dR : 0);

  if (secondOrder) {
    const Vector3 omega_cross2_t = omega.cross(omega.cross(n_t));
    dP(xi) -= (0.5 * dt2) * omega_cross2_t;
    dV(xi) -= dt * omega_cross2_t;
  }
  if (H) {
    H->setZero();
    const Matrix3 Omega = skewSymmetric(omega);
    const Matrix3 D_cross_state = Omega * R();
    H->setZero();
    D_R_R(H) << D_dP_R;
    D_t_v(H) << nRb.transpose() * (-dt2) * D_cross_state;
    D_v_v(H) << nRb.transpose() * (-2.0 * dt) * D_cross_state;
    D_t_R(H) << D_dP_dR;
    D_v_R(H) << D_dV_dR;

    if (secondOrder) {
      const Matrix3 D_cross2_state = Omega * D_cross_state;
      D_t_t(H) -= (0.5 * dt2) * D_cross2_state;
      D_v_t(H) -= dt * D_cross2_state;
    }
  }
  return xi;
}

//------------------------------------------------------------------------------
Vector9 NavState::correctPIM(const Vector9& pim, double dt,
    const Vector3& n_gravity, const std::optional<Vector3>& omegaCoriolis,
    bool use2ndOrderCoriolis, OptionalJacobian<9, 9> H1,
    OptionalJacobian<9, 9> H2) const {
  const Rot3& nRb = R_;
  const Velocity3& n_v = v_; // derivative is Ri !
  const Point3& n_t = t_; // derivative is Ri !
  const double dt22 = 0.5 * dt * dt;

  Vector9 xi;
  Matrix3 D_dP_Ri, D_dV_Ri, D_dP_Gx, D_dV_Gv;
  Matrix3 D_Gx_nt = Matrix::Zero(3, 3);
  Matrix3 D_Gv_nt = Matrix::Zero(3, 3);
  Matrix3 D_Gx_nv = dt * Matrix::Identity(3, 3);

  Vector3 Gamma_v, Gamma_p;
  if (omegaCoriolis) {
    Vector6 omega_n_gravity;
    omega_n_gravity << -(*omegaCoriolis)*dt, n_gravity*dt;
    Pose3 pose = Pose3::Expmap(omega_n_gravity);
    Rot3 Gamma_R = pose.rotation();
    Matrix3 GOmGt = skewSymmetric(Gamma_R.rotate(*omegaCoriolis));
    Vector3 n_v_prime = n_v + GOmGt * n_t;

    // add initial state contribution
    Gamma_v = pose.translation() + GOmGt * n_t;
    // add velocity and initial position contribution
    Matrix3 Omega = skewSymmetric(*omegaCoriolis);
    double phi = (*omegaCoriolis).norm();
    double c = phi * dt;
    double phi3 = phi * phi * phi;
    double phi4 = phi3 * phi;
    double a = (c * cos(c) - sin(c)) /(phi3);
    double b = (c*c/2 - cos(c) - c*sin(c) + 1) / (phi4);
    Matrix3 mat = dt22 * Matrix::Identity(3, 3);
    mat += a * Omega + b * Omega * Omega;
    Gamma_p =  dt * n_v_prime + mat * n_gravity;

    D_Gv_nt = GOmGt;
    D_Gx_nt = dt * GOmGt;
  } else {
    Gamma_v = dt * n_gravity;
    Gamma_p = dt * n_v  + dt22 * n_gravity;
  }
  dR(xi) = dR(pim);
  dP(xi) = dP(pim)+nRb.unrotate(Gamma_p, H1 ? &D_dP_Ri : 0, H1 ? &D_dP_Gx : 0);
  dV(xi) = dV(pim)+nRb.unrotate(Gamma_v, H1 ? &D_dV_Ri : 0, H1 ? &D_dV_Gv : 0);

  if (H1 || H2) {
    Matrix3 Ri = nRb.matrix();

    if (H1) {
      H1->setZero(); // if coriolis H1 is already initialized
      D_t_R(H1) += D_dP_Ri;
      D_t_v(H1) += D_dP_Gx * D_Gx_nv * Ri;
      D_v_R(H1) += D_dV_Ri;
      D_t_t(H1) += D_dP_Gx * D_Gx_nt * Ri;
      D_v_t(H1) += D_dV_Gv * D_Gv_nt * Ri;
    }
    if (H2) {
      H2->setIdentity();
    }
  }

  return xi;
}
//------------------------------------------------------------------------------

}/// namespace gtsam
