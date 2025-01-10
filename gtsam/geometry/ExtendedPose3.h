/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 *@file  ExtendedPose3.h
 *@brief 3D Extended Pose
 */

// \callgraph
#pragma once

#include <gtsam/config.h>

#include <gtsam/geometry/BearingRange.h>
#include <gtsam/geometry/Point3.h>
#include <gtsam/geometry/Rot3.h>
#include <gtsam/base/Lie.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <string>

namespace gtsam {

/**
 * A 3D extended pose (R,v,p) : (Rot3,Point3,Point3)
 * @addtogroup geometry
 * \nosubgrouping
 */
class GTSAM_EXPORT ExtendedPose3: public LieGroup<ExtendedPose3, 9> {
public:

  /** Pose Concept requirements */
  typedef Rot3 Rotation;
  typedef Point3 Velocity;
  typedef Point3 Position;
  typedef Position Translation;

private:

  Rot3 R_; ///< Rotation gRp, between global and pose frame
  Point3 v_; ///< Velocity gVp, from global origin to pose frame origin
  Point3 p_; ///< Position gPp, from global origin to pose frame origin

public:

  /// @name Standard Constructors
  /// @{

  /** Default constructor is origin */
 ExtendedPose3() : R_(traits<Rot3>::Identity()), v_(traits<Point3>::Identity()), p_(traits<Point3>::Identity()) {}

  /** Copy constructor */
  ExtendedPose3(const ExtendedPose3& pose) :
      R_(pose.R_), v_(pose.v_), p_(pose.p_) {
  }

  /** Construct from R,v,p */
  ExtendedPose3(const Rot3& R, const Point3& v, const Point3& p) :
      R_(R), v_(v), p_(p) {
  }

  //explicit Pose3(const Pose2& pose2);

  /** Constructor from 5*5 matrix */
  ExtendedPose3(const Matrix &T) :
      R_(T(0, 0), T(0, 1), T(0, 2), T(1, 0), T(1, 1), T(1, 2), T(2, 0), T(2, 1),
          T(2, 2)), v_(T(0, 3), T(1, 3), T(2, 3)), p_(T(0, 4), T(1, 4), T(2, 4)) {
  }

  /// Named constructor with derivatives
  static ExtendedPose3 Create(const Rot3& R, const Point3& v, const Point3& p,
                      OptionalJacobian<9, 3> HR = {},
                      OptionalJacobian<9, 3> Hv = {},
                      OptionalJacobian<9, 3> Hp = {}) {
    if (HR) *HR << I_3x3, Z_3x3, Z_3x3;
    if (Hv) *Hv << Z_3x3, R.transpose(), Z_3x3;
    if (Hp) *Hp << Z_3x3, Z_3x3, R.transpose();
    return ExtendedPose3(R, v, p);
  }

  /**
   *  Create Pose3 by aligning two point pairs
   *  A pose aTb is estimated between pairs (a_point, b_point) such that a_point = aTb * b_point
   *  Note this allows for noise on the points but in that case the mapping will not be exact.
   */
  static std::optional<ExtendedPose3> Align(const std::vector<Point3Pair>& abPointPairs) {
    const size_t n = abPointPairs.size();
    if (n < 3)
      return {};  // we need at least three pairs

    // calculate centroids
    Point3 aCentroid(0,0,0), bCentroid(0,0,0);
    for(const Point3Pair& abPair: abPointPairs) {
      aCentroid += abPair.first;
      bCentroid += abPair.second;
    }
    double f = 1.0 / n;
    aCentroid *= f;
    bCentroid *= f;

    // Add to form H matrix
    Matrix3 H = Z_3x3;
    for(const Point3Pair& abPair: abPointPairs) {
      Point3 da = abPair.first - aCentroid;
      Point3 db = abPair.second - bCentroid;
      H += da * db.transpose();
    }

    // ClosestTo finds rotation matrix closest to H in Frobenius sense
    Rot3 aRb = Rot3::ClosestTo(H);
    Point3 aTb = Point3(aCentroid) - aRb * Point3(bCentroid);
    Point3 v;
    return ExtendedPose3(aRb, v, aTb);
  }

  /// @}
  /// @name Testable
  /// @{

  /// print with optional string
  void print(const std::string& s = "") const {
    std::cout << s;
    R_.print("R:\n");
    std::cout << "v:" << v_ << ";" << std::endl;
    std::cout << "p:" << p_ << ";" << std::endl;
  }

  /// assert equality up to a tolerance
  bool equals(const ExtendedPose3& pose, double tol = 1e-9) const {
    return R_.equals(pose.R_, tol) && traits<Point3>::Equals(v_, pose.v_, tol) && traits<Point3>::Equals(p_, pose.p_, tol);
  }

  /// @}
  /// @name Group
  /// @{

  /// identity for group operation
  static ExtendedPose3 Identity() {
    return ExtendedPose3();
  }

  /// inverse transformation without derivative
  ExtendedPose3 inverse() const {
    Rot3 Rt = R_.inverse();
    return ExtendedPose3(Rt, Rt * (-v_), Rt * (-p_));
  }

  using LieGroup<ExtendedPose3, 9>::inverse; // version with derivative

  /// compose syntactic sugar
  ExtendedPose3 operator*(const ExtendedPose3& T) const {
    return ExtendedPose3(R_ * T.R_, v_ + R_ * T.v_, p_ + R_ * T.p_);
  }

  /// @}
  /// @name Lie Group
  /// @{

  /// Exponential map at identity - create a rotation from canonical coordinates \f$ [R_x,R_y,R_z,V_x,V_y,V_z,P_x,P_y,P_z] \f$
  static ExtendedPose3 Expmap(const Vector9& xi, OptionalJacobian<9, 9> Hxi = {}) {
    if (Hxi) *Hxi = ExpmapDerivative(xi);

    Vector3 phi(xi(0), xi(1), xi(2)), nu(xi(3), xi(4), xi(5)), rho(xi(6), xi(7), xi(8));

    Rot3 R = Rot3::Expmap(phi);
  
    double phi_norm = phi.norm();
    if (phi_norm < 1e-8)
      return ExtendedPose3(Rot3(), Point3(nu), Point3(rho));
    else {
      Matrix W = skewSymmetric(phi/phi_norm);
      Matrix A = I_3x3 + ((1 - cos(phi_norm)) / phi_norm) * W + ((phi_norm - sin(phi_norm)) / phi_norm) * (W * W);
      return ExtendedPose3(Rot3::Expmap(phi), Point3(A * nu), Point3(A * rho));
    }
  }

  /// Log map at identity - return the canonical coordinates \f$ [R_x,R_y,R_z,V_x,V_y,V_z,P_x,P_y,P_z] \f$ of this rotation
  static Vector9 Logmap(const ExtendedPose3& pose, OptionalJacobian<9, 9> Hpose = {}) {
    if (Hpose) *Hpose = LogmapDerivative(pose);
    const Vector3 phi = Rot3::Logmap(pose.rotation());
    const Vector3 v = pose.velocity();
    const Vector3 p = pose.position();
    const double t = phi.norm();
    if (t < 1e-8) {
      Vector9 log;
      log << phi, v, p;
      return log;
    } else {
      const Matrix3 W = skewSymmetric(phi / t);

      const double Tan = tan(0.5 * t);
      const Vector3 Wp = W * p;
      const Vector3 Wv = W * v;
      const Vector3 nu = v - (0.5 * t) * Wv + (1 - t / (2. * Tan)) * (W * Wv);
      const Vector3 rho = p - (0.5 * t) * Wp + (1 - t / (2. * Tan)) * (W * Wp);
      Vector9 log;
      log << phi, nu, rho;
      return log;
    }
  }

  /**
   * Calculate Adjoint map, transforming a twist in the this pose's (i.e, body) frame to the world spatial frame
   * Ad_pose is 9*9 matrix that when applied to twist xi \f$ [R_x,R_y,R_z,V_x,V_y,V_z,P_x,P_y,P_z] \f$, returns Ad_pose(xi)
   */
  Matrix9 AdjointMap() const {
    const Matrix3 R = R_.matrix();
    Matrix3 A = skewSymmetric(v_.x(), v_.y(), v_.z()) * R;
    Matrix3 B = skewSymmetric(p_.x(), p_.y(), p_.z()) * R;
    Matrix9 adj;
    adj << R, Z_3x3, Z_3x3, 
          A, R, Z_3x3,
          B, Z_3x3, R;
    return adj;
  }

  /**
   * Apply this pose's AdjointMap Ad_g to a twist \f$ \xi_b \f$, i.e. a body-fixed velocity, transforming it to the spatial frame
   * \f$ \xi^s = g*\xi^b*g^{-1} = Ad_g * \xi^b \f$
   */
  Vector9 Adjoint(const Vector9& xi_b) const {
    return AdjointMap() * xi_b;
  }

  /**
   * Compute the ad operator
   */
  static Matrix9 adjointMap(const Vector9 &xi) {
    Matrix3 w_hat = skewSymmetric(xi(0), xi(1), xi(2));
    Matrix3 v_hat = skewSymmetric(xi(3), xi(4), xi(5));
    Matrix3 a_hat = skewSymmetric(xi(6), xi(7), xi(8));
    Matrix9 adj;
    adj << w_hat, Z_3x3, Z_3x3,
          v_hat, w_hat, Z_3x3,
          a_hat, Z_3x3, w_hat;
    return adj;
  }

  /**
   * Action of the adjointMap on a Lie-algebra vector y, with optional derivatives
   */
  static Vector9 adjoint(const Vector9 &xi, const Vector9 &y,
      OptionalJacobian<9, 9> Hxi = {}) {
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
    return adjointMap(xi) * y;
  }

  /**
   * The dual version of adjoint action, acting on the dual space of the Lie-algebra vector space.
   */
  static Vector9 adjointTranspose(const Vector9& xi, const Vector9& y,
      OptionalJacobian<9, 9> Hxi = {}) {
    if (Hxi) {
      Hxi->setZero();
      for (int i = 0; i < 9; ++i) {
        Vector9 dxi;
        dxi.setZero();
        dxi(i) = 1.0;
        Matrix9 GTi = adjointMap(dxi).transpose();
        Hxi->col(i) = GTi * y;
      }
    }
    return adjointMap(xi).transpose() * y;
  }

  /// Derivative of Expmap
  static Matrix9 ExpmapDerivative(const Vector9& xi) {
    const Vector3 w = xi.head<3>();
    const Matrix3 Jw = Rot3::ExpmapDerivative(w);
    const Matrix63 Q = computeQforExpmapDerivative(xi);
    const Matrix3 Qv =  Q.topRows<3>();
    const Matrix3 Qp =  Q.bottomRows<3>();
    Matrix9 J;
    J << Jw, Z_3x3, Z_3x3,
      Qv, Jw, Z_3x3,
      Qp, Z_3x3, Jw;
    return J;
  }

  /// Derivative of Logmap
  static Matrix9 LogmapDerivative(const ExtendedPose3& pose) {
    const Vector9 xi = Logmap(pose);
    const Vector3 w = xi.head<3>();
    const Matrix3 Jw = Rot3::LogmapDerivative(w);
    const Matrix63 Q = computeQforExpmapDerivative(xi);
    const Matrix3 Qv =  Q.topRows<3>();
    const Matrix3 Qp =  Q.bottomRows<3>();
    const Matrix3 Qv2 = -Jw*Qv*Jw;
    const Matrix3 Qp2 = -Jw*Qp*Jw;
    Matrix9 J;

    J << Jw, Z_3x3, Z_3x3,
      Qv2, Jw, Z_3x3,
      Qp2, Z_3x3, Jw;
    return J;
  }

  // Chart at origin, depends on compile-time flag GTSAM_POSE3_EXPMAP
  struct ChartAtOrigin {
    static ExtendedPose3 Retract(const Vector9& xi, ChartJacobian Hxi = {}) {
#ifdef GTSAM_POSE3_EXPMAP
      return Expmap(xi, Hxi);
#else

#endif
    }
    static Vector9 Local(const ExtendedPose3& pose, ChartJacobian Hpose = {}) {
#ifdef GTSAM_POSE3_EXPMAP
      return Logmap(pose, Hpose);
#else

#endif
    }
  };

  Vector9 boxminus(const ExtendedPose3& g) const {
    // Matrix3 D_dR_R, D_dt_R, D_dv_R;
    // const Rot3 dR = R_.between(g.R_, H1 ? &D_dR_R : 0);
    // const Point3 dt = R_.unrotate(g.p_ - p_, H1 ? &D_dt_R : 0);
    // const Vector dv = R_.unrotate(g.v_ - v_, H1 ? &D_dv_R : 0);

    // Vector9 xi;
    // Matrix3 D_xi_R;
    // xi << Rot3::Logmap(dR, (H1 || H2) ? &D_xi_R : 0), dt, dv;
    // if (H1) {
    //   *H1 << D_xi_R * D_dR_R, Z_3x3, Z_3x3, //
    //   D_dt_R, -I_3x3, Z_3x3, //
    //   D_dv_R, Z_3x3, -I_3x3;
    // }
    // if (H2) {
    //   *H2 << D_xi_R, Z_3x3, Z_3x3, //
    //   Z_3x3, dR.matrix(), Z_3x3, //
    //   Z_3x3, Z_3x3, dR.matrix();
    // }
    Matrix3 D_dR_R, D_dt_R, D_dv_R;
    const Rot3 dR = g.R_;
    const Point3 dt = g.p_;
    const Vector dv = g.v_;

    Vector9 xi;
    xi << Rot3::Logmap(dR), dt, dv;
    return xi;
  }

  /**
   * wedge for ExtendedPose3:
   * @param xi 9-dim twist (omega,nu,rho) 
   * @return 5*5 element of Lie algebra
   */
  static Matrix wedge(double phix, double phiy, double phiz, double nux, double nuy, double nuz, double rhox, double rhoy, double rhoz) {
    return (Matrix(5, 5) << 0., -phiz, phiy, nux, rhox, phiz, 0., -phix, nuy, rhoy, -phiy, phix, 0., nuz, rhoz, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.).finished();
  }

  /**
   * Compute the 6x3 bottom-left block Qs of the SE_2(3) Expmap derivative matrix
   */
  static Matrix63 computeQforExpmapDerivative(const Vector9& xi) {
    const auto omega = xi.head<3>();
    const auto nu = xi.segment<3>(3);
    const auto rho = xi.tail<3>();
    const Matrix3 V = skewSymmetric(nu);
    const Matrix3 P = skewSymmetric(rho);
    const Matrix3 W = skewSymmetric(omega);

    Matrix3 Qv, Qp;
    Matrix63 Q;

#ifdef NUMERICAL_EXPMAP_DERIV

#else
    // The closed-form formula in Barfoot14tro eq. (102)
    double phi = omega.norm();
    if (std::abs(phi)>1e-5) {
      const double sinPhi = sin(phi), cosPhi = cos(phi);
      const double phi2 = phi * phi, phi3 = phi2 * phi, phi4 = phi3 * phi, phi5 = phi4 * phi;
      // Invert the sign of odd-order terms to have the right Jacobian
      Qv = -0.5*V + (phi-sinPhi)/phi3*(W*V + V*W - W*V*W)
              + (1-phi2/2-cosPhi)/phi4*(W*W*V + V*W*W - 3*W*V*W)
              - 0.5*((1-phi2/2-cosPhi)/phi4 - 3*(phi-sinPhi-phi3/6.)/phi5)*(W*V*W*W + W*W*V*W);
      Qp = -0.5*P + (phi-sinPhi)/phi3*(W*P + P*W - W*P*W)
              + (1-phi2/2-cosPhi)/phi4*(W*W*P + P*W*W - 3*W*P*W)
              - 0.5*((1-phi2/2-cosPhi)/phi4 - 3*(phi-sinPhi-phi3/6.)/phi5)*(W*P*W*W + W*W*P*W);
    }
    else {
      Qv = -0.5*V + 1./6.*(W*V + V*W - W*V*W)
          + 1./24.*(W*W*V + V*W*W - 3*W*V*W)
          - 0.5*(1./24. + 3./120.)*(W*V*W*W + W*W*V*W);
      Qp = -0.5*P + 1./6.*(W*P + P*W - W*P*W)
          + 1./24.*(W*W*P + P*W*W - 3*W*P*W)
          - 0.5*(1./24. + 3./120.)*(W*P*W*W + W*W*P*W);
    }
#endif
    Q << Qv, Qp;
    return Q;
  }

  /// @}
  /// @name Group Action on Point3
  /// @{

  /**
   * @brief takes point in Pose coordinates and transforms it to world coordinates
   * @param point point in Pose coordinates
   * @param Hself optional 3*6 Jacobian wrpt this pose
   * @param Hpoint optional 3*3 Jacobian wrpt point
   * @return point in world coordinates
   */
  Point3 transformFrom(const Point3& point,
                       OptionalJacobian<3, 9> Hself = {},
                       OptionalJacobian<3, 3> Hpoint = {}) const {
    // Only get matrix once, to avoid multiple allocations,
    // as well as multiple conversions in the Quaternion case
    const Matrix3 R = R_.matrix();
    if (Hself) {
      Hself->leftCols<3>() =
          R * skewSymmetric(-point.x(), -point.y(), -point.z());
      Hself->rightCols<3>() = R;
    }
    if (Hpoint) {
      *Hpoint = R;
    }
    return R_ * point + p_;
  }

  /** syntactic sugar for transformFrom */
  inline Point3 operator*(const Point3& point) const {
    return transformFrom(point);
  }

  /**
   * @brief takes point in world coordinates and transforms it to Pose coordinates
   * @param point point in world coordinates
   * @param Hself optional 3*6 Jacobian wrpt this pose
   * @param Hpoint optional 3*3 Jacobian wrpt point
   * @return point in Pose coordinates
   */
  Point3 transformTo(const Point3& point,
                     OptionalJacobian<3, 9> Hself = {},
                     OptionalJacobian<3, 3> Hpoint = {}) const {
    // Only get transpose once, to avoid multiple allocations,
    // as well as multiple conversions in the Quaternion case
    const Matrix3 Rt = R_.transpose();
    const Point3 q(Rt * (point - p_));
    if (Hself) {
      const double wx = q.x(), wy = q.y(), wz = q.z();
      (*Hself) << 0.0, -wz, +wy, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, +wz, 0.0, -wx,
          0.0, 0.0, 0.0, 0.0, -1.0, 0.0, -wy, +wx, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          -1.0;
    }
    if (Hpoint) {
      *Hpoint = Rt;
    }
    return q;
  }

  /// @}
  /// @name Standard Interface
  /// @{

  /// get rotation
  const Rot3& rotation(OptionalJacobian<3, 9> Hself = {}) const {
    if (Hself) {
      *Hself << I_3x3, Z_3x3, Z_3x3;
    }
    return R_;
  }

  /// get translation
  const Point3& velocity(OptionalJacobian<3, 9> Hself = {}) const {
    if (Hself) {
      *Hself << Z_3x3, rotation().matrix(), Z_3x3;
    }
    return v_;
  }

  /// get position
  const Point3& position(OptionalJacobian<3, 9> Hself = {}) const {
    if (Hself) *Hself << Z_3x3, Z_3x3, rotation().matrix();
    return p_;
  }

  const Point3& translation(OptionalJacobian<3, 9> Hself = {}) const {
    if (Hself) *Hself << Z_3x3, Z_3x3, rotation().matrix();
    return p_;
  }

  /// get x
  double x() const {
    return p_.x();
  }

  /// get y
  double y() const {
    return p_.y();
  }

  /// get z
  double z() const {
    return p_.z();
  }

  /** convert to 5*5 matrix */
  Matrix5 matrix() const {
    Matrix5 mat;
    mat << R_.matrix(), v_, p_, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        1.0;
    return mat;
  }

  /**
   * Assuming self == wTa, takes a pose aTb in local coordinates
   * and transforms it to world coordinates wTb = wTa * aTb.
   * This is identical to compose.
   */
  ExtendedPose3 transformPoseFrom(
      const ExtendedPose3& aTb, OptionalJacobian<9, 9> Hself = {},
      OptionalJacobian<9, 9> HaTb = {}) const {
    const ExtendedPose3& wTa = *this;
    return wTa.compose(aTb, Hself, HaTb);
  }

  /**
   *  Assuming self == wTa, takes a pose wTb in world coordinates
   * and transforms it to local coordinates aTb = inv(wTa) * wTb
   */
  ExtendedPose3 transformPoseTo(
      const ExtendedPose3& wTb, OptionalJacobian<9, 9> Hself = {},
      OptionalJacobian<9, 9> HwTb = {}) const {
    if (Hself) *Hself = -wTb.inverse().AdjointMap() * AdjointMap();
    if (HwTb) *HwTb = I_9x9;
    const ExtendedPose3& wTa = *this;
    return wTa.inverse() * wTb;
  }

  /**
   * Calculate range to a landmark
   * @param point 3D location of landmark
   * @return range (double)
   */
  double range(const Point3& point, OptionalJacobian<1, 9> Hself = {},
               OptionalJacobian<1, 3> Hpoint = {}) const {
    Matrix39 D_local_pose;
    Matrix3 D_local_point;
    Point3 local = transformTo(point, Hself ? &D_local_pose : 0,
                               Hpoint ? &D_local_point : 0);
    if (!Hself && !Hpoint) {
      return local.norm();
    } else {
      Matrix13 D_r_local;
      const double r = norm3(local, D_r_local);
      if (Hself) *Hself = D_r_local * D_local_pose;
      if (Hpoint) *Hpoint = D_r_local * D_local_point;
      return r;
    }
  }

  /**
   * Calculate range to another pose
   * @param pose Other SO(3) pose
   * @return range (double)
   */
  double range(const ExtendedPose3& pose,
               OptionalJacobian<1, 9> Hself = {},
               OptionalJacobian<1, 9> Hpose = {}) const {
    Matrix13 D_local_point;
    double r = range(pose.translation(), Hself, Hpose ? &D_local_point : 0);
    if (Hpose)
      *Hpose << Matrix13::Zero(), D_local_point * pose.rotation().matrix();
    return r;
  }

  /**
   * Calculate bearing to a landmark
   * @param point 3D location of landmark
   * @return bearing (Unit3)
   */
  Unit3 bearing(const Point3& point, OptionalJacobian<2, 9> Hself = {},
                OptionalJacobian<2, 3> Hpoint = {}) const {
    Matrix39 D_local_pose;
    Matrix3 D_local_point;
    Point3 local = transformTo(point, Hself ? &D_local_pose : 0,
                               Hpoint ? &D_local_point : 0);
    if (!Hself && !Hpoint) {
      return Unit3(local);
    } else {
      Matrix23 D_b_local;
      Unit3 b = Unit3::FromPoint3(local, D_b_local);
      if (Hself) *Hself = D_b_local * D_local_pose;
      if (Hpoint) *Hpoint = D_b_local * D_local_point;
      return b;
    }
  }

  /**
   * Calculate bearing to another pose
   * @param other 3D location and orientation of other body. The orientation
   * information is ignored.
   * @return bearing (Unit3)
   */
  Unit3 bearing(const ExtendedPose3& pose,
                OptionalJacobian<2, 9> Hself = {},
                OptionalJacobian<2, 9> Hpose = {}) const {
    if (Hpose) {
      Hpose->setZero();
      return bearing(pose.translation(), Hself, Hpose.cols<3>(3));
    }
    return bearing(pose.translation(), Hself, {});
  }

  /// @}
  /// @name Advanced Interface
  /// @{

  /**
   * Return the start and end indices (inclusive) of the translation component of the
   * exponential map parameterization
   * @return a pair of [start, end] indices into the tangent space vector
   */
  inline static std::pair<size_t, size_t> velocityInterval() {
    return std::make_pair(3, 5);
  }

  /**
   * Return the start and end indices (inclusive) of the translation component of the
   * exponential map parameterization
   * @return a pair of [start, end] indices into the tangent space vector
   */
  inline static std::pair<size_t, size_t> positionInterval() {
    return std::make_pair(6, 8);
  }
  inline static std::pair<size_t, size_t> translationInterval() {
    return std::make_pair(6, 8);
  }

  /**
   * Return the start and end indices (inclusive) of the rotation component of the
   * exponential map parameterization
   * @return a pair of [start, end] indices into the tangent space vector
   */
  static std::pair<size_t, size_t> rotationInterval() {
    return std::make_pair(0, 2);
  }

  /// Output stream operator
  GTSAM_EXPORT
  friend std::ostream& operator<<(std::ostream& os, const ExtendedPose3& pose) {
    os << pose.rotation();
    const Point3& v = pose.velocity();
    const Point3& p = pose.position();
    os << "v:[" << v.x() << ", " << v.y() << ", " << v.z() << "];\n";
    os << "p:[" << p.x() << ", " << p.y() << ", " << p.z() << "];\n";
    return os;
  }

 private:
  /** Serialization function */
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int /*version*/) {
    ar & BOOST_SERIALIZATION_NVP(R_);
    ar & BOOST_SERIALIZATION_NVP(v_);
    ar & BOOST_SERIALIZATION_NVP(p_);
  }
  /// @}

#ifdef GTSAM_USE_QUATERNIONS
  // Align if we are using Quaternions
  public:
    GTSAM_MAKE_ALIGNED_OPERATOR_NEW
#endif
};
// ExtendedPose3 class

/**
 * wedge for ExtendedPose3:
 * @param xi 9-dim twist (omega,nu,rho)
 * @return 5*5 element of Lie algebra
 */
template<>
inline Matrix wedge<ExtendedPose3>(const Vector& xi) {
  return ExtendedPose3::wedge(xi(0), xi(1), xi(2), xi(3), xi(4), xi(5), xi(6), xi(7), xi(8));
}



template <>
struct traits<ExtendedPose3> : public internal::LieGroup<ExtendedPose3> {};

template <>
struct traits<const ExtendedPose3> : public internal::LieGroup<ExtendedPose3> {};

// bearing and range traits, used in RangeFactor
template <>
struct Bearing<ExtendedPose3, Point3> : HasBearing<ExtendedPose3, Point3, Unit3> {};

template<>
struct Bearing<ExtendedPose3, ExtendedPose3> : HasBearing<ExtendedPose3, ExtendedPose3, Unit3> {};

template <typename T>
struct Range<ExtendedPose3, T> : HasRange<ExtendedPose3, T, double> {};
}  // namespace gtsam
