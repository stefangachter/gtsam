/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 *  @file  LieGroupPreintegration.h
 *  @author Luca Carlone
 *  @author Stephen Williams
 *  @author Richard Roberts
 *  @author Vadim Indelman
 *  @author David Jensen
 *  @author Frank Dellaert
 *  @author Martin Brossard
 **/

#pragma once

#include <gtsam/navigation/NavState.h>
#include <gtsam/navigation/PreintegrationBase.h>

namespace gtsam {

/**
 * IMU pre-integration on SE_2(3) Lie group.
 */
class GTSAM_EXPORT LieGroupPreintegration : public PreintegrationBase {
 protected:

  /**
   * Pre-integrated navigation state, from frame i to frame j
   * Note: relative position does not take into account velocity at time i, see deltap+, in [2]
   * Note: velocity is now also in frame i, as opposed to deltaVij in [2]
   */
  NavState deltaXij_;
  Matrix3 delRdelBiasOmega_; ///< Jacobian of preintegrated rotation w.r.t. angular rate bias
  Matrix3 delPdelBiasAcc_;   ///< Jacobian of preintegrated position w.r.t. acceleration bias
  Matrix3 delPdelBiasOmega_; ///< Jacobian of preintegrated position w.r.t. angular rate bias
  Matrix3 delVdelBiasAcc_;   ///< Jacobian of preintegrated velocity w.r.t. acceleration bias
  Matrix3 delVdelBiasOmega_; ///< Jacobian of preintegrated velocity w.r.t. angular rate bias

  /// Default constructor for serialization
  LieGroupPreintegration() {
    resetIntegration();
  }

public:
  /// @name Constructors
  /// @{

  /**
   *  Constructor, initializes the variables in the base class
   *  @param p    Parameters, typically fixed in a single application
   *  @param bias Current estimate of acceleration and rotation rate biases
   */
  LieGroupPreintegration(const boost::shared_ptr<Params>& p,
      const imuBias::ConstantBias& biasHat = imuBias::ConstantBias());

  /// @}

  /// @name Basic utilities
  /// @{
  /// Re-initialize PreintegratedMeasurements
  void resetIntegration() override;

  /// @}

  /// @name Instance variables access
  /// @{
  NavState deltaXij() const override { return deltaXij_; }
  Rot3     deltaRij() const override { return deltaXij_.attitude(); }
  Vector3  deltaPij() const override { return deltaXij_.position(); }
  Vector3  deltaVij() const override { return deltaXij_.velocity(); }

  Matrix3  delRdelBiasOmega() const { return delRdelBiasOmega_; }
  Matrix3  delPdelBiasAcc() const { return delPdelBiasAcc_; }
  Matrix3  delPdelBiasOmega() const { return delPdelBiasOmega_; }
  Matrix3  delVdelBiasAcc() const { return delVdelBiasAcc_; }
  Matrix3  delVdelBiasOmega() const { return delVdelBiasOmega_; }

  /// @name Testable
  /// @{
  bool equals(const LieGroupPreintegration& other, double tol) const;
  /// @}

  /// @name Main functionality
  /// @{

  /// Update preintegrated measurements and get derivatives
  /// It takes measured quantities in the j frame
  /// Modifies preintegrated quantities in place after correcting for bias and possibly sensor pose
  void update(const Vector3& measuredAcc, const Vector3& measuredOmega, const double dt,
              Matrix9* A, Matrix93* B, Matrix93* C) override;

  // Integrate forward in time given angular velocity and acceleration in body frame
  /// Uses second order integration for position, returns derivatives except dt.
  void updateFactor(const Vector3& b_acceleration, const Vector3& b_omega,
      const double dt, OptionalJacobian<9, 9> F, OptionalJacobian<9, 3> G1,
      OptionalJacobian<9, 3> G2);

  /// Given the estimate of the bias, return a NavState tangent vector
  /// summarizing the preintegrated IMU measurements so far
  Vector9 biasCorrectedDelta(const imuBias::ConstantBias& bias_i,
      OptionalJacobian<9, 6> H = boost::none) const override;

  /** Dummy clone for MATLAB */
  virtual boost::shared_ptr<LieGroupPreintegration> clone() const {
    return boost::shared_ptr<LieGroupPreintegration>();
  }

  /// @}

private:
  /** Serialization function */
  friend class boost::serialization::access;
  template<class ARCHIVE>
  void serialize(ARCHIVE & ar, const unsigned int /*version*/) {
    namespace bs = ::boost::serialization;
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(PreintegrationBase);
    ar & BOOST_SERIALIZATION_NVP(deltaXij_);
    ar & BOOST_SERIALIZATION_NVP(delRdelBiasOmega_);
    ar & BOOST_SERIALIZATION_NVP(delPdelBiasAcc_);
    ar & BOOST_SERIALIZATION_NVP(delPdelBiasOmega_);
    ar & BOOST_SERIALIZATION_NVP(delVdelBiasAcc_);
    ar & BOOST_SERIALIZATION_NVP(delVdelBiasOmega_);
  }
};

} /// namespace gtsam