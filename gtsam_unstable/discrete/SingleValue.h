/*
 * SingleValue.h
 * @brief domain constraint
 * @date Feb 6, 2012
 * @author Frank Dellaert
 */

#pragma once

#include <gtsam/discrete/DiscreteKey.h>
#include <gtsam_unstable/discrete/Constraint.h>

namespace gtsam {

/**
 * SingleValue constraint: ensures a variable takes on a certain value.
 * This could of course also be implemented by changing its `Domain`.
 */
class GTSAM_UNSTABLE_EXPORT SingleValue : public Constraint {
  size_t cardinality_;  /// < Number of values
  size_t value_;        ///<  allowed value

  DiscreteKey discreteKey() const {
    return DiscreteKey(keys_[0], cardinality_);
  }

 public:
  typedef std::shared_ptr<SingleValue> shared_ptr;

  /// Construct from key, cardinality, and given value.
  SingleValue(Key key, size_t n, size_t value)
      : Constraint(key), cardinality_(n), value_(value) {}

  /// Construct from DiscreteKey and given value.
  SingleValue(const DiscreteKey& dkey, size_t value)
      : Constraint(dkey.first), cardinality_(dkey.second), value_(value) {}

  // print
  void print(const std::string& s = "", const KeyFormatter& formatter =
                                            DefaultKeyFormatter) const override;

  /// equals
  bool equals(const DiscreteFactor& other, double tol) const override {
    if (!dynamic_cast<const SingleValue*>(&other))
      return false;
    else {
      const SingleValue& f(static_cast<const SingleValue&>(other));
      return (cardinality_ == f.cardinality_) && (value_ == f.value_);
    }
  }

  /// Compute error for each assignment and return as a tree
  AlgebraicDecisionTree<Key> errorTree() const override {
    throw std::runtime_error("SingleValue::error not implemented");
  }

  /// Calculate value
  double evaluate(const Assignment<Key>& values) const override;

  /// Convert into a decisiontree
  DecisionTreeFactor toDecisionTreeFactor() const override;

  /// Multiply into a decisiontree
  DecisionTreeFactor operator*(const DecisionTreeFactor& f) const override;

  /*
   * Ensure Arc-consistency: just sets domain[j] to {value_}.
   * @param j domain to be checked
   * @param (in/out) domains all domains, but only domains->at(j) will be checked.
   * @return true if domains->at(j) was changed, false otherwise.
   */
  bool ensureArcConsistency(Key j, Domains* domains) const override;

  /// Partially apply known values
  Constraint::shared_ptr partiallyApply(const DiscreteValues& values) const override;

  /// Partially apply known values, domain version
  Constraint::shared_ptr partiallyApply(
      const Domains& domains) const override;

  /// Get the number of non-zero values contained in this factor.
  uint64_t nrValues() const override { return 1; };

  DiscreteFactor::shared_ptr sum(size_t nrFrontals) const override {
    return toDecisionTreeFactor().sum(nrFrontals);
  }

  DiscreteFactor::shared_ptr sum(const Ordering& keys) const override {
    return toDecisionTreeFactor().sum(keys);
  }

  DiscreteFactor::shared_ptr max(size_t nrFrontals) const override {
    return toDecisionTreeFactor().max(nrFrontals);
  }

  DiscreteFactor::shared_ptr max(const Ordering& keys) const override {
    return toDecisionTreeFactor().max(keys);
  }
};

}  // namespace gtsam
