/*
 * concepts.h
 *
 * @data Dec 4, 2014
 * @author Mike Bosse
 * @author Frank Dellaert
 */

#pragma once

//#include "manifold.h"
//#include "chart.h"
#include <gtsam/base/Matrix.h>
#include <gtsam/base/Testable.h>
#include <boost/concept_check.hpp>
#include <boost/concept/requires.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_base_of.hpp>

namespace gtsam {

namespace traits {

/**
 * @name Algebraic Structure Traits
 * @brief Associate a unique tag with each of the main GTSAM concepts
 */
//@{
template<typename T> struct structure_category;
//@}

/**
 * @name Algebraic Structure Tags
 * @brief Possible values for traits::structure_category<T>::type
 */
//@{
struct manifold_tag {};
struct group_tag {};
struct lie_group_tag: public manifold_tag, public group_tag {};
struct vector_space_tag: public lie_group_tag {};
//@}

}// namespace traits

namespace manifold {

/** @name Free functions any Manifold needs to define */
//@{
//@}

namespace traits {

/** @name Manifold Traits */
//@{
template<typename Manifold> struct dimension;
template<typename Manifold> struct TangentVector;
template<typename Manifold> struct DefaultChart;
//@}

}// \ namespace traits

/// Check invariants for Manifold type
template<typename T>
BOOST_CONCEPT_REQUIRES(((Testable<T>)),(bool)) //
check_invariants(const T& a, const T& b) {
  typedef typename traits::DefaultChart<T>::type Chart;
  return true;
}

/**
 * Base class for Charts
 * Derived has to implement local and retract as static methods
 */
template <class T, class Derived>
struct Chart {
  typedef T ManifoldType;
  typedef typename traits::TangentVector<T>::type TangentVector;

  // TODO, maybe we need Retract and Local to be unary, or both
  // TOOD, also, this indirection mechanism does not seem to help
  static TangentVector Local(const ManifoldType& p, const ManifoldType& q) {
    return Derived::local(p, q);
  }
  static ManifoldType Retract(const ManifoldType& p, const TangentVector& v) {
    return Derived::retract(p, v);
  }
protected:
  Chart() {
    (void) &Local;
    (void) &Retract;
  } // enforce early instantiation. TODO does not seem to work
};

} // \ namespace manifold

template<typename T>
class IsManifold {
public:
  typedef typename traits::structure_category<T>::type structure_category_tag;
  static const size_t dim = manifold::traits::dimension<T>::value;
  typedef typename manifold::traits::TangentVector<T>::type TangentVector;
  typedef typename manifold::traits::DefaultChart<T>::type DefaultChart;

  BOOST_CONCEPT_USAGE(IsManifold) {
    BOOST_STATIC_ASSERT_MSG(
        (boost::is_base_of<traits::manifold_tag, structure_category_tag>::value),
        "This type's structure_category trait does not assert it as a manifold (or derived)");
    BOOST_STATIC_ASSERT(TangentVector::SizeAtCompileTime == dim);
    BOOST_STATIC_ASSERT_MSG(
        (boost::is_base_of<manifold::Chart<T,DefaultChart>, DefaultChart>::value),
        "This type's DefaultChart does not derive from manifold::Chart, as required");
    // make sure Derived methods in Chart are defined
    v = DefaultChart::local(p,q);
    q = DefaultChart::retract(p,v);
  }
private:
  T p,q;
  TangentVector v;
};

namespace group {

/** @name Free functions any Group needs to define */
//@{
template<typename T> T compose(const T&g, const T& h);
template<typename T> T between(const T&g, const T& h);
template<typename T> T inverse(const T&g);
//@}

namespace traits {

/** @name Group Traits */
//@{
template<typename T> struct identity;
template<typename T> struct flavor;
//@}

/** @name Group Flavor Tags */
//@{
struct additive_tag {
};
struct multiplicative_tag {
};
//@}

}// \ namespace traits

/// Check invariants
template<typename T>
BOOST_CONCEPT_REQUIRES(((Testable<T>)),(bool)) //
check_invariants(const T& a, const T& b, double tol = 1e-9) {
  T e = traits::identity<T>::value;
  return compose(a, inverse(a)).equals(e, tol)
      && between(a, b).equals(compose(inverse(a), b), tol)
      && compose(a, between(a, b)).equals<T>(b, tol);
}
} // \ namespace group

/**
 * Group Concept
 */
template<typename T>
class IsGroup {
public:

  typedef typename traits::structure_category<T>::type structure_category_tag;
  typedef typename group::traits::identity<T>::value_type identity_value_type;
  typedef typename group::traits::flavor<T>::type flavor_tag;

  void operator_usage(group::traits::multiplicative_tag) {
    g = g * h;
  }
  void operator_usage(group::traits::additive_tag) {
    g = g + h;
    g = h - g;
    g = -g;
  }

  BOOST_CONCEPT_USAGE(IsGroup) {
    using group::compose;
    using group::between;
    using group::inverse;
    BOOST_STATIC_ASSERT_MSG(
        (boost::is_base_of<traits::group_tag, structure_category_tag>::value),
        "This type's structure_category trait does not assert it as a group (or derived)");
    e = group::traits::identity<T>::value;
    g = compose(g, h);
    g = between(g, h);
    g = inverse(g);
    operator_usage(flavor);
  }

private:
  flavor_tag flavor;
  T e, g, h;
};

namespace lie_group {

/** @name Free functions any Group needs to define */
//@{
// TODO need Jacobians
//template<typename T> T compose(const T&g, const T& h);
//template<typename T> T between(const T&g, const T& h);
//template<typename T> T inverse(const T&g);
//@}

namespace traits {

/** @name Lie Group Traits */
//@{
//@}

}// \ namespace traits

/// Check invariants
//template<typename T>
//BOOST_CONCEPT_REQUIRES(((Testable<T>)),(bool)) check_invariants(const T& a,
//    const T& b) {
//  bool check_invariants(const V& a, const V& b) {
//    return equal(Chart::retract(a, b), a + b)
//        && equal(Chart::local(a, b), b - a);
//  }
//}
}// \ namespace lie_group

/**
 * Lie Group Concept
 */
template<typename T>
class IsLieGroup: public IsGroup<T>, public IsManifold<T> {
public:

  typedef typename traits::structure_category<T>::type structure_category_tag;

  BOOST_CONCEPT_USAGE(IsLieGroup) {
    BOOST_STATIC_ASSERT_MSG(
        (boost::is_base_of<traits::lie_group_tag, structure_category_tag>::value),
        "This type's trait does not assert it as a Lie group (or derived)");
    // TODO Check with Jacobian
//    using lie_group::compose;
//    using lie_group::between;
//    using lie_group::inverse;
//    g = compose(g, h);
//    g = between(g, h);
//    g = inverse(g);
  }
private:

  T g, h;
};

template<typename T>
class IsVectorSpace: public IsLieGroup<T> {
public:

  typedef typename traits::structure_category<T>::type structure_category_tag;

  BOOST_CONCEPT_USAGE(IsVectorSpace) {
    BOOST_STATIC_ASSERT_MSG(
        (boost::is_base_of<traits::vector_space_tag, structure_category_tag>::value),
        "This type's trait does not assert it as a vector space (or derived)");
    r = p + q;
    r = -p;
    r = p - q;
  }

private:
  T p, q, r;
};

} // namespace gtsam

