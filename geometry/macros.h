#pragma once

namespace nlls {
namespace internal {

#define NLLS_INHERIT_GROUP_PROPERTIES \
    using BaseType::DIM; \
    using BaseType::DOF; \
    using Scalar = typename BaseType::Scalar; \
    using LieGroup = typename BaseType::LieGroup; \
    using Tangent = typename BaseType::Tangent; \
    using Adjoint = typename BaseType::Adjoint; \
    using Vector = typename BaseType::Vector; \
    using Rotation = typename BaseType::Rotation; \
    using Transform = typename BaseType::Transform; \
    using Dvdl = typename BaseType::Dvdl;

#define NLLS_INHERIT_TANGENT_PROPERTIES \
    using BaseType::DIM; \
    using BaseType::DOF; \
    using Scalar = typename BaseType::Scalar; \
    using LieGroup = typename BaseType::LieGroup; \
    using Tangent = typename BaseType::Tangent; \
    using LieAlg = typename BaseType::LieAlg; \
    using Vector = typename BaseType::Vector; \
    using Storage = typename BaseType::Storage;

#define NLLS_EXTRACT_GROUP_STORAGE(type) \
    using Storage = typename internal::LieGroupTraits<type>::Storage;

} // namespace internal
} // namespace nlls