#pragma once

namespace hitnlls {
namespace geometry {

#define HITNLLS_INHERIT_GROUP_PROPERTIES \
    using BaseType::DIM; \
    using BaseType::DOF; \
    using Scalar = typename BaseType::Scalar; \
    using Storage = typename BaseType::Storage; \
    using Vector = typename BaseType::Vector; \
    using Adjoint = typename BaseType::Adjoint; \
    using LieGroup = typename BaseType::LieGroup; \
    using Tangent = typename BaseType::Tangent;

#define HITNLLS_INHERIT_TANGENT_PROPERTIES \
    using BaseType::DIM; \
    using BaseType::DOF; \
    using Scalar = typename BaseType::Scalar; \
    using Storage = typename BaseType::Storage; \
    using LieGroup = typename BaseType::LieGroup; \
    using Tangent = typename BaseType::Tangent;

#define PI (3.14159265358979323846)

} // namespace geometry
} // namespace hitnlls