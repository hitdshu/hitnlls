#pragma once

#include <list>

#include "../matrix/dense.h"
#include "../ils/vertex.h"

namespace nlls {

class CameraBase : public VertexBase {
public:
    virtual int Dof() const = 0;
    virtual void Reset() = 0;
    virtual void Update(const float *inc) = 0;
    virtual void Push() = 0;
    virtual void Pop() = 0;

    virtual Vector2f Ray2Pixel(const Vector3f &pt3d, Matrix23f *dp = nullptr) const = 0;
    virtual Vector3f Pixel2Ray(const Vector2f &pt2d) const = 0;
    virtual int Width() const = 0;
    virtual int Height() const = 0;
};

NLLS_REGISTER_REGISTER(CameraBase)
#define NLLS_REGISTER_CAMERA(name) \
    NLLS_REGISTER_VERTEX(name) \
    NLLS_REGISTER_CLASS(CameraBase, name)

template <int N, class D>
class Camera : public CameraBase {
public:
    static constexpr int DOF = N + D::DOF;
    using CameraParamVector = Matrix<float, N, 1>;
    using TotalParamVector = Matrix<float, DOF, 1>;

    explicit Camera() { width_ = -1; height_ = -1; data_.SetZero(); }
    explicit Camera(int width, int height) : width_(width), height_(height) { data_.SetZero(); }
    explicit Camera(int width, int height, const TotalParamVector &tp) : width_(width), height_(height), data_(tp) {}

    NLLS_INLINE void Init(int width, int height, const TotalParamVector &tp) { width_ = width; height_ = height; data_ = tp; }
    NLLS_INLINE void Init(int width, int height) { width_ = width; height_ = height; }
    NLLS_INLINE void Init(const TotalParamVector &tp) { data_ = tp; }
    NLLS_INLINE TotalParamVector &ToVector() { return data_; }
    NLLS_INLINE const TotalParamVector &ToVector() const { return data_; }
    NLLS_INLINE float *Data() { return data_.Data(); }
    NLLS_INLINE const float *Data() const { return data_.Data(); }
    NLLS_INLINE float *CameraData() { return Data(); }
    NLLS_INLINE const float *CameraData() const { return Data(); }
    NLLS_INLINE float *DistortionData() { return Data() + N; }
    NLLS_INLINE const float *DistortionData() const { return Data() + N; }

    virtual int Dof() const override final { return DOF; }
    virtual void Reset() override final { data_.SetZero(); }
    virtual void Update(const float *inc) override final { for (int i = 0; i < DOF; ++i) { data_[i] += inc[i]; } }
    virtual void Push() override final { stack_.push_back(data_); }
    virtual void Pop() override final { data_ = stack_.back(); stack_.pop_back(); }

    virtual int Width() const override final { return width_; }
    virtual int Height() const override final { return height_; }

    virtual Vector2f Ray2Pixel(const Vector3f &pt3d, Matrix23f *dp = nullptr) const { return Project(pt3d, dp, nullptr); }
    virtual Vector3f Pixel2Ray(const Vector2f &pt2d) const { return UnProject(pt2d).Normalized(); }

    virtual Vector2f Project(const Vector3f &pt3d, Matrix23f *dp = nullptr, Matrix<float, 2, DOF> *dc = nullptr) const = 0;
    virtual Vector3f UnProject(const Vector2f &pt2d) const = 0;

private:
    int width_;
    int height_;
    TotalParamVector data_;
    std::list<TotalParamVector> stack_;
};

} // namespace nlls