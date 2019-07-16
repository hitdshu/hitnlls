#pragma once

#include "matrix/matrixx.h"
#include "common/register.h"

namespace hitnlls {
namespace node {

class BaseNode {
public:
    BaseNode(int ndim = -1) { fixed_ = false; id_ = -1; ndim_ = ndim; marginalized_ = false; }
    virtual ~BaseNode() = default;

    virtual void SetId(int id) final { id_ = id; }
    virtual void SetFixed(bool fixed) final { fixed_ = fixed; }
    virtual void SetMarginalized(bool marg) final { marginalized_ = marg; }
    virtual int GetId() const final { return id_; }
    virtual int GetDim() const final { return ndim_; }
    virtual bool GetFixed() const final { return fixed_; }
    virtual bool GetMarginalized() const final { return marginalized_; }

    virtual void UpdateImp(const ::hitnlls::matrix::Matrixxf &inc) = 0;
    virtual void SaveToCache() = 0;
    virtual void RestoreCache() = 0;

    virtual void UpdatePlus(const ::hitnlls::matrix::Matrixxf &inc) final {
        if (CheckUpdate(inc)) {
            UpdateImp(inc); 
        } else {
            ::std::cout << "Update dimension is not consistent for node " << id_ << ::std::endl;
        }
    }
    virtual bool CheckUpdate(const ::hitnlls::matrix::Matrixxf &inc) final {
        std::cout << "Inc rows " << inc.Rows() << " inc cols " << inc.Cols() << std::endl;
        if(inc.Rows() == ndim_ && inc.Cols() == 1) {
            return true; 
        } else {
            return false;
        }
    }

    BaseNode(const BaseNode &) = delete;
    BaseNode &operator=(const BaseNode &) = delete;

protected:
    bool fixed_;
    bool marginalized_;
    int id_;
    int ndim_;
};

HITNLLS_REGISTER_REGISTERER(BaseNode);
#define HITNLLS_REGISTER_NODE(name) \
    HITNLLS_REGISTER_CLASS(BaseNode, name)

} // namespace node
} // namespace hitnlls