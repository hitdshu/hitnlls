#include "factor/factor_point3f_se3.h"
#include "node/node_point3f.h"
#include "node/node_pose3.h"

namespace hitnlls {
namespace factor {

::hitnlls::matrix::Matrix<float, 2, 1> FactorPoint3fSE3::Evaluate() {
    ::hitnlls::node::NodePoint3f *node_pw = dynamic_cast<::hitnlls::node::NodePoint3f *>(GetNode(0));
    ::hitnlls::node::NodePose3 *node_tcw = dynamic_cast<::hitnlls::node::NodePose3 *>(GetNode(1));

    ::hitnlls::matrix::Vector3f pw = node_pw->GetEstimate();
    ::hitnlls::node::SE3 se3 = node_tcw->GetEstimate();
    ::hitnlls::matrix::Vector3f pc = se3.MapPoint(pw);
    return cam_->Project(pc);
}

::hitnlls::matrix::Matrixxf FactorPoint3fSE3::Jacobian(int nidx) {
    ::hitnlls::node::NodePoint3f *node_pw = dynamic_cast<::hitnlls::node::NodePoint3f *>(GetNode(0));
    ::hitnlls::node::NodePose3 *node_tcw = dynamic_cast<::hitnlls::node::NodePose3 *>(GetNode(1));

    ::hitnlls::matrix::Vector3f pw = node_pw->GetEstimate();
    ::hitnlls::node::SE3 tcw = node_tcw->GetEstimate();
    ::hitnlls::matrix::Vector3f pc = tcw.MapPoint(pw);
    ::hitnlls::matrix::Matrix22f rcw = tcw.GetRotMat33f();

    ::hitnlls::matrix::Matrixxf dpw = rcw;
    ::hitnlls::matrix::Matrixxf dse3(3, 6);
    dse3.SetBlock(0, 3, 3, 3, ::hitnlls::matrix::Matrix33f::Identity());
    ::hitnlls::matrix::Matrix33f dso3;
    dso3(0, 1) = - pw[2];
    dso3(1, 0) = pw[2];
    dso3(0, 2) = pw[1];
    dso3(2, 0) = - pw[1];
    dso3(1, 2) = - pw[0];
    dso3(2, 1) = pw[0];
    dso3 = - rcw * dso3;
    dse3.SetBlock(0, 0, 3, 3, dso3);

    if (0 == nidx) {
        return cam_->Derivative(pc) * dpw;
    } else {
        return cam_->Derivative(pc) * dse3;
    }
}

HITNLLS_REGISTER_FACTOR(FactorPoint3fSE3);

} // namespace factor
} // namespace hinlls