#include "factor/factor_point2f_se2.h"
#include "node/node_point2f.h"
#include "node/node_pose2.h"

namespace hitnlls {
namespace factor {

::hitnlls::matrix::Matrix<float, 2, 1> FactorPoint2fSE2::Evaluate() {
    ::hitnlls::node::NodePoint2f *node_pw = dynamic_cast<::hitnlls::node::NodePoint2f *>(GetNode(0));
    ::hitnlls::node::NodePose2 *node_tcw = dynamic_cast<::hitnlls::node::NodePose2 *>(GetNode(1));

    ::hitnlls::matrix::Vector2f pw = node_pw->GetEstimate();
    ::hitnlls::node::SE2 se2 = node_tcw->GetEstimate();
    return se2.MapPoint(pw);
}

::hitnlls::matrix::Matrixxf FactorPoint2fSE2::Jacobian(int nidx) {
    ::hitnlls::node::NodePoint2f *node_pw = dynamic_cast<::hitnlls::node::NodePoint2f *>(GetNode(0));
    ::hitnlls::node::NodePose2 *node_tcw = dynamic_cast<::hitnlls::node::NodePose2 *>(GetNode(1));

    ::hitnlls::matrix::Vector2f pw = node_pw->GetEstimate();
    ::hitnlls::node::SE2 tcw = node_tcw->GetEstimate();
    ::hitnlls::matrix::Matrix22f rcw = tcw.GetRotMat22f();

    ::hitnlls::matrix::Matrixxf dpw = rcw;
    ::hitnlls::matrix::Matrixxf dse2(2, 3);
    dse2(0, 0) = rcw(0, 1) * pw[0] - rcw(0, 0) * pw[1];
    dse2(1, 0) = rcw(0, 0) * pw[0] + rcw(0, 1) * pw[1];
    dse2(0, 1) = 1.0;
    dse2(1, 2) = 1.0;

    if (0 == nidx) {
        return dpw;
    } else {
        return dse2;
    }
}

HITNLLS_REGISTER_FACTOR(FactorPoint2fSE2);

} // namespace factor
} // namespace hinlls