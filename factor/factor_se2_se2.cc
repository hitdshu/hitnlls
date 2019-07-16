#include "factor/factor_se2_se2.h"
#include "node/node_pose2.h"

namespace hitnlls {
namespace factor {

::hitnlls::matrix::Matrix<float, 3, 1> FactorSE2SE2::Evaluate() {
    ::hitnlls::node::NodePose2 *node_tc1w = dynamic_cast<::hitnlls::node::NodePose2 *>(GetNode(0));
    ::hitnlls::node::NodePose2 *node_tc2w = dynamic_cast<::hitnlls::node::NodePose2 *>(GetNode(1));

    ::hitnlls::node::SE2 tc1w = node_tc1w->GetEstimate();
    ::hitnlls::node::SE2 tc2w = node_tc2w->GetEstimate();
    ::hitnlls::node::SE2 twc2 = tc2w.Inverse();
    ::hitnlls::node::SE2 tc1c2 = tc1w * twc2;

    return tc1c2.ToVector3f();
}

::hitnlls::matrix::Matrixxf FactorSE2SE2::Jacobian(int nidx) {
    ::hitnlls::node::NodePose2 *node_tc1w = dynamic_cast<::hitnlls::node::NodePose2 *>(GetNode(0));
    ::hitnlls::node::NodePose2 *node_tc2w = dynamic_cast<::hitnlls::node::NodePose2 *>(GetNode(1));

    ::hitnlls::node::SE2 tc1w = node_tc1w->GetEstimate();
    ::hitnlls::matrix::Matrix22f Rc1w = tc1w.GetRotMat22f();
    ::hitnlls::node::SE2 twc1 = tc1w.Inverse();
    ::hitnlls::node::SE2 tc2w = node_tc2w->GetEstimate();
    ::hitnlls::node::SE2 twc2 = tc2w.Inverse();
    ::hitnlls::matrix::Matrix22f Rwc2 = twc2.GetRotMat22f();
    ::hitnlls::matrix::Matrix22f Rc1c2 = Rc1w * Rwc2;

    ::hitnlls::matrix::Matrixxf dtc1w(3, 3);
    ::hitnlls::matrix::Matrixxf dtc2w(3, 3);
    dtc1w(0, 0) = 1;
    dtc1w(1, 0) = Rc1c2(1, 0) * tc2w.GetPosition(0) + Rc1c2(0, 0) * tc2w.GetPosition(1);
    dtc1w(1, 1) = 1;
    dtc1w(2, 0) = - Rc1c2(0, 0) * tc2w.GetPosition(0) - Rc1c2(0, 1) * tc2w.GetPosition(1);
    dtc1w(2, 2) = 1;

    dtc2w(0, 0) = - 1;
    dtc2w(1, 0) = - dtc1w(1, 0);
    dtc2w(1, 1) = - Rc1c2(0, 0);
    dtc2w(1, 2) = - Rc1c2(0, 1);
    dtc2w(2, 0) = - dtc1w(2, 0);
    dtc2w(2, 1) = - Rc1c2(1, 0);
    dtc2w(2, 2) = - Rc1c2(1, 1);

    if (0 == nidx) {
        return dtc1w;
    } else {
        return dtc2w;
    }
}

HITNLLS_REGISTER_FACTOR(FactorSE2SE2);

} // namespace factor
} // namespace hitnlls