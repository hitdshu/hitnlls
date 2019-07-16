#include "common/register.h"

namespace hitnlls {
namespace common {

BaseClassMap &GlobalFactoryMap() {
    static BaseClassMap factory_map;
    return factory_map;
}

bool GetRegisteredClasses(const std::string &bcls_name, std::vector<std::string> &dcls_names) {
    BaseClassMap &map = GlobalFactoryMap();
    auto iter = map.find(bcls_name);
    if (iter == map.end()) {
        ::std::cout << "class not registered:" << bcls_name << std::endl;
        return false;
    }
    for (auto pair : iter->second) {
        dcls_names.push_back(pair.first);
    }
    return true;
}


} // namespace common
} // namespace hitnlls