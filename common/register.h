#pragma once

#include <memory>
#include <map>
#include <string>

#define HITNLLS_REGISTER_REGISTER(base_class)                                                                   \
class Register##base_class {                                                                                    \
public:                                                                                                         \
    typedef std::shared_ptr<base_class> (*base_class##Creator)();                                               \
    typedef std::map<std::string, base_class##Creator> base_class##CreatorRegistry;                             \
    static base_class##CreatorRegistry &Registry() {                                                            \
        static base_class##CreatorRegistry registry;                                                            \
        return registry;                                                                                        \
    }                                                                                                           \
    static void Add##base_class##Creator(const std::string &type, base_class##Creator creator) {                \
        base_class##CreatorRegistry &registry = Registry();                                                     \
        registry[type] = creator;                                                                               \
    }                                                                                                           \
    static std::shared_ptr<base_class> Create##base_class(const std::string &type) {                            \
        base_class##CreatorRegistry &registry = Registry();                                                     \
        return registry[type]();                                                                                \
    }                                                                                                           \
};                                                                                                              \
class RegisterRegister##base_class {                                                                            \
public:                                                                                                         \
    RegisterRegister##base_class(const std::string &type, std::shared_ptr<base_class> (*creator)()) {           \
        Register##base_class::Add##base_class##Creator(type, creator);                                          \
    }                                                                                                           \
};

#define HITNLLS_REGISTER_CLASS(base_class, name)                                                                \
namespace {                                                                                                     \
std::shared_ptr<base_class> name##Creator() {                                                                   \
    return std::shared_ptr<base_class>(new name());                                                             \
}                                                                                                               \
RegisterRegister##base_class name##_register(#name, name##Creator);                                             \
}
