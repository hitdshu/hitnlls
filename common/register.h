#pragma once

#include <map>
#include <string>
#include <vector>
#include <iostream>

namespace hitnlls {
namespace common {

class Any {
public:
    template <typename ValueType>
    explicit Any(const ValueType &value) : content_(new Holder<ValueType>(value)) {}
    
    Any() : content_(nullptr) {}
    Any(const Any &other) : content_(other.content_ ? other.content_->Clone() : nullptr) {}
    ~Any() { 
        if (content_)
            delete content_;
    }

    template <typename ValueType>
    ValueType *AnyCast() {
        return content_ ? &(static_cast<Holder<ValueType> *>(content_)->held_) : nullptr;
    }

private:
    class PlaceHolder {
    public:
        virtual ~PlaceHolder() {}
        virtual PlaceHolder *Clone() const = 0;
    };
    template<typename ValueType>
    class Holder : public PlaceHolder {
    public:
        explicit Holder(const ValueType &value) : held_(value) {}
        virtual ~Holder() {}
        virtual PlaceHolder *Clone() const { return new Holder(held_); }
        ValueType held_;
    };

    PlaceHolder *content_;
};

class ObjectFactory {
public:
    ObjectFactory() {}
    virtual ~ObjectFactory() {}
    virtual Any NewInstance() { return Any(); }
    ObjectFactory(const ObjectFactory &) = delete;
    ObjectFactory &operator=(const ObjectFactory &) = delete;
};

typedef std::map<std::string, ObjectFactory *> FactoryMap;
typedef std::map<std::string, FactoryMap> BaseClassMap;
BaseClassMap &GlobalFactoryMap();

bool GetRegisteredClasses(const std::string &bcls_name, std::vector<std::string> &dcls_names);

} // namespace common
} // namespace hitnlls

#define HITNLLS_REGISTER_REGISTERER(base_class)                                                 \
class base_class##Registerer {                                                                  \
typedef ::hitnlls::common::Any Any;                                                             \
typedef ::hitnlls::common::FactoryMap FactoryMap;                                               \
public:                                                                                         \
    static base_class *GetInstanceByName(const ::std::string &name) {                           \
        FactoryMap &map = ::hitnlls::common::GlobalFactoryMap()[#base_class];                   \
        FactoryMap::iterator iter = map.find(name);                                             \
        if (iter == map.end()) {                                                                \
            ::std::cout << "Get instance " << name << " failed.";                               \
            return nullptr;                                                                     \
        }                                                                                       \
        Any object = iter->second->NewInstance();                                               \
        return *(object.AnyCast<base_class *>());                                               \
    }                                                                                           \
    static std::vector<base_class *> GetAllInstances() {                                        \
        std::vector<base_class *> instances;                                                    \
        FactoryMap &map = ::hitnlls::common::GlobalFactoryMap()[#base_class];                   \
        instances.reserve(map.size());                                                          \
        for (auto item : map) {                                                                 \
            Any object = item.second->NewInstance();                                            \
            instances.push_back(*(object.AnyCast<base_class *>()));                             \
        }                                                                                       \
        return instances;                                                                       \
    }                                                                                           \
    static const ::std::string GetUniqInstanceName() {                                          \
        FactoryMap &map = ::hitnlls::common::GlobalFactoryMap()[#base_class];                   \
        return map.begin()->first;                                                              \
    }                                                                                           \
    static base_class *GetUniqInstance() {                                                      \
        FactoryMap &map = ::hitnlls::common::GlobalFactoryMap()[#base_class];                   \
        Any object = map.begin()->second->NewInstance();                                        \
        return *(object.AnyCast<base_class *>());                                               \
    }                                                                                           \
    static bool IsValid(const ::std::string &name) {                                            \
        FactoryMap &map = ::hitnlls::common::GlobalFactoryMap()[#base_class];                   \
        return map.find(name) != map.end();                                                     \
    }                                                                                           \
};

#define HITNLLS_REGISTER_CLASS(base_class, name)                                                \
namespace {                                                                                     \
class ObjectFactory##name : public hitnlls::common::ObjectFactory {                             \
public:                                                                                         \
    virtual ~ObjectFactory##name() {}                                                           \
    virtual ::hitnlls::common::Any NewInstance() {                                              \
        return ::hitnlls::common::Any(new name());                                              \
    }                                                                                           \
};                                                                                              \
__attribute__((constructor)) void RegisterFactory##name() {                                     \
    ::hitnlls::common::FactoryMap &map = ::hitnlls::common::GlobalFactoryMap()[#base_class];    \
    if (map.find(#name) == map.end())                                                           \
        map[#name] = new ObjectFactory##name();                                                 \
}                                                                                               \
}
