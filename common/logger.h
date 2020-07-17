#pragma once

#include <string>
#include <limits>
#include <mutex>
#include <iostream>
#include <functional>
#include <unordered_map>
#include <fstream>

#include "macros.h"

namespace nlls {

struct LogLevel {
    int priority;
    std::string name;
    LogLevel(int priority, const std::string &name) : priority(priority), name(name) { this->name.resize(5, ' '); }
    bool operator< (const LogLevel& rhs) const { return priority <  rhs.priority; }
    bool operator<=(const LogLevel& rhs) const { return priority <= rhs.priority; }
    bool operator> (const LogLevel& rhs) const { return priority >  rhs.priority; }
    bool operator>=(const LogLevel& rhs) const { return priority >= rhs.priority; }
};
std::ostream &operator<<(std::ostream &out, const LogLevel &ll) { out << ll.name; return out; }
const LogLevel LOG_ALL(std::numeric_limits<int>::min(), "ALL");
const LogLevel LOG_NONE(std::numeric_limits<int>::max(), "NONE");
const LogLevel LOG_DEBUG(0, "DEBUG");
const LogLevel LOG_INFO(100, "INFO");
const LogLevel LOG_WARN(200, "WARN");
const LogLevel LOG_ERROR(300, "ERROR");
const LogLevel LOG_FATAL(400, "FATAL");

using PrintHeaderFunc = std::function<void(std::ostream &, const std::string &, const LogLevel &, int)>;

int ThreadId();

class Logger final {
public:
    NLLS_NONCOPYABLE(Logger)
    explicit Logger(const std::string &name = "", const LogLevel &ll = LOG_ALL);
    void SetLogLevel(const LogLevel &ll);
    void SetPrintHeader(const PrintHeaderFunc &pf);
    void SetLoggerFile(const std::string &file_name);
    const LogLevel &GetLogLevel() const { return ll_; }
    const PrintHeaderFunc &GetPrintHeader() const { return pf_; }
private:
    struct GlobalLoggerHelper {
        std::mutex m;
        std::unordered_map<std::string, std::ofstream> name2ofs;
    };
    static GlobalLoggerHelper &GetGlobalHelper();
    class LoggerStream {
    public:
        explicit LoggerStream(const LogLevel &ll, Logger &logger) : ll_(ll), logger_(logger) {
            enabled_ = ll > logger_.ll_;
            if (enabled_) {
                GetGlobalHelper().m.lock();
                logger_.pf_(logger_.os_, logger_.name_, ll_, ThreadId());
            }
        }
        ~LoggerStream() {
            if (enabled_) {
                logger_.os_ << std::endl;
                GetGlobalHelper().m.unlock();
            }
        }
        template <typename T>
        NLLS_INLINE LoggerStream &operator<<(const T &v) {
            if (!enabled_) {
                return *this;
            } else {
                logger_.os_ << v;
                return *this;
            }
        }
    private:
        const LogLevel ll_;
        Logger &logger_;
        bool enabled_;
    };
    friend class LoggerStream;
public:
    static Logger &GetGlobalLogger();
    LoggerStream operator<<(const LogLevel &ll) {
        return LoggerStream(ll, *this);
    }

private:
    const std::string name_;
    PrintHeaderFunc pf_;
    LogLevel ll_;
    std::ostream os_;
};

#define ALL LOG_ALL
#define NONE LOG_NONE
#define DEBUG LOG_DEBUG
#define INFO LOG_INFO
#define WARN LOG_WARN
#define ERROR LOG_ERROR
#define FATAL LOG_FATAL

#define LOG(type) \
    nlls::Logger::GetGlobalLogger() << type
#define LOG_IF(type, cond, exp) \
    if (cond) { LOG(type) << exp; }
#define LOG_DETAILED(type) \
    nlls::Logger::GetGlobalLogger() << type << __FILE__ << " " << __LINE__ << " "
#define LOG_DETAILED_IF(type, cond, exp) \
    if (cond) { LOG_DETAILED(type) << exp; }
#define LOG_LEVEL(type) \
    nlls::Logger::GetGlobalLogger().SetLogLevel(type);
#define LOG2FILE(file_name) \
    nlls::Logger::GetGlobalLogger().SetLoggerFile(file_name);

} // namespace nlls