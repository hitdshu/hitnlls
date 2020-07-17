#include <iomanip>
#include <atomic>
#include <ctime>

#include "logger.h"

namespace {
void DefaultPrintHeaderFunc(std::ostream &os, const std::string &ln, const nlls::LogLevel &ll, int tid) {
    std::string str_date;
    std::string str_time;
    {
        time_t tp = time(0);
	    struct tm *tmp = localtime(&tp);
        char ca_date[11];
	    char ca_time[9];
	    strftime(ca_date, 11, "%Y.%m.%d", tmp);
        strftime(ca_time, 9, "%H:%M:%S", tmp);
        str_date = std::string(ca_date);
        str_time = std::string(ca_time);
    }
    if (ln.empty()) {
        os << "[" << std::setw(3) << std::right << std::setfill('0') << tid << "] " << ll << " " << str_date << " " << str_time << ": ";
    } else {
        os << "[" << std::setw(3) << std::right << std::setfill('0') << tid << "] " << ln << " " << ll << " " << str_date << " " << str_time << ": ";
    }
}
} // namespace 

namespace nlls {

int ThreadId() {
    static std::atomic<int> thread_cnt(0);
    thread_local int tid = ++thread_cnt;
    return tid;
}

Logger::Logger(const std::string &name, const LogLevel &ll) : name_(name), ll_(ll), os_(std::cout.rdbuf()) {
    pf_ = &DefaultPrintHeaderFunc;
}

void Logger::SetLogLevel(const LogLevel &ll) {
    std::unique_lock<std::mutex> lock(GetGlobalHelper().m);
    ll_ = ll;
}

void Logger::SetPrintHeader(const PrintHeaderFunc &pf) {
    std::unique_lock<std::mutex> lock(GetGlobalHelper().m);
    pf_ = pf;
}

void Logger::SetLoggerFile(const std::string &file_name) {
    std::unique_lock<std::mutex> lock(GetGlobalHelper().m);
    std::ofstream tmp(file_name);
    if (!tmp) {
        return;
    }
    GetGlobalHelper().name2ofs[name_] = std::move(tmp);
    os_.rdbuf(GetGlobalHelper().name2ofs[name_].rdbuf());
}

Logger::GlobalLoggerHelper &Logger::GetGlobalHelper() {
    static GlobalLoggerHelper glh;
    return glh;
}

Logger &Logger::GetGlobalLogger() {
    static Logger logger;
    return logger;
}

} // namespace nlls