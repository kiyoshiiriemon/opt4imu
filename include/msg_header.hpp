/****************************************************************************
 * Copyright (C) 2023 Kiyoshi Irie
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this file,
 * You can obtain one at https://mozilla.org/MPL/2.0/.
 ****************************************************************************/

#pragma once
#include <iomanip>
#ifdef _MSC_VER
#include "vc_gettimeofday.hpp"
#else
#include <sys/time.h>
#endif
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/c_local_time_adjustor.hpp>
#include <msgpack.hpp>

namespace opt4imu {

struct Time
{
    int32_t sec;
    int32_t nsec;
    Time() : sec(0), nsec(0) {}
    Time(double sec_since_epoch)
    {
        double t_sec;
        nsec = (int32_t) (std::modf(sec_since_epoch, &t_sec) * 1e9);
        sec = (int32_t)t_sec;
    }
    Time(int32_t in_sec, int32_t in_nsec) : sec(in_sec), nsec(in_nsec)
    {
    }
    static Time now() {
#ifdef _MSC_VER
        auto time_since_epoch = std::chrono::system_clock::now().time_since_epoch();
        auto sec = std::chrono::duration_cast<std::chrono::seconds>(time_since_epoch).count();
        auto nsec = std::chrono::duration_cast<std::chrono::nanoseconds>(time_since_epoch).count();
        return Time(sec, nsec - sec * 1000000000);
#else
        struct timeval t;
        gettimeofday(&t, NULL);
        //std::cout << "current time: " << t.tv_sec << std::endl;
        return Time(static_cast<int32_t>(t.tv_sec), static_cast<int32_t>(t.tv_usec)*1000);
#endif
    }
    bool is_newer_than(const Time &t) const
    {
        //std::cout << "this: " << sec << ", " << nsec << std::endl;
        //std::cout << "in  : " << t.sec << ", " << t.nsec << std::endl;
        return (t.sec < sec ) || (t.sec == sec && t.nsec < nsec);
    }
    double secondsSinceEpoch()
    {
        return sec + nsec * 1e-9;
    }
    Time add(double duration)
    {
        return Time(secondsSinceEpoch() + duration);
    }
    static double duration(const Time &t_old, const Time &t_new)
    {
        double v = t_new.sec - t_old.sec;
        v += (t_new.nsec - t_old.nsec)/1e9;
        return v;
    }
};

inline std::ostream& operator<<(std::ostream &os, const Time &obj)
{
    os << obj.sec << "." << std::setfill('0') << std::setw(9) << obj.nsec;
    return os;
}

struct Header
{
    uint32_t seq = 0;
    Time stamp;
    std::string frame_id;
};

}

