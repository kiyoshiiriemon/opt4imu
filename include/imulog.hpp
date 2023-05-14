/****************************************************************************
 * Copyright (C) 2023 Kiyoshi Irie
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this file,
 * You can obtain one at https://mozilla.org/MPL/2.0/.
 ****************************************************************************/

#pragma once
#include <fstream>
#include <iostream>
#include <list>
#include <vector>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include "imu.hpp"

namespace opt4imu {

class IMULog
{
public:
	std::vector<Imu> data;
	IMULog() {}
	size_t size() const { return data.size(); }
	bool loadImuFile(const char *filename, double start_time = -1, double end_time = -1)
    {
        std::string delim(", \t\r\n");
        bool header_processed = false;
        std::ifstream ifs(filename);
        if (!ifs) return false;

        std::string buf;
        int time_sec_ind = 0;
        int time_nsec_ind = 1;
        int acc_x_ind = 2, acc_y_ind = 3, acc_z_ind = 4;
        int gyr_x_ind = 5, gyr_y_ind = 6, gyr_z_ind = 7;
        int q0_ind = 14, q1_ind = 15, q2_ind = 16, q3_ind = 17;

        while (ifs && std::getline(ifs, buf)) {
            std::vector<std::string> items;
            boost::split(items, buf, boost::is_any_of(delim));
            //std::cout << buf << " : split: " << items.size() << std::endl;
            if (items[q0_ind].size() == 0) {
                continue;
            }
            Imu imu;
            int sec = boost::lexical_cast<int>(items[time_sec_ind]);
            int nsec = boost::lexical_cast<int>(items[time_nsec_ind]);
            //std::cerr << "time: " << sec << " " << nsec << std::endl;
            imu.header.stamp = opt4imu::Time(sec, nsec);
            if (start_time > 0 && opt4imu::Time::duration(opt4imu::Time(0, 0), imu.header.stamp) < start_time) continue;
            if (end_time > 0 && opt4imu::Time::duration(opt4imu::Time(0, 0), imu.header.stamp) > end_time) continue;
            if (q0_ind > 0) {
                imu.orientation.w() = std::stof(items[q0_ind]);
                imu.orientation.x() = std::stof(items[q1_ind]);
                imu.orientation.y() = std::stof(items[q2_ind]);
                imu.orientation.z() = std::stof(items[q3_ind]);
//            imu.orientation.w() = boost::lexical_cast<double>(items[q0_ind]);
//            imu.orientation.x() = boost::lexical_cast<double>(items[q1_ind]);
//            imu.orientation.y() = boost::lexical_cast<double>(items[q2_ind]);
//            imu.orientation.z() = boost::lexical_cast<double>(items[q3_ind]);
            }
            //std::cerr << "q: " << imu.orientation << std::endl;
            imu.linear_acceleration << boost::lexical_cast<double>(items[acc_x_ind]),
                    boost::lexical_cast<double>(items[acc_y_ind]),
                    boost::lexical_cast<double>(items[acc_z_ind]);
            imu.angular_velocity << boost::lexical_cast<double>(items[gyr_x_ind]),
                    boost::lexical_cast<double>(items[gyr_y_ind]),
                    boost::lexical_cast<double>(items[gyr_z_ind]);
            data.push_back(imu);
        }
        return true;
    }
	bool loadMTFile(const char *filename) {
        std::string delim(", \t");
        bool header_processed = false;
		std::ifstream ifs(filename);
		if (!ifs) return false;

		std::string buf;
        int sample_time_ind = -1;
		int acc_x_ind = 0, acc_y_ind = 0, acc_z_ind = 0;
		int gyr_x_ind = 0, gyr_y_ind = 0, gyr_z_ind = 0;
		int q0_ind = 0, q1_ind = 0, q2_ind = 0, q3_ind = 0;

		while (ifs && std::getline(ifs, buf)) {
			if (boost::starts_with( buf, "//" )) continue;

			if (!header_processed && boost::starts_with( buf, "PacketCounter")) {
		        std::vector<std::string> items;
		        boost::split(items, buf, boost::is_any_of(delim));
				//std::cout << items.size() << std::endl;
				int ind = 0;
				for(const std::string &s : items) {
					//std::cout << s << std::endl;
                    if (s == std::string("SampleTimeFine")) {
                        sample_time_ind = ind;
                    }
					if (s == std::string("Acc_X")) {
						acc_x_ind = ind;
					}
					if (s == std::string("Acc_Y")) {
						acc_y_ind = ind;
					}
					if (s == std::string("Acc_Z")) {
						acc_z_ind = ind;
					}
					if (s == std::string("Gyr_X")) {
						gyr_x_ind = ind;
					}
					if (s == std::string("Gyr_Y")) {
						gyr_y_ind = ind;
					}
					if (s == std::string("Gyr_Z")) {
						gyr_z_ind = ind;
					}
					if (s == std::string("Quat_q0")) {
						q0_ind = ind;
					}
					if (s == std::string("Quat_q1")) {
						q1_ind = ind;
					}
					if (s == std::string("Quat_q2")) {
						q2_ind = ind;
					}
					if (s == std::string("Quat_q3")) {
						q3_ind = ind;
					}
					++ind;
				}
				header_processed = true;
			}
			else {
		        std::vector<std::string> items;
		        boost::split(items, buf, boost::is_any_of(delim));
				if (items[q0_ind].size() == 0) {
					continue;
				}
				Imu imu;
				if (sample_time_ind > 0) {
                    double usec = boost::lexical_cast<double>(items[1]) * 100;
                    int sec = (int)(usec / 1e6);
                    //std::cerr << "time: " << sec << " " << usec << std::endl;
                    imu.header.stamp = opt4imu::Time(sec, int((usec - sec * 1e6) * 1000));
				}
				if (q0_ind > 0) {
                    imu.orientation.w() = boost::lexical_cast<double>(items[q0_ind]);
                    imu.orientation.x() = boost::lexical_cast<double>(items[q1_ind]);
                    imu.orientation.y() = boost::lexical_cast<double>(items[q2_ind]);
                    imu.orientation.z() = boost::lexical_cast<double>(items[q3_ind]);
                }
				imu.linear_acceleration << boost::lexical_cast<double>(items[acc_x_ind]),
					boost::lexical_cast<double>(items[acc_y_ind]),
					boost::lexical_cast<double>(items[acc_z_ind]);
				imu.angular_velocity << boost::lexical_cast<double>(items[gyr_x_ind]),
					boost::lexical_cast<double>(items[gyr_y_ind]),
					boost::lexical_cast<double>(items[gyr_z_ind]);
				data.push_back(imu);
				if (boost::lexical_cast<int>(items[0]) % 10 == 0) {
					std::cout << 
						items[q0_ind] << " " <<
						items[q1_ind] << " " <<
						items[q2_ind] << " " <<
						items[q3_ind] << " " <<
						items[acc_x_ind] << " " <<
						items[acc_y_ind] << " " <<
						items[acc_z_ind] << " " <<
						items[gyr_x_ind] << " " <<
						items[gyr_y_ind] << " " <<
						items[gyr_z_ind] << " " <<
						std::endl;
				}
			}
			//std::cout << buf << std::endl;
		}
		return true;
	}
	void subtractAccBias(double ax0, double ay0, double az0)
	{
		opt4imu::Vec3D offset(ax0, ay0, az0);
		for (auto &d : data) {
			d.linear_acceleration = d.linear_acceleration - offset;
		}
	}
	void setUniformTimeStamp(double step_usec) {
        double sec = 0;
        for(auto &v : data) {
            double usec = (sec - (int)sec) * 1e6;
            //std::cout << sec << std::endl;
            v.header.stamp = opt4imu::Time((int)sec, (int)usec * 1000);
            sec += step_usec * 1e-6;
	    }
    }
    Eigen::Vector3d getInitialMeanGyro(double dur_sec) {
        Eigen::Vector3d gyro;
        gyro << 0, 0, 0;
        double cnt = 0;
        for(auto &v : data) {
            double t = opt4imu::Time::duration(data[0].header.stamp, v.header.stamp);
            if (t >= dur_sec) break;
            gyro += v.angular_velocity;
            ++cnt;
        }
        return gyro / cnt;
    }
    Eigen::Vector3d getInitialMeanAcc(double dur_sec) {
        Eigen::Vector3d acc;
        acc << 0, 0, 0;
        double cnt = 0;
        for(auto &v : data) {
            double t = opt4imu::Time::duration(data[0].header.stamp, v.header.stamp);
            if (t >= dur_sec) break;
            acc += v.linear_acceleration;
            ++cnt;
        }
        return acc / cnt;
	}
	void applyCalibrationParams(const opt4imu::Vec3D &accbias, const opt4imu::Vec3D &gyrobias, const opt4imu::Vec3D &accgain, const opt4imu::Vec3D &gyrogain)
    {
        for(auto &v : data) {
            v.linear_acceleration -= accbias;
            v.linear_acceleration.array() *= accgain.array();
            v.angular_velocity -= gyrobias;
            v.angular_velocity.array() *= gyrogain.array();
        }
    }
    double duration(int i, int j) const { return opt4imu::Time::duration(data[i].header.stamp, data[j].header.stamp); }
    const opt4imu::Vec3D &acc(int i) const { return data[i].linear_acceleration; }
    const opt4imu::Vec3D &gyro(int i) const { return data[i].angular_velocity; }
};
} // opt4imu
