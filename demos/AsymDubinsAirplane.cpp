/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2024, Metron, Inc.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the Metron, Inc. nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *********************************************************************/

/* Author: Mark Moll */

#include <ompl/base/spaces/OwenStateSpace.h>
#include <ompl/base/spaces/VanaStateSpace.h>
#include <ompl/base/spaces/VanaOwenStateSpace.h>
#include <ompl/base/spaces/TrochoidAirplaneStateSpace.h>
#include <ompl/base/objectives/StateCostIntegralObjective.h>
#include <ompl/geometric/SimpleSetup.h>
#include <ompl/tools/debug/Profiler.h>
#include <boost/program_options.hpp>
#include <cmath>
#include <fstream>

namespace ob = ompl::base;
namespace og = ompl::geometric;
namespace po = boost::program_options;

std::string toString(ob::State const *state, unsigned int numDims)
{
    const auto &st = *state->as<ob::VanaStateSpace::StateType>();
    std::stringstream s;
    for (unsigned int i = 0; i < numDims; ++i)
        s << st[i] << ' ';
    s << st.yaw() << '\n';
    return s.str();
}
std::string toString(ob::ScopedState<> const &state)
{
    return toString(state.get(), state.getSpace()->getDimension() - 1);
}

template <class Space>
typename Space::PathType getPath(ob::ScopedState<> const &start, ob::ScopedState<> const &goal)
{
    auto path = start.getSpace()->as<Space>()->getPath(start.get(), goal.get());
    if (!path)
    {
#if ENABLE_PROFILING == 0
        std::cout << "start: " << toString(start);
        std::cout << "goal:  " << toString(goal);
#endif
        throw std::runtime_error("Could not find a valid path");
    }
    return *path;
}

template <class Space>
void saveStatistic(std::ostream &logfile, ob::ScopedState<> const &start, ob::ScopedState<> const &goal)
{
    double length;
    unsigned int success;
    char type = '?';
    const auto& name = start.getSpace()->getName();
    try
    {
        ompl::tools::Profiler::ScopedBlock _(name);
        auto path = getPath<Space>(start, goal);
        success = 1;
        length = path.length();
        if constexpr (!std::is_same_v<Space,ob::VanaStateSpace>)
            type = (char)path.category();
    }
    catch (std::runtime_error &e)
    {
        success = 0;
        length = std::numeric_limits<double>::quiet_NaN();
    }
    logfile << ',' << success << ',' << length << ',' << type;
}

template <class Space>
void savePath(ob::ScopedState<> const &start, ob::ScopedState<> const &goal, std::string const &pathName)
{
    auto path = getPath<Space>(start, goal);
    auto space = start.getSpace()->as<Space>();
    ob::ScopedState<Space> state(start);
    // double dist = path.length(), d1, d2;
    // for (unsigned i=1; i<100; ++i)
    // {
    //     space->interpolate(start.get(), goal.get(), (double)i/50., path, st);
    //     d1 = space->distance(start.get(), state.get());
    //     d2 = space->distance(state.get(), goal.get());
    //     if (std::abs(d1 + d2 - dist) > .01)
    //     {
    //         std::cout << i << "," << dist << "," << d1 << "," << d2 << "," << d1+d2-dist << "," << toString(state) <<
    //         '\n';
    //     }
    // }

    if (!pathName.empty())
    {
        std::ofstream outfile(pathName);
        for (double t = 0.; t <= 1.001; t += .001)
        {
            space->interpolate(start.get(), goal.get(), t, path, state.get());
            outfile << t << ' ' << toString(state);
        }
    }
    std::cout << "start: " << toString(start) << "goal: " << toString(goal) << "path: " << path << '\n';
}

int main(int argc, char *argv[])
{
    try
    {
        std::string pathName;
        double radius, maxPitch, minPitch, windRatio, windHeading;
        unsigned numSamples;
        po::options_description desc("Options");
        // clang-format off
        desc.add_options()
            ("help", "show help message")
            ("trochoidairplane", "generate a trochoid airplane path starting from (x,y,z,yaw)=(0,0,0,0) to a random pose")
            ("savepath", po::value<std::string>(&pathName), "save an (approximate) solution path to file") 
            ("start", po::value<std::vector<double>>()->multitoken(),
                "use (x,y,z,yaw) as the start")
            ("goal", po::value<std::vector<double>>()->multitoken(),
                "use (x,y,z,yaw) as the goal instead of a random state")
            ("radius", po::value<double>(&radius)->default_value(1.), "turn radius")
            ("maxpitch", po::value<double>(&maxPitch)->default_value(.5), "maximum pitch angle")
            ("minpitch", po::value<double>(&minPitch)->default_value(-0.5), "minimum pitch angle")
            ("windratio", po::value<double>(&windRatio)->default_value(.5), "wind ratio")
            ("windheading", po::value<double>(&windHeading)->default_value(.5), "wind heading");
            // clang-format on

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc,
                                         po::command_line_style::unix_style ^ po::command_line_style::allow_short),
                  vm);
        po::notify(vm);

        if ((vm.count("help") != 0u) || argc == 1)
        {
            std::cout << desc << "\n";
            return 1;
        }

        ob::StateSpacePtr space = [&]() -> ob::StateSpacePtr
        {
            // set the bounds for the R^3 part of the space
            ob::RealVectorBounds bounds(3);
            bounds.setLow(-10);
            bounds.setHigh(10);
            if (vm.count("trochoidairplane") != 0)
            {
                auto space = std::make_shared<ob::TrochoidAirplaneStateSpace>(radius, windRatio, windHeading, maxPitch);
                // space->setTolerance(1e-16);
                space->setBounds(bounds);
                return space;
            }
            else
            {
                auto space = std::make_shared<ob::OwenStateSpace>(radius, maxPitch, minPitch);
                space->setBounds(bounds);
                return space;
            }
        }();
        ob::ScopedState<> start(space);
        ob::ScopedState<> goal(space);

        if (vm.count("start") == 0u)
        {
            for (unsigned int i = 0; i < space->getDimension(); ++i)
                start[i] = 0.;
        }
        else
        {
            std::cout << "Start dimensions: " << space->getDimension() << std::endl;
            auto startVec = vm["start"].as<std::vector<double>>();
            for (unsigned int i = 0; i < space->getDimension(); ++i)
                start[i] = startVec[i];
        }

        if (vm.count("goal") == 0u)
        {
            goal.random();
        }
        else
        {
            auto goalVec = vm["goal"].as<std::vector<double>>();
            std::cout << "Goal dimensions: " << space->getDimension() << std::endl;
            for (unsigned int i = 0; i < space->getDimension(); ++i)
                goal[i] = goalVec[i];
        }

        savePath<ob::OwenStateSpace>(start, goal, pathName);
    }
    catch (std::exception &e)
    {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch (...)
    {
        std::cerr << "Exception of unknown type!\n";
    }

    return 0;
}
