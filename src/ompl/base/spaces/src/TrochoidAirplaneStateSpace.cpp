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

#include "ompl/base/spaces/TrochoidAirplaneStateSpace.h"
#include "ompl/base/spaces/Dubins3DMotionValidator.h"
#include "ompl/base/SpaceInformation.h"
#include "ompl/tools/config/MagicConstants.h"
#include "ompl/util/Exception.h"
#include <boost/math/tools/toms748_solve.hpp>

using namespace ompl::base;

namespace
{
    constexpr double twopi = 2. * boost::math::constants::pi<double>();

    // tolerance for boost root finding algorithm
    const boost::math::tools::eps_tolerance<double> TOLERANCE(20);

    // max # iterations for searching for roots
    constexpr std::uintmax_t MAX_ITER = 32;
}  // namespace

TrochoidAirplaneStateSpace::TrochoidAirplaneStateSpace(double turningRadius, double windRatio, double windHeading, double maxPitch, double verticalWindRatio)
  : rho_(turningRadius), eta_(windRatio), psi_(windHeading), tanMaxPitch_(std::tan(maxPitch)), eta_z_(verticalWindRatio), trochoidSpace_(turningRadius, windRatio, windHeading)
{
    setName("TrochoidAirplane" + getName());
    type_ = STATE_SPACE_TROCHOID_AIRPLANE;
    addSubspace(std::make_shared<RealVectorStateSpace>(3), 1.0);
    addSubspace(std::make_shared<SO2StateSpace>(), 0.5);
    lock();
}

State *TrochoidAirplaneStateSpace::allocState() const
{
    auto *state = new StateType();
    allocStateComponents(state);
    return state;
}

void TrochoidAirplaneStateSpace::registerProjections()
{
    class OwenDefaultProjection : public ProjectionEvaluator
    {
    public:
        OwenDefaultProjection(const StateSpace *space) : ProjectionEvaluator(space)
        {
        }

        unsigned int getDimension() const override
        {
            return 3;
        }

        void defaultCellSizes() override
        {
            cellSizes_.resize(3);
            bounds_ = space_->as<TrochoidAirplaneStateSpace>()->getBounds();
            cellSizes_[0] = (bounds_.high[0] - bounds_.low[0]) / magic::PROJECTION_DIMENSION_SPLITS;
            cellSizes_[1] = (bounds_.high[1] - bounds_.low[1]) / magic::PROJECTION_DIMENSION_SPLITS;
            cellSizes_[2] = (bounds_.high[2] - bounds_.low[2]) / magic::PROJECTION_DIMENSION_SPLITS;
        }

        void project(const State *state, Eigen::Ref<Eigen::VectorXd> projection) const override
        {
            projection = Eigen::Map<const Eigen::VectorXd>(
                state->as<TrochoidAirplaneStateSpace::StateType>()->as<RealVectorStateSpace::StateType>(0)->values, 3);
        }
    };

    registerDefaultProjection(std::make_shared<OwenDefaultProjection>(this));
}

std::optional<TrochoidAirplaneStateSpace::PathType> TrochoidAirplaneStateSpace::getPath(const State *state1, const State *state2) const
{
    auto s1 = state1->as<StateType>();
    auto s2 = state2->as<StateType>();
    auto path = trochoidSpace_.getPath(state1, state2, rho_, eta_, psi_);
    double dz = (*s2)[2] - (*s1)[2], len = path.length();
    double dz_sign = dz > 0 ? 1 : -1;
    if (std::abs(dz) <= len * tanMaxPitch_ * std::abs(eta_z_ + dz_sign))
    {
        // low altitude path
        return PathType{path, rho_, eta_, psi_, dz, eta_z_};
    }
    auto periodic_path = trochoidSpace_.getPath(state1, state1, rho_, eta_, psi_, true);
    double lengthPeriodicPath = periodic_path.length();
    if (std::abs(dz) > (len + lengthPeriodicPath) * tanMaxPitch_ * std::abs(eta_z_ + dz_sign))
    {
        // high altitude path
        unsigned int k = std::floor((std::abs(dz) / (tanMaxPitch_ * std::abs(eta_z_ + dz_sign)) - len) / lengthPeriodicPath);
        auto radius = rho_;
        auto radiusFun = [&, this](double r)
        { return (trochoidSpace_.getPath(state1, state2, r, eta_, psi_).length() + trochoidSpace_.getPath(state1, state1, r, eta_, psi_, true).length() * k) * tanMaxPitch_ * std::abs(eta_z_ + dz_sign) - std::abs(dz); };
        std::uintmax_t iter = MAX_ITER;
        auto result = boost::math::tools::bracket_and_solve_root(radiusFun, radius, 2., true, TOLERANCE, iter);
        radius = .5 * (result.first + result.second);
        path = trochoidSpace_.getPath(state1, state2, radius, eta_, psi_);
        periodic_path = trochoidSpace_.getPath(state1, state1, radius, eta_, psi_, true);
        return PathType{path, radius, eta_, psi_, dz, eta_z_, k, periodic_path};
    } else {
        
        // medium altitude path
        auto zi = trochoidSpace_.allocState()->as<TrochoidStateSpace::StateType>();
        auto phi = 0.1;
        auto phiFun = [&, this](double phi)
        {
            turn(state1, rho_, eta_, psi_, phi, zi);
            return (std::abs(phi) + trochoidSpace_.getPath(zi, state2, rho_, eta_, psi_).length()) * tanMaxPitch_ - std::abs(dz);
        };

        try
        {
            std::uintmax_t iter = MAX_ITER;
            auto result = boost::math::tools::bracket_and_solve_root(phiFun, phi, 2., true, TOLERANCE, iter);
            phi = .5 * (result.first + result.second);
            if (std::abs(phiFun(phi)) > 2e-5)
                throw std::domain_error("fail");
        }
        catch (...)
        {
            try
            {
                std::uintmax_t iter = MAX_ITER;
                phi = -.1;
                auto result = boost::math::tools::bracket_and_solve_root(phiFun, phi, 2., true, TOLERANCE, iter);
                phi = .5 * (result.first + result.second);
            }
            catch (...)
            {
                // OMPL_ERROR("this shouldn't be happening!");
                return {};
            }
        }
        assert(std::abs(phiFun(phi)) < 1e-5);
        turn(state1, rho_, eta_, psi_, phi, zi);
        path = trochoidSpace_.getPath(zi, state2, rho_, eta_, psi_);
        trochoidSpace_.freeState(zi);
        return PathType{path, rho_, eta_, psi_, dz, eta_z_, phi};
    }
}

void TrochoidAirplaneStateSpace::turn(const State *from, double turnRadius, double windRatio, double windHeading, double angle, State *state) const
{
    auto s0 = from->as<TrochoidStateSpace::StateType>();
    auto s1 = state->as<TrochoidStateSpace::StateType>();
    double phi = s0->getYaw(), v = angle;
    double delta = angle > 0 ? 1 : -1;
    double x_t10 = s0->getX() - (turnRadius* delta) * sin(phi);
    double y_t10 = s0->getY() + (turnRadius * delta) * cos(phi);
    s1->setXY(x_t10 + (turnRadius * delta) * sin(delta *(1.0/turnRadius) * v + phi)  + windRatio * v, \
    y_t10 - (turnRadius * delta) * cos(delta *(1.0/turnRadius) * v + phi));
    s1->setYaw(phi + delta *(1.0/turnRadius) *v);
}

double TrochoidAirplaneStateSpace::distance(const State *state1, const State *state2) const
{
    auto path = trochoidSpace_.getPath(state1, state2, rho_, eta_, psi_);
    double dz = state1->as<ompl::base::TrochoidAirplaneStateSpace::StateType>()->as<ompl::base::RealVectorStateSpace::StateType>(0)->values[2] - state2->as<ompl::base::TrochoidAirplaneStateSpace::StateType>()->as<ompl::base::RealVectorStateSpace::StateType>(0)->values[2];
    if (std::abs(dz) <= path.length() * tanMaxPitch_) {
        return path.length();
    } else {
        return std::abs(dz) / tanMaxPitch_;
    }
    return getMaximumExtent();
}

unsigned int TrochoidAirplaneStateSpace::validSegmentCount(const State *state1, const State *state2) const
{
    return StateSpace::validSegmentCount(state1, state2);
}

void TrochoidAirplaneStateSpace::interpolate(const State *from, const State *to, double t, State *state) const
{
    if (auto path = getPath(from, to))
        interpolate(from, to, t, *path, state);
    else if (from != state)
        copyState(state, from);
}

void TrochoidAirplaneStateSpace::interpolate(const State *from, const State *to, double t, PathType &path, State *state) const
{
    if (t >= 1.)
    {
        if (to != state)
            copyState(state, to);
        return;
    }
    if (t <= 0.)
    {
        if (from != state)
            copyState(state, from);
        return;
    }

    auto f = from->as<StateType>();
    auto s = state->as<StateType>();
    (*s)[2] = (*f)[2] + t * path.deltaZ_;
    if (path.phi_ == 0.)
    {
        if (path.numTurns_ == 0)
        {
            // low altitude path
            trochoidSpace_.interpolate(from, path.path_, t, state, path.turnRadius_, path.windRatio_, path.windHeading_);
        }
        else
        {
            // high altitude path
            
            // Parse Trochoidal periodic paths
            double lengthPeriodicPath = path.periodic_path_.length();
            auto lengthSpiral = lengthPeriodicPath * path.numTurns_;
            
            auto lengthPath = path.path_.length();
            auto length = lengthSpiral + lengthPath, dist = t * length;
            if (dist > lengthSpiral) {
                trochoidSpace_.interpolate(from, path.path_, (dist - lengthSpiral) / lengthPath, state, path.turnRadius_, path.windRatio_, path.windHeading_);
            } else {
                double periodic_t = (dist - lengthPeriodicPath * std::floor(dist/lengthPeriodicPath))/lengthPeriodicPath;
                trochoidSpace_.interpolate(from, path.periodic_path_, periodic_t, state, path.turnRadius_, path.windRatio_, path.windHeading_);
            }
        }
    }
    else
    {
        // medium altitude path
        auto lengthTurn = std::abs(path.phi_);
        auto lengthPath = path.path_.length();
        auto length = lengthTurn + lengthPath, dist = t * length;
        if (dist > lengthTurn)
        {
            State *s = (state == to) ? trochoidSpace_.allocState() : state;
            turn(from, path.turnRadius_, path.windRatio_, path.windHeading_, path.phi_, s);
            trochoidSpace_.interpolate(s, path.path_, (dist - lengthTurn) / lengthPath, state, path.turnRadius_, path.windRatio_, path.windHeading_);
            if (state == to)
                trochoidSpace_.freeState(s);
        }
        else
        {
            auto angle = dist;
            if (path.phi_ < 0)
                angle = -angle;
            turn(from, path.turnRadius_, path.windRatio_, path.windHeading_, angle, state);
        }
    }

    getSubspace(1)->enforceBounds(state->as<CompoundStateSpace::StateType>()->as<SO2StateSpace::StateType>(1));
}

double TrochoidAirplaneStateSpace::PathType::length() const
{
    double hlen = turnRadius_ * (path_.length() + twopi * numTurns_ + phi_);
    return std::sqrt(hlen * hlen + deltaZ_ * deltaZ_);
}

TrochoidAirplaneStateSpace::PathCategory TrochoidAirplaneStateSpace::PathType::category() const
{
    if (phi_ == 0.)
    {
        if (numTurns_ == 0)
            return PathCategory::LOW_ALTITUDE;
        else
            return PathCategory::HIGH_ALTITUDE;
    }
    else
    {
        if (numTurns_ == 0)
            return PathCategory::MEDIUM_ALTITUDE;
        else
            return PathCategory::UNKNOWN;
    }
}

namespace ompl::base
{
    std::ostream &operator<<(std::ostream &os, const TrochoidAirplaneStateSpace::PathType &path)
    {
        static const TrochoidStateSpace trochoidStateSpace;

        os << "TrochoidAirplanePath[ category = " << (char)path.category() << "\n\tlength = " << path.length()
           << "\n\tturnRadius=" << path.turnRadius_ << "\n\twindRatio=" << path.windRatio_
           <<  "\n\twindHeading=" << path.windHeading_ << "\n\tverticalWindRatio=" << path.verticalWindRatio_
           << "\n\tdeltaZ=" << path.deltaZ_ << "\n\tphi=" << path.phi_
           << "\n\tnumTurns=" << path.numTurns_ << "\n\tpath=" << path.path_;
        os << "]";
        return os;
    }
}  // namespace ompl::base
