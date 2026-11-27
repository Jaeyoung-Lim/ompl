/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2025, Autonomous Systems Laboratory, ETH Zurich
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
 *   * Neither the name of the ETH Zurich nor the names of its
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

/* Author: Jaeyoung Lim */

#include "ompl/base/spaces/TrochoidStateSpace.h"
#include "ompl/base/SpaceInformation.h"
#include "ompl/util/Exception.h"
#include <array>
#include <cmath>
#include <queue>
#include <algorithm>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/roots.hpp>

using namespace ompl::base;

namespace
{
    constexpr double twopi = 2. * boost::math::constants::pi<double>();
    const double TROCHOID_EPS = 1e-6;

    double trochoid_delx(double omega, double dir, double delt, double phi0, double wind_ratio){
        return 1.0/(dir*omega)*(std::sin(dir*omega*delt+phi0) - std::sin(phi0)) + wind_ratio*delt;
    }

    double straight_delx(double wind_ratio, double delt, double hInitRad ) {
      return ( std::cos(hInitRad) + wind_ratio )*delt;
    }

    double straight_dely(double delt, double hInitRad ) {
      return std::sin(hInitRad)*delt;
    }

    double trochoid_dely(double omega, double dir, double delt, double hInitRad ){
    return 1.0/(dir*omega)*(-std::cos(dir*omega*delt+hInitRad) + std::cos(hInitRad));
    }

    double trochoid_delh(double omega, double dir, double delt ) {
        return dir*omega*delt;
    }

    /**
     * @brief Compute end points of the trochoid path
     * 
     * @param delta_1 Turning direction of the starting trochoid
     * @param delta_2 Turning direction of the ending trochoid
     * @param tA      Time value of the starting trochoid
     * @param tBeta   Duration of first arc and line segment
     * @param T       Total duration of path
     * @param phi0    Initial heading
     * @param radius  Minimum turn radius in air-relative frame
     * @param wind_ratio Wind ratio
     * @param xf      Relative target position
     * @param yf      Relative target position
     * @param phif    Target heading
     */
    void computeEndpointBSB(double delta_1, double delta_2, double tA, double tBeta, double T, double phi0, double radius, double wind_ratio, double &xf, double &yf, double &phif) {
      // state after first trochoid straight
        double x0 = 0.0, y0 = 0.0;
        double omega = (1.0/radius);
        double xBinit = x0 + trochoid_delx(omega, delta_1, tA, phi0, wind_ratio);
        double yBinit = y0 + trochoid_dely(omega, delta_1, tA, phi0);
        double hBinit = phi0 + trochoid_delh(omega, delta_1, tA ); 

        double xBfinal, yBfinal, hBfinal;
        xBfinal = xBinit + straight_delx(wind_ratio, tBeta - tA, hBinit);
        yBfinal = yBinit + straight_dely(tBeta - tA, hBinit );
        hBfinal = hBinit;

        // final state 
        xf = xBfinal + trochoid_delx(omega, delta_2, T-tBeta , hBfinal, wind_ratio);
        yf = yBfinal + trochoid_dely(omega, delta_2, T-tBeta , hBfinal );
        phif = hBfinal + trochoid_delh(omega, delta_2, T - tBeta);
        phif = fmod(phif, twopi);
        if ( phif < 0 ){
            phif = phif + twopi;
        }
    }
 
    /**
     * @brief Check conditions for the BSB path type. trochoid_path is updated in case the path is shorter 
     * 
     * @param delta_1 Turning direction of the starting trochoid
     * @param delta_2 Turning direction of the ending trochoid
     * @param tA      Time value of the starting trochoid
     * @param tB      Time value of the ending trochoid
     * @param xt10    Constant defined in Eq. (20) in Techy (2009)
     * @param yt10    Constant defined in Eq. (21) in Techy (2009)
     * @param phi0    Initial heading
     * @param phit1   Initial heading in trochoid frame
     * @param xt20    Constaint defined in Eq. (22) in Techy (2009)
     * @param yt20    Constant defined in Eq. (23) in Techy (2009)
     * @param phit2   Target heading in trochoid frame 
     * @param xf      Relative target position
     * @param yf      Relative target position
     * @param radius  Minimum turn radius in air-relative frame
     * @param wind_ratio Wind ratio
     * @param periodic Whether the BSB path is periodic
     * @param trochoid_path Calculated trochoid path object
     * @return true   Returns true if the condiitons are conistent
     * @return false  Reutrns false if the computed path are degenerate
     */
    bool checkConditionsBSB(double delta_1, double delta_2, double tA, double tB, double xt10, double yt10, double phi0, double phit1, double xt20, double yt20, double phit2, double xf, double yf, double radius, double wind_ratio, bool periodic, TrochoidStateSpace::PathType &trochoid_path) {
        // check that tA, tB, and T are valid numbers 
        if ( std::isnan(tA) || std::isnan(tB) ){
            return false;
        }

        double xtAdot = cos(delta_1*(1./radius)*tA + phit1) + wind_ratio;
        double ytAdot = sin(delta_1*(1./radius)*tA + phit1);
        double xtBdot = cos(delta_2*(1./radius)*tB + phit2) + wind_ratio;
        double ytBdot = sin(delta_2*(1./radius)*tB + phit2);

        double xtA = (radius/delta_1)*sin(delta_1*(1./radius)*tA + phit1) + wind_ratio*tA + xt10;
        double ytA = (-radius/delta_1)*cos(delta_1*(1./radius)*tA + phit1) + yt10;
        double xtB = (radius/delta_2) *sin(delta_2*(1./radius)*tB + phit2) + wind_ratio*tB + xt20;
        double ytB = (-radius/delta_2)*cos(delta_2*(1./radius)*tB + phit2) + yt20;

        double eps = 0.001;
        double t2pi = twopi * radius;

        if ( tB < -t2pi || tB > t2pi ){
            return false;
        } 
        if ( tA < 0 || tA > 2*t2pi ){
            return false;
        }
        if ( (xtB-xtA)*xtAdot < -eps ){
            return false;
        }  
        if ( (ytB-ytA)*ytAdot < -eps ){
            return false;
        } 
        if ( (xtB-xtA)*xtBdot < -eps ){
            return false;
        }  
        if ( (ytB-ytA)*ytBdot < -eps ){
            return false;
        } 
        // check that equal 
        if (fabs(ytB-ytA) > eps && fabs(ytAdot) > eps && fabs( (xtB - xtA)/(ytB-ytA) - xtAdot/ytAdot ) >= eps ){
            return false;
        }  
        if (fabs(ytB-ytA) > eps && fabs(ytBdot) > eps && fabs( (xtB - xtA)/(ytB-ytA) - xtBdot/ytBdot ) >= eps ){
            return false;
        }  
        // check endpoint
        double p = sqrt( (xtB-xtA)*(xtB-xtA) + (ytB-ytA)*(ytB-ytA) )/sqrt(xtBdot*xtBdot + ytBdot*ytBdot);
        double T = tA + p + (t2pi - tB);

        // check if the candidate satisfies the desired final (x,y) position 
        // double xFinal, yFinal, hFinalRad;
        double tBeta = tA + p;
        double endpoint_x, endpoint_y, endoint_phi;
        computeEndpointBSB(delta_1, delta_2, tA, tBeta, T, phi0, radius, wind_ratio, endpoint_x, endpoint_y, endoint_phi);
        if ( fabs(endpoint_x - xf) > eps || fabs(endpoint_y - yf) > eps ){
            return false;
        }

        if (!periodic) {
            if (T < trochoid_path.length()) {//Update path if it is shorter
                trochoid_path.length_[0] = tA;
                trochoid_path.length_[1] = p;
                trochoid_path.length_[2] = t2pi - tB;
            }
            return true;
        } else {
            if (T < trochoid_path.length() && T > 0.1) {//Update path if it is shorter
                trochoid_path.length_[0] = tA;
                trochoid_path.length_[1] = p;
                trochoid_path.length_[2] = t2pi - tB;
                return true;
            } else {
                return false;
            }
        }
    }



    /**
     * @brief Solve for trochoid RSL / LSR paths (delta_2 == -delta_1).
     *
     * The problem is solved in the air-relative frame, where the trochoids become circular arcs of radius r and the
     * target drifts with velocity -wind_ratio along the x axis (Wu and Chen, "Minimum-Time Paths for Dubins Airplane
     * in Steady Wind", JGCD 2025, Lemma 3). Let beta be the heading of the straight segment, d its length, alpha and
     * gamma the angles of the first and last arc, and c1, c3(T) the centers of the first and last circle, where
     * c3(T) - c1 = (Kx - w T, Ky). Tangency of the straight segment to both circles gives
     *
     *   c3(T) - c1 = d (cos beta, sin beta) + (delta_1 - delta_2) r (sin beta, -cos beta).
     *
     * The arc angles alpha = delta_1 (beta - phi0) and gamma = delta_2 (phif - beta) (both mod 2 pi) depend on beta
     * only, and T = r (alpha + gamma) + d. Projecting onto the segment direction yields
     *
     *   T(beta) = (r (alpha + gamma) + Kx cos beta + Ky sin beta) / (1 + w cos beta),
     *
     * and projecting onto its normal leaves the scalar equation
     *
     *   F(beta) = -(Kx - w T(beta)) sin beta + Ky cos beta + (delta_1 - delta_2) r = 0.
     *
     * F is smooth except where alpha or gamma wrap around (beta = phi0, beta = phif). All sign changes on each smooth
     * piece are bracketed on a fixed grid and refined with TOMS 748, so the cost is constant.
     *
     * @param x0 Start position x
     * @param y0 Start position y
     * @param phi0 Start heading
     * @param xf Target position x
     * @param yf Target position y
     * @param phif End heading
     * @param delta_1 Turning direction of starting arc
     * @param delta_2 Turning direction of ending arc, must be -delta_1
     * @param radius Minimum air-relative turning radius
     * @param wind_ratio Wind ratio
     * @param periodic Whether the path is periodic
     * @param path Solved trochoid path, updated if a shorter path is found
     */
    void trochoidBSBMixed(double x0, double y0, double phi0, double xf, double yf, double phif, double delta_1,
                          double delta_2, double radius, double wind_ratio, bool periodic,
                          TrochoidStateSpace::PathType &path)
    {
        const double r = radius, w = wind_ratio, s = delta_1 - delta_2;
        const double Kx = (xf - delta_2 * r * std::sin(phif)) - (x0 - delta_1 * r * std::sin(phi0));
        const double Ky = (yf + delta_2 * r * std::cos(phif)) - (y0 + delta_1 * r * std::cos(phi0));
        auto wrap = [](double a)
        {
            a = std::fmod(a, twopi);
            return a < 0. ? a + twopi : a;
        };

        // Smooth pieces of F are separated by the headings where the first or last arc vanishes
        std::array<double, 4> breaks{0., wrap(phi0), wrap(phif), twopi};
        std::sort(breaks.begin(), breaks.end());

        constexpr int gridPerTwoPi = 32;
        // In wind, an extra full loop on the last arc lets the target drift into reach, which can be the only
        // feasible RSL / LSR path (e.g. for periodic loiter paths), so gamma in [0, 4 pi) is searched as well.
        for (int extraLoops = 0; extraLoops <= 1; ++extraLoops)
        for (std::size_t piece = 0; piece + 1 < breaks.size(); ++piece)
        {
            const double lo = breaks[piece], hi = breaks[piece + 1];
            if (hi - lo < 1e-12)
                continue;
            // alpha and gamma are linear in beta on this piece
            const double mid = .5 * (lo + hi);
            const double alphaMid = wrap(delta_1 * (mid - phi0));
            const double gammaMid = wrap(delta_2 * (phif - mid)) + twopi * extraLoops;
            auto arcLength = [&](double beta)
            { return r * (alphaMid + gammaMid + (delta_1 - delta_2) * (beta - mid)); };
            auto totalTime = [&](double beta, double sb, double cb)
            { return (arcLength(beta) + Kx * cb + Ky * sb) / (1. + w * cb); };
            auto F = [&](double beta)
            {
                const double sb = std::sin(beta), cb = std::cos(beta);
                return -(Kx - w * totalTime(beta, sb, cb)) * sb + Ky * cb + s * r;
            };

            auto accept = [&](double beta)
            {
                const double sb = std::sin(beta), cb = std::cos(beta);
                const double T = totalTime(beta, sb, cb);
                const double arcs = arcLength(beta);
                const double d = T - arcs;
                if (d < -TROCHOID_EPS || T <= 0.)
                    return;
                const double tA = r * (alphaMid + delta_1 * (beta - mid));
                const double tC = r * (gammaMid - delta_2 * (beta - mid));
                if (tA < -TROCHOID_EPS || tC < -TROCHOID_EPS)
                    return;
                if (periodic ? (T <= 0.1 || T >= path.length()) : T >= path.length())
                    return;
                // Guard against spurious roots at the discontinuities by checking the endpoint in the ground frame
                double ex, ey, eh;
                computeEndpointBSB(delta_1, delta_2, std::max(tA, 0.), std::max(tA, 0.) + std::max(d, 0.), T, phi0,
                                   radius, wind_ratio, ex, ey, eh);
                if (std::abs(ex + x0 - xf) > 1e-3 || std::abs(ey + y0 - yf) > 1e-3)
                    return;
                path.length_[0] = std::max(tA, 0.);
                path.length_[1] = std::max(d, 0.);
                path.length_[2] = std::max(tC, 0.);
            };
            auto solveBracket = [&](double a, double b, double fa, double fb)
            {
                std::uintmax_t maxIter = 50;
                auto res = boost::math::tools::toms748_solve(F, a, b, fa, fb,
                                                             boost::math::tools::eps_tolerance<double>(48), maxIter);
                accept(.5 * (res.first + res.second));
            };

            constexpr int maxGrid = gridPerTwoPi + 4;
            const int n = std::max(4, static_cast<int>(std::ceil(gridPerTwoPi * (hi - lo) / twopi)));
            std::array<double, maxGrid + 1> bs, fs;
            for (int i = 0; i <= n; ++i)
            {
                bs[i] = lo + (hi - lo) * i / n;
                fs[i] = F(bs[i]);
            }
            for (int i = 1; i <= n; ++i)
            {
                if (fs[i - 1] == 0.)
                    accept(bs[i - 1]);
                else if (fs[i - 1] * fs[i] < 0.)
                    solveBracket(bs[i - 1], bs[i], fs[i - 1], fs[i]);
            }
            if (fs[n] == 0.)
                accept(bs[n]);

            // Two close roots can hide between grid points, where |F| has a local minimum that does not change sign
            // on the grid. Refine such minima if a parabola through the neighbouring samples suggests a crossing.
            for (int i = 1; i < n; ++i)
            {
                const double f0 = fs[i - 1], f1 = fs[i], f2 = fs[i + 1];
                if (f0 * f1 <= 0. || f1 * f2 <= 0. || std::abs(f1) > std::abs(f0) || std::abs(f1) > std::abs(f2))
                    continue;
                const double curvature = f0 - 2. * f1 + f2;  // second difference, in units of the grid step
                const double vertex = f1 - (f2 - f0) * (f2 - f0) / (8. * curvature);
                if (vertex * f1 > 0. && std::abs(vertex) > .25 * std::abs(f1))
                    continue;
                const double sign = f1 > 0. ? 1. : -1.;
                auto res = boost::math::tools::brent_find_minima([&](double b) { return sign * F(b); }, bs[i - 1],
                                                                 bs[i + 1], 40);
                const double bMin = res.first, fMin = sign * res.second;
                if (fMin * f1 >= 0.)
                    continue;
                solveBracket(bs[i - 1], bMin, f0, fMin);
                solveBracket(bMin, bs[i + 1], fMin, f2);
            }
        }
    }

    /**
     * @brief Solve for trochoid BSB paths
     *
     * @param x0 Start position x
     * @param y0 Start position y
     * @param phi0 Start heading
     * @param xf Target position x
     * @param yf Target position y
     * @param phif End heading
     * @param delta_1 Turning direction of starting arc
     * @param delta_2 Turning direction of ending arc
     * @param radius Minimum air-relative turning radius
     * @param wind_ratio Wind ratio
     * @param periodic Whether the path is periodic
     * @param path Solved trochoid path
     */
    void trochoidBSB(double x0, double y0, double phi0, double xf, double yf, double phif, double delta_1, double delta_2, double radius, double wind_ratio, bool periodic, TrochoidStateSpace::PathType &path)
    {

        double V_w = wind_ratio;
        double omega = 1.0/radius;
        double t2pi = twopi * radius;

        double phit1 = fmod(phi0, twopi);
        double xt10 = x0 - delta_1 * radius * std::sin(phit1);
        double yt10 = y0 + delta_1 * radius * std::cos(phit1);

        double phit2 = fmod(phif - delta_2 * omega * t2pi, twopi);
        double xt20 = xf - delta_2 * radius * std::sin(delta_2 * omega * t2pi + phit2) - V_w * t2pi;
        double yt20 = yf + delta_2 * radius * cos(delta_2 * omega * t2pi + phit2);

        if (delta_1 == delta_2) { // RSR/LSL: analytic solution
            for (int k_int = -2; k_int < 2; k_int++){
                double k = (double)k_int;
                double alpha = std::atan2(
                    yt20 - yt10, xt20 - xt10 + V_w * (fmod(phit1 - phit2, twopi) - k * twopi) / (delta_1 * omega));
                double tA = t2pi / (delta_1 * twopi) * (std::asin(V_w * std::sin(alpha)) + alpha - phit1);
                if ( tA > t2pi || tA < 0){
                    tA = tA - t2pi*floor(tA/t2pi);
                }
                double tB = tA + (fmod(phit1 - phit2, twopi) - k * twopi) / (delta_1 * twopi) * t2pi;
                checkConditionsBSB(delta_1, delta_2, tA, tB, xt10, yt10, phi0, phit1, xt20, yt20, phit2, xf, yf, radius, wind_ratio, periodic, path);
            }
        } else {  // RSL/LSR: scalar root finding in the air-relative frame
            trochoidBSBMixed(x0, y0, phi0, xf, yf, phif, delta_1, delta_2, radius, wind_ratio, periodic, path);
        }

        return;
    }

    /**
     * @brief Check conditions for BBB Path types
     * 
     * @param delta_1 Turning direction of starting trochoid
     * @param tA      Time value of the starting trochoid
     * @param tB      Time value of the ending trochoid
     * @param T       Total duration of path
     * @param xt10    Constant defined in Eq. (20) in Techy (2009)
     * @param yt10    Constant defined in Eq. (21) in Techy (2009)
     * @param phit1   Initial heading in trochoid frame
     * @param xf      Relative target position
     * @param yf      Relative target position
     * @param phif    Target heading
     * @param radius  Minimum turn radius in air-relative frame
     * @param wind_ratio Wind ratio
     * @param periodic Whether the BSB path is periodic
     * @param trochoid_path Calculated trochoid path object
     * @return true   Returns true if the condiitons are conistent
     * @return false  Reutrns false if the computed path are degenerate
     */
    bool checkConditionsBBB( double delta_1, double tA, double tB, double T, double xt10, double yt10, double phit1, double xf, double yf, double phif, double radius, double wind_ratio, bool periodic, TrochoidStateSpace::PathType &trochoid_path){

        if ( std::isnan(tA) || std::isnan(tB) || std::isnan(T) ){
            return false;
        }
        double omega = 1./radius;
        double t2pi = twopi * radius;
        double d2 = -delta_1;
        double d3 = delta_1;

        double phit3 = fmod(phif - d3 * omega * T, twopi);
        double xt30 = xf - 1.0/(d3*omega)*sin( phif ) - wind_ratio*T;
        double yt30 = yf + 1.0/(d3*omega)*cos( phif );

        double phit2 = 2*delta_1*omega*tA + phit1;
        double xt20 = xt30 - 2*1.0/(d2*omega)*sin(d2*omega*tB + phit2);
        double yt20 = yt30 + 2*1.0/(d2*omega)*cos(d2*omega*tB + phit2);


        double xt1tA = 1.0/(delta_1*omega)*sin(delta_1*omega*tA + phit1)
                                    + wind_ratio*tA + xt10;
        double yt1tA = -1.0/(delta_1*omega)*cos(delta_1*omega*tA + phit1) + yt10;

        double xt2tA = 1.0/(d2*omega)*sin(d2*omega*tA + phit2)
                                    + wind_ratio*tA + xt20;
        double yt2tA = -1.0/(d2*omega)*cos(d2*omega*tA + phit2) + yt20;

        double xt2tB = 1.0/(d2*omega)*sin(d2*omega*tB + phit2)
                                    + wind_ratio*tB + xt20;
        double yt2tB = -1.0/(d2*omega)*cos(d2*omega*tB + phit2) + yt20;
        double xt3tB = 1.0/(d3*omega)*sin(d3*omega*tB + phit3)
                                    + wind_ratio*tB + xt30;
        double yt3tB = -1.0/(d3*omega)*cos(d3*omega*tB + phit3) + yt30;  

        double xt1dottA = 1.0*cos(delta_1*omega*tA + phit1) + wind_ratio;
        double yt1dottA = 1.0*sin(delta_1*omega*tA + phit1);
        double xt2dottA = 1.0*cos(d2*omega*tA + phit2) + wind_ratio;
        double yt2dottA = 1.0*sin(d2*omega*tA + phit2);

        double xt2dottB= 1.0*cos(d2*omega*tB + phit2) + wind_ratio;
        double yt2dottB = 1.0*sin(d2*omega*tB + phit2);
        double xt3dottB = 1.0*cos(d3*omega*tB + phit3) + wind_ratio;
        double yt3dottB = 1.0*sin(d3*omega*tB + phit3);

        double eps = 0.001;


        // check that: 0 <= tA <= tB <= T <= 4.0*t2pi   
        if ( tA < 0 || tB < tA || T < tB || T > 4.0*t2pi ){
            return false;
        }
        if ( fabs(xt1tA - xt2tA) >= eps ){
            return false; 
        }
        if ( fabs(yt1tA - yt2tA) >= eps ){
            return false; 
        }
        if ( fabs(xt2tB - xt3tB) >= eps ){
            return false; 
        }
        if ( fabs(yt2tB - yt3tB) >= eps ){
            return false; 
        }
        if ( fabs(xt1dottA - xt2dottA) >= eps ){
            return false; 
        }
        if ( fabs(yt1dottA - yt2dottA) >= eps ){
            return false; 
        }
        if ( fabs(xt2dottB - xt3dottB) >= eps ){
            return false; 
        }
        if ( fabs(yt2dottB - yt3dottB) >= eps ){
            return false; 
        }

        if (!periodic) {
            if (T < trochoid_path.length()) {//Update path if it is shorter
                trochoid_path.length_[0] = tA;
                trochoid_path.length_[1] = tB - tA;
                trochoid_path.length_[2] = T - tB;
            }
            return true;
        } else {
            if (T < trochoid_path.length() && T > 0.1) {//Update path if it is shorter
                trochoid_path.length_[0] = tA;
                trochoid_path.length_[1] = tB - tA;
                trochoid_path.length_[2] = T - tB;
                return true;
            } else {
                return false;
            }
        }

    }

    /**
     * @brief Find all roots of a scalar function on [lo, hi] given the points that split it into intervals on which
     * the function is monotonic. A root is reported where the function changes sign over an interval (refined by
     * bisection) or vanishes at a break point.
     *
     * @param f       Function to find the roots of
     * @param breaks  Sorted break points, breaks[0] = lo and breaks[numBreaks - 1] = hi
     * @param numBreaks Number of break points
     * @param roots   Output array of roots, sorted in increasing order
     * @return Number of roots found
     */
    template <typename Function, std::size_t NB, std::size_t NR>
    std::size_t monotonicIntervalRoots(const Function &f, const std::array<double, NB> &breaks, std::size_t numBreaks,
                                       std::array<double, NR> &roots)
    {
        constexpr double zeroTol = 1e-12;
        constexpr int maxBisections = 100;
        std::size_t numRoots = 0;
        double fLo = f(breaks[0]);
        for (std::size_t i = 0; i + 1 < numBreaks && numRoots < NR; ++i)
        {
            double lo = breaks[i], hi = breaks[i + 1];
            double fHi = f(hi);
            if (std::abs(fLo) < zeroTol)
            {
                if (numRoots == 0 || roots[numRoots - 1] != lo)
                    roots[numRoots++] = lo;
            }
            else if (std::abs(fHi) >= zeroTol && (fLo < 0) != (fHi < 0))
            {
                double a = lo, b = hi, fa = fLo;
                for (int j = 0; j < maxBisections && b - a > zeroTol; ++j)
                {
                    double mid = 0.5 * (a + b), fMid = f(mid);
                    if ((fMid < 0) == (fa < 0))
                    {
                        a = mid;
                        fa = fMid;
                    }
                    else
                    {
                        b = mid;
                    }
                }
                roots[numRoots++] = 0.5 * (a + b);
            }
            fLo = fHi;
        }
        if (numRoots < NR && numBreaks > 1 && std::abs(fLo) < zeroTol)
            roots[numRoots++] = breaks[numBreaks - 1];
        return numRoots;
    }

    /**
     * @brief Solve for trochoid BBB paths.
     *
     * The problem is solved in the air-relative frame, where the trochoids become circular arcs and the target drifts
     * with velocity -wind_ratio along the x axis (Wu and Chen, "Closed-Form Solutions for Minimum-Time Paths of Dubins
     * Airplane in Steady Wind", 2024). Let alpha, beta, gamma be the angles of the three arcs and c1, c3(T) the centers
     * of the first and last circle. Since the middle circle is tangent to both, |c3(T) - c1| = 4 r sin(beta / 2), and
     * the heading constraint gives T = r (c + 2 beta + 2 pi k). Squaring yields a single equation in beta,
     *
     *   G(beta) = a1 beta^2 + a2 beta + a3 cos(beta) + a4 = 0,
     *
     * whose roots on [0, 2 pi] are all found with the improved bisection method (Algorithm 1 in Wu and Chen). G'' only
     * vanishes at cos(beta) = wind_ratio^2, so the intervals on which G' is monotonic are known in closed form. alpha
     * and gamma then follow from the direction of c3 - c1 and the heading constraint.
     *
     * @param x0 Start position x
     * @param y0 Start position y
     * @param phit1 Start heading
     * @param xf Target position x
     * @param yf Target position y
     * @param phif End heading
     * @param delta_1 Turn direction of starting arc
     * @param radius minimum air-relative turning radius
     * @param wind_ratio wind ratio
     * @param periodic Whether the path is periodic
     * @param path Computed path object
     */
    void trochoidBBB(double x0, double y0, double phit1, double xf, double yf, double phif, double delta_1,
                     double radius, double wind_ratio, bool periodic, TrochoidStateSpace::PathType &path)
    {
        const double w = wind_ratio;
        const double xt10 = x0 - delta_1 * radius * sin(phit1);
        const double yt10 = y0 + delta_1 * radius * cos(phit1);

        // vector from the first circle center to the last circle center at T = 0, normalized by radius
        const double xc = (xf - delta_1 * radius * sin(phif) - xt10) / radius;
        const double yc = (yf + delta_1 * radius * cos(phif) - yt10) / radius;

        // total turning, alpha + gamma = c + beta + 2 pi k
        double c = fmod(delta_1 * (phif - phit1), twopi);
        if (c < 0)
            c += twopi;

        // G'' (beta) = 8 (w^2 - cos(beta)) vanishes at most twice on [0, 2 pi]
        std::array<double, 4> secondDerivBreaks;
        std::size_t numSecondDerivBreaks = 0;
        secondDerivBreaks[numSecondDerivBreaks++] = 0.0;
        if (w * w < 1.0)
        {
            double b = std::acos(w * w);
            secondDerivBreaks[numSecondDerivBreaks++] = b;
            secondDerivBreaks[numSecondDerivBreaks++] = twopi - b;
        }
        secondDerivBreaks[numSecondDerivBreaks++] = twopi;

        // alpha + gamma in [0, 4 pi) and c + beta in [0, 4 pi] restrict k to {-1, 0, 1}
        for (int k = -1; k <= 1; ++k)
        {
            // horizontal center offset at T = r (c + 2 pi k), the remaining drift is -2 w beta
            const double e = xc - w * (c + twopi * k);
            auto G = [&](double beta)
            { return 4. * w * w * beta * beta - 4. * w * e * beta + e * e + yc * yc - 8. + 8. * cos(beta); };
            auto dG = [&](double beta) { return 8. * w * w * beta - 4. * w * e - 8. * sin(beta); };

            // G' is monotonic between consecutive roots of G'', hence has at most 3 roots (+1 at 2 pi)
            std::array<double, 4> dGRoots;
            std::size_t numDGRoots = monotonicIntervalRoots(dG, secondDerivBreaks, numSecondDerivBreaks, dGRoots);

            // G is monotonic between consecutive roots of G'
            std::array<double, 6> firstDerivBreaks;
            std::size_t numFirstDerivBreaks = 0;
            firstDerivBreaks[numFirstDerivBreaks++] = 0.0;
            for (std::size_t i = 0; i < numDGRoots; ++i)
                if (dGRoots[i] > firstDerivBreaks[numFirstDerivBreaks - 1] && dGRoots[i] < twopi)
                    firstDerivBreaks[numFirstDerivBreaks++] = dGRoots[i];
            firstDerivBreaks[numFirstDerivBreaks++] = twopi;

            std::array<double, 6> betaRoots;
            std::size_t numBetaRoots = monotonicIntervalRoots(G, firstDerivBreaks, numFirstDerivBreaks, betaRoots);

            for (std::size_t i = 0; i < numBetaRoots; ++i)
            {
                const double beta = betaRoots[i];
                const double T = radius * (c + 2. * beta + twopi * k);
                if (T < 0)
                    continue;
                // c3 - c1 = 4 r sin(beta / 2) (cos(psi), sin(psi)) with psi = phi1 - delta_1 beta / 2
                const double psi = atan2(yc, xc - w * T / radius);
                const double phi1 = psi + delta_1 * beta / 2.;
                const double phi2 = phi1 - delta_1 * beta;
                double alpha = fmod(delta_1 * (phi1 - phit1), twopi);
                if (alpha < 0)
                    alpha += twopi;
                double gamma = fmod(delta_1 * (phif - phi2), twopi);
                if (gamma < 0)
                    gamma += twopi;
                // the winding of alpha and gamma must match k; this rejects roots of G belonging to another k
                if (std::abs(alpha + gamma - (c + beta + twopi * k)) > 1e-6)
                    continue;
                const double tA = radius * alpha;
                const double tB = tA + radius * beta;
                checkConditionsBBB(delta_1, tA, tB, T, xt10, yt10, phit1, xf, yf, phif, radius, wind_ratio, periodic,
                                   path);
            }
        }
    }

    bool isLongPathCase(double x0, double y0, double phi0, double xf, double yf, double phif, double radius, double /* wind_ratio */)
    {   
        double dx = (xf - x0)/radius, dy = (yf - y0)/radius;
        double d = sqrt(dx*dx + dy*dy), theta = atan2(dy, dx);
        double alpha = phi0 - theta, beta = phif - theta;
        bool is_geometric_long_path = (std::abs(std::sin(alpha)) + std::abs(std::sin(beta)) +
            std::sqrt(4 - std::pow(std::cos(alpha) + std::cos(beta), 2)) - d) < 0;
        bool is_wind_blowingaway = dx < 0;
        return is_geometric_long_path && is_wind_blowingaway;
    }

    TrochoidStateSpace::PathType trochoidRSR(double x0, double y0, double phi0, double xf, double yf, double phif, double radius, double wind_ratio, bool periodic)
    {
        TrochoidStateSpace::PathType path(TrochoidStateSpace::dubinsPathType()[1]);
        trochoidBSB(x0, y0, phi0, xf, yf, phif, -1, -1, radius, wind_ratio, periodic, path);
        return path;
    }

    TrochoidStateSpace::PathType trochoidLSL(double x0, double y0, double phi0, double xf, double yf, double phif, double radius, double wind_ratio, bool periodic)
    {
        TrochoidStateSpace::PathType path(TrochoidStateSpace::dubinsPathType()[0]);
        trochoidBSB(x0, y0, phi0, xf, yf, phif, 1, 1, radius, wind_ratio, periodic, path);
        return path; 
    }

    TrochoidStateSpace::PathType trochoidRSL(double x0, double y0, double phi0, double xf, double yf, double phif, double radius, double wind_ratio, bool periodic)
    {
        TrochoidStateSpace::PathType path(TrochoidStateSpace::dubinsPathType()[2]);
        trochoidBSB(x0, y0, phi0, xf, yf, phif, -1, 1, radius, wind_ratio, periodic, path);
        return path;
    }

    TrochoidStateSpace::PathType trochoidLSR(double x0, double y0, double phi0, double xf, double yf, double phif, double radius, double wind_ratio, bool periodic)
    {
        TrochoidStateSpace::PathType path(TrochoidStateSpace::dubinsPathType()[3]);
        trochoidBSB(x0, y0, phi0, xf, yf, phif, 1, -1, radius, wind_ratio, periodic, path);
        return path;
    }

    TrochoidStateSpace::PathType trochoidLRL(double x0, double y0, double phi0, double xf, double yf, double phif, double radius, double wind_ratio, bool periodic)
    {
        TrochoidStateSpace::PathType path(TrochoidStateSpace::dubinsPathType()[5]);
        trochoidBBB(x0, y0, phi0, xf, yf, phif, 1, radius, wind_ratio, periodic, path);
        return path;
    }

    TrochoidStateSpace::PathType trochoidRLR(double x0, double y0, double phi0, double xf, double yf, double phif, double radius, double wind_ratio, bool periodic)
    {
        TrochoidStateSpace::PathType path(TrochoidStateSpace::dubinsPathType()[4]);
        trochoidBBB(x0, y0, phi0, xf, yf, phif, -1, radius, wind_ratio, periodic, path);
        return path;
    }

    TrochoidStateSpace::PathType getPath(const double x0, const double y0, const double phi0, const double xf, const double yf, const double phif, double radius, double wind_ratio, bool periodic)
    {
        if (fabs(x0 - xf) < TROCHOID_EPS && fabs(y0 - yf) < TROCHOID_EPS && fabs(phi0 - phif) < TROCHOID_EPS && !periodic)
            return {TrochoidStateSpace::dubinsPathType()[0], 0.0, 0.0};

        TrochoidStateSpace::PathType path(trochoidLSL(x0, y0, phi0, xf, yf, phif, radius, wind_ratio, periodic)), tmp(trochoidRSR(x0, y0, phi0, xf, yf, phif, radius, wind_ratio, periodic));
        double len, minLength = path.length();

        if ((len = tmp.length()) < minLength)
        {
            minLength = len;
            path = tmp;
        }
        tmp = trochoidRSL(x0, y0, phi0, xf, yf, phif, radius, wind_ratio, periodic);
        if ((len = tmp.length()) < minLength)
        {
            minLength = len;
            path = tmp;
        }
        tmp = trochoidLSR(x0, y0, phi0, xf, yf, phif, radius, wind_ratio, periodic);
        if ((len = tmp.length()) < minLength)
        {
            minLength = len;
            path = tmp;
        }
        if (!isLongPathCase(x0, y0, phi0, xf, yf, phif, radius, wind_ratio)) {
            tmp = trochoidRLR(x0, y0, phi0, xf, yf, phif, radius, wind_ratio, periodic);
            if ((len = tmp.length()) < minLength)
            {
                minLength = len;
                path = tmp;
            }
            tmp = trochoidLRL(x0, y0, phi0, xf, yf, phif, radius, wind_ratio, periodic);
            if ((len = tmp.length()) < minLength)
                path = tmp;
        }
        return path;
    }

    TrochoidStateSpace::PathType getPeriodicPath(const double phi0, double radius, double wind_ratio, double direction)
    {
        const bool periodic = true;
        TrochoidStateSpace::PathType path, tmp;
        double len, minLength;
        if (direction > 0.0)
        {
            path = trochoidLSL(0.0, 0.0, phi0, 0.0, 0.0, phi0, radius, wind_ratio, periodic);
            minLength = path.length();
            tmp = trochoidLSR(0.0, 0.0, phi0, 0.0, 0.0, phi0, radius, wind_ratio, periodic);
            if ((len = tmp.length()) < minLength)
            {
                minLength = len;
                path = tmp;
            }
            if (!isLongPathCase(0.0, 0.0, phi0, 0.0, 0.0, phi0, radius, wind_ratio))
            {
                tmp = trochoidLRL(0.0, 0.0, phi0, 0.0, 0.0, phi0, radius, wind_ratio, periodic);
                if (tmp.length() < minLength)
                    path = tmp;
            }
        }
        else
        {
            path = trochoidRSR(0.0, 0.0, phi0, 0.0, 0.0, phi0, radius, wind_ratio, periodic);
            minLength = path.length();
            tmp = trochoidRSL(0.0, 0.0, phi0, 0.0, 0.0, phi0, radius, wind_ratio, periodic);
            if ((len = tmp.length()) < minLength)
            {
                minLength = len;
                path = tmp;
            }
            if (!isLongPathCase(0.0, 0.0, phi0, 0.0, 0.0, phi0, radius, wind_ratio))
            {
                tmp = trochoidRLR(0.0, 0.0, phi0, 0.0, 0.0, phi0, radius, wind_ratio, periodic);
                if (tmp.length() < minLength)
                    path = tmp;
            }
        }
        return path;
    }

}  // namespace

namespace ompl::base
{
    std::ostream &operator<<(std::ostream &os, const TrochoidStateSpace::PathType &path)
    {
        os << "TrochoidPath[ type=";
        for (unsigned i = 0; i < 3; ++i)
            if (path.type_->at(i) == TrochoidStateSpace::TROCHOID_LEFT)
                os << "L";
            else if (path.type_->at(i) == TrochoidStateSpace::TROCHOID_STRAIGHT)
                os << "S";
            else
                os << "R";
        os << ", length=" << path.length_[0] << '+' << path.length_[1] << '+' << path.length_[2] << '=' << path.length()
           << ", reverse=" << path.reverse_ << " ]";
        return os;
    }
}  // namespace ompl::base

TrochoidStateSpace::PathType TrochoidStateSpace::trochoidRSR(double x0, double y0, double phi0, double xf, double yf, double phif, double radius, double wind_ratio, bool periodic)
{
    return ::trochoidRSR(x0, y0, phi0, xf, yf, phif, radius, wind_ratio, periodic);
}

TrochoidStateSpace::PathType TrochoidStateSpace::trochoidLSL(double x0, double y0, double phi0, double xf, double yf, double phif, double radius, double wind_ratio, bool periodic)
{
    return ::trochoidLSL(x0, y0, phi0, xf, yf, phif, radius, wind_ratio, periodic);
}

TrochoidStateSpace::PathType TrochoidStateSpace::trochoidRSL(double x0, double y0, double phi0, double xf, double yf, double phif, double radius, double wind_ratio, bool periodic)
{
    return ::trochoidRSL(x0, y0, phi0, xf, yf, phif, radius, wind_ratio, periodic);
}

TrochoidStateSpace::PathType TrochoidStateSpace::trochoidLSR(double x0, double y0, double phi0, double xf, double yf, double phif, double radius, double wind_ratio, bool periodic)
{
    return ::trochoidLSR(x0, y0, phi0, xf, yf, phif, radius, wind_ratio, periodic);
}

TrochoidStateSpace::PathType TrochoidStateSpace::trochoidLRL(double x0, double y0, double phi0, double xf, double yf, double phif, double radius, double wind_ratio, bool periodic)
{
    return ::trochoidLRL(x0, y0, phi0, xf, yf, phif, radius, wind_ratio, periodic);
}

TrochoidStateSpace::PathType TrochoidStateSpace::trochoidRLR(double x0, double y0, double phi0, double xf, double yf, double phif, double radius, double wind_ratio, bool periodic)
{
    return ::trochoidRLR(x0, y0, phi0, xf, yf, phif, radius, wind_ratio, periodic);
}

TrochoidStateSpace::PathType TrochoidStateSpace::periodicTrochoidLRL(double x0, double y0, double phi0, double radius, double wind_ratio)
{
    return trochoidLRL(x0, y0, phi0, x0, y0, phi0, radius, wind_ratio, true);
}

TrochoidStateSpace::PathType TrochoidStateSpace::periodicTrochoidRLR(double x0, double y0, double phi0, double radius, double wind_ratio)
{
    return trochoidRLR(x0, y0, phi0, x0, y0, phi0, radius, wind_ratio, true);
}

const std::vector<std::vector<TrochoidStateSpace::TrochoidPathSegmentType> >& TrochoidStateSpace::dubinsPathType() {
    static const std::vector<std::vector<TrochoidStateSpace::TrochoidPathSegmentType>> *pathType =
        new std::vector<std::vector<TrochoidStateSpace::TrochoidPathSegmentType>>(
            {{{TROCHOID_LEFT, TROCHOID_STRAIGHT, TROCHOID_LEFT},
              {TROCHOID_RIGHT, TROCHOID_STRAIGHT, TROCHOID_RIGHT},
              {TROCHOID_RIGHT, TROCHOID_STRAIGHT, TROCHOID_LEFT},
              {TROCHOID_LEFT, TROCHOID_STRAIGHT, TROCHOID_RIGHT},
              {TROCHOID_RIGHT, TROCHOID_LEFT, TROCHOID_RIGHT},
              {TROCHOID_LEFT, TROCHOID_RIGHT, TROCHOID_LEFT}}});
    return *pathType;
    }

double TrochoidStateSpace::distance(const State *state1, const State *state2) const
{
    return isSymmetric_ ? symmetricDistance(state1, state2, rho_, eta_, psi_w_) : distance(state1, state2, rho_, eta_, psi_w_);
}
double TrochoidStateSpace::distance(const State *state1, const State *state2, double radius, double wind_ratio, double wind_heading)
{
    return getPath(state1, state2, radius, wind_ratio, wind_heading).length();
}
double TrochoidStateSpace::symmetricDistance(const State *state1, const State *state2, double radius, double wind_ratio, double wind_heading)
{
    return std::min(getPath(state1, state2, radius, wind_ratio, wind_heading).length(), getPath(state2, state1, radius, wind_ratio, wind_heading).length());
}

void TrochoidStateSpace::interpolate(const State *from, const State *to, const double t, State *state) const
{
    bool firstTime = true;
    PathType path;
    interpolate(from, to, t, firstTime, path, state);
}

void TrochoidStateSpace::interpolate(const State *from, const State *to, const double t, bool &firstTime,
                                   PathType &path, State *state) const
{
    if (firstTime)
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
        if (equalStates(from, to)) {
            path = getPath(from, to, rho_, eta_, psi_w_, true);
        }
        else {
            path = getPath(from, to, rho_, eta_, psi_w_);
        }
        if (isSymmetric_)
        {
            PathType path2(getPath(to, from));
            if (path2.length() < path.length())
            {
                path2.reverse_ = true;
                path = path2;
            }
        }
        firstTime = false;
    }
    interpolate(from, path, t, state, rho_, eta_, psi_w_);
}

void TrochoidStateSpace::interpolate(const State *from, const PathType &path, double t, State *state,
                                   double radius, double wind_ratio, double wind_heading) const
{
    auto *s = allocState()->as<StateType>();
    double seg = t * path.length(), phi, v, x_t10, y_t10, delta;
    s->setXY(0., 0.);
    s->setYaw(from->as<StateType>()->getYaw() - wind_heading);
    if (!path.reverse_) {
        for (unsigned int i = 0; i < 3 && seg > 0; ++i)
        {
            v = std::min(seg, path.length_[i]);
            phi = s->getYaw();
            seg -= v;
            switch (path.type_->at(i))
            {
                case TROCHOID_LEFT:
                    delta = 1;
                    x_t10 = s->getX() - (radius* delta) * sin(phi);
                    y_t10 = s->getY() + (radius * delta) * cos(phi);
                    s->setXY(x_t10 + (radius * delta) * sin(delta *(1.0/radius) * v + phi)  + wind_ratio * v, \
                    y_t10 - (radius * delta) * cos(delta *(1.0/radius) * v + phi));
                    s->setYaw(phi + delta *(1.0/radius) *v);
                    break;
                case TROCHOID_RIGHT:
                    delta = -1;
                    x_t10 = s->getX() - (radius* delta) * sin(phi);
                    y_t10 = s->getY() + (radius * delta) * cos(phi);
                    s->setXY(x_t10 + (radius * delta) * sin(delta *(1.0/radius) * v + phi)  + wind_ratio * v, \
                    y_t10 - (radius * delta) * cos(delta *(1.0/radius) * v + phi));
                    s->setYaw(phi + delta *(1.0/radius) *v);
                    break;
                case TROCHOID_STRAIGHT:
                    s->setXY(s->getX() + v * cos(phi) + v * wind_ratio, s->getY() + v * sin(phi));
                    break;
            }
        }
    }

    state->as<StateType>()->setX(s->getX()* std::cos(wind_heading) - s->getY()* std::sin(wind_heading) + from->as<StateType>()->getX());
    state->as<StateType>()->setY(s->getX()* std::sin(wind_heading) + s->getY()* std::cos(wind_heading) + from->as<StateType>()->getY());
    getSubspace(1)->enforceBounds(s->as<SO2StateSpace::StateType>(1));
    state->as<StateType>()->setYaw(s->getYaw() + wind_heading);
    freeState(s);
}

unsigned int TrochoidStateSpace::validSegmentCount(const State *state1, const State *state2) const
{
    return StateSpace::validSegmentCount(state1, state2);
}

TrochoidStateSpace::PathType TrochoidStateSpace::getPath(const State *state1, const State *state2, bool periodic) const
{
    return getPath(state1, state2, rho_, eta_, psi_w_, periodic);
}

TrochoidStateSpace::PathType TrochoidStateSpace::getPath(const State *state1, const State *state2, double radius, double wind_ratio, double wind_heading, bool periodic)
{
    const auto *s1 = static_cast<const TrochoidStateSpace::StateType *>(state1);
    const auto *s2 = static_cast<const TrochoidStateSpace::StateType *>(state2);
    double x0 = s1->getX(), y0 = s1->getY(), psi_0 = s1->getYaw();
    double xf = s2->getX(), yf = s2->getY(), psi_f = s2->getYaw();

    double phi0 = psi_0 - wind_heading;
    double phif = psi_f - wind_heading;

    // Transform into trochoid frame
    double x_trochoid = (xf- x0) * std::cos(wind_heading) + (yf- y0) * std::sin(wind_heading);
    double y_trochoid = -(xf-x0) * std::sin(wind_heading) + (yf- y0) * std::cos(wind_heading);

    return ::getPath(0.0, 0.0, phi0, x_trochoid, y_trochoid, phif, radius, wind_ratio, periodic);
}

TrochoidStateSpace::PathType TrochoidStateSpace::getPeriodicPath(const State *state1, double direction) const
{
    return getPeriodicPath(state1, direction, rho_, eta_, psi_w_);
}

TrochoidStateSpace::PathType TrochoidStateSpace::getPeriodicPath(const State *state1, double direction, double radius, double wind_ratio, double wind_heading)
{
    const auto *s1 = static_cast<const TrochoidStateSpace::StateType *>(state1);
    double phi0 = s1->getYaw() - wind_heading;
    return ::getPeriodicPath(phi0, radius, wind_ratio, direction);
}
