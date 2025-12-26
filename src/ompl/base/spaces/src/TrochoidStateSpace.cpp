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
#include <queue>
#include <boost/math/constants/constants.hpp>

using namespace ompl::base;

namespace
{
    constexpr double twopi = 2. * boost::math::constants::pi<double>();
    const double TROCHOID_EPS = 1e-6;

    double trochoid_delx(double omega, double dir, double delt, double phi0, double wind_ratio){
        return .0/(dir*omega)*(std::sin(dir*omega*delt+phi0) - std::sin(phi0)) + wind_ratio*delt;
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
        if ( fabs(endpoint_x - xf) > eps && fabs(endpoint_y - yf) > eps ){
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
     * @brief Solve BSB Trochoid path. Closed form solution for delta_1==delta_2, or numerically solved with Eq. 39 in
     * Techy (2009)
     *
     * @param p0      Guess value of t_A for root solving
     * @param delta_1 Turning direction of the starting trochoid
     * @param delta_2 Turning direction of the ending trochoid
     * @param k       Number of turns
     * @param xt10    Constant defined in Eq. (20) in Techy (2009)
     * @param yt10    Constant defined in Eq. (21) in Techy (2009)
     * @param phit1   Initial heading in trochoid frame
     * @param xt20    Constaint defined in Eq. (22) in Techy (2009)
     * @param yt20    Constant defined in Eq. (23) in Techy (2009)
     * @param phit2   Target heading in trochoid frame
     * @param radius  Minimum turn radius in air-relative frame
     * @param wind_ratio Wind ratio
     * @return double
     */
    double fixedPointBSB( double p0, double delta_1, double delta_2, double k, double xt10, double yt10, double phit1, double xt20, double yt20, double phit2, double radius, double wind_ratio )
    {
        constexpr double tol = 0.0001;
        constexpr int N = 50;
        std::array<double, N> pvec;
        pvec[0] = p0;
        int i = 1;
        double omega = 1./radius;
        while (i < N)
        {
            double p = pvec[i-1];
            double E = (delta_1-delta_2) / ( delta_2*delta_1*omega ) - (yt20-yt10)/wind_ratio;
            double F = (xt20 - xt10) / wind_ratio + (delta_1 / delta_2 - 1) * p +
                       (phit1 - phit2 + k * twopi) / (delta_2 * omega);
            double G = (yt20 - yt10) + 1.0/wind_ratio*(delta_2-delta_1)/(delta_1*delta_2*omega);
            double f = E*cos(delta_1*omega*p + phit1) + F*sin(delta_1*omega*p + phit1) - G;
            double fBar = - E*sin(delta_1*omega*p + phit1)*delta_1*omega 
                        + F*cos(delta_1*omega*p + phit1)*delta_1*omega
                        + sin(delta_1*omega*p + phit1)*(delta_1/delta_2-1);
            pvec[i] = p - f/fBar;
            if ( fabs( pvec[i] - pvec[i-1]) < tol ){
            return pvec[i];
            }
            i++;
        }
        return std::numeric_limits<double>::quiet_NaN();
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
        } else {  // RSL/LSR: numerical solution
            for (int k_int = -2; k_int < 2; k_int++){
                double k = (double)k_int;
                constexpr int numTestPts = 10;
                std::array<double, numTestPts> rootsComputed;
                constexpr double sameRootEpsilon = 0.001;
                for ( int l = 0; l < numTestPts; l++){
                    bool rootAlreadyfound = false;
                    // initial guess
                    double p0 = (double)(l+1) * t2pi / (double)numTestPts;
                    // root solving 
                    double tA = fixedPointBSB( p0, delta_1, delta_2, k, xt10, yt10, phit1, xt20, yt20, phit2, radius, wind_ratio);
                    // check bounds
                    if (tA < 0 && fabs(tA) < TROCHOID_EPS) {
                        tA = 0.0;
                    }
                    if ( tA > t2pi || tA < 0){
                        tA = tA - t2pi*floor(tA/t2pi);
                    }
                    // check if the root has already been computed 
                    if ( l > 0 ){
                    for (int i = 0; i < l; i++){
                        if ( fabs( rootsComputed[i] - tA ) <= sameRootEpsilon ){
                        rootAlreadyfound = true;
                        } 
                    }
                    }
                    // store the computed root
                    rootsComputed[l] = tA;
                    // if the root is unique
                    if ( !rootAlreadyfound ){
                        double tB = delta_1 / delta_2 * tA + (phit1 - phit2 + k * twopi) / (delta_2 * omega);
                        checkConditionsBSB(delta_1, delta_2, tA, tB, xt10, yt10, phi0, phit1, xt20, yt20, phit2, xf, yf, radius, wind_ratio, periodic, path);
                    }
                }
            }
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
     * @brief Original 2D Newton-Raphson solver for BBB path type (kept for benchmarking).
     *
     * This is the original exhaustive 2D grid search implementation. It is kept for
     * performance comparison with the optimized 1D solver.
     *
     * @param p0 Guess value for tA
     * @param p1 Guess value for T
     * @param delta_1 Turning direction of start trochoid
     * @param k Number of turns
     * @param xt10 Constant defined in Eq. (A1)
     * @param yt10 Constant defined in Eq. (A2)
     * @param phit1 Constant defined in Eq. (A3)
     * @param xf Relative Target position x
     * @param yf Relative Target position y
     * @param phif Target heading
     * @param radius Minimum turn radius in air-relative frame
     * @param wind_ratio Wind ratio
     * @return Result of roots solved as a pair of doubles
     */
    std::pair<double, double> fixedpointBBB(double p0, double p1, double delta_1, double k, double xt10, double yt10,
                                            double phit1, double xf, double yf, double phif, double radius,
                                            double wind_ratio)
    {
        constexpr double tol = 0.0001;
        constexpr int N = 50;
        double d2 = -delta_1;
        double d3 = delta_1;
        double omega = 1.0/radius;
        std::array<double, N> pvec1;
        std::array<double, N> pvec2;
        pvec1[0] = p0;
        pvec2[0] = p1;
        int i = 1;
        while (i < N)
        {
            double tA = pvec1[i-1];
            double T = pvec2[i-1];
            // vector F = [f1; f2];
            double f1 =
                2. / (delta_1 * omega) * sin(delta_1 * omega * tA + phit1) + wind_ratio * T + xt10 - xf +
                1.0 / (d3 * omega) * sin(phif) +
                2 / (d2 * omega) * sin(d2 * omega * T / 2 + (phif + phit1 + k * twopi) / 2 + delta_1 * omega * tA);
            double f2 = -2 * 1.0 / (delta_1 * omega) * cos(delta_1 * omega * tA + phit1) + yt10 - yf -
                        1.0 / (d3 * omega) * cos(phif) -
                        2 * 1.0 / (d2 * omega) *
                            cos(d2 * omega * T / 2 + (phif + phit1 + k * twopi) / 2 + delta_1 * omega * tA);
            // matrix FBar = [ a b ; c d ];
            double a = 2 * 1.0 *
                       (cos(delta_1 * omega * tA + phit1) -
                        cos(d2 * omega * T / 2 + (phif + phit1 + k * twopi) / 2 + delta_1 * omega * tA));
            double b =
                wind_ratio + 1.0 * cos(d2 * omega * T / 2 + (phif + phit1 + k * twopi) / 2 + delta_1 * omega * tA);
            double c = 2 * 1.0 *
                       (sin(delta_1 * omega * tA + phit1) -
                        sin(d2 * omega * T / 2 + (phif + phit1 + k * twopi) / 2 + delta_1 * omega * tA));
            double d = 1.0 * sin(d2 * omega * T / 2 + (phif + phit1 + k * twopi) / 2 + delta_1 * omega * tA);
            double det = a*d - b*c;
            // pvec = x - FBar*FBarInv;
            pvec1[i] = pvec1[i-1] - (d/det*f1 - b/det*f2);
            pvec2[i] = pvec2[i-1] - (-c/det*f1 + a/det*f2);
            // convergence criteria
            if ( (pvec1[i]-pvec1[i-1])*(pvec1[i]-pvec1[i-1]) 
                + (pvec2[i]-pvec2[i-1])*(pvec2[i]-pvec2[i-1]) < tol ){
                return {pvec1[i], pvec2[i]};
            }
            i++;
        }
        return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    }

    /**
     * @brief Original BBB path solver using exhaustive 2D grid search (kept for benchmarking).
     *
     * This is the original O(n²) implementation that tests a grid of initial conditions
     * and uses 2D Newton-Raphson. Kept for performance comparison.
     *
     * @param x0 Start position x
     * @param y0 Start position y
     * @param phit1 Start heading
     * @param xf Target position x
     * @param yf Target position y
     * @param phif End heading
     * @param delta_1 Turn direction of starting arc
     * @param radius Minimum air-relative turning radius
     * @param wind_ratio Wind ratio
     * @param periodic Whether the path is periodic
     * @param path Computed path object
     */
    void trochoidBBB_original(double x0, double y0, double phit1, double xf, double yf, double phif, double delta_1, double radius, double wind_ratio, bool periodic, TrochoidStateSpace::PathType &path) {
        double omega = (1/radius);
        double t2pi = twopi * radius;

        double xt10 = x0 - 1.0/(delta_1*omega)*sin(phit1);
        double yt10 = y0 + 1.0/(delta_1*omega)*cos(phit1);

        double delta_2 = - delta_1;

        // test a grid of initial conditions, grid resolution
        constexpr int numTestPts = 10;
        std::array<std::array<double, 3>, numTestPts * numTestPts> rootsComputed;
        double sameRootEpsilon = 0.001;
        for (int kint = -2; kint < 3; kint++){
            double k = (double)kint;
            for ( int l = 0; l < numTestPts; l++){
            for ( int n = 0; n < numTestPts; n++){
                bool rootAlreadyfound = false;
                // initial guess
                double tA_test = (double)(l+1)*t2pi/numTestPts;
                double T_test = (double)(n+1)*3*t2pi/numTestPts;
                auto [tA, T] =
                    fixedpointBBB(tA_test, T_test, delta_1, k, xt10, yt10, phit1, xf, yf, phif, radius, wind_ratio);
                double tB = tA + T / 2 + (phif - phit1 + k * twopi) / (2 * delta_2 * omega);

                int curRow = (l)*numTestPts+(n);
                // check if the root has already been computed 
                if ( l > 0 || n > 0 ){
                    for (int i = 0; i < curRow; i++){
                        if (    fabs( rootsComputed[i][0] - tA ) <= sameRootEpsilon 
                            && fabs( rootsComputed[i][1] - tB ) <= sameRootEpsilon 
                            && fabs( rootsComputed[i][2] - T ) <= sameRootEpsilon ){
                        rootAlreadyfound = true;
                        } 
                    }
                }
                // store the computed roots
                rootsComputed[curRow] = {tA, tB, T};

                // if the root is unique
                if ( !rootAlreadyfound ){
                    checkConditionsBBB(delta_1, tA, tB, T, xt10, yt10, phit1, xf, yf, phif, radius, wind_ratio, periodic, path);
                }
            }
            }
        }  
    }

    /**
     * @brief 1D Newton-Raphson solver for BBB paths.
     * 
     * Solves for tA using the constraint that both f1=0 and f2=0 must be satisfied.
     * T is computed from tA using the closed-form solution derived from f2=0.
     * 
     * The solver handles the multi-valued nature of arccos by tracking which branch
     * (arccos or -arccos, plus 2π offsets) to use based on the sign parameter.
     * 
     * @param tA_init Initial guess for tA
     * @param sign Sign for arccos branch (+1 for arccos, -1 for -arccos)
     * @param theta_shift Additional 2π shift for theta_mid
     * @param k The k parameter for phase offset
     * @param C_y Precomputed y-constraint constant
     * @param delta_1 Direction parameter
     * @param omega 1/radius
     * @param phit1 Initial heading
     * @param phif Final heading
     * @param xt10 x-coordinate of first turning center
     * @param xf Final x position
     * @param wind_ratio Wind speed ratio
     * @return Pair of (tA, T) if converged, (NaN, NaN) otherwise
     */
    std::pair<double, double> solveBBB_1D(double tA_init, int sign, double theta_shift, double k, double C_y,
                                           double delta_1, double omega, double phit1, double phif,
                                           double xt10, double xf, double wind_ratio)
    {
        constexpr double tol = 1e-8;
        constexpr int maxIter = 30;
        
        double d2 = -delta_1;
        double d3 = delta_1;
        double phi_offset = (phif + phit1 + k * twopi) / 2.0;
        double signD = static_cast<double>(sign);
        
        double tA = tA_init;
        
        for (int iter = 0; iter < maxIter; iter++) {
            double theta1 = delta_1 * omega * tA + phit1;
            double sin_theta1 = sin(theta1);
            double cos_theta1 = cos(theta1);
            
            // Compute cos(θ_mid) from f2=0 constraint
            double cos_theta_mid = (C_y - 2.0 / (delta_1 * omega) * cos_theta1) * (d2 * omega / 2.0);
            
            // Check if cos_theta_mid is in valid range
            if (fabs(cos_theta_mid) > 1.0) {
                return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
            }
            
            // theta_mid = sign * arccos(cos_theta_mid) + theta_shift
            // where sign ∈ {-1, +1} and theta_shift ∈ {-2π, 0, 2π}
            double arccos_val = acos(cos_theta_mid);
            double theta_mid = signD * arccos_val + theta_shift;
            double sin_theta_mid = sin(theta_mid);
            double abs_sin_theta_mid = fabs(sin_theta_mid);
            
            // Compute T from closed-form
            double T = 2.0 * (theta_mid - phi_offset - delta_1 * omega * tA) / (d2 * omega);
            
            // Check T bounds
            if (T < 0 || T > 4.0 * twopi / omega) {
                return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
            }
            
            // Compute residual g = f1(tA, T(tA))
            double g = 2.0 / (delta_1 * omega) * sin_theta1 + wind_ratio * T + xt10 - xf +
                       1.0 / (d3 * omega) * sin(phif) + 2.0 / (d2 * omega) * sin_theta_mid;
            
            // Check convergence
            if (fabs(g) < tol) {
                return {tA, T};
            }
            
            // Avoid division by zero
            if (abs_sin_theta_mid < 1e-10) {
                return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
            }
            
            // Compute derivative dg/dtA using chain rule
            // theta_mid = sign * arccos(cos_theta_mid) + theta_shift
            // d(arccos(u))/du = -1/sqrt(1-u²) = -1/|sin(arccos(u))|
            // d(cos_theta_mid)/dtA = d2 * omega * sin(theta1)
            // d(theta_mid)/dtA = sign * (-1/|sin(arccos_val)|) * d(cos_theta_mid)/dtA
            //                  = -sign / |sin(arccos_val)| * d2 * omega * sin(theta1)
            // Note: |sin(arccos_val)| = sqrt(1 - cos²_theta_mid), but we can use sin_theta_mid if careful with sign
            double sin_arccos_val = sqrt(1.0 - cos_theta_mid * cos_theta_mid);  // Always positive
            if (sin_arccos_val < 1e-10) {
                return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
            }
            double dtheta_mid_dtA = -signD / sin_arccos_val * d2 * omega * sin_theta1;
            
            // dT/dtA = 2 * (dtheta_mid/dtA - delta_1 * omega) / (d2 * omega)
            double dT_dtA = 2.0 * (dtheta_mid_dtA - delta_1 * omega) / (d2 * omega);
            
            // dg/dtA = d(f1)/d(tA) where f1 depends on tA directly and through T(tA) and theta_mid
            // = 2*cos(theta1) + (2/(d2*omega))*cos(theta_mid)*dtheta_mid/dtA + wind_ratio*dT/dtA
            double dg_dtA = 2.0 * cos_theta1 + 
                           (2.0 / (d2 * omega)) * cos_theta_mid * dtheta_mid_dtA + 
                           wind_ratio * dT_dtA;
            
            // Avoid division by zero or near-zero derivative
            if (fabs(dg_dtA) < 1e-10) {
                return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
            }
            
            // Newton update
            tA = tA - g / dg_dtA;
            
            // Keep tA in reasonable bounds
            if (tA < 0) tA = 0;
            if (tA > 4.0 * twopi / omega) tA = 4.0 * twopi / omega;
        }
        
        // Did not converge
        return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    }

    /**
     * @brief Optimized BBB path solver using true 1D Newton-Raphson.
     * 
     * This approach:
     * 1. Eliminates T analytically from f2=0
     * 2. Solves the resulting 1D equation in tA using Newton-Raphson
     * 
     * The search space is: 
     * - tA initial guesses from 0 to t2pi (10 points)
     * - k from -2 to 2 (5 values)
     * - arccos sign: +1 or -1 (2 values)
     * - theta_shift: -2π, 0, 2π (3 values)
     * 
     * Total: 10 * 5 * 2 * 3 = 300 initial guesses, compared to 10*10*5 = 500 for original
     */
    void trochoidBBB(double x0, double y0, double phit1, double xf, double yf, double phif, double delta_1, double radius, double wind_ratio, bool periodic, TrochoidStateSpace::PathType &path) {
        double omega = 1.0 / radius;
        double t2pi = twopi * radius;

        double xt10 = x0 - 1.0 / (delta_1 * omega) * sin(phit1);
        double yt10 = y0 + 1.0 / (delta_1 * omega) * cos(phit1);

        double d2 = -delta_1;
        double d3 = delta_1;
        
        // Precompute y-constraint constant: C_y = yt10 - yf - (1/(d3*omega))*cos(phif)
        double C_y = yt10 - yf - 1.0 / (d3 * omega) * cos(phif);

        // Grid resolution for tA search
        constexpr int numTestPts = 10;
        
        // theta_shift values to try
        constexpr double shifts[] = {-twopi, 0.0, twopi};
        constexpr int numShifts = 3;
        
        // Track unique roots to avoid duplicates  
        std::array<std::array<double, 3>, numTestPts * 5 * 2 * numShifts> rootsComputed;
        int numRootsFound = 0;
        constexpr double sameRootEpsilon = 0.001;
        
        for (int kint = -2; kint < 3; kint++) {
            double k = static_cast<double>(kint);
            
            for (int l = 0; l < numTestPts; l++) {
                double tA_init = static_cast<double>(l + 1) * t2pi / numTestPts;
                
                // Try both arccos signs
                for (int sign = -1; sign <= 1; sign += 2) {
                    // Try different theta_shift values
                    for (int si = 0; si < numShifts; si++) {
                        double theta_shift = shifts[si];
                        
                        auto [tA, T] = solveBBB_1D(tA_init, sign, theta_shift, k, C_y, delta_1, omega, 
                                                   phit1, phif, xt10, xf, wind_ratio);
                        
                        if (std::isnan(tA) || std::isnan(T)) {
                            continue;
                        }
                        
                        double tB = tA + T / 2 + (phif - phit1 + k * twopi) / (2 * d2 * omega);
                        
                        // Check if this root was already found
                        bool rootAlreadyFound = false;
                        for (int i = 0; i < numRootsFound; i++) {
                            if (fabs(rootsComputed[i][0] - tA) <= sameRootEpsilon &&
                                fabs(rootsComputed[i][1] - tB) <= sameRootEpsilon &&
                                fabs(rootsComputed[i][2] - T) <= sameRootEpsilon) {
                                rootAlreadyFound = true;
                                break;
                            }
                        }
                        
                        if (!rootAlreadyFound && numRootsFound < static_cast<int>(rootsComputed.size())) {
                            rootsComputed[numRootsFound] = {tA, tB, T};
                            numRootsFound++;
                            
                            checkConditionsBBB(delta_1, tA, tB, T, xt10, yt10, phit1, 
                                              xf, yf, phif, radius, wind_ratio, periodic, path);
                        }
                    }
                }
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

    TrochoidStateSpace::PathType getPath(const double x0, const double y0, const double phi0, const double xf, const double yf, const double phif, double radius, double wind_ratio, bool periodic)
    {
        if (fabs(x0 - xf) < TROCHOID_EPS && fabs(y0 - yf) < TROCHOID_EPS && fabs(phi0 - phif) < TROCHOID_EPS && !periodic)
            return {TrochoidStateSpace::dubinsPathType()[0], 0.0, 0.0};

        TrochoidStateSpace::PathType path(TrochoidStateSpace::trochoidLSL(x0, y0, phi0, xf, yf, phif, radius, wind_ratio, periodic)), tmp(TrochoidStateSpace::trochoidRSR(x0, y0, phi0, xf, yf, phif, radius, wind_ratio, periodic));
        double len, minLength = path.length();

        if ((len = tmp.length()) < minLength)
        {
            minLength = len;
            path = tmp;
        }
        tmp = TrochoidStateSpace::trochoidRSL(x0, y0, phi0, xf, yf, phif, radius, wind_ratio, periodic);
        if ((len = tmp.length()) < minLength)
        {
            minLength = len;
            path = tmp;
        }
        tmp = TrochoidStateSpace::trochoidLSR(x0, y0, phi0, xf, yf, phif, radius, wind_ratio, periodic);
        if ((len = tmp.length()) < minLength)
        {
            minLength = len;
            path = tmp;
        }
        // if (!isLongPathCase(x0, y0, phi0, xf, yf, phif, radius, wind_ratio)) {
        //     tmp = trochoidRLR(x0, y0, phi0, xf, yf, phif, radius, wind_ratio, periodic);
        //     if ((len = tmp.length()) < minLength)
        //     {
        //         minLength = len;
        //         path = tmp;
        //     }
        //     tmp = trochoidLRL(x0, y0, phi0, xf, yf, phif, radius, wind_ratio, periodic);
        //     if ((len = tmp.length()) < minLength)
        //         path = tmp;
        // }
        return path;
    }

    TrochoidStateSpace::PathType getPeriodicPath(const double phi0, double radius, double wind_ratio, double direction)
    {
        bool periodic = true;
        TrochoidStateSpace::PathType path;
        if (direction > 0.0) {
            path = TrochoidStateSpace::trochoidLSL(0.0, 0.0, phi0, 0.0, 0.0, phi0, radius, wind_ratio, periodic);
            double len, minLength = path.length();

            TrochoidStateSpace::PathType tmp = TrochoidStateSpace::trochoidLSR(0.0, 0.0, phi0, 0.0, 0.0, phi0, radius, wind_ratio, periodic);
            if ((len = tmp.length()) < minLength)
            {
                minLength = len;
                path = tmp;
            }
            if (!isLongPathCase(0.0, 0.0, phi0, 0.0, 0.0, phi0, radius, wind_ratio)) {
                tmp = TrochoidStateSpace::trochoidLRL(0.0, 0.0, phi0, 0.0, 0.0, phi0, radius, wind_ratio, periodic);
                if ((len = tmp.length()) < minLength)
                    minLength = len;
                    path = tmp;
            }
        } else {
            path = TrochoidStateSpace::trochoidRSR(0.0, 0.0, phi0, 0.0, 0.0, phi0, radius, wind_ratio, periodic);
            double len, minLength = path.length();

            TrochoidStateSpace::PathType tmp = TrochoidStateSpace::trochoidRSL(0.0, 0.0, phi0, 0.0, 0.0, phi0, radius, wind_ratio, periodic);
            if ((len = tmp.length()) < minLength)
            {
                minLength = len;
                path = tmp;
            }
            if (!isLongPathCase(0.0, 0.0, phi0, 0.0, 0.0, phi0, radius, wind_ratio)) {
                tmp = TrochoidStateSpace::trochoidRLR(0.0, 0.0, phi0, 0.0, 0.0, phi0, radius, wind_ratio, periodic);
                if ((len = tmp.length()) < minLength)
                {
                    minLength = len;
                    path = tmp;
                }
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
    TrochoidStateSpace::PathType path(TrochoidStateSpace::dubinsPathType()[1]);
    trochoidBSB(x0, y0, phi0, xf, yf, phif, -1, -1, radius, wind_ratio, periodic, path);
    return path;
}

TrochoidStateSpace::PathType TrochoidStateSpace::trochoidLSL(double x0, double y0, double phi0, double xf, double yf, double phif, double radius, double wind_ratio, bool periodic)
{
    TrochoidStateSpace::PathType path(TrochoidStateSpace::dubinsPathType()[0]);
    trochoidBSB(x0, y0, phi0, xf, yf, phif, 1, 1, radius, wind_ratio, periodic, path);
    return path; 
}

TrochoidStateSpace::PathType TrochoidStateSpace::trochoidRSL(double x0, double y0, double phi0, double xf, double yf, double phif, double radius, double wind_ratio, bool periodic)
{
    TrochoidStateSpace::PathType path(TrochoidStateSpace::dubinsPathType()[2]);
    trochoidBSB(x0, y0, phi0, xf, yf, phif, -1, 1, radius, wind_ratio, periodic, path);
    return path;
}

TrochoidStateSpace::PathType TrochoidStateSpace::trochoidLSR(double x0, double y0, double phi0, double xf, double yf, double phif, double radius, double wind_ratio, bool periodic)
{
    TrochoidStateSpace::PathType path(TrochoidStateSpace::dubinsPathType()[3]);
    trochoidBSB(x0, y0, phi0, xf, yf, phif, 1, -1, radius, wind_ratio, periodic, path);
    return path;
}

TrochoidStateSpace::PathType TrochoidStateSpace::trochoidLRL(double x0, double y0, double phi0, double xf, double yf, double phif, double radius, double wind_ratio, bool periodic)
{
    TrochoidStateSpace::PathType path(TrochoidStateSpace::dubinsPathType()[5]);
    trochoidBBB(x0, y0, phi0, xf, yf, phif, 1, radius, wind_ratio, periodic, path);
    return path;
}

TrochoidStateSpace::PathType TrochoidStateSpace::trochoidRLR(double x0, double y0, double phi0, double xf, double yf, double phif, double radius, double wind_ratio, bool periodic)
{
    TrochoidStateSpace::PathType path(TrochoidStateSpace::dubinsPathType()[4]);
    trochoidBBB(x0, y0, phi0, xf, yf, phif, -1, radius, wind_ratio, periodic, path);
    return path;
}

TrochoidStateSpace::PathType TrochoidStateSpace::trochoidLRL_original(double x0, double y0, double phi0, double xf, double yf, double phif, double radius, double wind_ratio, bool periodic)
{
    TrochoidStateSpace::PathType path(TrochoidStateSpace::dubinsPathType()[5]);
    trochoidBBB_original(x0, y0, phi0, xf, yf, phif, 1, radius, wind_ratio, periodic, path);
    return path;
}

TrochoidStateSpace::PathType TrochoidStateSpace::trochoidRLR_original(double x0, double y0, double phi0, double xf, double yf, double phif, double radius, double wind_ratio, bool periodic)
{
    TrochoidStateSpace::PathType path(TrochoidStateSpace::dubinsPathType()[4]);
    trochoidBBB_original(x0, y0, phi0, xf, yf, phif, -1, radius, wind_ratio, periodic, path);
    return path;
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

TrochoidStateSpace::PathType TrochoidStateSpace::periodicTrochoidLRL(double x0, double y0, double phi0, double radius, double wind_ratio) {
    return trochoidLRL(x0, y0, phi0, x0, y0, phi0, radius, wind_ratio, true);
}

TrochoidStateSpace::PathType TrochoidStateSpace::periodicTrochoidRLR(double x0, double y0, double phi0, double radius, double wind_ratio) {
    return trochoidRLR(x0, y0, phi0, x0, y0, phi0, radius, wind_ratio, true);
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

TrochoidStateSpace::PathType TrochoidStateSpace::getPeriodicPath(const State *state1, double direction) const
{
    return getPeriodicPath(state1, rho_, eta_, psi_w_, direction);
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

TrochoidStateSpace::PathType TrochoidStateSpace::getPeriodicPath(const State *state1, double direction, double radius, double wind_ratio, double wind_heading)
{
    const auto *s1 = static_cast<const TrochoidStateSpace::StateType *>(state1);
    double phi0 = s1->getYaw() - wind_heading;

    return ::getPeriodicPath(phi0, radius, wind_ratio, direction);
}
