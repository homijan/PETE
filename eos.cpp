// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
// reserved. See files LICENSE and NOTICE for details.
//
// This file is part of CEED, a collection of benchmarks, miniapps, software
// libraries and APIs for efficient high-order finite element and spectral
// element discretizations for exascale applications. For more information and
// source code availability see http://github.com/ceed.
//
// The CEED research is supported by the Exascale Computing Project 17-SC-20-SC,
// a collaborative effort of two U.S. Department of Energy organizations (Office
// of Science and the National Nuclear Security Administration) responsible for
// the planning and preparation of a capable exascale ecosystem, including
// software, applications, hardware, advanced system engineering and early
// testbed platforms, in support of the nation's exascale computing imperative.

#include "eos.hpp"

#ifdef MFEM_USE_MPI

using namespace std;

namespace mfem
{

namespace nth
{

double IGEOS::GetGamma_eps(double index, double rho, double eps)
{
   // Materials are indexed by their gamma in the case of ideal gas.
   double gamma = index;
   return gamma;
}

double IGEOS::GetP_eps(double index, double rho, double eps)
{ 
   // Ideal gas provides a simple isentrop based sound speed.
   double gamma = GetGamma_eps(index, rho, eps);
   return (gamma - 1.0) * rho * eps;
}

double IGEOS::GetSoundSpeed_eps(double index, double rho, double eps)
{
   // Ideal gas provides a simple isentrop based sound speed.
   double gamma = GetGamma_eps(index, rho, eps);
   return sqrt(gamma * (gamma-1.0) * eps);
}

double HEREOS::GetT_eps(double index, double rho, double eps)
{ 
   // HotFix
   rho = max(1e-7, rho);
   // Index of material within hereos tables.
   int hindex = 1; 
   // TODO works only for electrons
   return her_eos->Get_Te_from_epse(hindex, rho, eps);
   //return her_eos->Get_T_from_eps(hindex, rho, eps);
}

double HEREOS::GetP_eps(double index, double rho, double eps)
{ 
   // HotFix
   rho = max(1e-7, rho);
   // Index of material within hereos tables.
   int hindex = 1; 
   // TODO works only for electrons
   double Te = her_eos->Get_Te_from_epse(hindex, rho, eps);
   return her_eos->Get_pe(hindex, rho, Te);  
   //double T = her_eos->Get_T_from_eps(hindex, rho, eps);
   //return her_eos->Get_p(hindex, rho, T);
}

double HEREOS::GetGamma_eps(double index, double rho, double eps)
{
   // HotFix
   rho = max(1e-7, rho);
   // Index of material within hereos tables.
   int hindex = 1; 
   // TODO works only for electrons
   double Te = her_eos->Get_Te_from_epse(hindex, rho, eps);
   double pe = her_eos->Get_pe(hindex, rho, Te);
   double gamma = 1.0 + pe / rho / eps;  
   //double T = her_eos->Get_T_from_eps(hindex, rho, eps);
   //double p = her_eos->Get_p(hindex, rho, T);
   //double gamma = 1.0 + p / rho / eps;
   return gamma;
}

double HEREOS::GetSoundSpeed_eps(double index, double rho, double eps)
{
   // HotFix
   rho = max(1e-7, rho);
   // Sound speed approximation based on ideal gas isentrop.
   double gamma = GetGamma_eps(index, rho, eps);
   return sqrt(gamma * (gamma-1.0) * eps);
}

} // namespace nth

} // namespace mfem

#endif // MFEM_USE_MPI
