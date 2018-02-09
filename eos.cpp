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

double IGEOS::GetP_eps(double index, double rho, double eps)
{ 
   // Materials are indexed by their gamma in the case of ideal gas.
   double gamma = index;
   return (gamma - 1.0) * rho * eps;
}

double IGEOS::GetSoundSpeed_eps(double index, double rho, double eps)
{
   // Materials are indexed by their gamma in the case of ideal gas.
   double gamma = index;
   return sqrt(gamma * (gamma-1.0) * eps);
}

} // namespace nth

} // namespace mfem

#endif // MFEM_USE_MPI
