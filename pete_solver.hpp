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

#ifndef MFEM_NTH_PETE_SOLVER
#define MFEM_NTH_PETE_SOLVER

#include "mfem.hpp"
#include "laghos_solver.hpp"
#include "eos.hpp"

#ifdef MFEM_USE_MPI

#include <memory>
#include <iostream>
#include <fstream>

namespace mfem
{

namespace nth
{

// Given a solutions state (x, v, e), this class performs all necessary
// computations to evaluate the new slopes (dx_dt, dv_dt, de_dt).
class PETELagHydroOperator : public hydrodynamics::LagrangianHydroOperator
{
protected:
   // General multi-material equation of state (materials indexed by gamma).
   nth::EOS *eos;
    
   virtual void ComputeMaterialProperties(int nvalues, const double gamma[],
                                          const double rho[], const double e[],
                                          double p[], double cs[]) const									  
   {
      for (int v = 0; v < nvalues; v++)
      {
         //p[v]  = (gamma[v] - 1.0) * rho[v] * e[v];
         //cs[v] = sqrt(gamma[v] * (gamma[v]-1.0) * e[v]);
         // The materials are indexed by gamma in eos.
		 p[v]  = eos->GetP_eps(gamma[v], rho[v], e[v]);
         cs[v] = eos->GetSoundSpeed_eps(gamma[v], rho[v], e[v]);
      }
   }

public:
   PETELagHydroOperator(int size, ParFiniteElementSpace &h1_fes,
                        ParFiniteElementSpace &l2_fes,
                        Array<int> &essential_tdofs, ParGridFunction &rho0,
                        int source_type_, double cfl_,
                        Coefficient *material_, bool visc, bool pa,
                        double cgt, int cgiter, nth::EOS *eos_)
   : LagrangianHydroOperator(size, h1_fes, l2_fes, essential_tdofs, rho0,
                             source_type_, cfl_, material_, visc, pa, cgt, 
                             cgiter), eos(eos_) { }

   ~PETELagHydroOperator() {}
};

} // namespace nth

} // namespace mfem

#endif // MFEM_USE_MPI

#endif // MFEM_NTH_PETE_SOLVER
