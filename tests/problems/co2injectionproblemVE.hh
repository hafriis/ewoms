// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Ewoms::Co2InjectionVEProblem
 */
#ifndef EWOMS_CO2_INJECTIONVE_PROBLEM_HH
#define EWOMS_CO2_INJECTIONVE_PROBLEM_HH

#include <ewoms/models/immiscible/immisciblemodel.hh>
#include <ewoms/linear/parallelamgbackend.hh>

#include <opm/material/fluidsystems/H2ON2FluidSystem.hpp>
#include <opm/material/fluidsystems/BrineCO2FluidSystem.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/RegularizedBrooksCoreyVE.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/thermal/SomertonThermalConductionLaw.hpp>
#include <opm/material/thermal/ConstantSolidHeatCapLaw.hpp>
#include <opm/material/binarycoefficients/Brine_CO2.hpp>
#include <opm/material/common/UniformTabulated2DFunction.hpp>
#include <opm/material/common/Unused.hpp>

//#include <dune/grid/yaspgrid.hh>
//#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <opm/grid/polyhedralgrid/dgfparser.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <sstream>
#include <iostream>
#include <string>

namespace Ewoms {
//! \cond SKIP_THIS
template <class TypeTag>
class Co2InjectionVEProblem;

namespace Co2Injection {
#include <opm/material/components/co2tables.inc>
}
//! \endcond
}

BEGIN_PROPERTIES

NEW_TYPE_TAG(Co2InjectionVEBaseProblem);

// declare the CO2 injection problem specific property tags
NEW_PROP_TAG(FluidSystemPressureLow);
NEW_PROP_TAG(FluidSystemPressureHigh);
NEW_PROP_TAG(FluidSystemNumPressure);
NEW_PROP_TAG(FluidSystemTemperatureLow);
NEW_PROP_TAG(FluidSystemTemperatureHigh);
NEW_PROP_TAG(FluidSystemNumTemperature);

NEW_PROP_TAG(MaxDepth);
NEW_PROP_TAG(Temperature);
NEW_PROP_TAG(SimulationName);

// Set the grid type
//SET_TYPE_PROP(Co2InjectionVEBaseProblem, Grid, Dune::YaspGrid<2>);
SET_TYPE_PROP(Co2InjectionVEBaseProblem, Grid, Dune::PolyhedralGrid< 2, 2 >);

// Set the problem property
SET_TYPE_PROP(Co2InjectionVEBaseProblem, Problem,
              Ewoms::Co2InjectionVEProblem<TypeTag>);

// Set fluid configuration
SET_PROP(Co2InjectionVEBaseProblem, FluidSystem)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Ewoms::Co2Injection::CO2Tables CO2Tables;

public:
    typedef Opm::BrineCO2FluidSystem<Scalar, CO2Tables> type;
    //typedef Opm::H2ON2FluidSystem<Scalar, /*useComplexRelations=*/false> type;
};

// Set the material Law
SET_PROP(Co2InjectionVEBaseProblem, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    enum { liquidPhaseIdx = FluidSystem::liquidPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Opm::TwoPhaseMaterialTraits<Scalar,
                                        /*wettingPhaseIdx=*/FluidSystem::liquidPhaseIdx,
                                        /*nonWettingPhaseIdx=*/FluidSystem::gasPhaseIdx> Traits;

public:
    // define the material law which is parameterized by effective
    // saturations
    typedef Opm::RegularizedBrooksCoreyVE<Traits> type;
};

// Set the thermal conduction law
SET_PROP(Co2InjectionVEBaseProblem, ThermalConductionLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    // define the material law parameterized by absolute saturations
    typedef Opm::SomertonThermalConductionLaw<FluidSystem, Scalar> type;
};

// set the energy storage law for the solid phase
SET_TYPE_PROP(Co2InjectionVEBaseProblem, SolidEnergyLaw,
              Opm::ConstantSolidHeatCapLaw<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Use the algebraic multi-grid linear solver for this problem
SET_TAG_PROP(Co2InjectionVEBaseProblem, LinearSolverSplice, ParallelAmgLinearSolver);

// Write the Newton convergence behavior to disk?
SET_BOOL_PROP(Co2InjectionVEBaseProblem, NewtonWriteConvergence, false);

// Enable gravity
SET_BOOL_PROP(Co2InjectionVEBaseProblem, EnableGravity, false);

// set the defaults for the problem specific properties
SET_SCALAR_PROP(Co2InjectionVEBaseProblem, FluidSystemPressureLow, 3e7);
SET_SCALAR_PROP(Co2InjectionVEBaseProblem, FluidSystemPressureHigh, 4e7);
SET_INT_PROP(Co2InjectionVEBaseProblem, FluidSystemNumPressure, 100);
SET_SCALAR_PROP(Co2InjectionVEBaseProblem, FluidSystemTemperatureLow, 290);
SET_SCALAR_PROP(Co2InjectionVEBaseProblem, FluidSystemTemperatureHigh, 500);
SET_INT_PROP(Co2InjectionVEBaseProblem, FluidSystemNumTemperature, 100);

SET_SCALAR_PROP(Co2InjectionVEBaseProblem, MaxDepth, 2500);
SET_SCALAR_PROP(Co2InjectionVEBaseProblem, Temperature, 293.15);
SET_STRING_PROP(Co2InjectionVEBaseProblem, SimulationName, "co2injectionVE");

// The default for the end time of the simulation
SET_SCALAR_PROP(Co2InjectionVEBaseProblem, EndTime, 1e4);

// The default for the initial time step size of the simulation
SET_SCALAR_PROP(Co2InjectionVEBaseProblem, InitialTimeStepSize, 250);

// The default DGF file to load
SET_STRING_PROP(Co2InjectionVEBaseProblem, GridFile, "data/co2injectionVE.dgf");

END_PROPERTIES

namespace Ewoms {
/*!
 * \ingroup TestProblems
 *
 * \brief Problem where \f$CO_2\f$ is injected under a low permeable
 *        layer at a depth of 2700m.
 *
 * The domain is sized 60m times 40m and consists of two layers, one
 * which is moderately permeable (\f$K = 10^{-12}\;m^2\f$) for \f$ y >
 * 22\; m\f$ and one with a lower intrinsic permeablility (\f$
 * K=10^{-13}\;m^2\f$) in the rest of the domain.
 *
 * \f$CO_2\f$ gets injected by means of a forced-flow boundary
 * condition into water-filled aquifer, which is situated 2700m below
 * sea level, at the lower-right boundary (\f$5m<y<15m\f$) and
 * migrates upwards due to buoyancy. It accumulates and eventually
 * enters the lower permeable aquitard.
 *
 * The boundary conditions applied by this problem are no-flow
 * conditions on the top bottom and right boundaries and a free-flow
 * boundary condition on the left. For the free-flow condition,
 * hydrostatic pressure is assumed.
 */
template <class TypeTag>
class Co2InjectionVEProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
        typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { numPhases = FluidSystem::numPhases };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { liquidPhaseIdx = FluidSystem::liquidPhaseIdx };
    enum { CO2Idx = FluidSystem::CO2Idx };
    enum { BrineIdx = FluidSystem::BrineIdx };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { contiCO2EqIdx = conti0EqIdx + CO2Idx };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, ThermalConductionLaw) ThermalConductionLaw;
    typedef typename GET_PROP_TYPE(TypeTag, SolidEnergyLawParams) SolidEnergyLawParams;
    typedef typename ThermalConductionLaw::Params ThermalConductionLawParams;

    //**********************HAF**************************************
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef Opm::MathToolbox<Evaluation> Toolbox;
    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    Co2InjectionVEProblem(Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        eps_ = 1e-6;

        temperatureLow_ = EWOMS_GET_PARAM(TypeTag, Scalar, FluidSystemTemperatureLow);
        temperatureHigh_ = EWOMS_GET_PARAM(TypeTag, Scalar, FluidSystemTemperatureHigh);
        nTemperature_ = EWOMS_GET_PARAM(TypeTag, unsigned, FluidSystemNumTemperature);

        pressureLow_ = EWOMS_GET_PARAM(TypeTag, Scalar, FluidSystemPressureLow);
        pressureHigh_ = EWOMS_GET_PARAM(TypeTag, Scalar, FluidSystemPressureHigh);
        nPressure_ = EWOMS_GET_PARAM(TypeTag, unsigned, FluidSystemNumPressure);

        maxDepth_ = EWOMS_GET_PARAM(TypeTag, Scalar, MaxDepth);
        temperature_ = EWOMS_GET_PARAM(TypeTag, Scalar, Temperature);

        // initialize the tables of the fluid system
        // FluidSystem::init();
        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);

        fineLayerBottom_ = 22.0;

        // intrinsic permeabilities
        fineK_ = this->toDimMatrix_(1e-13);
        //carseK_ = this->toDimMatrix_(1e-12);
        coarseK_ = this->toDimMatrix_(1e-13); //NOTE: No heterogeneity yet!!!
        computeIntegratedPermeabilities(); //NOTE: For testing without VE influence!!!

        // porosities
        finePorosity_ = 0.15;
        coarsePorosity_ = 0.15; //NOTE: No heterogeneity yet!!!

        // parameters for the Brooks-Corey law
        fineMaterialParams_.setEntryPressure(1e4);
        coarseMaterialParams_.setEntryPressure(5e3);
        fineMaterialParams_.setLambda(2.0);
        coarseMaterialParams_.setLambda(2.0);

        fineMaterialParams_.finalizePlain();
        coarseMaterialParams_.finalizePlain();

        typedef typename MaterialLaw::RegularizedBrooksCorey PlainLaw;

#warning TODO: calculate the VE parameters! --- Seems OK now..
        // residual saturations
        //Scalar srw = 0.27;
        Scalar srw = 0.2;
        Scalar srn = 0.2;
        //Scalar srw = 0.0;
        //Scalar srn = 0.0;
        fineMaterialParams_.setResidualSaturation(liquidPhaseIdx, srw); //OK???????
        fineMaterialParams_.setResidualSaturation(gasPhaseIdx, srn); //OK???????
        coarseMaterialParams_.setResidualSaturation(liquidPhaseIdx, srw); //OK???????
        coarseMaterialParams_.setResidualSaturation(gasPhaseIdx, srn); //OK???????

        //End points for the relative permeabilities:
        Scalar krnFine = PlainLaw::twoPhaseSatKrn(fineMaterialParams_, srw);//Should be changed???
        fineMaterialParams_.setKrnEndPoint(krnFine);
        Scalar krwFine = PlainLaw::twoPhaseSatKrw(fineMaterialParams_, 1.0-srn);//Should be changed???
        fineMaterialParams_.setKrwEndPoint(krwFine);
        Scalar krnC = PlainLaw::twoPhaseSatKrn(coarseMaterialParams_, srw);//Should be changed???
        coarseMaterialParams_.setKrnEndPoint(krnC);
        Scalar krwC = PlainLaw::twoPhaseSatKrw(coarseMaterialParams_, 1.0-srn);//Should be changed???
        coarseMaterialParams_.setKrwEndPoint(krwC);

        fineMaterialParams_.finalize();
        coarseMaterialParams_.finalize();
        
        // parameters for the somerton law thermal conduction
        computeThermalCondParams_(fineThermalCondParams_, finePorosity_);
        computeThermalCondParams_(coarseThermalCondParams_, coarsePorosity_);

        // assume constant heat capacity and granite
        solidEnergyLawParams_.setSolidHeatCapacity(790.0 // specific heat capacity of granite [J / (kg K)]
                                                   * 2700.0); // density of granite [kg/m^3]
        solidEnergyLawParams_.finalize();

        //**********************HAF**************************************
        initializeSmax_VE(); //NOTE: For testing without VE influence!!!
        initializeH_VE(); //NOTE: For testing without VE influence!!!
    }


    //**********************HAF**************************************
    void initializeH_VE()
    {
#warning TODO: Vil griddet fra topSurf og 2D griddet her ha samme nummerering? Hvis ikke har vi vel et problem...f.eks. i statementet H_VE_[idx] = verteqUtil.get_H_VE(compressedDofIdx) nedenfor.

        const auto& simulator = this->simulator();

        ElementContext elemCtx(simulator);
        const auto& vanguard = simulator.vanguard();
        const auto& gridView = vanguard.gridView();
        int numElements = gridView.size(/*codim=*/0);
        std::cout << " numElements= " <<  numElements << std::endl;
        H_VE_.resize(numElements, 1.0);

        // /* //NOTE: For testing without VE influence!!!   
        Dune::VerteqColumnUtility<Grid> verteqUtil (elemCtx.problem().simulator().vanguard().grid());
        auto elemIt = vanguard.gridView().template begin<0>();
        const auto& elemEndIt = vanguard.gridView().template end<0>();
        int idx = 0;
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& elem = *elemIt;

            elemCtx.updatePrimaryStencil(elem);
            const auto& stencil = elemCtx.stencil(0); //OK med 0 ???
            unsigned compressedDofIdx = elemCtx.globalSpaceIndex(0, 0); //OK med 0 i argumentene???
            //const auto& entity = stencil.entity(spaceIdx); //PRØV DENNE!!! spaceIdx er vanligvis 0 (tilsv. for timeIdx)
            //std::cout << "compressedDofIdx= " <<  compressedDofIdx << std::endl;
            //const auto& entity = stencil.entity(compressedDofIdx); //OK???NEI, blir feil her!!!
            //H_VE_.push_back(verteqUtil.get_H_VE(entity));
            //H_VE_[idx] = verteqUtil.get_H_VE(entity); //@HAF: Does NOT work. Seg. error!!!
            H_VE_[idx] = verteqUtil.get_H_VE(compressedDofIdx);
            //std::cout << " compressedDofIdx= " << compressedDofIdx << " H_VE_[idx]= " <<  H_VE_[idx] << std::endl;
            
            idx++;
        }
        // */

        //exit(0);
    }


    template <class Context>
    void setH_VE(const Context& context, unsigned dofIdx, unsigned timeIdx, Scalar value)
    {
        unsigned compressedDofIdx = context.globalSpaceIndex(dofIdx, timeIdx);
        H_VE_[compressedDofIdx] = value;
    }

    template <class Context>
    Scalar getH_VE(const Context& context, unsigned dofIdx, unsigned timeIdx) const
    {
        unsigned compressedDofIdx = context.globalSpaceIndex(dofIdx, timeIdx);
        return H_VE_[compressedDofIdx];
    }
    
    
    void initializeSmax_VE()
    {
        const auto& simulator = this->simulator();

        const auto& vanguard = simulator.vanguard();
        const auto& gridView = vanguard.gridView();
        int numElements = gridView.size(/*codim=*/0);
        SmaxVE_.resize(numElements, 0.0);

        //ElementContext elemCtx(simulator);
        //const auto& vanguard = simulator.vanguard();
        //auto elemIt = vanguard.gridView().template begin</*codim=*/0>();
        //const auto& elemEndIt = vanguard.gridView().template end</*codim=*/0>();
        //for (; elemIt != elemEndIt; ++elemIt) {
            //const Element& elem = *elemIt;
            
            //elemCtx.updatePrimaryStencil(elem);
            //elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
            
            //unsigned compressedDofIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            
            //const auto& iq = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
            //const auto& fs = iq.fluidState();
            //Scalar So = Opm::decay<Scalar>(fs.saturation(oilPhaseIdx));
            //maxOilSaturation_[compressedDofIdx] = std::max(maxOilSaturation_[compressedDofIdx], So);

            //Scalar initialValue = 0.0; //Should be changed...
            //SmaxVE_.push_back(initialValue);
        //}
        
    }
    

    template <class Context>
    void setSmax_VE(const Context& context, unsigned dofIdx, unsigned timeIdx, Scalar value)
    {
        unsigned compressedDofIdx = context.globalSpaceIndex(dofIdx, timeIdx);
        SmaxVE_[compressedDofIdx] = value;
    }

    template <class Context>
    Scalar getSmax_VE(const Context& context, unsigned dofIdx, unsigned timeIdx) const
    {
        unsigned compressedDofIdx = context.globalSpaceIndex(dofIdx, timeIdx);
        return SmaxVE_[compressedDofIdx];
    }

    
    /*!
     * \copydoc FvBaseMultiPhaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemTemperatureLow,
                             "The lower temperature [K] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemTemperatureHigh,
                             "The upper temperature [K] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, FluidSystemNumTemperature,
                             "The number of intervals between the lower and "
                             "upper temperature");

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemPressureLow,
                             "The lower pressure [Pa] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemPressureHigh,
                             "The upper pressure [Pa] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, FluidSystemNumPressure,
                             "The number of intervals between the lower and "
                             "upper pressure");

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, Temperature,
                             "The temperature [K] in the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, MaxDepth,
                             "The maximum depth [m] of the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, SimulationName,
                             "The name of the simulation used for the output "
                             "files");
    }

    /*!
     * \name Problem parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    {
        std::ostringstream oss;
        oss << EWOMS_GET_PARAM(TypeTag, std::string, SimulationName)
            << "_" << Model::name();
        if (GET_PROP_VALUE(TypeTag, EnableEnergy))
            oss << "_ni";
        oss << "_" << Model::discretizationName();
        return oss.str();
    }

    /*!
     * \copydoc FvBaseProblem::endTimeStep
     */
    void endTimeStep()
    {
#ifndef NDEBUG
        Scalar tol = this->model().newtonMethod().tolerance()*1e5;
        this->model().checkConservativeness(tol);

        // Calculate storage terms
        PrimaryVariables storageL, storageG;
        this->model().globalPhaseStorage(storageL, /*phaseIdx=*/0);
        this->model().globalPhaseStorage(storageG, /*phaseIdx=*/1);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout << "Storage: liquid=[" << storageL << "]"
                      << " gas=[" << storageG << "]\n" << std::flush;
        }
#endif // NDEBUG
        //***************************HAF******START***********************************
        //@HAF: Should this code snippet be here???
        const auto& simulator = this->simulator();
        ElementContext elemCtx(simulator);
        const auto& vanguard = simulator.vanguard();
        auto elemIt = vanguard.gridView().template begin</*codim=*/0>();
        const auto& elemEndIt = vanguard.gridView().template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            const Element& elem = *elemIt;

            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            unsigned compressedDofIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& iq = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& fs = iq.fluidState();

            Scalar Sn = Opm::decay<Scalar>(fs.saturation(gasPhaseIdx)); //@HAF: Is this the index we want?

            //auto& iqTEST = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
            //auto& fsTEST = iqTEST.fluidState();
            //if (Sn < 0.0)
            //{
            //  fsTEST.setSaturation(gasPhaseIdx, 0.0);
            //  fsTEST.setSaturation(liquidPhaseIdx, 1.0);
            //}
            
            SmaxVE_[compressedDofIdx] = std::max(SmaxVE_[compressedDofIdx], Sn); //NOTE: For testing without VE influence!!!
            std::cout << " compressedDofIdx= " << compressedDofIdx << " SmaxVE= " << SmaxVE_[compressedDofIdx] << " Sn= " << Sn << std::endl;
        }
        //***************************HAF******END*************************************
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        const auto& pos = context.pos(spaceIdx, timeIdx);
        if (inHighTemperatureRegion_(pos))
            return temperature_ + 100;
        return temperature_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context, unsigned spaceIdx,
                                           unsigned timeIdx) const
    {
        //const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        //if (isFineMaterial_(pos))
        //    return fineK_;
        //return coarseK_; //NOTE: For testing without VE influence!!!

        //////unsigned compressedDofIdx = context.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
        // /* //NOTE: For testing without VE influence!!!
        unsigned compressedDofIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return integratedPermA_[compressedDofIdx];
        // */
    }


    /*!
     * \copydoc ???????
     */
    //*****HAF---START*********************************
    void computeIntegratedPermeabilities()
    {
        //******************************************************************
        //NOTE: We assume (for now) isotropic permeability!!!
        //******************************************************************
        //@HAF

        const auto& simulator = this->simulator();
        ElementContext elemCtx(simulator);

        const auto& vanguard = simulator.vanguard();
        const auto& gridView = vanguard.gridView();
        int numElements = gridView.size(/*codim=*/0);
        integratedPermA_.resize(numElements);
                
        auto elemIt = vanguard.gridView().template begin</*codim=*/0>();
        const auto& elemEndIt = vanguard.gridView().template end</*codim=*/0>();
        int idx = 0;

        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& elem = *elemIt;

            elemCtx.updatePrimaryStencil(elem);
            computeSpecificIntegratedPermeability(elemCtx, /*spaceIdx*/ 0, /*timeIdx*/ 0, idx);
            idx++;
        }
    }


    template <class Context>
        void computeSpecificIntegratedPermeability(const Context& context, unsigned spaceIdx, unsigned timeIdx, const int& idx)
    {
        //******************************************************************
        //NOTE: We assume (for now) isotropic permeability!!!
        //******************************************************************
        //@HAF

        // /* //NOTE: For testing without VE influence!!!   
        Dune::VerteqColumnUtility< Grid > verteqUtil ( context.problem().simulator().vanguard().grid() );
        const auto& stencil = context.stencil(timeIdx);
        const auto& entity = stencil.entity(spaceIdx);
        //std::cout << "Start column for entity " << entity.impl().index() << std::endl;
        const auto endCol = verteqUtil.end( entity );
        Scalar perm = 0.0;
        Scalar fineIsotropicPerm = fineK_[0][0];
        Scalar coarseIsotropicPerm = coarseK_[0][0];
        
        for( auto col = verteqUtil.begin( entity ); col != endCol; ++col )
        {
            const auto& colCell = *col;

            if (isFineMaterial_(colCell.h()))
            {
                perm += fineIsotropicPerm;
            }
            else
            {
                perm += coarseIsotropicPerm;
            }
            perm *= colCell.dz();
        }
        
        //std::cout << "permeability= " << perm << std::endl; exit(0);

        integratedPermA_[idx] = this->toDimMatrix_(perm);
        // */
    }
    //*****HAF---END***********************************

    
    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
       // /* //NOTE: For testing without VE influence!!!   
      Dune::VerteqColumnUtility< Grid > verteqUtil ( context.problem().simulator().vanguard().grid() );
      const auto& stencil = context.stencil(timeIdx);
      const auto& entity = stencil.entity(spaceIdx);
      //std::cout << "Start column for entity " << entity.impl().index() << std::endl;
      const auto endCol = verteqUtil.end( entity );
      Scalar porosity = 0.0;
      for( auto col = verteqUtil.begin( entity ); col != endCol; ++col )
      {
        const auto& colCell = *col;
        
        //std::cout << "Column cell [ " << colCell.index()
        //        << " ]: h = " << colCell.h()
        //    << " fine cell idx " << colCell.fineCellIndex() 
        //        << " dz = "   << colCell.dz() << std::endl;
        
      
      
        //if (isFineMaterial_(colCell.h()))
        //   porosity += finePorosity_;

        if (isFineMaterial_(colCell.h()))
        {
            porosity += colCell.dz()*finePorosity_;
        }
        else
        {
            porosity += colCell.dz()*coarsePorosity_;
        }
      }
       // */
      //std::cout << "porosity= " << finePorosity_ << std::endl; //exit(0);
      //std::cout << "porosity= " << coarsePorosity_ << std::endl; //exit(0);
      //std::cout << "porosity= " << porosity << std::endl; exit(0);

      return porosity;
      // return 0.15; //NOTE: For testing without VE influence!!!
    }

    

    

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& context,
                                               unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineMaterialParams_;
        return coarseMaterialParams_;
    }

    /*!
     * \brief Return the parameters for the heat storage law of the rock
     *
     * In this case, we assume the rock-matrix to be granite.
     */
    template <class Context>
    const SolidEnergyLawParams&
    solidEnergyLawParams(const Context& context OPM_UNUSED,
                         unsigned spaceIdx OPM_UNUSED,
                         unsigned timeIdx OPM_UNUSED) const
    { return solidEnergyLawParams_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::thermalConductionParams
     */
    template <class Context>
    const ThermalConductionLawParams &
    thermalConductionLawParams(const Context& context,
                            unsigned spaceIdx,
                            unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineThermalCondParams_;
        return coarseThermalCondParams_;
    }

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::boundary
     */
    template <class Context>
    void boundary(BoundaryRateVector& values, const Context& context,
                  unsigned spaceIdx, unsigned timeIdx) const
    {
#warning TODO: make this meaningful: We let (at least for the time being) setFreeFlow at all the boundaries.
        const auto& pos = context.pos(spaceIdx, timeIdx);

        // Temporary BC's from HAF ****START*************************************

        if (onLeftBoundary_(pos))
        {
            Opm::CompositionalFluidState<Scalar, FluidSystem> fs;

            initialFluidState_(fs, context, spaceIdx, timeIdx);
            fs.checkDefined();
            
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
            //values.setInFlow(context, spaceIdx, timeIdx, fs);
            //values.setOutFlow(context, spaceIdx, timeIdx, fs);
            /*
            Scalar HVE = fs.getH_VE();
            Scalar v0 = values[0];
            values[0] = v0*HVE;
            Scalar v1 = values[1];
            values[1] = v1*HVE;
            */
        }
        else if (onRightBoundary_(pos))
        //if ((onLeftBoundary_(pos))  || (onRightBoundary_(pos)))
        {
            //values.setFreeFlow(context, spaceIdx, timeIdx, fs);
            //values.setInFlow(context, spaceIdx, timeIdx, fs);

            Scalar pTop = 1e5; //@HAF: New.
            const auto& insideFluidState = context.intensiveQuantities(spaceIdx, timeIdx).fluidState();
            
            //Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
            
            Opm::ImmiscibleFluidState<Scalar, FluidSystem> fs;
            fs.assign(insideFluidState);
            // /*
            Scalar gasPressure = pTop - fs.pressure(FluidSystem::liquidPhaseIdx) + fs.pressure(FluidSystem::gasPhaseIdx);
            fs.setPressure(FluidSystem::liquidPhaseIdx, pTop);
            fs.setPressure(FluidSystem::gasPhaseIdx, gasPressure);
            //std::cout << " SMax= " << fs.getSmax() << " satur= " << fs.saturation(FluidSystem::gasPhaseIdx) << std::endl;
            std::cout << " gasPressure= " << gasPressure << " fgasPressure= " << fs.pressure(FluidSystem::gasPhaseIdx) << " fsLiquidPressure= " << fs.pressure(FluidSystem::liquidPhaseIdx) << std::endl;
            // */
            
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
             /*
            Scalar HVE = fs.getH_VE();
            Scalar v0 = values[0];
            values[0] = v0*HVE;
            Scalar v1 = values[1];
            values[1] = v1*HVE;
             */
            std::cout << " values[0]= " << values[0] << " values[1]= " << values[1] << std::endl;
        }
        else
        {
            // impose a no flow boundary conditon ---
            values.setNoFlow(); 
        }
        // Temporary BC's from HAF ****END***************************************

        /*
        if (onLeftBoundary_(pos)) {
            Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
            initialFluidState_(fs, context, spaceIdx, timeIdx);
            fs.checkDefined();

            // impose an freeflow boundary condition
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else if (onInlet_(pos)) {
            RateVector massRate(0.0);
            massRate[contiCO2EqIdx] = -1e-3; // [kg/(m^3 s)] ... arealet istedet??? eller per meter i 2D??????

            typedef Opm::ImmiscibleFluidState<Scalar, FluidSystem> FluidState;
            FluidState fs;
            fs.setSaturation(gasPhaseIdx, 1.0);
            const auto& pg =
                context.intensiveQuantities(spaceIdx, timeIdx).fluidState().pressure(gasPhaseIdx);
            fs.setPressure(gasPhaseIdx, Toolbox::value(pg));
            fs.setTemperature(temperature(context, spaceIdx, timeIdx));

            typename FluidSystem::template ParameterCache<Scalar> paramCache;
            paramCache.updatePhase(fs, gasPhaseIdx);
            Scalar h = FluidSystem::template enthalpy<FluidState, Scalar>(fs, paramCache, gasPhaseIdx);

            // impose an forced inflow boundary condition for pure CO2
            values.setMassRate(massRate);
            values.setEnthalpyRate(massRate[contiCO2EqIdx] * h);
        }
        else
            // no flow on top and bottom
            values.setNoFlow();
        */
        
    }

    // \}

    /*!
     * \name Volumetric terms
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables& values, const Context& context, unsigned spaceIdx,
                 unsigned timeIdx) const
    {
        Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
        initialFluidState_(fs, context, spaceIdx, timeIdx);

        // const auto& matParams = this->materialLawParams(context, spaceIdx,
        // timeIdx);
        // values.assignMassConservative(fs, matParams, /*inEquilibrium=*/true);
        values.assignNaive(fs);
    }

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all components is 0
     * everywhere.
     */
    template <class Context>
    void source(RateVector& rate,
                const Context& context OPM_UNUSED,
                unsigned spaceIdx OPM_UNUSED,
                unsigned timeIdx OPM_UNUSED) const
    {
        rate = Scalar(0.0);
#warning do something more sensible here
        //Do we have one or more point sources here?
        //In the case of a point source, it must be placed in the 2D grid anyway,
        //since vertical integration does not matter in this test case when dz=1.......
        //But if the source
        //spreads out over several column cells, we must perform an integration
        //using Dune::VerteqColumnUtility< Grid >
        //In any case it seems like the position(s) of the source(s) must be hardcoded in this routine.
        //I guess such kind of information can be found from the grid...
        //@haf: how to find the geometry from the grid???
        //@haf: NOTE: For the time being we just skip the integration (which generally MUST
        //be done even if we have a point source since we need dz), and just assume that
        //dz = 1 !!!
        //@haf: Moreover, we "cheat" for the time being and (for test purposes) just put in a
        //point source in a valid point in the (2D) grid. We just choose the index 100.
        //Later on this should be coded properly!!!
        unsigned globalDofIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        // /*
        int HAFTestIndex = 192 + 11;
        int HAFTestIndexp = 192 + 11 + 1;
        int HAFTestIndexm = 192 + 11 - 1;
        // */

        //int HAFTestIndex = 4; //For 3x3 grid
        //int HAFTestIndex = 5; //For 4x4 grid (litt spesiellt...)
        //int HAFTestIndex = 12; //For 5x5 grid

        //int BMC1_TestIndexL = 2*800 + 400 - 2;
        //int BMC1_TestIndexU = 2*800 + 400 + 1;

        //int BMC2_TestIndexL = 2*800 + 200 - 2;
        //int BMC2_TestIndexU = 2*800 + 200 + 1;
        int BMC2_TestIndexL = 200 - 2;
        int BMC2_TestIndexU = 200 + 1;
        int BMC2_TestIndex = 200;

        int BMC3_TestIndexL = 75 - 2;
        int BMC3_TestIndexU = 75 + 1;
        int BMC3_TestIndex = 75;

        const auto& simulator = this->simulator();
        Scalar actTime = simulator.time();

        //if (globalDofIdx == HAFTestIndex)
        //if ((globalDofIdx >= BMC1_TestIndexL) && (globalDofIdx <= BMC1_TestIndexU))
        if ((globalDofIdx >= BMC2_TestIndexL) && (globalDofIdx <= BMC2_TestIndexU))
        //if (globalDofIdx == BMC2_TestIndex)
        //if ((globalDofIdx >= BMC3_TestIndexL) && (globalDofIdx <= BMC3_TestIndexU))
        {
            //rate[Indices::conti0EqIdx + CO2Idx] = 1e-4; // just a demo, [kg / m^3 / s]
            
            // /*
            //if (actTime <= 315360000.0)
            if (actTime <= 630720000.0)
                //if (actTime <= 1e+10)
            {
                rate[Indices::conti0EqIdx + CO2Idx] = 1e-4; // just a demo, [kg / m^3 / s]
                //rate[Indices::conti0EqIdx + CO2Idx] = 1e-2; // just a demo, [kg / m^3 / s]
                //rate[Indices::conti0EqIdx + CO2Idx] = 5e-4; // just a demo, [kg / m^3 / s]
                //rate[Indices::conti0EqIdx + CO2Idx] = 9.3e-4; // just a demo, [kg / m^3 / s]
                //rate[Indices::conti0EqIdx + CO2Idx] = 2.537e-6; //NB: TH!!! // just a demo, [kg / m^3 / s]
            }
            else
            {
                //std::cout << " SLUTT!!! " << std::endl;
                rate[Indices::conti0EqIdx + CO2Idx] = 0.0; // [kg / m^3 / s]
            }
            // */
        }
        else
        {
            rate[Indices::conti0EqIdx + CO2Idx] = 0.0; // [kg / m^3 / s]
        }
    }

    //! \}

private:
    template <class Context, class FluidState>
    void initialFluidState_(FluidState& fs,
                            const Context& context,
                            unsigned spaceIdx,
                            unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        //std::cout << " pos[0]= " << pos[0] << " pos[1]= " << pos[1] << std::endl;
        
        //////
        // set temperature
        //////
        fs.setTemperature(temperature(context, spaceIdx, timeIdx));

        //////
        // set saturations
        //////
        fs.setSaturation(FluidSystem::liquidPhaseIdx, 1.0);
        fs.setSaturation(FluidSystem::gasPhaseIdx, 0.0);

        //////
        // set pressures
        /////

        ////Scalar densityL = FluidSystem::Brine::liquidDensity(temperature_, Scalar(1e0));//@HAF: OK for VE? (see also line above)
        //Scalar depth = maxDepth_ - pos[dim - 1]; //@HAF: Not relevant for VE???
        //Scalar pl = 1e5 - densityL * this->gravity()[dim - 1] * depth; //@HAF: Not relevant for VE???
        Scalar pl = 2.6e7; //@HAF: New.
        //Scalar pl = 1.0e8; //@HAF: Newest...
        Scalar pTop = 1e5; //@HAF: New.

        Scalar densityL = FluidSystem::Brine::liquidDensity(temperature_, pl);
        Scalar densityCO2 = FluidSystem::CO2::gasDensity(temperature_, pl);
        //Scalar densityCO2 = FluidSystem::Brine::gasDensity(temperature_, pl);

        Scalar pC[numPhases];
        const auto& matParams = this->materialLawParams(context, spaceIdx, timeIdx);
        //MaterialLaw::capillaryPressures(pC, matParams, fs);
        //MaterialLaw::RegularizedBrooksCoreyVE::capillaryPressures(pC, matParams, fs); //@HAF: New.

        //fs.setPressure(liquidPhaseIdx, pl + (pC[liquidPhaseIdx] - pC[liquidPhaseIdx]));
        //fs.setPressure(gasPhaseIdx, pl + (pC[gasPhaseIdx] - pC[liquidPhaseIdx]));
        
        //NOTE: HAF: Mangler "vanlig capillary pressure"!!!!!!!!!!!!!!!! VE
        //******************************************************
        //HAF: Special coding for "BENCHMARK: CASE 2": VE 
        Scalar dipIgrader = 1.0;
        Scalar dipAngle = (3.14159265359/180)*dipIgrader;
        Scalar gravity = 9.80665;
        //Scalar densityL = 1016.96;
        //Scalar densityCO2 = 523.384;
        //std::cout << " HAFdensityL= " << densityL << " HAFdensityCO2= " << densityCO2 << std::endl;
        //******************************************************
        //fs.setPressure(liquidPhaseIdx, pl);
        //fs.setPressure(gasPhaseIdx, pl);
        //fs.setPressure(liquidPhaseIdx, pl - densityL*gravity*sin(dipAngle)*pos[0]);
        ////fs.setPressure(gasPhaseIdx, pl - densityL*gravity*sin(dipAngle)*pos[0]);
        //fs.setPressure(gasPhaseIdx, pl - densityCO2*gravity*sin(dipAngle)*pos[0]);

        Scalar xMax = this->boundingBoxMax()[0];
        fs.setPressure(liquidPhaseIdx, pTop + densityL*gravity*sin(dipAngle)*(xMax-pos[0]));
        fs.setPressure(gasPhaseIdx, pTop + densityL*gravity*sin(dipAngle)*(xMax-pos[0]));
        //fs.setPressure(gasPhaseIdx, pTop + densityCO2*gravity*sin(dipAngle)*(xMax-pos[0]));
        
        //////
        // set composition of the liquid phase
        //////
        //fs.setMoleFraction(liquidPhaseIdx, CO2Idx, 0.005);
        //fs.setMoleFraction(liquidPhaseIdx, BrineIdx, 1.0 - fs.moleFraction(liquidPhaseIdx, CO2Idx));

        // /*
        fs.setMoleFraction(liquidPhaseIdx, CO2Idx, 0.0); //@HAF: OK???
        fs.setMoleFraction(liquidPhaseIdx, BrineIdx, 1.0 - fs.moleFraction(liquidPhaseIdx, CO2Idx));
        fs.setMoleFraction(gasPhaseIdx, CO2Idx, 1.0); //@HAF: OK???
        fs.setMoleFraction(gasPhaseIdx, BrineIdx, 0.0); //@HAF: OK???
        // */
        
        //Fjerne dette???
        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        typedef Opm::ComputeFromReferencePhase<Scalar, FluidSystem> CFRP;
        CFRP::solve(fs, paramCache,
                    /*refPhaseIdx=*/liquidPhaseIdx,
                    /*setViscosity=*/true,
                    /*setEnthalpy=*/true);
    }

    bool onLeftBoundary_(const GlobalPosition& pos) const
    { return pos[0] < eps_; }

    bool onRightBoundary_(const GlobalPosition& pos) const
    { return pos[0] > this->boundingBoxMax()[0] - eps_; }

    bool onInlet_(const GlobalPosition& pos) const
    { return onRightBoundary_(pos) && (5 < pos[1]) && (pos[1] < 15); }

    bool inHighTemperatureRegion_(const GlobalPosition& pos) const
    { //return (pos[0] > 20) && (pos[0] < 30) && (pos[1] > 5) && (pos[1] < 35);
        return true;
    }

    void computeThermalCondParams_(ThermalConductionLawParams& params, Scalar poro)
    {
        Scalar lambdaWater = 0.6;
        Scalar lambdaGranite = 2.8;

        Scalar lambdaWet = std::pow(lambdaGranite, (1 - poro))
                           * std::pow(lambdaWater, poro);
        Scalar lambdaDry = std::pow(lambdaGranite, (1 - poro));

        params.setFullySaturatedLambda(gasPhaseIdx, lambdaDry);
        params.setFullySaturatedLambda(liquidPhaseIdx, lambdaWet);
        params.setVacuumLambda(lambdaDry);
    }

    bool isFineMaterial_(const GlobalPosition& pos) const
    { //return pos[dim - 1] > fineLayerBottom_;
        return true;
    }

    bool isFineMaterial_(const Scalar h) const
    { //return h > fineLayerBottom_;
        return true;
    }

    DimMatrix fineK_;
    DimMatrix coarseK_;
    Scalar fineLayerBottom_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;

    ThermalConductionLawParams fineThermalCondParams_;
    ThermalConductionLawParams coarseThermalCondParams_;
    SolidEnergyLawParams solidEnergyLawParams_;

    Scalar temperature_;
    Scalar maxDepth_;
    Scalar eps_;

    unsigned nTemperature_;
    unsigned nPressure_;

    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;

    //**********************HAF**************************************
    std::vector<Scalar> SmaxVE_;
    std::vector<Scalar> H_VE_;
    std::vector<DimMatrix> integratedPermA_; //HAF
};
} // namespace Ewoms

#endif
