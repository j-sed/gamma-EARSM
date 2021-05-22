/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::RASModels::gammaEARSM

Group
    grpRASTurbulence

Description
   Implementation of the gammaEARSM turbulence model.

    This turbulence model is described in:
    \verbatim
        Link to thesis, once it is released
    \endverbatim

    Turbulence model based on existing models:
    \verbatim
        Hellsten, Antti K. (2005). 
        "New Advanced K-W Turbulence Model for High-Lift Aerodynamics." 
        AIAA Journal 43 (9): 1857–69. doi:10.2514/1.13754.
    \endverbatim
    and
    \verbatim
        MENTER, Florian R., Pavel E. SMIRNOV, Tao LIU a Ravikanth AVANCHA. 
        A one-equation local correlation-based transition model. 
        Flow, Turbulence and Combustion. 2015, 95(4), 583–619. ISSN 1386-6184. 
        doi:10.1007/s10494-015-9622-4
    \endverbatim
	Implementation of gammaSST and EARSM is taken over from Jiri Furst

    The default model coefficients are
    \verbatim
        gammaEARSMCoeffs
        {
        betaStar       0.09
	    gamma1         0.518
        gamma2         0.44
	    beta1          0.0747
	    beta2          0.0828
	    alphaK1        1.1
	    alphaK2        1.1
	    alphaOmega1    0.53
	    alphaOmega2    1.0
	    alphaD1        1.0
	    alphaD2        0.4
        kInf           1.e-10
        A0             -0.72
        Flength        28;
        ca2            0.06;
        ce2            50;
        sigmaGamma     1;
        CPG1           100.0;
        CPG1lim        1.5;
        CPG2           -7.34;
        CPG3           10.0;
        CPG2lim        3.0;
        CTU1           110.0;
        CTU2           2000.0;
        CTU3           0.85;
        curvatureCorrection false
        }
    \endverbatim

Author
    Jan Sedlacek <jan.sedlacek95@gmail.com> and Jiri Furst <Jiri.Furst@gmail.com>

SourceFiles
    gammaEARSM.C

\*---------------------------------------------------------------------------*/

#ifndef gammaEARSM_H
#define gammaEARSM_H

#include "RASModel.H"
#include "nonlinearEddyViscosity.H"
#include "fvOptions.H"
#include "kOmegaSST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class gammaEARSM Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class gammaEARSM
:
    public nonlinearEddyViscosity<RASModel<BasicTurbulenceModel> >
    //,public kOmegaSST<BasicTurbulenceModel>
{
    // Private memmber functions
    
    // Disallow default bitwise copy construct and assignment
    gammaEARSM(const gammaEARSM&);
    gammaEARSM& operator=(const gammaEARSM&);
    
    
protected:
    
    // Protected data
    
    // Model coefficients

    dimensionedScalar betaStar_;

    dimensionedScalar gamma1_;
    dimensionedScalar gamma2_;

    dimensionedScalar beta1_;
    dimensionedScalar beta2_;

    dimensionedScalar alphaK1_;
    dimensionedScalar alphaK2_;
    
    dimensionedScalar alphaOmega1_;
    dimensionedScalar alphaOmega2_;
    
    dimensionedScalar alphaD1_;
    dimensionedScalar alphaD2_;

    dimensionedScalar kInf_;

    Switch curvatureCorrection_;
    
    dimensionedScalar A0_;
    
    dimensionedScalar Flength_;
    
    dimensionedScalar ca2_;
    dimensionedScalar ce2_;
    dimensionedScalar sigmaGamma_;
    dimensionedScalar CPG1_;
    dimensionedScalar CPG1lim_;
    dimensionedScalar CPG2_;
    dimensionedScalar CPG3_;
    dimensionedScalar CPG2lim_;
    dimensionedScalar CTU1_;
    dimensionedScalar CTU2_;
    dimensionedScalar CTU3_;
    dimensionedScalar ReThetacLim_;
    dimensionedScalar Ck_;
    dimensionedScalar CSEP_;
    // Fields
    
    volScalarField k_;
    volScalarField omega_;
    
    volScalarField gammaInt_;
    //- Wall distance
    //  Note: different to wall distance in parent RASModel
    //  which is for near-wall cells only
    const volScalarField& y_;


    // Protected Member Functions

    virtual tmp<volScalarField> fMix(const volScalarField& gradKgradOmegaByOmega) const;
    
    tmp<volScalarField> blend
    (
        const volScalarField& fMix,
        const dimensionedScalar& psi1,
        const dimensionedScalar& psi2
    ) const
    {
        return fMix*(psi1 - psi2) + psi2;
    }

    tmp<volScalarField> alphaK(const volScalarField& fMix) const
    {
        return blend(fMix, alphaK1_, alphaK2_);
    }

    tmp<volScalarField> alphaOmega(const volScalarField& fMix) const
    {
        return blend(fMix, alphaOmega1_, alphaOmega2_);
    }

    tmp<volScalarField> alphaD(const volScalarField& fMix) const
    {
        return blend(fMix, alphaD1_, alphaD2_);
    }

    tmp<volScalarField> beta(const volScalarField& fMix) const
    {
        return blend(fMix, beta1_, beta2_);
    }

    tmp<volScalarField> gamma(const volScalarField& fMix) const
    {
        return blend(fMix, gamma1_, gamma2_);
    }

    virtual void correctNut();

    virtual volScalarField N(
        const volScalarField& A3p,
        const volScalarField& P1,
        const volScalarField& P2
    ) const;
    
    // Model functions of gammaSST model
    tmp<volScalarField> ReThetac() const;
    tmp<volScalarField> Fonset(const volScalarField& S) const;
    tmp<volScalarField> Fonset1(const volScalarField& S) const;
    tmp<volScalarField> Fturb() const;
    tmp<volScalarField> TuL() const;
    tmp<volScalarField> FPG() const;
    // tmp<volScalarField> LambdaThetaTest() const; // this is just for testing purposes
    tmp<volScalarField> F3() const; 
    
    virtual void correctNonlinearStress(const volTensorField& gradU);

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("gammaEARSM");


    // Constructors

        //- Construct from components
       gammaEARSM
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~gammaEARSM()
    {}


    // Member Functions
    
    //- Re-read model coefficients if they have changed
    virtual bool read();
    
    //- Return the effective diffusivity for k
    tmp<volScalarField> DkEff(const volScalarField& fMix) const
        {
            return tmp<volScalarField>
                (
                    new volScalarField("DkEff", alphaK(fMix)*this->nut_ + this->nu())
                );
        }
    
    //- Return the effective diffusivity for omega
    tmp<volScalarField> DomegaEff(const volScalarField& fMix) const
        {
            return tmp<volScalarField>
                (
                    new volScalarField("DomegaEff", alphaOmega(fMix)*this->nut_ + this->nu())
                );
        }
    
    //- Return the total fluctuation kinetic energy
    virtual tmp<volScalarField> k() const
        {
            return k_;
        }
    
    
    //- Return the turbulence specific dissipation rate
    virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }
    
    //- Return the total fluctuation kinetic energy dissipation rate
    virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
                (
                    new volScalarField
                    (
                    IOobject
                    (
                        "epsilon",
                        this->mesh_.time().timeName(),
                        this->mesh_
                    ),
                    betaStar_*k_*omega_,
                    omega_.boundaryField().types()
                    )
                );
            
        }
    
    
    //- Return the effective diffusivity for gamma
    tmp<volScalarField> DgammaEff() const
        {
            return tmp<volScalarField>
                (
                    new volScalarField("DgammaEff", (this->nut_/sigmaGamma_) + this->nu())
                );
        }
    
    //- Return the intermittency
    virtual tmp<volScalarField> gammaInt() const
        {
            return gammaInt_;
        }
    
    virtual tmp<volScalarField> Rt() const
        {
            return tmp<volScalarField>
                (
                    new volScalarField("Rt", this->k_/(this->nu()*this->omega_)) 
                    // LANGTRY's definition
                );
        }
    
    //- Validate the turbulence fields after construction
    //  Update turbulence viscosity and other derived fields as requires
    virtual void validate();
    
    //- Solve the turbulence equations and correct the turbulence viscosity
    virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#   include "gammaEARSM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
