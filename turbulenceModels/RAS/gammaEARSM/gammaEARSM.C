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


\*---------------------------------------------------------------------------*/

#include "gammaEARSM.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"
#include "kOmegaSST.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> gammaEARSM<BasicTurbulenceModel>::fMix
(
    const volScalarField& gradKgradOmegaByOmega
) const
{
    tmp<volScalarField> Gamma = min
        (
            max
            (
                sqrt(k_) / (betaStar_ * omega_ * y_),
                500.0 * this->nu() / (omega_ * sqr(y_))
            ),
            20.0 * k_ /
            max( sqr(y_) * gradKgradOmegaByOmega, 200.0 * kInf_)
        );
    // Modification of the f_mix function based on gamma-SST
    return max
            (
                tanh(1.5 * pow4(Gamma)),
                F3()
            );
}


template<class BasicTurbulenceModel>
void gammaEARSM<BasicTurbulenceModel>::correctNut()
{
  correctNonlinearStress(fvc::grad(this->U_));
}


template<class BasicTurbulenceModel>
volScalarField gammaEARSM<BasicTurbulenceModel>::N
(
    const volScalarField& A3p,
    const volScalarField& P1,
    const volScalarField& P2
) const
{
    volScalarField N = A3p / 3.0;

    forAll(N, i)
    {
        if (P2[i] < 0)
        {
            N[i] += 2*pow(sqr(P1[i]) - P2[i], 1./6.)
                * cos( 1./3.*acos( P1[i]/sqrt(sqr(P1[i]) - P2[i]))) ;
        }
        else
        {
            scalar a = max(P1[i] + sqrt(P2[i]), 0.0);
            scalar b = P1[i] - sqrt(P2[i]);
            N[i] += pow(a, 1./3.) + sign(b) * pow(fabs(b), 1./3.);
        }
    };

    forAll(N.boundaryField(), patchi)
    {
        fvPatchScalarField& pN = N.boundaryFieldRef()[patchi];
        const fvPatchScalarField& pP1 = P1.boundaryField()[patchi];
        const fvPatchScalarField& pP2 = P2.boundaryField()[patchi];

        forAll(pN, i)
        {
            if (pP2[i] < 0)
            {
                pN[i] += 2*pow(sqr(pP1[i]) - pP2[i], 1./6.)
                    * cos( 1./3.*acos( pP1[i]/sqrt(sqr(pP1[i]) - pP2[i]))) ;
            }
            else
            {
                scalar a = max(pP1[i] + sqrt(pP2[i]), 0.0);
                scalar b = pP1[i] - sqrt(pP2[i]);
                pN[i] += pow(a, 1./3.) + sign(b) * pow(fabs(b), 1./3.);
            }
        };
        
    };

    return N;
}


template<class BasicTurbulenceModel>
void gammaEARSM<BasicTurbulenceModel>::correctNonlinearStress(const volTensorField& gradU)
{
    scalar Ctau = 6.0;
    volScalarField tau(
        max
        (
            1.0 / (this->betaStar_ * this->omega_),
            Ctau * sqrt(this->nu() / (this->betaStar_ * max(this->k_, this->kMin_) * this->omega_))
        ));
    
    volSymmTensorField S("S", tau * dev(symm(gradU)));
    volTensorField     W("W", -tau * skew(gradU));
    // NOTE: Wij = 1/2(dui/dxj - duj/dxi) = - skew(grad(U))

    volScalarField IIS  = tr(S & S);

    if (this->curvatureCorrection_)
    {
        const volVectorField& U = this->U_;
        const surfaceScalarField& phi = this->phi_;
        const rhoField& rho = this->rho_;

        volSymmTensorField DSDt = tau*dev(symm(fvc::grad(fvc::ddt(U)))) 
		+ fvc::div(phi/fvc::interpolate(rho),S,"div(phiv,S)");
        volVectorField SDeps = *(skew(S & DSDt))*2;
        volScalarField IIIS = tr(S & S & S);
        
        volTensorField B = (pow(IIS,2)*I + 12*IIIS*S + 6*IIS*(S&S))/
            max(2*pow(IIS,3) - 12*pow(IIIS,2), 1.e-10);

        volVectorField BSDeps = B & SDeps;
        volTensorField Wr = *(BSDeps);
        
        W -= (tau/this->A0_)*Wr;
    }
    
    volScalarField IIW  = tr(W & W);
    volScalarField IV   = tr(S & W & W);

    scalar Neq = 81.0 / 20.0;
    scalar CDiff = 2.2;
    volScalarField beta1eq = - 6.0/5.0 * Neq / (sqr(Neq) - 2*IIW);
    volScalarField A3p = 9.0/5.0 + 9.0/4.0 * CDiff * max(1 + beta1eq*IIS, 0.0);
    volScalarField P1 = (sqr(A3p)/27 + (9.0/20.0)*IIS - (2.0/3.0)*IIW) * A3p;
    volScalarField P2 = sqr(P1) - pow3(sqr(A3p)/9 + 0.9*IIS + (2.0/3.0)*IIW);
    
    volScalarField N = this->N(A3p, P1, P2);

    volScalarField Q = 5.0/6.0*(sqr(N) - 2*IIW)*(2*sqr(N)-IIW);

    volScalarField beta1 = -N*(2.0*sqr(N) - 7.0*IIW) / Q;
    volScalarField beta3 = -12.0 * IV / (N * Q);
    volScalarField beta4 = -2.0 * (sqr(N)  - 2.0*IIW) / Q;
    volScalarField beta6 = -6.0 * N / Q;
    volScalarField beta9 =  6.0 / Q;

    
    volScalarField Cmu = - 0.5 * (beta1 + IIW * beta6);

    this->nut_ = Cmu * this->k_ * tau;
    this->nut_.correctBoundaryConditions();

    
    this->nonlinearStress_ = this->k_ * symm(
        beta3 * ( (W & W) - (1.0/3.0) * IIW * I )
        + beta4 * ( (S & W) - (W & S) )
        + beta6 * ( (S & W & W) + (W & W & S) - IIW * S - (2.0/3.0) * IV * I)
        + beta9 * ( (W & S & W & W) - (W & W & S & W) )
    );

    this->nonlinearStress_.correctBoundaryConditions();

    BasicTurbulenceModel::correctNut();

}

template<class BasicTurbulenceModel>
tmp<volScalarField> gammaEARSM<BasicTurbulenceModel>::ReThetac() const
{
  return CTU1_ + CTU2_*exp(-CTU3_*TuL()*FPG() );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> gammaEARSM<BasicTurbulenceModel>::Fonset(const volScalarField& S) const
{
    return tmp<volScalarField>
	(
	    new volScalarField
     	    (
                IOobject
       	        (
	            "Fonset",
    	            this->runTime_.timeName(),
	            this->mesh_,
	            IOobject::NO_READ,
	            IOobject::NO_WRITE
	        ),
		max
		(
                    min(Fonset1(S), 2.0) - max(1.0-pow3(Rt()/3.5),0.0),
                    0.0
		)
   	    )
	);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> gammaEARSM<BasicTurbulenceModel>::Fonset1(const volScalarField& S) const
{
    return sqr(this->y_)*S/this->nu() / (2.2*ReThetac());
}

template<class BasicTurbulenceModel>
tmp<volScalarField> gammaEARSM<BasicTurbulenceModel>::Fturb() const
{
    return exp(-pow4(Rt()/2));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> gammaEARSM<BasicTurbulenceModel>::TuL() const
{
    return min(100 * sqrt(2.0/3.0*this->k_) / (this->omega_ * this->y_), 100.0);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> gammaEARSM<BasicTurbulenceModel>::LambdaThetaTest() const
{
    volVectorField n = fvc::grad(this->y_);
    volScalarField lambdaThetaLTest = 
        min( 1.0, max( -1.0, 
        -7.57e-3 * ( fvc::grad(this->U_ & n) & n) * sqr(this->y_) / this->nu() + 0.0128));

    tmp<volScalarField> tLambdaThetaTest(new volScalarField("LambdaThetaTest", lambdaThetaLTest));
    volScalarField& LambdaThetaTest_ = tLambdaThetaTest.ref();
    
    forAll(LambdaThetaTest_,i){
        LambdaThetaTest_[i] = lambdaThetaLTest[i];
    }
    return tLambdaThetaTest;
}
template<class BasicTurbulenceModel>
tmp<volScalarField> gammaEARSM<BasicTurbulenceModel>::FPG() const
{
    volVectorField n = fvc::grad(this->y_);
    volScalarField lambdaThetaL = 
        min( 1.0, max( -1.0, 
        -7.57e-3 * ( fvc::grad(this->U_ & n) & n) * sqr(this->y_) / this->nu() + 0.0128));

    tmp<volScalarField> tFPG(new volScalarField("FPG", lambdaThetaL));
 
    volScalarField& FPG_ = tFPG.ref();
    forAll(FPG_, i) {
        if (lambdaThetaL[i]>=0) 
            FPG_[i] = min(1 + CPG1_.value()*lambdaThetaL[i], CPG1lim_.value());
        else
            FPG_[i] = min(1 + CPG2_.value()*lambdaThetaL[i] + 
            CPG3_.value()*min(lambdaThetaL[i]+0.0681,0), 
            CPG2lim_.value());
        FPG_[i] = max(FPG_[i], 0.0);
    }

    return tFPG;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> gammaEARSM<BasicTurbulenceModel>::F3() const
{
    return exp(-sqr(pow4(this->y_*sqrt(this->k_)/(scalar(120)*this->nu()))));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
gammaEARSM<BasicTurbulenceModel>::gammaEARSM
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
    :
    nonlinearEddyViscosity<RASModel<BasicTurbulenceModel> >
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    // gammaEARSM coefficients
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),

    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            0.518
        )
    ),

    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),

    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.0747
        )
    ),

    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            1.1
        )
    ),

    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.1
        )
    ),

    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.53
        )
    ),

    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            1.0
        )
    ),

    alphaD1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaD1",
            this->coeffDict_,
            1.0
        )
    ),

    alphaD2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaD2",
            this->coeffDict_,
            0.4
        )
    ),

// gammaSST coefficients
    kInf_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kInf_",
            this->coeffDict_,
            sqr(dimVelocity),
            1.e-10
        )
    ),

    curvatureCorrection_
    (
        Switch::lookupOrAddToDict
        (
            "curvatureCorrection",
            this->coeffDict_,
            false
        )
    ),

    A0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A0",
            this->coeffDict_,
            -0.72
        )
    ),
            Flength_
            (
               dimensioned<scalar>::lookupOrAddToDict
              (
            "Flength",
            this->coeffDict_,
            28
        )
    ),
    ca2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ca2",
            this->coeffDict_,
            0.06
        )
    ),
    ce2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ce2",
            this->coeffDict_,
            50.0
        )
    ),
    sigmaGamma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaGamma",
            this->coeffDict_,
            1.0
        )
    ),
    CPG1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CPG1",
            this->coeffDict_,
            100.0
        )
    ),
    CPG1lim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CPG1lim",
            this->coeffDict_,
            1.5
        )
    ),
    CPG2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CPG2",
            this->coeffDict_,
            -7.34
        )
    ),
    CPG3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CPG3",
            this->coeffDict_,
            0.0
        )
    ),
    CPG2lim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CPG2lim",
            this->coeffDict_,
            10
        )
    ),
    CTU1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CTU1",
            this->coeffDict_,
            110.0
        )
    ),
    CTU2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CTU2",
            this->coeffDict_,
            2000.0
        )
    ),
    CTU3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CTU3",
            this->coeffDict_,
            0.85
        )
    ),
    ReThetacLim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ReThetacLim",
            this->coeffDict_,
            1100.
        )
    ),
    Ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ck",
            this->coeffDict_,
            1.
        )
    ),
    CSEP_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CSEP",
            this->coeffDict_,
            1.
        )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    gammaInt_ // gammaSST dictionary object
    (
        IOobject
        (
            "gamma",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
	this->mesh_
    ),
    y_(wallDist::New(this->mesh_).y())
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->correctNut();
        this->printCoeffs(type);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool gammaEARSM<BasicTurbulenceModel>::read()
{
    if (nonlinearEddyViscosity<RASModel<BasicTurbulenceModel> >::read())
    {    
        betaStar_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        alphaD1_.readIfPresent(this->coeffDict());
        alphaD2_.readIfPresent(this->coeffDict());
        kInf_.readIfPresent(this->coeffDict());

        Flength_.readIfPresent(this->coeffDict());
        ca2_.readIfPresent(this->coeffDict());
        ce2_.readIfPresent(this->coeffDict());
        sigmaGamma_.readIfPresent(this->coeffDict());
        CPG1_.readIfPresent(this->coeffDict());
        CPG1lim_.readIfPresent(this->coeffDict());
        CPG2_.readIfPresent(this->coeffDict());
        CPG3_.readIfPresent(this->coeffDict());
        CPG2lim_.readIfPresent(this->coeffDict());
        CTU1_.readIfPresent(this->coeffDict());
        CTU2_.readIfPresent(this->coeffDict());
        CTU3_.readIfPresent(this->coeffDict());
        ReThetacLim_.readIfPresent(this->coeffDict());
        Ck_.readIfPresent(this->coeffDict());
        CSEP_.readIfPresent(this->coeffDict());curvatureCorrection_.readIfPresent("curvatureCorrection", this->coeffDict());
        A0_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }    
    
}

template<class BasicTurbulenceModel>
void gammaEARSM<BasicTurbulenceModel>::validate()
{
    this->correctNut();
}


template<class BasicTurbulenceModel>
void gammaEARSM<BasicTurbulenceModel>::correct()
{

    if (!this->turbulence_)
    {
        return;
    }

    nonlinearEddyViscosity<RASModel<BasicTurbulenceModel> >::correct();

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));
    
    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField>   tgradU = fvc::grad(U);
    const volScalarField S2(2*magSqr(symm(tgradU())));
    const volScalarField S("S", sqrt(S2));
    const volScalarField W("Omega", sqrt(2*magSqr(skew(tgradU()))));

    volScalarField G
    (
        this->GName(),
        (nut * dev(twoSymm(tgradU())) - this->nonlinearStress_) && tgradU()
    );
    
    

    omega_.boundaryFieldRef().updateCoeffs();


    volScalarField gradKgradOmegaByOmega
    (
        (fvc::grad(k_) & fvc::grad(omega_)) / omega_
    );

    volScalarField fMix( this->fMix(gradKgradOmegaByOmega) );
    // Formulation that is non-temporary. This is useful for writing out this field
    // const volScalarField CDOmega
    // ( "CDOmega",
    //     this->alphaD(fMix) * alpha * rho * max( gradKgradOmegaByOmega, dimensionedScalar("zero", inv(sqr(dimTime)), 0.0))
    // );
    {
        volScalarField gamma( this->gamma(fMix) );
        volScalarField beta( this->beta(fMix) ); // theoretically could be tmp
        volScalarField alphaD( this->alphaD(fMix) );
        tmp<volScalarField> CDOmega = alphaD * alpha * rho *
            max( gradKgradOmegaByOmega, dimensionedScalar("zero", inv(sqr(dimTime)), 0.0));
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*DomegaEff(fMix), omega_)
         ==
            alpha*rho*gamma * omega_ / max(k_, this->kMin_) * G
            - fvm::SuSp((2.0/3.0)*alpha*rho*gamma*divU, omega_)
            - fvm::Sp(alpha*rho*beta*omega_, omega_)
            + CDOmega
            + fvOptions(alpha, rho, omega_)            
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(omega_);
        bound(omega_, this->omegaMin_);
    }
    
    
    // Turbulent kinetic energy equation
    const volScalarField FonLim(
        "FonLim",
        min( max(sqr(this->y_)*S/this->nu() / (
            2.2*ReThetacLim_) - 1., 0.), 3.)
    );
    const volScalarField PkLim(
        "PkLim",
        5*Ck_ * max(gammaInt()-0.2,0.) * (1-gammaInt()) * FonLim * 
        max(3*CSEP_*this->nu() - this->nut_, 0.*this->nut_) * S * W
    );
    
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(fMix), k_)
     ==
        alpha*rho*(G*gammaInt() +PkLim)
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU*gammaInt(), k_)
        - fvm::Sp(max(gammaInt(), scalar(0.1))*betaStar_*alpha*rho*omega_, k_)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNonlinearStress(tgradU());
    
    
    
    // Intermittency equation (2)

    volScalarField Pgamma1 = Flength_ * S * gammaInt_ * Fonset(S);
    volScalarField Pgamma2 = ca2_ * W * gammaInt_ * Fturb();
    
    tmp<fvScalarMatrix> gammaEqn
    (
        fvm::ddt(alpha, rho, gammaInt_)
        + fvm::div(alphaRhoPhi, gammaInt_)
        - fvm::laplacian(alpha*rho*this->DgammaEff(), gammaInt_)
        ==
        alpha*rho*Pgamma1 - fvm::Sp(alpha*rho*Pgamma1, gammaInt_) +
        alpha*rho*Pgamma2 - fvm::Sp(alpha*rho*ce2_*Pgamma2, gammaInt_)
    );
        
    gammaEqn.ref().relax();
    solve(gammaEqn);

    bound(gammaInt_,scalar(0));

    if (debug && this->runTime_.outputTime()) {
        //Examples of possible debug fields
        
        //S.write();
        //W.write();
        // gradKgradOmegaByOmega.write();
        // fMix.write();
        // CDOmega.write();
        //const volScalarField Pgamma("Pgamma", Pgamma1*(scalar(1)-gammaInt_));
        //Pgamma.write();
        //const volScalarField Egamma("Egamma", Pgamma2*(scalar(1)-ce2_*gammaInt_));
        //Egamma.write();
        // FonLim.write();
        //PkLim.write();
        // Fonset(S)().write();
        // F3()().write();
        // FPG()().write();
        // TuL()().write();
        // LambdaThetaTest()().write();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
