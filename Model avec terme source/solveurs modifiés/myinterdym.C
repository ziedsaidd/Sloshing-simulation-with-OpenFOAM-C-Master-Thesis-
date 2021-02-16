/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    interDyMFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "MULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "twoPhaseMixture.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "mathematicalConstants.H"
using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"

    pimpleControl pimple(mesh);

#   include "readGravitationalAcceleration.H"
#   include "initContinuityErrs.H"
#   include "createFields.H"
#   include "createControls.H"
#   include "correctPhi.H"
#   include "CourantNo.H"
#   include "setInitialDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readControls.H"
#       include "CourantNo.H"

        // Make the fluxes absolute
        fvc::makeAbsolute(phi, U);

#       include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        bool meshChanged = mesh.update();
        reduce(meshChanged, orOp<bool>());

#       include "volContinuity.H"

        volScalarField gh("gh", g & mesh.C());
        surfaceScalarField ghf("ghf", g & mesh.Cf());

        if (correctPhi && meshChanged)
        {
#           include "correctPhi.H"
        }

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi, U);

        if (checkMeshCourantNo)
        {
#           include "meshCourantNo.H"
        }
//


//

//
//
// Physical time value at current iteration
    dimensionedScalar Tn_
    (
        "Tn",
        dimensionSet(0,0,1,0,0,0,0),
        runTime.value()
    );
// Dummy unit scale acceleration term
    dimensionedVector dummyAcc_
    (
        "dummyAcc",
        dimensionSet(0,0,0,0,0,0,0),
        vector(0, 1, 0)
    );
// Dimensioned acceleration term that is push to UEqn
    dimensionedVector refFrameAcc_
    (
       "refFrameAcc_",
       dimensionSet(1,-2,-2,0,0,0,0),
       vector::zero
    );

//Period
 dimensionedScalar period_
 (
        "period",
        dimensionSet(0,0,1,0,0,0,0),
        1.0372058
    );
//amplitude
dimensionedScalar amplitude_
 (
        "amplitude",
        dimensionSet(0,1,0,0,0,0,0),
        0.005
    );

		

//acceleration.dimensions();
//acceleration.val()[1]=rho*(-amplitude_*(sqr(2*pi/period_)) * (sin(2*pi*(Tn_)/period_)));
forAll(acceleration, cellI)
{
acceleration[cellI][1]= (rho[cellI]*(amplitude_*(sqr(2*pi/period_)) * (sin(2*pi*(Tn_)/period_)))).value();
}


//amplitude
//dimensionedScalar roo_
// (
//        "amplitude",
//        dimensionSet(1,-3,0,0,0,0,0),
//        1000
//    );


// Generate a scalar of the current reference frame acceleration magnitude
 //dimensionedScalar currAcc_ =
 //(-amplitude_*(sqr(2*pi/period_))*(sin(2*pi*(Tn_)/period_)));
//rho*(-amplitude_*(sqr(2*pi/period_)) * (sin(2*pi*(Tn_)/period_)));


// Scalar multiplication to obtain acceleration vector and recast
 //  acceleration[1] = currAcc_*rho ;
       // currAcc_ * dummyAcc_;

//Info << "Current Frame Acceleration: " << acceleration << endl;
//
//
//

        // Pressure-velocity corrector
        while (pimple.loop())
        {
            twoPhaseProperties.correct();

#           include "alphaEqnSubCycle.H"

#           include "UEqn.H"

            // --- PISO loop
            while (pimple.correct())
            {
                    #           include "pEqn.H"
            }

            p = pd + rho*gh;

            if (pd.needReference())
            {
                p += dimensionedScalar
                (
                    "p",
                    p.dimensions(),
                    pRefValue - getRefCellValue(p, pdRefCell)
                );
            }

            turbulence->correct();
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
