//    Evaluate pressure gradient for channel flow
      dimensionedScalar dpdx
      (
       "dpdx",
       dimensionSet(1,-2,-2,0,0,0,0),
       scalar(1.)
      );
//
      scalar dpdx_temp = 0;
      scalar drhoUx    = 0;
      scalar Tbulk     = 0;
      scalar vol       = 0;
      scalar rhob      = 0;
      volScalarField rhoUx = rhoU.component(0) ;
      volVectorField fluxU   = fvc::div(momFlux) ;
      volScalarField fluxUx  = fluxU.component(0)  ;
//      forAll (flUx, cellI)
      forAll (fluxUx, cellI)
      {
       dpdx_temp += momFl[cellI].component(0)*mesh.V()[cellI];//fluxUx[cellI]*mesh.V()[cellI];
       drhoUx    += rhoUx[cellI]*mesh.V()[cellI];
       Tbulk     += rhoUx[cellI]*T[cellI]*mesh.V()[cellI];
       rhob      += rho[cellI]*mesh.V()[cellI];
       vol       += mesh.V()[cellI];
      }
//
      reduce(dpdx_temp,sumOp<scalar>());//MPI all reduce
      reduce(drhoUx,sumOp<scalar>())   ;//MPI all reduce
      reduce(Tbulk,sumOp<scalar>())    ;//MPI all reduce
      reduce(vol      ,sumOp<scalar>());
      dpdx       = dpdx_temp/vol*dpdx;
      dpdx_prt   = dpdx;
      rhob       = rhob/vol;
      drhoUx     = drhoUx/rhob/vol;
      Tbulk      = Tbulk/vol/drhoUx;
      drhoUx_prt = drhoUx;
      Tbulk_prt  = Tbulk;
      rhob_prt   = rhob;

