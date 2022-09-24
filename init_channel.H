// Initialize with Poiseulle velocity profile
   forAll(U,cellI)
    {
	     float rand1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	     float rand2 = 0.5*(2*rand1-1);
//    U[cellI].x()=3./2.*Foam::sqrt(1.4)*rm*(1.-Foam::pow(U.mesh().C()[cellI].y(),2));
//    U[cellI].x()=3./2.*Foam::sqrt(1.4)*rm*(1.-Foam::pow(U.mesh().C()[cellI].y(),4));
      U[cellI].x()=Foam::sqrt(1.4)*rm*(1.-Foam::pow(U.mesh().C()[cellI].y(),4)) + rand2 ;
//	if (U.mesh().C()[cellI].y()<0.)
//	{
 //        U[cellI].x()= Foam::tanh(0.5*(U.mesh().C()[cellI].y()+1.)/0.4);
//	}
//	 else if (U.mesh().C()[cellI].y()>0.)
//	{
 //        U[cellI].x()= Foam::tanh((1.-0.5*(U.mesh().C()[cellI].y()+1.))/0.4);
//	}
//
        U[cellI].y()= 0.;//05*Foam::sqrt(1.4)*rm*Foam::sin(0.5*M_PI*U.mesh().C()[cellI].y())*Foam::cos(2.*M_PI*U.mesh().C()[cellI].z()); 
        U[cellI].z()= 0.;//05*Foam::sqrt(1.4)*rm*Foam::sin(0.5*M_PI*U.mesh().C()[cellI].y())*Foam::sin(2.*M_PI*U.mesh().C()[cellI].z()); 
    }


/*      scalar dUx    = 0;
      scalar vol       = 0;
      volScalarField Ux = U.component(0) ;
      forAll (Ux, cellI)
      {
       dUx    += U[cellI].component(0)*mesh.V()[cellI];
       vol    += mesh.V()[cellI];
      }
//
      reduce(dUx,sumOp<scalar>())   ;//MPI all reduce
      reduce(vol      ,sumOp<scalar>());
      dUx     = dUx/vol;

      forAll (Ux, cellI)
      {
       U[cellI].component(0) = U[cellI].component(0) - dUx + Foam::sqrt(1.4)*rm;
      }
//
*/
