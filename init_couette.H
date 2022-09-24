// Initialize with Poiseulle velocity profile
   forAll(U,cellI)
    {
    //U[cellI].x()=Foam::sqrt(1.4)*rm*0.5*(1.+U.mesh().C()[cellI].y());
    //U[cellI].x()=Foam::sqrt(1.4)*rm*(U.mesh().C()[cellI].y());
      U[cellI].x()=Foam::sqrt(1.4)*rm*(U.mesh().C()[cellI].y()+Foam::sinh(3.*U.mesh().C()[cellI].y()))/(1.+4.*Foam::sinh(3.*U.mesh().C()[cellI].y()));
//      U[cellI].x()=3./2.*Foam::sqrt(1.4)*rm*(1.-Foam::pow(U.mesh().C()[cellI].y(),8));
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
