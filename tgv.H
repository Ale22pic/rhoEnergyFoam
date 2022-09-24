// Initial velocity and pressure field with Taylor-Green vortex
   forAll(U,cellI)
    {
        U[cellI].x()=Foam::sqrt(1.4)*0.01*Foam::sin(U.mesh().C()[cellI].x())
        *Foam::cos(U.mesh().C()[cellI].y())
        *Foam::cos(U.mesh().C()[cellI].z());
//
        U[cellI].y()= -Foam::sqrt(1.4)*0.01*Foam::cos(U.mesh().C()[cellI].x())
        *Foam::sin(U.mesh().C()[cellI].y())
        *Foam::cos(U.mesh().C()[cellI].z());

        U[cellI].z()=0.;

        p[cellI]=1.+pow((Foam::sqrt(1.4)*0.01),2)*
        (
         (Foam::cos(2*U.mesh().C()[cellI].z())+2.)*  
         (Foam::cos(2*U.mesh().C()[cellI].x())
        +Foam::cos(2*U.mesh().C()[cellI].y()))-2.)/16;
    }

//        thermo.correct();


