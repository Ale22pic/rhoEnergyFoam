     double tke = 0. ;
     double vol = 0. ;
     double ens = 0. ; //To compute enstrophy during simulation

     forAll (U,icell) {
     scalar uu = U[icell].component(0);
     scalar vv = U[icell].component(1);
     scalar ww = U[icell].component(2);
     tke += 0.5*mesh.V()[icell]*(uu*uu+vv*vv+ww*ww) ;
     vol += mesh.V()[icell];
     }
     tke = tke/vol;

     volVectorField Vorticity=Foam::fvc::curl(U);

     forAll (Vorticity,icell) {
     scalar v1 = Vorticity[icell].component(0);
     scalar v2 = Vorticity[icell].component(1);
     scalar v3 = Vorticity[icell].component(2);
     ens += 0.5*mesh.V()[icell]*(v1*v1+v2*v2+v3*v3) ;
     vol += mesh.V()[icell];
     }
     ens = ens/vol;
//
     if (Pstream::master())//master process
     {
       std::ofstream diagn;
       diagn.open ("diagnostics.dat", std::ofstream::out | std::ofstream::app);
       diagn << runTime.timeName() << " " <<  tke  << " " << ens << " " << dpdx_prt.value() << " " << drhoUx_prt.value() << " " << Tbulk_prt.value() << " " << rhob_prt.value() << std :: endl ;
       diagn.close();
}
