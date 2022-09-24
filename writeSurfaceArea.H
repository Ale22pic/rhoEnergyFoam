       forAll( mesh.boundaryMesh(), iPatch )
       {
        const polyPatch& patch = mesh.boundaryMesh()[iPatch] ;
        label patchID = mesh.boundaryMesh()[iPatch];
        if ((patch.type()=="wall"))
        {
         surfaceArea.open (label, std::ofstream::out | std::ofstream::app);
         forAll( patch , iface )
         {
          std::ofstream surfaceArea;
          surfaceArea <<  mesh.Sf().boundaryField()[iPatch][iface] << std::endl;
          surfaceArea.close();
         }
        }
       }
