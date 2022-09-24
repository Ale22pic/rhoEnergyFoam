// Evaluate Ducors sensor

// First are initialized the elements neccessary
// in the Ducors sensor formula. Here is 
// calculated the divergence of U, the 
// square divergence of and the magnitude of 
// the square rotor of U.


   volScalarField divU  = fvc :: div(U) ;

   volScalarField divU2 = divU * divU ;

   volScalarField rotU2 = magSqr( fvc :: curl(U) ) ;


// Here the Ducors sensor is applied to all
// the internal field cells in order to detect
// where there is a compression due to a shock.


   forAll( U, icell )

   {
	   ducSensor[icell] = max(-divU[icell]/Foam::sqrt(divU2[icell] + rotU2[icell] + epsilon),0.) ; 
   }


// Here the Ducors sensor is applied to all the
// boundary field cells.

   forAll( mesh.boundaryMesh(), iPatch )
   {
	   const polyPatch& patch = mesh.boundaryMesh()[iPatch] ;
    	
	   if ((patch.type()=="processor") )
     
	   {
      
		   forAll( patch , iface )
      
		   {
			   ducSensor.boundaryFieldRef()[iPatch][iface] = 
				   
				   max(-divU.boundaryField()[iPatch][iface] /
						   
				   Foam::sqrt(divU2.boundaryField()[iPatch][iface] + 
						  
			           rotU2.boundaryField()[iPatch][iface] + epsilon),0.) ; 
 
      	           }
           
	   }
   
    }

