//     	This file is included in rhoEnergyFoam.C if
//     	the flag for the pressure artificial diffusion
//     	is activated.
//     	Add AUSM dissipation on pressure term.

        //c = Foam :: sqrt ( thermo.Cp() / thermo.Cv() / psi ) ;
        //c = Foam :: sqrt ( 1.4 / psi ) ;


//     	Reconstruct p_L; p_R; rho_L; rho_R; U_L; U_R; 
     	
//      _L = _j is relative to the positive flux while _R = _j+1 
//      to the negative flux. The positive flux is directed from 
//      the owner cell to the neighbour cell and the negative flux 
//      is directed from the neighbour cell to the owner. 
//      Here the function interpolate is always used to 
//      compute the values at the intercell. Here the first 
//      argument is the field on which do the interpolation. "pos"
//      and "neg" are dimensionless scalar field difined in
//      createFields.H. "pos" has value 1 and "neg" has value -1.

//	In OF, surface fields have just one single value which is 
//	related to the face, but when you apply Godunov-like scheme 
//	you actually want to compute the flux through a face solving 
//	(maybe approxiamtely) the Riemann problem defined by the 
//	TWO different states of the fluid on the opposite sides 
//	of each face. Thus, for each unknown (say p, U and T), you 
//	need to evaluate two surface fields with the values of the 
//	variable on the two sides of the face. The fvc::interpolate 
//	function normally use the sign of the velocity flux to "decide" 
//	how to upwind the cell values in order to interpolate a face value. 
//	So the idea was to "fool" the interpolation routine calling it twice: 
//	once with a "positive" flux (e.g. +1) and one with a "negative" flux 
//	(e.g. -1). Depending on the selected interpolation scheme, you 
//	will get two (most likely) different estimates of the face 
//	values from which the Godunov-like scheme (in this case a 
//	central-upwind scheme) can be applied. The nice (and somewhat unusual) 
//	thing in OF is that given an internal face it always has an "owner" 
//	and a "neighbour": in OF notation the flux is positive if the flow 
//	is going through the face from the owner cell to the neighbour cell, 
//	and negative otherwise (i.e. the face normal is positive from 
//	owner to neighbour). Thus, because of the unstructured grid 
//	used in OF you cannot easily think interms of i-1,j and i,j but 
//	you better think in terms of owner/neighbour relation. So you actually 
//	solve a Riemann problem at each face between an "owner-side" and 
//	a "neighbour-side" with the convention that the flux will be positive 
//	in the "owner-to-neighbour" direction. It doesn't matter if 
//	(in a more conventional way of thinking) face normals are outward or 
//	inward with respect to any given cell. At the end, 
//	everything will be consistent.

       	p_L = fvc :: interpolate(p, pos, "reconstruct(p)") ;
       	
	p_R = fvc :: interpolate(p, neg, "reconstruct(p)") ;

       	rho_L = fvc :: interpolate(rho, pos, "reconstruct(rho)") ;
       	
	rho_R = fvc :: interpolate(rho, neg, "reconstruct(rho)") ;


//	On the velocity field U (which is a volVectorField) is
//	performed a interpolation in order to have a surfaceVectorField
//	and then is multiplied, with a scalar product, by the unit normal
//	of the intercell obtained by (mesh.Sf() / (mesh.magSf())). In 
//	this way we have now a surfaceScalarField at the left and right
//	side of the intercell.


	U_L  = fvc :: interpolate(U, pos, "reconstruct(U)") & (mesh.Sf() / (mesh.magSf())) ;
       	
	U_R  = fvc :: interpolate(U, neg, "reconstruct(U)") & (mesh.Sf() / (mesh.magSf())) ;

	
	const labelUList& own = mesh.owner();
       	
	surfaceScalarField duc = fvc :: interpolate(ducSensor) ;
       	
//	volScalarField gamma = thermo.Cp() / thermo.Cv() ;
       
//	surfaceScalarField gamm = fvc :: interpolate(gamma) ;
       
//	volScalarField cts = Foam :: sqrt(2. * (gamma-1.) / (gamma+1.)*H) ;
       
//	surfaceScalarField css = fvc :: interpolate(cts) ;
       
//	surfaceScalarField cs_R = css * css / max( css , -U_R ) ;
       
//	surfaceScalarField cs_L = css * css / max( css , U_L) ;
       
	c12 = fvc :: interpolate(c) ;
      
	M_L = U_L / c12 ;
       
	M_R = U_R / c12 ;
       
	m2 = ( U_L * U_L +U_R * U_R ) / ( 2. * c12 * c12 ) ;
       
	m0 = Foam :: sqrt ( min( pos , max( m2 , minfty ) ) ) ;
 	
//      Loop on all cells

//	"forAll" is a very quick and simple way of iterating 
//	over all items of a list. You can read it as: for all items in 
//	"own", look at each one as "iface". So "iface" will
//	be now the index of the fields in this for cycle. "iface"
//	is actually "own". Here in this code own is = mesh.owner() 
//	which means that "own" contains the list of cell which are
//	owner in the mesh. This is how it looks the mesh.owner file
//	(the owner file it can be found in the polymesh folder):
//
//	(
//	0 // this is face 0
//	0 // this is face 1
//	0 // this is face 2
//	1 // this is face 3
//	1 // this is face 4
//	... etc.
//
//	The owner and neighbour files each contain a list of cell numbers. 
//	The position in that list is the face number. 
//	So, face 2's "owner" is cell 0, face 4's "owner" is cell 1, and so on.
//	Note that the owner file is longer than the neighbour file, 
//	that's because the owner file also lists boundary faces, 
//	whereas the neighbour file does not.

      
       	


	forAll(own,iface)
       	
	{
        
		if ( duc[iface] > ducLevelPress )
        
		{
           
//      		Left and Right state
         
			scalar pl = p_L[iface] ;
        
			scalar pr = p_R[iface] ;
        
			scalar ml = M_L[iface] ;
        
			scalar mr = M_R[iface] ;
        
			scalar ul = U_L[iface] ;
        
			scalar ur = U_R[iface] ;
        
			scalar dl = rho_L[iface] ;
        
			scalar dr = rho_R[iface] ;

//
        
			scalar fa = m0[iface] * ( 2. - m0[iface] ) ;
        
			scalar alpha = 3. / 16. * ( -4. + 5. * fa * fa ) ; 

//
        
			scalar p5p = p5 ( ml , 1 , alpha ) ;
        
			scalar p5m = p5 ( mr , -1 , alpha) ;

//
        
			scalar dpr = p5 (mr , 1 , alpha) - p5m ;
        
			scalar dpl = p5p  - p5(ml, -1, alpha)  ;

//       
        
			scalar pu = -ku * p5p * p5m * ( dl + dr ) * c12[iface] * fa * ( ur - ul ) ;

//
        
			scalar dp12 = pr * dpr - pl * dpl ; 

// 
         

//      		Update p
         
			
			pave[iface] += duc[iface]*(- 0.5*(dp12) + pu)  ; // Pressure dissipation proportional to Ducros sensor
      
        	}
       	
	}


//                        Info << "pave in AUSM.H at "<< runTime.timeName() << "\n"  << pave << endl ;


//     	Add artificial dissipation on the processor boundaries
       
	forAll( mesh.boundaryMesh() , iPatch )
       	
	{
        	const polyPatch& patch = mesh.boundaryMesh()[iPatch] ;
        	
		if ((patch.type()=="processor"))
        	
		{
         	
			forAll( patch , iface )
         	
			{
          		
				if ( duc.boundaryField()[iPatch][iface] > ducLevelPress )
          	
				{
           	
//					Left and Right state
           
					scalar pl = p_L.boundaryField()[iPatch][iface] ;
           
					scalar pr = p_R.boundaryField()[iPatch][iface] ;
           
					scalar ml = M_L.boundaryField()[iPatch][iface] ;
           
					scalar mr = M_R.boundaryField()[iPatch][iface] ;
           
					scalar ul = U_L.boundaryField()[iPatch][iface] ;
           
					scalar ur = U_R.boundaryField()[iPatch][iface] ;
           
					scalar dl = rho_L.boundaryField()[iPatch][iface] ;
           
					scalar dr = rho_R.boundaryField()[iPatch][iface] ;

//
           
					scalar fa = m0.boundaryField()[iPatch][iface] * ( 2. - m0.boundaryField()[iPatch][iface] ) ;
           
					scalar alpha = 3. / 16. * ( -4. + 5. * fa * fa ) ; 

// 

			     		scalar p5p = p5 ( ml , 1 , alpha ) ;
           
					scalar p5m = p5 ( mr , -1 , alpha ) ;

//
           
					scalar dpr = p5 ( mr , 1 , alpha ) - p5m ;
           
					scalar dpl = p5p  - p5( ml , -1 , alpha ) ;

//       
           
					scalar pu = -ku * p5p * p5m * ( dl + dr ) * c12.boundaryField()[iPatch][iface] * fa * ( ur - ul ) ;

//
           
					scalar dp12 = pr * dpr - pl * dpl ; 
 
            
//	   				Update p

			     		pave.boundaryFieldRef()[iPatch][iface] += duc.boundaryField()[iPatch][iface] * ( -0.5 * ( dp12 ) + pu ) ;
        
          
				}
         
			}
        
		}
       
	}
