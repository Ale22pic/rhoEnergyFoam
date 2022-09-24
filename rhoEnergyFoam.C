//  =========                 |
//  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
//   \\    /   O peration     |
//    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
//     \\/     M anipulation  |
//-------------------------------------------------------------------------------
//License
//
//    This file is part of OpenFOAM.
//
//    OpenFOAM is free software: you can redistribute it and/or modify it
//    under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
//    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//    for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
//
//Application
//
//    rhoEnergyFoam 
//
//Description
//
//    Numerical solver for the solution of compressible shock-free flows.
//    The convective terms are discretized using Pirozzoli's scheme(JCP 2010),
//    density based splitting. The scheme allows conservation of the kinetic
//    energy in the inviscid incompressible limit.
//    midPoint interpolation must be selected in fvScheme dictionary. 
//    Viscous terms(Laplacians only) are evaluated directly, computing
//    the face normal gradients.
//    A third-order low-storage RK scheme is used for time integration.
//    The OpenFOAM labrary for turbulence models is included.
//      
//    Author: Davide Modesti (davide.modesti@uniroma1.it)
//    
//    Last update 06/04/2017
//
//    Reference
//
//    D. Modesti and S. Pirozzoli A high-fidelity solver for 
//    turbulent compressible flows on
//    unstructured meshes. Comput. & Fluids (2017)
//
//---------------------------------------------------------------------------

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "fixedRhoFvPatchScalarField.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include <fstream>      // std::ofstream

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//  Here is included StopWatch.H which is a file that let to count
//  the seconds required by the code for each line you want. In order
//  to use it, after "StopWatch.H" write a line with "StopWatch <Name of the line want analyze>"
//  and then go to the desired line and write "<Name of the line want analyze>.start()" before
//  the line and after the line "<Name of the line want analyze>.stop()". At the end
//  of the code there are Info messeges who tell the total second required. 
//  The profiling lines are indicated with //P. Check the following lines in 
//  the code to understand better how it works.


//P
    #include "StopWatch.H"
    StopWatch totalTime;
    StopWatch cycleTime;
    StopWatch turbvalTime;
    StopWatch oldvarTime;
    StopWatch flagIfTime;
    StopWatch soundspeedTime;
    StopWatch rhoaveTime;
    StopWatch UaveTime;
    StopWatch phiTime;
    StopWatch phitTime; 
    StopWatch HTime;
    StopWatch HaveTime;
    StopWatch paveTime;
    StopWatch sensorTime;
    StopWatch AUSMPTime;
    StopWatch vecDivUTime;
    StopWatch muEffTime;
    StopWatch tauMCTime;
    StopWatch kTime;
    StopWatch kaveTime;
    StopWatch momVisFluxTime;
    StopWatch heatFluxTime;
    StopWatch visWorkTime;
    StopWatch enVisFluxTime;
    StopWatch TotFluTime;
    StopWatch AUSMCTime;
    StopWatch TotDivFluTime;
    StopWatch RKSubTime;
    StopWatch correctTime;
    StopWatch turbcorrTime;
//P


// The formulas in this code are all written in the document 
// "RhoEnergyFoam Formulas"

// AUSM useful functions

// It is important to remember that, the inputs inside
// the function are just a way to construct the function
// itself. As in all the programming languages, once I call
// the function I can rename the inputs. The function will 
// know that the first input is, for exemple, rm and the second
// input sgn.

// m1 is a function that calculates the split Mach number as
// a polynomial function of 1st degree, that's what 1 means in
// "m1". Related formula (1). It returns M_1. Here and on in
// the formulas "rm" is the Mach number, while "sgn" is the
// sign functione (which is equal to -1 if the value is negative;
// equal to 0 if the value is zero and equal to 1 if the 
// value is positive. "mag" is the function magnitude.


//P
//   m1Time.start();
//P

   float m1 (float rm , float sgn ) // Return M_1
   
   {
   
   float r ;  

   r = 0.5*(rm + sgn * mag( rm )) ;
   
   return r ; 
   
   }

//P
//  m1Time.stop();
//P


// m2 is a function that calculates the split Mach number as
// a polynomial function of 2nd degree, that's what 2 means in
// "m2". Related formula (2). It returns M_2.  


//P
//  m2Time.start();
//P

   float m2 (float rm , float sgn ) // Return M_2
   
   {
   
   float r ; 
   
   r = sgn * 0.25 * ( rm + sgn ) * ( rm + sgn ) ;
   
   return r ; 
   
   }

//P
//  m2Time.stop();
//P


// p5 is a function in terms of split Mach number, usefull
// in some functions of the AUSM method. 
// It relays on the value of rm. So if rm is less than 1,
// p5 is a polynomial function of 5th order ("5" in p5
// means the order). If rm is more than 1, the value of 
// p5 is just 1 (due to the fact that m1 in considered in 
// the formula. Related formula (3). The "alpha" in the 
// function is defined in the file AUSM.H. 


//P
//  p5Time.start();
//P

   float p5 (float rm , float sgn , float alpha ) // Return p_5
   
   {
   
   float r ; 
   
   if (abs ( rm ) < 1. )
   
   {
    
    r = m2 (rm , sgn) * ( (sgn * 2 - rm ) - sgn * 16 * alpha * rm * m2( rm , -sgn ) ) ;
   
   }
   
   else
   
   {
   
    r = 1. / rm * m1 ( rm , sgn) ;
   
   }
   
   return r ; 
   
   }

//P
//  p5Time.stop();
//P


// m4 is a function that calculates the split Mach number as
// a polynomial function of 4th degree, that's what 4 means in
// "m4". Related formula (4). It returns M_4. Even here the value
// of m4 depends on rm value. So if rm is > 1, m4 is equal to m1.
// If rm < 1, m4 is 4th degree polynomials. The "beta" in the 
// function is a constant that is defined in the file AUSM.H.


//P
//  m4Time.start();
//P
 
   float m4 (float rm , float sgn , float beta) // Return M_4
   {

   float r ; 
   
   if ( abs ( rm ) < 1 )
   {

    r = m2 ( rm , sgn ) * ( 1 - sgn * 16 * beta * m2 ( rm , -sgn ) ) ;
   
   }

   else

   {

    r = m1 ( rm , sgn ) ;
   
   }

   return r ; 

   }

//P
//  m4Time.stop();
//P

// Main


int main(int argc, char *argv[])

// Here we start to include the files necessary for the
// solver. Descriptions for each one is given below.
 

{

    #define NO_CONTROL
    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createTimeControls.H"
    #include "readThermophysicalProperties.H"
    #include "variables.H"


// The "postProcess.H" file execute application functionObjects 
// to post-process existing results. If the "dict" argument is 
// specified the functionObjectList is constructed 
// from that dictionary otherwise the functionObjectList 
// is constructed from the "functions" sub-dictionary 
// of "system/controlDict". 

// The "addCheckCaseOptions.H" file gives the istruction to 
// check the set-up only using a single time step. 

// The "setRootCaseLists.H" file declare some "standard" list options.

// The "createTime.H" file create the time function as "runTime".

// The "createMesh.H" file create a fvMesh (specified region or defaultRegion) 
// with additional handling of -dry-run and -dry-run-write options (which
// mean to run for one step the simulation).

// The "createFields.H" file creates the necessary field for
// the simulation (such as volScalarField, surfaceScalarField, etc.) for
// exemple rho or U fields.

// The "createFieldRefs.H" file sets field such as p, T, psi and mu.

// The "createTimeControls.H" file reads the control parameters 
// used by setDeltaT. It searchs for parameters such as 
// adjustTimeStep, maxCo and maxDeltaT in the system/controlDict.

// The "readThermophysicalProperties.H" file looks for the Pr number
// in order to calculate "k" which is not accessible through "thermo".

// The "variables.H" file sets variables as c_L; M_R; rhoR etc. which
// are useful in the AUSM methods, which is called including AUSM.H.	


//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//P
    totalTime.start();
//P



//  Check on the turbulence model. The function "validate()" just calls 
//  the "correctNut" function (which calculate the new
//  value of the turbulent viscosity). This just makes sure that 
//  the turbulent viscosity is computed correctly before starting the
//  loop. 


//P
    turbvalTime.start();
//P

    turbulence->validate() ;

//P
    turbvalTime.stop();
//p
    
//  the object named turbulence is a pointer, since the functions of the object
//  are called using "->".The operation of those member functions 
//  depend on which turbulence model we are using.


//  Info message about the beginning of the loop

    
    Info << nl << "Starting time loop" << endl ;

    
//  Info message about the starting time


    Info << "Start Timing = " << runTime.clockTimeIncrement() << " s" << nl << endl ;


//  Initialize the counter for the iterations


    int Iter = 0 ;


    //Info << "Initialize Channel" << endl ;

    //#include "init_channel.H"

    //rhoU = rho*U;
    //rhoE = rho*(e + 0.5*magSqr(U)); 


//  Start time loop

    while (runTime.run()) 
    

    {
     

// 	Iteration counter
	

        Iter++ ; // Add 1 to the variable Iter
        

	Info << "Iter No. = " << Iter << endl ; // Print the variable Iter


//	Give info message abount the time of the current cycle using the 
//	function "timeName()" of the class "runTime"


        Info << "Iteration starting time = " << runTime.timeName() <<  endl ;
        

	Info << "Current iteration deltaT = " <<  runTime.deltaTValue() << endl;

	
	runTime++ ;

//      Saving quantities at preavious time step, those are resolved 
//      in the Runge-Kutta (RK) of 4-th order cycle.


//P
  	oldvarTime.start();
//P


	rhoOld = rho ; 
       
	
	rhoUOld = rhoU ; 
       
	
	rhoEOld = rhoE ; 


//P
 	oldvarTime.stop();
//P



//	Evaluation of the activation flags about AUSM diffusion terms
//	in order to optimize the calling of the file "sensor.H". If the
//	flags are set on true in the controlDict, a value of 1 it will
//	assigned, else 0. Before the if cycle, the flags are initialized.


//P
  	flagIfTime.start();
//P


	int flagP = 0 ;

	int flagC = 0 ;

	if ( pressArtDiff )

                {
                        flagP = 1 ;
                }

        if ( convArtDiff )

                {
                        flagC = 1 ;
                }


//P
  	flagIfTime.stop();
//P


//      RK Time step. The indices here is "cycle" and it goes till is
//      less than rkCoeff size, which is 4. It goes till 3 because in 
//      C++ the indices start from 0, so since rkCoeff in RK4 are 4 
//      coefficients, the rkCoeff is a vector of 4 scalars, and 
//      the indices of this vector goes from 0 to 3.
       

//P
        cycleTime.start(); //Start profiling Runge-Kutte cycle
//P


	for (int cycle =0 ; cycle < rkCoeff.size() ; cycle++)


	{


//      	Speed of sound and Mach number using thermo.Cp, thermo.Cv and
//      	psi which are defined in createFields.H and createFieldRef.H, 
//      	so they are a way to express Cp, Cv and Psi in OpenFOAM. The
//      	Mach number is simply calculated taking the velocity field
//      	"U" and dividing each element by "c", so even Mach will be
//      	a volVectorField, a file that contain the vector component
//      	(which are 3) of Mach for each cells. So if cells in the 
//      	mesh are 12, then we will have 12 vectors in the file. 
//      	The related formulas are (5) and (6).
        

//		c = Foam :: sqrt ( thermo.Cp() / thermo.Cv() / psi ) ;


//P	
  		soundspeedTime.start();
//P


		c = Foam :: sqrt ( 1.4 / psi ) ;


//P
  		soundspeedTime.stop();
//p


//		Mach = U / c ;


//      	Interpolated quantities at cell faces. A surfaceScalarField is 
//      	defined on the faces of the mesh, so, for exemple rhoave[1] 
//      	gives you the flux in the second face in the mesh (remember 
//      	that indices start from 0). 
//      	Face ordering is defined by the faces file in polyMesh. 
//      	For each of the faces, you can get the owner/neigbour cell 
//      	info by asking the mesh, using:
//
//      	label faceOwnerCell = mesh.owner()[1];
//		label faceNeighbourCell = mesh.neighbour()[1];
//
//		In surfaceScalarFields, the internal faces belong to the 
//		internalField() and patch fields are done in the usual way. 
//		Note the limitations on the patch field types: you 
//		cannot have, say, a zeroGradient patch field 
//		on the surfaceScalarField. 
//		By doing rhoave[1], you are really getting 
//		rhoave.internalField()[1] because the Geometric 
//		field is derived from the internal field.
//		Interpolation schemes are used to transform 
//		cell-centre quantities to face centres. The operation 
//		is used in many of the finite volume calculations, 
//		e.g. for the calculation of gradient, divergence, 
//		and Laplacian terms. In this case we start from a 
//		volVectorField which is "U" and go to a surface
//		vector field by interpolating. So basically we 
//		are bringing the value from the cell center to 
//		the cell faces, useful for finite volume method.
//		The interpolation scheme is obtained from 
//		system/controlDict, under the subDict InterpolationSchemes
//		and in the case of rhoEnergyFoam is always set as 
//		default the midPoint scheme, as it can be seen in
//		the tutorials file.
     
//P
  		rhoaveTime.start();
//P

        	surfaceScalarField rhoave = fvc :: interpolate ( rho ) ;


//P
  		rhoaveTime.stop();
//P


//P
  		UaveTime.start();
//P 


		surfaceVectorField Uave = fvc :: interpolate ( U ) ;


//P
  		UaveTime.stop();
//P


//      	Calculation of the flux at the intercell	


//		"phi" is the flux. The flux is for compressible solvers 
//		rho*u*A ("phit" in this case) or otherwise u*A ("phi") 
//		on the faces. Or more precisely the flux between cells. 
//		Your velocity field U in your timestep folders is saved 
//		at the center of each cell. The flux on the other hand 
//		is the value on the faces between cells, the flow from 
//		one cell to the next. Therefore it is a 
//		surfaceScalarField and not a volScalarField like rho. 
//		You can calculate your flux with different methods 
//		(upwind, linear, etc).


//P
  		phiTime.start();
//P

  
//        	phi = fvc :: interpolate ( U ) & mesh.Sf() ;     // cambiare phi con flux
        	phi = Uave & mesh.Sf() ;     // cambiare phi con flux

//P
  		phiTime.stop();
//P


//P
 		phitTime.start();
//P

        	phit = fvc :: interpolate ( rhoU ) & mesh.Sf() ; // Flux in turbulent model     //cambiare phit con fluxt

//P
  		phitTime.stop();
//P


//		"fvc ::" is finite volume calculus (it is a class)
//		and is used for calculations on those fields. 
//		"fvm ::" on the other hand
//		is used for equation systems. So interpolate(U)
//		(which is a function of the class fvc) 
//		gives you the velocity on the patch, mesh.Sf is 
//		the surface vector (area and normal of the face, so
//		is a vector which is normal to the surface and its
//		magnitude is equal to the surface area). mesh.Sf()
//		contains a vector for each cell surface.
//		And "&" is the scalar product. Hence the scalar 
//		product of "fvc::interpolate(U)" (which is a 
//		surfaceVectorField) and mesh.Sf() (which is also
//		a surfaceVectorField) and represents the area, gives
//		you the flux. While U is m/s your phi (or phit) file 
//		tells you the flux is m^3/s in the case of incompressible
//		flow, while kg/s for compressible flow.


//      	Enthalpy calculation. It is a volVectorField, so a file
//      	that containt a value of enthalpy for each cell.


//P
  		HTime.start();
//P  


		H = ( rhoE + p ) / rho ;


//P
  		HTime.stop();
//P


//      	Enthalpy at the intercell, in order to bring the
//      	enthalpy at the intercell and create a 
//      	surfaceScalarField.


//P
  		HaveTime.start();
//P


        	surfaceScalarField Have = fvc :: interpolate ( H ) ;


//P
  		HaveTime.stop();
//P


//      	Pressure at the intercell, so from a volVectorField
//      	we pass to a surfaceScalarField, having the pressure
//      	at the intercell.


//P
  		paveTime.start();
//P


        	surfaceScalarField pave = fvc :: interpolate ( p ) ;    


//P
  		paveTime.stop();
//P

//                Info << "pave before AUSM.H"<< runTime.timeName() << "\n"  << pave << endl ;


//      	Here are used the flags about the activation of the
//      	dissipation. If at least one flag is activated the 
//      	file "sensor.H" it will be included, else not. This is 
//      	done in order to optimize the code.      	 


                if ( flagP + flagC > 0 )

                {
//P
		      sensorTime.start();
//P
 
		      #include "sensor.H"

//P
		      sensorTime.stop();
//P

                }


//              Activation of the dissipation on pressure terms if the 
//              flag pressArtDiff is found set on true inside the controlDict


//P
  		AUSMPTime.start();
//P


        	if ( pressArtDiff )

        	{
                

		      #include "AUSM.H"      // AUSM+up dissipation on 
					     // pressure term, 
					     // add dissipation on "pave".
        
		}


//P
  		AUSMPTime.stop();
//P


//      	Evaluate viscous terms

//      	Divergence of the velocity is evaluated. In this part from a volScalarField
//      	it is brougth to the intercell by using interpolate and now we
//      	have a surfaceScalarField. It is constructed in three same 
//      	component since in the N-S equation, in the viscous terms
//      	the divergence of the velocity U is multiplied then for an
//      	identity matrix. We need vecDivU in this form below in order 
//      	to subtract it to the snGrad(U). This operations are done
//      	next in the calculation of viscous momentum flux.        


//P
  		vecDivUTime.start();
//P


		vecDivU.component(0) = fvc :: interpolate ( fvc :: div(U) ) ;
        
		vecDivU.component(1) = vecDivU.component(0) ;

		vecDivU.component(2) = vecDivU.component(0) ;
			
		vecDivU.setOriented(true) ;


//P
  		vecDivUTime.stop();
//P 


// 		The setOriented function is only included in the ESI version and
// 		not in the Foundation. Here is reported the explaination:
// 		We introduced the oriented flag for surface fields 
// 		to discriminate between flux fields that are oriented, 
// 		i.e. where the sign depends on owner->neighbour connectivity, 
// 		and plain interpolated fields, e.g. pressure on faces.
//		If you apply a field operation that includes the area vector 
//		field Sf() the resultant field becomes oriented; similarly, 
//		taking, e.g. the mag() of an oriented field makes it 'unoriented'. 
//		The rules applied when dealing with [un]oriented fields also 
//		gives us an extra layer of protection/verification 
//		that the mathematics is consistent.
//		If the field is constructed without using field operations 
//		and it is a flux field, you should set the oriented type to true.
//
//		If the line vecDivU.setoriented is not call there will be a problem
//		in the next calculations since there are operator like "-" (minus) 
//		which require an oriented or unoriented field since "-" in the next calculation
//		is used in the definition of a surfaceScalarField, which is 
//		oriented. If there isn't setOriented line in the vecDivU file
//		it won't be specified if it is oriented or unoriented.
        	
		
//		Here from the turbulence model is taken the effective dinamic
//		viscosity. It is constructed in "turbulenceModel.H" and it is
//		a volScalarField, so the value for each center cell. The muEff
//		is muEff = mu_laminar + mu_turbulent. Here we can't write muEff
//		as "volScalarField turbulence -> muEff() ;" or 
//		"volScalarField muEff() ;" or "volScalarField muEff ;" or
//		"turbulence -> muEff();" but we 
//		have to write as it is written below in order to declare a 
//		volScalarField named muEff. 


//P
  		muEffTime.start();
//P

		volScalarField muEff( turbulence -> muEff() ) ;


//P
  		muEffTime.stop();
//P


//		This is declaring a volTensorField object tauMC (which is
//		the stress tensor),  
//		and initialising it with its name and its contents, 
//		where the latter is calculated from the expression 
//		after the comma. "fvc::" is the finite volume calculus
//		class; grad() is the gradient operator that acts on U
//		and gives a rank tensor. "Foam::T()" is the class
//		Foam which contains the function "T" which is simply
//		the function that give the transpose a tensor. While "dev2" is 
//		the function that give deviatoric part of a tensor.
//		Remember that The stress tensor can be separated into two 
//		components. One component is a hydrostatic or dilatational 
//		stress that acts to change the volume of the material only; 
//		the other is the deviatoric stress that acts to change the shape only. 
//		All this tensor is then multiplied by muEff. Sometimes in the
//		expression we can also find "(!inviscid)" time the expression,
//		which is an expression that is either zero if the flow 
//		is inviscid or 1 if the flow is not inviscid.
//		So this line calculates the shear stress tensor.


//P
  		tauMCTime.start();
//P


                volTensorField tauMC( "tauMC" , muEff * dev2 (Foam :: T ( fvc :: grad (U) ) ) ) ;


//P
  		tauMCTime.stop();
//P


//		Here the muEff, from the center cell is brought to the
//		cell faces (or intercell) by simply interpolating, so 
//		from a volScalarField we go to a surfaceScalarField.

		surfaceScalarField muave = fvc :: interpolate ( muEff ) ; 

//		Here is defined the thermal conductivity "k" with the
//		name "k" and the expression Cp times muEff and divided
//		by the Prandtl number Pr. We remember that the
//		effective thermal diffusivity alphaEff is 
//		equal to muEff/Prt where Prt is the turbulent Prandtl
//		number.


//P
  		kTime.start();
//P


                volScalarField k( "k" , thermo.Cp() * muEff / Pr ) ; 


//P
  		kTime.stop();
//P


//		Then here from the center cell k is interpolated to
//		the intercell.


//P
  		kaveTime.start();
//P


        	surfaceScalarField kave = fvc :: interpolate ( k ) ;


//P
  		kaveTime.stop();
//P


//        	Momentum viscous flux is calculated here. Let's see in 
//        	details what functions there are. "fvc::snGrad()" 
//        	calculates the surface normal gradient. Since the surface 
//        	normal is only one direction, you will get one component 
//        	for every component of velocity. However if you 
//        	input a scalar field, you will get only a scalar. 
//        	So the output really depends on the input. So the 
//        	full gradient of U at a face can be interpolated from
//        	the gradient at the center cell. Different schemes are
//        	available based on the angle between the vector that
//        	joins the two cell centers and the normal of the surface,
//        	and this angle represents the degree of non-orthogonality.
//        	For exemple, in the "ForwardStep" tutorial, an "uncorrected"
//        	snGradScheme is used. "mesh.magSf()" represents the face
//        	area. 


//P
  		momVisFluxTime.start();
//P


		surfaceVectorField momVisFlux = muave * ( fvc :: snGrad(U) * mesh.magSf() ) ; 
			       	

//P
  		momVisFluxTime.stop();
//P


//		Energy viscous flux. The functions here are the same used upper,
//		here the snGrad is applied on the temperature T, which is a scalar.
  

//P
  		heatFluxTime.start();
//P


		surfaceScalarField heatFlux =  kave * fvc :: snGrad(T) * mesh.magSf() ;


//P
  		heatFluxTime.stop();
//P


//		The momVisFlux multiplied for Uave (which is a surfaceVectorField) 
//		doing a scalar product. So the scalar product between two 
//		surfaceVectorField give a scalar, so it will give a surfaceScalarField,
//		which is the viscous work. 


//P
  		visWorkTime.start();
//P


	        surfaceScalarField visWork = ( momVisFlux + 
				
				fvc :: dotInterpolate(mesh.Sf(), tauMC) ) & Uave ;


//P
  		visWorkTime.stop();
//P


//		Here is calculated the energy by summing the heat flux and
//		the viscous work.


//P
  		enVisFluxTime.start();
//P


		enVisFlux = heatFlux + visWork ;


//P
  		enVisFluxTime.stop();
//P


//		Total fluxes, Eulerian - Viscous. First is calculated
//		the flux of rho by multipling the rho at the intercell
//		by phi, which is the flux at the intercell (is given
//		by the multiplication between mesh.Sf() and Uave 
//		which is the velocity vector at the interface.  


//P
  		TotFluTime.start();
//P


		surfaceScalarField rhoFlux = rhoave * phi ; 
        	

		momFlux = rhoave * Uave * phi + pave * mesh.Sf() - momVisFlux ;     


        	enFlux = rhoave * Have * phi - enVisFlux ; 


//P
  		TotFluTime.stop();
//P


//		Activation of the dissipation on convective terms if the 
//		flag conArtDiff is found set on true inside the controlDict


//P
  		AUSMCTime.start();
//P


        	if ( convArtDiff ) 
        
		{
         
		
			#include "AUSM_conv.H"    // AUSM+up dissipation on convective 
						  // terms, add dissipation on 
						  // rhoFlux,
						  // momFlux and enFlux
        	}


//P
  		AUSMCTime.stop();
//P

  
//		After the calculation of the flux, the divergence is
//		applied on them in order to have the elements for
//		the resolution of the RK sub-step. This quantities
//		are the approximation of the sum of the integral of
//		the fluxes in the formula (2) in the Modesti-Pirozzoli
//		paper. Also from surfaceScalarField we are now back to
//		volScalarField and volVectorField,
//		so we bring the value from the interface of the 
//		cell to the center cell, since we are interested to it.		


//P
  		TotDivFluTime.start();
//P


	        volScalarField rhoFl = fvc :: div ( rhoFlux ) ;

        	volVectorField momFl = fvc :: div ( momFlux ) - fvc :: div ( tauMC ) ;

                volScalarField enFl = fvc :: div ( enFlux ) ;


//P
  		TotDivFluTime.stop();
//P


	      	//#include "pressureGrad.H" // forcing term to keep a constant 
      	                                  // mass flow rate in channel flow.
					  // This is a file useful in 
					  // turbulent channel flow.

//         	RK sub-step are performed. The way how it works is that for
//         	each iteration in time, the code resolve 4 sub-step of the RK
//         	method and at end of this 4 sub-step it has the value of 
//         	rhoU; rho and rhoE at that time step. The sub-step are four
//         	since we are dealing with a 4th order RK method. So to better
//         	explain how it works, At the first sub-step it calculates
//         	a value of rho, for exemple, with the first rkCoeff in 
//         	the list and using the current deltaT. Then save the value
//         	of rho and goes to the second sub-step where recall the 
//         	value of rho as rhoOld, and goes to the end of the sub-step
//         	calculating the new value of rho, using now the second
//         	rkCoeff of the list and the same deltaT as before. This is
//         	done till the end of the rkCoeff, which in our case are 4.
//         	Then the RK for this time iteration is finished, it can 
//         	goes on end evaluate a new deltaT and start a new time
//         	iteration.
  

//P
  		RKSubTime.start();
//P


		rho = rhoOld + rkCoeff[cycle] * runTime.deltaT() * ( -rhoFl ) ;         

 	        rhoU = rhoUOld + rkCoeff[cycle] * runTime.deltaT() * ( -momFl ) ;

	        rhoE = rhoEOld + rkCoeff[cycle] * runTime.deltaT() * ( -enFl ) ;


//P
  		RKSubTime.stop();
//P


//        	Update primitive variables and boundary conditions. 
//        	To understand better the meaning of' ()' after a variable
//              let's keep rho, which is a 'volScalarField' type object 
//              for fluid density. This object includes 'internalField' 
//              and 'boundaryField' info. rho() in the source terms is 
//              'volScalarField::Internal' type object wherein the 
//              boundary information is absent. This is fairly new functionality, 
//              and it is useful to reduce computational costs for parallel 
//              computations by reducing parallel communications which are 
//              mostly needed for boundaryFields rather than internalFields.
//              For a given 'vol*Field', say 'Object', '()' operator is 
//              defined, and you can call it by appending it to the given 
//              field, like 'Object()'. This 'Internal' field is only 
//              usable for sources on the right hand side of a constructed 
//              equation. Constructed equation are, for exemple, of the type: 
//
//		tmp<fvScalarMatrix> epsEqn
//		(
//			fvm::ddt(alpha, rho, epsilon_)
//			==
//			C1_*alpha()*rho()*G*epsilon_()/k_()
//		);
//
//		The equations below are not the case of constructed 
//		equation. In the calculation of U.ref() it is not
//		necessasy to declare U.ref, since U is a field
//		declared at the beginning of the code. Where with
//		".ref" we are defining a function of the object U
//		(which is a volVectorField type object). And so
//		U.ref() is equal to the devision of rhoU in the 
//		internal field by rho in the internal field. 


//P
  		correctTime.start();
//P

	
		U.ref() = rhoU() / rho() ;


//              The fields boundary condition are corrected


        	U.correctBoundaryConditions() ;

//		Once the U at boundary field is corrected we can
//		use it tho evaluate others fields at the 
//		bounadaries.

		rhoU.boundaryFieldRef() = rho.boundaryField() * U.boundaryField() ;

        	e = rhoE / rho - 0.5 * magSqr (U) ;
        
		e.correctBoundaryConditions();
        
//		Thermodinamic library is called. Thermo->correct updates 
//		thermodynamic properties. Another interesting info is
//		that thermo->rho calculates density, based on the 
//		thermodynamic model. The field "e" is required in the
//		thermo.correct, so this is the reason why it is
//		evaluated.
        
		thermo.correct() ;

//		"rhoE" at boundary field is evaluated using the value of
//		"rho", "e" and "U" at the boundary field, using the same
//		formula to evaluate "rhoE" in the internal field.
        
        	rhoE.boundaryFieldRef() = rho.boundaryField() * ( e.boundaryField() + 
				0.5 * magSqr ( U.boundaryField() ) ) ;

        
		p.ref() = rho() / psi() ;

        	p.correctBoundaryConditions() ;

	        rho.boundaryFieldRef() = psi.boundaryField() * p.boundaryField() ; // psi=1/(R*T)


//P
  		correctTime.stop();
//P


		//if ( Tbulk_target > 0. ) 
		//{
		
		//#include "tbforce.H"
		
		//}


        } 
	
//      End of RK time integration
	

//	"turbulence->correct()" updates the turbulence model 
//	based upon the now solved for U. What it does 
//	depends on whether you have a RANS or LES 
//	turbulence model enabled.


//P
  	turbcorrTime.start();
//P


        turbulence -> correct() ; 
  

//P
  	turbcorrTime.stop();
//P


	runTime.write() ;
   
	
	#include "diagnostics.H" // print tke and enstrophy on diagnostics.dat


//P	
	cycleTime.intervall();

	cycleTime.stop();
//P


        #include "cycletime.H"   // To printe the number of iteration and time
				 // required for each iteraration and each
				 // deltaT used

        #include "step.H"        // Evaluate Courant Number

	#include "setDeltaT.H"   // Adjust time step and it will be used at next iteration


    }



//  End of the main

    runTime.write() ;

//P
    totalTime.stop();
//P


    Info << "Start Timing = " << runTime.clockTimeIncrement() << " s" << nl << endl ;

    Info << "\nTotal time for the simulation (s) = " << totalTime.getTotalTime() << nl <<  endl ;

    Info << "\nCycle time (s) = " << cycleTime.getTotalTime() << nl << endl ;

    Info << "\nTurbulence validation time (s) = " << turbvalTime.getTotalTime() << nl << endl ;
    
    Info << "\nInitialize old variable time (s) = " << oldvarTime.getTotalTime() << nl << endl ;
    
    Info << "\nFlag checking time (s) = " << flagIfTime.getTotalTime() << nl << endl ;

    Info << "\nSound speed calculation time (s) = " << soundspeedTime.getTotalTime() << nl << endl ;

    Info << "\nRho interpolation time (s) = " << rhoaveTime.getTotalTime() << nl << endl ;

    Info << "\nU interpolation time (s) = " << UaveTime.getTotalTime() << nl << endl ;

    Info << "\nphi calculation time (s) = " << phiTime.getTotalTime() << nl << endl ;

    Info << "\nphit calculation time (s) = " << phitTime.getTotalTime() << nl << endl ;

    Info << "\nH calculation time (s) = " << HTime.getTotalTime() << nl << endl ;

    Info << "\nHave calculation time (s) = " << HaveTime.getTotalTime() << nl << endl ;

    Info << "\npave calculation time (s) = " << paveTime.getTotalTime() << nl << endl ;

    Info << "\nsensor.H inclusion time (s) = " << sensorTime.getTotalTime() << nl << endl ;

    Info << "\nAUSM.H inclusion time (s) = " << AUSMPTime.getTotalTime() << nl << endl ;

    Info << "\nvecDivU definition time (s) = " << vecDivUTime.getTotalTime() << nl << endl ;

    Info << "\nmuEff calculation time (s) = " << muEffTime.getTotalTime() << nl << endl ;

    Info << "\ntuaMC calculation time (s) = " << tauMCTime.getTotalTime() << nl << endl ;

    Info << "\nk definition time (s) = " << kTime.getTotalTime() << nl << endl ;

    Info << "\nkave definition time (s) = " << kaveTime.getTotalTime() << nl << endl ;

    Info << "\nmomVisFlux calculation time (s) = " << momVisFluxTime.getTotalTime() << nl << endl ;

    Info << "\nheatFlux calculation time (s) = " << heatFluxTime.getTotalTime() << nl << endl ;

    Info << "\nvisWork calculation (s) = " << visWorkTime.getTotalTime() << nl << endl ;

    Info << "\nenVisFlux calculation time (s) = " << enVisFluxTime.getTotalTime() << nl << endl ;

    Info << "\nTotal fluxes calculation time (s) = " << TotFluTime.getTotalTime() << nl << endl ;

    Info << "\nAUSM_conv.H activation time (s) = " << AUSMCTime.getTotalTime() << nl << endl ;

    Info << "\nDivergence of total fluxes calculation time (s) = " << TotDivFluTime.getTotalTime() << nl << endl ;

    Info << "\nRK sub step time (s) = " << RKSubTime.getTotalTime() << nl << endl ;

    Info << "\nField correction time (s) = " << correctTime.getTotalTime() << nl << endl ;

    Info << "\nTurbulence correction time (s) = " << turbcorrTime.getTotalTime() << nl << endl ;



    Info << "End\n" << endl ;

    return 0 ;

}
// ************************************************************************* //
