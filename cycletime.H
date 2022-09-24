	if ( Pstream :: master() )
	{
	
		std :: ofstream cycle ;

		cycle.open ( "cycletime.dat" , std :: ofstream :: out | std :: ofstream :: app ) ;

		cycle << Iter << " " << cycleTime.getIntervallTime() << " " << 
			runTime.deltaTValue() << std :: endl ;

		cycle.close() ;

	}	


