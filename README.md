I. Synopsis
------------
This program simulates a single-phase DC/AC inverter controlled by a Hybrid Predictive Control algorithm. The switching configuration of the H-bridge in the inverter cirucit is controlled so as to maximize the time in between consecutive switches by means of the prediction of trajectories. 
This code allows to simulate standard 3 levels as well as more complex single-phase multilevel inverters.

II. File List
------------
	
	run_inverter.m
	C_inverter.m
	D_inverter.m
	f_inverter.m
	g_inverter.m
	hybridsolver.m
	postprocessing.m
	example.eps
	f_reference.m
	fftTests.m
	miscPlots.m
	nfftAlgorithm.m
	odeJumpEvent.m
	patchell.m
	PlotFrame.m
	plotPredTrajectories.m
	VDCnoise.m
	VideoScript.m
	

III. To Run
------------
	Using Matlab 2015a or higher, run the script run_inverter.m. 
	