/*

Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TESTCREATETRACES_HPP_
#define TESTCREATETRACES_HPP_

#include <cxxtest/TestSuite.h>
#include "DetectAfterDepolarisations.hpp"
#include "ThresholdIntervention.hpp"
#include "ColumnDataReader.hpp"
#include "DrugDataReader.hpp"
#include <sstream>
#include <iomanip>
#include "AbstractCardiacCell.hpp"
#include "AbstractCvodeCell.hpp"
#include "CellProperties.hpp"
#include "OutputFileHandler.hpp"
#include "shannon_wang_puglisi_weber_bers_2004Cvode.hpp"
#include "fink_noble_giles_model_2008Cvode.hpp"
#include "ten_tusscher_model_2006_epiCvode.hpp"
#include "ten_tusscher_model_2006_endoCvode.hpp"
#include "ten_tusscher_model_2006_MCvode.hpp"
#include "iyer_2004Cvode.hpp"
#include "noble_model_1998Cvode.hpp"
#include "ohara_rudy_2011_endoCvode.hpp"
#include "AbstractParameterisedSystem.hpp"
#include "AbstractCvodeSystem.hpp"
#include <boost/lexical_cast.hpp>

#include "SteadyStateRunner.hpp"

#include "Debug.hpp"

class TestCreateTraces : public CxxTest::TestSuite
{
public:
	void TestMakeTraces() throw(Exception)
	{
		/* Make sure this code is only run if CVODE is installed/enabled on the computer */
#ifdef CHASTE_CVODE

		// check that command line arguments have been supplied
		CommandLineArguments* p_args = CommandLineArguments::Instance();
		unsigned argc = *(p_args->p_argc); // has the number of arguments.
		std::cout << (argc-1)/2 << " arguments supplied.\n" << std::flush;
		if (argc == 1)
		{
			std::cerr << "Please input an argument\n"
					"--model\n"
					"	1 = shannon_wang_puglisi_weber_bers_2004\n"
					"	2 = fink_noble_giles_model_2008\n"
					"	3 = ten_tusscher_model_2006_epi\n"
					"	4 = ten_tusscher_model_2006_endo\n"
					"	5 = ten_tusscher_model_2006_M\n"
					"	6 = ohara_rudy_2011_endo\n"
					"--intervention\n"
					"	1 = membrane_L_type_calcium_current_conductance\n"
					"	2 = membrane_fast_sodium_current_shift_inactivation\n"
					"	3 = membrane_rapid_delayed_rectifier_potassium_current_conductance\n"
					"--value <intervention level>\n";
			return;
		}


		//		CommandLineArguments* p_args = CommandLineArguments::Instance();
		//		unsigned argc = *(p_args->p_argc); // has the number of arguments.
		char* model = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--model");
		unsigned model_index = atoi(model);
		char* inter = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--intervention");
		unsigned inter_index = atoi(inter);
		double value = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("--value");

		// for running simulations
		double start_time = 0.0;
		double end_time = 10000.0;
		double dt = 0.1;

		// stimulus and solver for all models
		double stim_magnitude = -25.5;
		double stim_duration = 3;
		double stim_start = 100;
		double stim_period = 3000;
		boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(stim_magnitude,stim_duration,stim_period,stim_start));
		boost::shared_ptr<AbstractIvpOdeSolver> p_solver;

		// make new model
		boost::shared_ptr<AbstractCvodeCell> p_model;
		if (model_index ==1)
		{
			boost::shared_ptr<AbstractCvodeCell> p_new_model(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLCvode(p_solver, p_stimulus));
			p_model = p_new_model;
		}
		else if (model_index == 2)
		{
			boost::shared_ptr<AbstractCvodeCell> p_new_model(new Cellfink_noble_giles_model_2008FromCellMLCvode(p_solver, p_stimulus));
			p_model = p_new_model;
		}
		else if (model_index == 3)
		{
			boost::shared_ptr<AbstractCvodeCell> p_new_model(new Cellten_tusscher_model_2006_epiFromCellMLCvode(p_solver, p_stimulus));
			p_model = p_new_model;
		}
		else if (model_index == 4)
		{
			boost::shared_ptr<AbstractCvodeCell> p_new_model(new Cellten_tusscher_model_2006_endoFromCellMLCvode(p_solver, p_stimulus));
			p_model = p_new_model;
		}
		else if (model_index == 5)
		{
			boost::shared_ptr<AbstractCvodeCell> p_new_model(new Cellten_tusscher_model_2006_MFromCellMLCvode(p_solver, p_stimulus));
			p_model = p_new_model;
		}
		else if (model_index == 6)
		{
			boost::shared_ptr<AbstractCvodeCell> p_new_model(new Cellohara_rudy_2011_endoFromCellMLCvode(p_solver, p_stimulus));
			p_model = p_new_model;
		}
		else
		{
			boost::shared_ptr<AbstractCvodeCell> p_new_model(new Cellten_tusscher_model_2006_MFromCellMLCvode(p_solver, p_stimulus));
			p_model = p_new_model;
		}

		std::string modelname = p_model->GetSystemName();

		std::string intervention;
		double starting_value;
		OdeSolution solution;
		if (inter_index ==1)
		{
			intervention = "membrane_L_type_calcium_current_conductance";
			starting_value = p_model->GetParameter(intervention);
		}
		if (inter_index == 2)
		{
			intervention = "membrane_fast_sodium_current_shift_inactivation";
			starting_value = 1;
		}
		if (inter_index == 3)
		{
			intervention = "membrane_rapid_delayed_rectifier_potassium_current_conductance";
			starting_value = p_model->GetParameter(intervention);
		}
		std::cout << intervention << " begins at: " << starting_value << "\n";

		// run to steady state
		SteadyStateRunner steady_runner(p_model);
		// set parameter as constant
		p_model->SetParameter(intervention,starting_value);
		bool result = steady_runner.RunToSteadyState();
		std::cout << "Running to steady state: " << result << "\n";

		p_model->SetParameter(intervention,starting_value*value);

		std::cout << "Set to " << starting_value*value << "\n";

		// // for nai
		//				p_model->SetStateVariable("cytosolic_sodium_concentration",starting_value);
		//				p_stimulus->SetPeriod(1000.0);
		//				solution = p_model->Compute(start_time,30000,dt);
		//				p_stimulus->SetPeriod(4000.0);
		//p_stimulus->SetPeriod(4000.0);

		std::string foldername = modelname + intervention + boost::lexical_cast<std::string>(value);

		solution = p_model->Compute(start_time,end_time,dt);
		solution.WriteToFile(("SingleRuns/" + foldername),foldername,"ms");

#else
		/* CVODE is not enable or installed*/
		std::cout << "Cvode is not enabled.\n";
#endif
	}
};

#endif
