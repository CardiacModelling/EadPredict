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

/*

   Copyright (c) 2005-2015, University of Oxford.
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

#ifndef TESTOUTPUTTRACES_HPP_
#define TESTOUTPUTTRACES_HPP_

#include <cxxtest/TestSuite.h>


#include "DrugDataReader.hpp"
#include "ColumnDataReader.hpp"
#include "Trace.hpp"
#include <sstream>
#include <iomanip>
#include "AbstractCardiacCell.hpp"
#include "AbstractCvodeCell.hpp"
#include "CellProperties.hpp"
#include "OutputFileHandler.hpp"
#include "ohara_rudy_2011_endoCvode.hpp"
#include "AbstractParameterisedSystem.hpp"
#include "AbstractCvodeSystem.hpp"
#include "SteadyStateRunner.hpp"
#include "Debug.hpp"
#include "ColumnDataWriter.hpp"

class TestOutputTraces : public CxxTest::TestSuite
{
	public:
		void TestRunSimulation() throw(Exception)
		{
			/* Make sure this code is only run if CVODE is installed/enabled on the computer */
#ifdef CHASTE_CVODE

			// check that command line arguments have been supplied
			CommandLineArguments* p_args = CommandLineArguments::Instance();
			unsigned argc = *(p_args->p_argc); // has the number of arguments.
			std::cout << (argc-1)/2 << " arguments supplied.\n" << std::flush;
			if (argc == 1)
			{
				std::cerr << "Please input arguments\n"
					"--drug [0-40]\n"
					"--intervention\n"
					"	1: I_CaL increase\n"
					"--value\n"
					"--multiplier\n";
				return;
			}

			std::string file_location = getenv("CHASTE_TEST_OUTPUT");
			//Take command line arguments
			unsigned drug_index = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("--drug");
			unsigned intervention = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("--intervention");
			double multiplier = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("--multiplier");
			double value = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("--value");
			// set up simulation
			double dt = 0.1;
			double start_time = 0;
			double end_time = 3099;
			double stim_magnitude = -25.5;
			double stim_duration = 3;
			double stim_start = 100;
			double stim_period = 3000;
			boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(stim_magnitude,stim_duration,stim_period,stim_start));
			boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
			boost::shared_ptr<AbstractCvodeCell> p_model(new Cellohara_rudy_2011_endoFromCellMLCvode(p_solver, p_stimulus));

			// take in drug data
			DrugDataReader drug_data("projects/EadPredict/test/curated_dataset.dat");
			// make sure there are enough drugs
			TS_ASSERT_DIFFERS(drug_data.GetNumDrugs(),0u);
			// write down all the channel names
			std::vector<std::string> ap_predict_channels;
			ap_predict_channels.push_back("membrane_fast_sodium_current_conductance");
			ap_predict_channels.push_back("membrane_L_type_calcium_current_conductance");
			ap_predict_channels.push_back("membrane_rapid_delayed_rectifier_potassium_current_conductance");
			ap_predict_channels.push_back("membrane_slow_delayed_rectifier_potassium_current_conductance");
			ap_predict_channels.push_back("membrane_persistent_sodium_current_conductance");
			ap_predict_channels.push_back("membrane_transient_outward_current_conductance");
			ap_predict_channels.push_back("membrane_inward_rectifier_potassium_current_conductance");


			// get original values for each channel
			std::vector<double> original_values;
			for (unsigned int i= 0; i<ap_predict_channels.size(); i++)
			{
				original_values.push_back(p_model->GetParameter(ap_predict_channels[i]));
			}
			// make sure there are enough drugs
			TS_ASSERT_DIFFERS(drug_data.GetNumDrugs(),0u);

			// get out the name of the drug
			std::string drug_name = drug_data.GetDrugName(drug_index) ;
			TS_ASSERT_DIFFERS(drug_name.length(), 0u);
			std::cout << "\n" << drug_name << "\n";
			double original_drug_conc = drug_data.GetClinicalDoseRange(drug_index,1u);
			// loop through 26 concentrations
			double drug_conc = original_drug_conc * multiplier;

			for (unsigned channel_idx = 0; channel_idx < ap_predict_channels.size(); channel_idx++)
			{
				double conductance_factor = DrugDataReader::CalculateConductanceFactor(drug_conc,drug_data.GetIC50Value(drug_index,channel_idx));
				double ic50 = drug_data.GetIC50Value(drug_index,channel_idx);
				p_model->SetParameter(ap_predict_channels[channel_idx],original_values[channel_idx]*conductance_factor);

				double new_value = p_model->GetParameter(ap_predict_channels[channel_idx]);
				double new_multiplier = new_value / original_values[channel_idx];
				std::cout << ap_predict_channels[channel_idx] << " (IC50 = " << ic50 << ")";
				std::cout << " is set to ";
				std::cout << new_value << ", which is ";
				std::cout << new_multiplier << "x the original value\n";
								
			}



				
			// run to steady state
			SteadyStateRunner steady_runner(p_model);
			bool result = steady_runner.RunToSteadyState();
			std::cout << "Running to steady state: " << result << "\n";
			// do the intervention

			if (intervention==1)
			{
				double starting_ical = p_model->GetParameter(ap_predict_channels[1]);
				std::cout << "Setting " << ap_predict_channels[1] << " to " << value << "x original\n";
				p_model->SetParameter(ap_predict_channels[1],starting_ical*value);
			}
			OdeSolution solution = p_model->Compute(start_time,end_time,dt);

			std::vector<double> late_sodium_trace, calcium_trace, herg_trace, ito_trace, ik1_trace, fast_sodium_trace, iks_trace;
			std::vector<double> voltage_trace;
			
			solution.CalculateDerivedQuantitiesAndParameters(p_model.get());
			voltage_trace = solution.GetAnyVariable("membrane_voltage");
			calcium_trace = solution.GetAnyVariable("membrane_L_type_calcium_current");
			late_sodium_trace = solution.GetAnyVariable("membrane_persistent_sodium_current");
			herg_trace = solution.GetAnyVariable("membrane_rapid_delayed_rectifier_potassium_current");
			ito_trace = solution.GetAnyVariable("membrane_transient_outward_current");
			ik1_trace = solution.GetAnyVariable("membrane_inward_rectifier_potassium_current");
			fast_sodium_trace = solution.GetAnyVariable("membrane_fast_sodium_current");
			iks_trace = solution.GetAnyVariable("membrane_slow_delayed_rectifier_potassium_current_conductance");
			
			std::ofstream trace_file;

			std::stringstream file_name_stream;
			file_name_stream << file_location << "Tox_Res_Paper/" << intervention << "_" << "drug_effect_traces_" << drug_name << "_" << multiplier << "_" << value;
			std::string file_name = file_name_stream.str();
			trace_file.open(file_name.c_str());
			trace_file << "voltage\t late_sodium\t calcium\t herg\t ito\t ik1\t fast_sodium\t iks\n";
			trace_file << std::setprecision(20);
			for (unsigned int d=0; d < late_sodium_trace.size(); d++)
			{
				trace_file << voltage_trace[d] << "\t" << late_sodium_trace[d]<< "\t" << calcium_trace[d]<< "\t";
				trace_file << herg_trace[d] << "\t" << ito_trace[d] << "\t" << ik1_trace[d] << "\t";
				trace_file << fast_sodium_trace[d] << "\t" << iks_trace[d] << "\n";
			}
			trace_file.close();



#else
			/* CVODE is not enable or installed*/
			std::cout << "Cvode is not enabled.\n";
#endif
		}

};
#endif
