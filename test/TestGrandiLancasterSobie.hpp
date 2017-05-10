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

#ifndef TESTLANCASTERSOBIE_HPP_
#define TESTLANCASTERSOBIE_HPP_

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
#include <sys/stat.h>
#include "SteadyStateRunner.hpp"
#include "ohara_rudy_2011_endoCvode.hpp"
#include "grandi_pasqualini_bers_2010_ssCvode.hpp"
#include "projects/ApPredict/src/single_cell/SingleActionPotentialPrediction.hpp"

#include "Debug.hpp"

class TestLancasterSobie : public CxxTest::TestSuite
{
	public:
		void TestGetAPD50AndCai()
		{
			/* Make sure this code is only run if CVODE is installed/enabled on the computer */
#ifdef CHASTE_CVODE

			// take in drug data
			DrugDataReader drug_data("projects/BethM/test/Tox_Res_Paper/curated_dataset.dat");
			// make sure there are enough drugs
			TS_ASSERT_DIFFERS(drug_data.GetNumDrugs(),0);

			// for running simulations
			double start_time = 0.0;
			double end_time = 1000.0;
			double dt = 0.1;

			// set up a model so I can steal the default stimulus
			boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
			boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
			boost::shared_ptr<AbstractCvodeCell> p_model(new Cellgrandi_pasqualini_bers_2010_ssFromCellMLCvode(p_solver, p_stimulus));
			p_stimulus = p_model->UseCellMLDefaultStimulus();
			//output file for all drug results
			std::ofstream apd50_file;
			std::string file_location = getenv("CHASTE_TEST_OUTPUT");
			std::string apd50_filename = file_location + "Tox_Res_Paper/curated_dataset_grandi_lancaster_sobie.dat";
			apd50_file.open(apd50_filename.c_str());
			// check that the file exists
			//TS_ASSERT_EQUALS(stat(apd50_filename.c_str(), &sb),0);

			// loop through all the drugs
			for (unsigned int drug_index=0; drug_index < drug_data.GetNumDrugs(); drug_index++)
			{
				// get out the name of the drug
				std::string drug_name = drug_data.GetDrugName(drug_index) ;
				TS_ASSERT_DIFFERS(drug_name.length(), 0);

				double drug_conc = drug_data.GetClinicalDoseRange(drug_index,1u);
				// make a new model
				boost::shared_ptr<AbstractCvodeCell> p_model(new Cellgrandi_pasqualini_bers_2010_ssFromCellMLCvode(p_solver, p_stimulus));
				apd50_file << drug_name;

				std::vector<double> conductance_factors;

				for (int b=0; b<7; b++)
				{
					conductance_factors.push_back(DrugDataReader::CalculateConductanceFactor(drug_conc,drug_data.GetIC50Value(drug_index,b)));
				}

				// take conductance factors and apply them to the model
				std::vector<double> original_values;
				std::vector<std::string> block_channel_names;
				block_channel_names.push_back("membrane_fast_sodium_current_conductance");
				block_channel_names.push_back("membrane_L_type_calcium_current_conductance");
				block_channel_names.push_back("membrane_rapid_delayed_rectifier_potassium_current_conductance");
				block_channel_names.push_back("membrane_slow_delayed_rectifier_potassium_current_conductance");
				block_channel_names.push_back("membrane_persistent_sodium_current_conductance");
				block_channel_names.push_back("membrane_transient_outward_current_conductance");
				block_channel_names.push_back("membrane_inward_rectifier_potassium_current_conductance");
				for (unsigned int i= 0; i<block_channel_names.size(); i++)
				{
					try
					{
						original_values.push_back(p_model->GetParameter(block_channel_names[i]));
					}
					catch(Exception &e)
					{
						original_values.push_back(-2);
						std::cout << "No " << block_channel_names[i] << "\n";
					}
				}

				for (unsigned int j= 0; j<block_channel_names.size(); j++)
				{
					try
					{
						p_model->SetParameter(block_channel_names[j],original_values[j]*conductance_factors[j]);
					}
					catch(Exception &e)
					{
						std::cout << "\n";
					}
				}

				// run to steady state
				SteadyStateRunner steady_runner(p_model);
				bool result = steady_runner.RunToSteadyState();
				std::cout << "Running to steady state: " << result << "\n";
				SingleActionPotentialPrediction ap_runner(p_model);
				ap_runner.SetMaxNumPaces(1000u);
				ap_runner.RunSteadyPacingExperiment();
				double apd50;
				if (ap_runner.DidErrorOccur())
				{
					std::cout << "An error occurred\n";
				}
				else
				{
					apd50 = ap_runner.GetApd50();
				}

				apd50_file << "\t" << apd50 << "\t";
				// now run once to get out the calcium concentration
				OdeSolution solution = p_model->Compute(start_time,end_time,dt);
				std::vector<double> calcium_conc = solution.GetAnyVariable("cytosolic_calcium_concentration");

				apd50_file << "\t" << calcium_conc[calcium_conc.size()-10] << "\n";


			}
			apd50_file.close();
#else
			/* CVODE is not enable or installed*/
			std::cout << "Cvode is not enabled.\n";
#endif

		}

};
#endif
