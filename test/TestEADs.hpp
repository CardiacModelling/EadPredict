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

#ifndef TESTEADS_HPP_
#define TESTEADS_HPP_

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
#include "CreateModel.hpp"
#include <sys/stat.h>

#include <boost/assign.hpp>

#include "SteadyStateRunner.hpp"

#include "Debug.hpp"

class TestEads : public CxxTest::TestSuite
{
	public:
		void TestCreateThresholdTables() throw(Exception)
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
					"	6 = ohara_rudy_2011\n"
					"	7 = grandi_pasqualini_bers_2010_ss\n"
					"--intervention\n"
					"	1 = cytosolic_sodium_concentration\n"
					"	2 = membrane_rapid_delayed_rectifier_potassium_current_conductance\n"
					"	3 = membrane_fast_sodium_current_shift_inactivation\n"
					"	4 = membrane_L_type_calcium_current_conductance\n"
					"	5 = membrane_L_type_calcium_current_conductance_double_extracellular_potassium_concentration\n"
					"	6 = membrane_persistent_sodium_current_conductance\n"
					"	7 = membrane_rapid_delayed_rectifier_potassium_current_conductance_double_extracellular_potassium_concentration\n"
					"	8 =  membrane_fast_sodium_current_shift_inactivation_double_extracellular_potassium_concentration\n"
					"	9 = membrane_L_type_calcium_current_conductance_three_quarters_extracellular_potassium_concentration\n"
					"	10 = membrane_rapid_delayed_rectifier_potassium_current_conductance_three_quarters_extracellular_potassium_concentration\n"
					"	11 = membrane_fast_sodium_current_shift_inactivation_three_quarters_extracellular_potassium_concentration\n"
					"	12 = membrane_fast_sodium_current_reduced_inactivation\n"
					"--sample \n"
					"--ll [lower limit]\n"
					"--hl [higher limit]\n"
					"--IP\n"
					"	0 = do not create intervention profile \n"
					"	1 = do create intervention profile \n";
				return;
			}

			//Take command line arguments
			unsigned model_index = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("--model");
			unsigned inter_index = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("--intervention");
			int conc_index = 2;
			int sample_index = CommandLineArguments::Instance()->GetIntCorrespondingToOption("--sample");
			bool intervention_profile = CommandLineArguments::Instance()->GetBoolCorrespondingToOption("--IP");

			// put together vector of intervention values for testing
			double low_limit = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("--ll");
			double high_limit = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("--hl");

			double step = (high_limit - low_limit)/99;
			std::vector<double> intervention_values;
			intervention_values.push_back(low_limit);
			for (int a=1; a<100; a++)
			{
				intervention_values.push_back(intervention_values[a-1]+step);
			}

			// check that the end values line up
			TS_ASSERT_DELTA(intervention_values[intervention_values.size()-1], high_limit, high_limit/100);


			// take in drug data
			DrugDataReader drug_data("projects/EadPredict/test/curated_dataset.dat");
			// make sure there are enough drugs
			TS_ASSERT_DIFFERS(drug_data.GetNumDrugs(),0u);
			
				std::vector<std::string> ap_predict_channels;
				ap_predict_channels.push_back("membrane_fast_sodium_current_conductance");
				ap_predict_channels.push_back("membrane_L_type_calcium_current_conductance");
				ap_predict_channels.push_back("membrane_rapid_delayed_rectifier_potassium_current_conductance");
				ap_predict_channels.push_back("membrane_slow_delayed_rectifier_potassium_current_conductance");
				ap_predict_channels.push_back("membrane_persistent_sodium_current_conductance");
				ap_predict_channels.push_back("membrane_transient_outward_current_conductance");
				ap_predict_channels.push_back("membrane_inward_rectifier_potassium_current_conductance");



			// for output, check that folders exist
			std::string file_location = getenv("CHASTE_TEST_OUTPUT");
			std::vector<std::string> pathname;
			std::string threshold_folder = "Tox_Res_Paper/DruggedSteadyStateThresholds";
			std::string ip_folder = "Tox_Res_Paper/ThreeColumnInterventionProfile";
			std::string trace_folder = "Tox_Res_Paper/ADCheck";
			pathname.push_back(file_location + threshold_folder);
			pathname.push_back(file_location + ip_folder);
			pathname.push_back(file_location + trace_folder);
			pathname.push_back(file_location + trace_folder + "/Bisection");
			pathname.push_back(file_location + trace_folder + "/Sequential");
			struct stat sb;
			for (unsigned e=0; e< pathname.size(); e++)
			{
				//and if they don't, create them
				if (stat(pathname[e].c_str(), &sb) != 0)
				{
					boost::filesystem::create_directories(pathname[e]);
				}
				TS_ASSERT(stat(pathname[e].c_str(), &sb) == 0);
			}


			// for running simulations
			double start_time = 0.0;
			double dt = 0.1;

			// stimulus info
			double stim_magnitude = -25.5;
			double stim_duration = 3;
			double stim_start = 100;
			double stim_period = 3000;
			double end_time = 4*stim_period;
			boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(stim_magnitude,stim_duration,stim_period,stim_start));

			//output file for all drug results
			std::ofstream fine_file;
			std::stringstream fine_filestream;
			std::string fine_filename;

			// make sure there are enough drugs
			TS_ASSERT_DIFFERS(drug_data.GetNumDrugs(),0u);

			// loop through all the drugs
			for (unsigned int drug_index=0; drug_index < drug_data.GetNumDrugs(); drug_index++)
			{
				// get out the name of the drug
				std::string drug_name = drug_data.GetDrugName(drug_index) ;
				TS_ASSERT_DIFFERS(drug_name.length(), 0u);

				double drug_conc = drug_data.GetClinicalDoseRange(drug_index,1u);
				std::cout << drug_name << "\n";
				// make new model and give model_creator some more info
				CreateModel model_creator(inter_index, model_index, p_stimulus, intervention_values, start_time, end_time, dt, drug_name, conc_index, trace_folder);
				// create the file for recording the threshold results
				if (drug_index == 0)
				{
					fine_filestream << file_location << threshold_folder << "/" << model_creator.GetIdentifier() << "_" << (sample_index);
					fine_filename = fine_filestream.str();
					fine_file.open(fine_filename.c_str());
				}
				// check that the file exists
				TS_ASSERT_EQUALS(stat(fine_filename.c_str(), &sb),0);

				fine_file << drug_name;




				// calculate the proportion of the different channels which are still active
				std::vector<double> conductance_factors;

				for (unsigned channel_idx = 0; channel_idx < ap_predict_channels.size(); channel_idx++)
				{
                                        conductance_factors.push_back(DrugDataReader::CalculateConductanceFactor(drug_conc,drug_data.GetIC50Value(drug_index,channel_idx)));
                                }
							// apply drug to model
							model_creator.Drug(conductance_factors, ap_predict_channels);
							// run to steady state
							model_creator.GetToSteadyState();

							//for storing threshold values
							std::vector<double> coarse_threshold;

							// Create an intervention profile if required
							if (intervention_profile)
							{

								// set up file for outputting results of different levels of intervention
								std::ofstream outputfile;
								std::stringstream filestream;
								std::string filename;
								filestream << file_location << ip_folder << "/" << drug_name << (conc_index) << model_creator.GetIdentifier();
								filename = filestream.str();
								outputfile.open(filename.c_str());
								std::cout << "\nSaving output to: " << filename << "\n";
								// check that the file exists
								TS_ASSERT_EQUALS(stat(filename.c_str(), &sb),0);
								bool found_threshold = false;
								bool expected_change = model_creator.GetExpectedChange();

								////loop over intervention levels to create profile and find coarse threshold
								for( int c=0; c < 100; c++)
								{
									model_creator.ResetAndCheck();
									// apply intervention and see if there's an AD
									int causes_afterdepolarisation = model_creator.Intervene(c, conc_index, drug_name);
									std::cout << c << " " << intervention_values[c] << " " << causes_afterdepolarisation << "\n";
									outputfile << c << " " << intervention_values[c] << " "<< causes_afterdepolarisation << "\n";


									if ((causes_afterdepolarisation == expected_change) && !found_threshold)
									{
										// if this is the first time an AD is detected or the first time the AD disappears
										coarse_threshold.push_back(intervention_values[c]);
										coarse_threshold.push_back(intervention_values[c-1]);
										found_threshold = true;
										std::cout << "\nThe coarse threshold is: " << coarse_threshold[0] << coarse_threshold[1] << "\n\n";
									}
								}
								if (!found_threshold)
								{
									coarse_threshold.push_back(low_limit);
									coarse_threshold.push_back(high_limit);
									std::cout << "Threshold is not between higher and lower limits\n";
								}

								outputfile.close();
							}
							else // or assume that the threshold is between the higher and lower limits
							{

								coarse_threshold.push_back(low_limit);
								coarse_threshold.push_back(high_limit);
							}


							// find the exact value of the threshold
							model_creator.SetLimits(coarse_threshold);
							bool found_fine_threshold = false;
							for (unsigned bisection_iterations = 0; bisection_iterations<100; bisection_iterations++)
							{
								model_creator.ResetAndCheck();

								// intervention_result could be 1 (AD) 0 (no AD) or -2 (parameter value out of range error)
								int intervention_result = model_creator.InterveneAndCheck();
								bool result;
								if (intervention_result != -2)
								{
									result = intervention_result;
								}

								// move to next interval to check
								found_fine_threshold = model_creator.IntervalBisection(result);

								if (found_fine_threshold)
								{
									std::cout << "The threshold is: " << model_creator.GetThreshold() << "\n";
									fine_file << " " << model_creator.GetThreshold() << "\n";
									break;
								}
							}
							if (!found_fine_threshold)
							{
								fine_file << " -2\n";
								std::cout << "Failed to find threshold.\n";
							}

			}
			fine_file.close();

		}
#else
		/* CVODE is not enable or installed*/
		std::cout << "Cvode is not enabled.\n";
#endif

};

#endif
