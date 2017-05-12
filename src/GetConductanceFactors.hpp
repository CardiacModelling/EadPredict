/*

Copyright (c) 2005-2017, University of Oxford.
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

#ifndef GETCONDUCTANCEFACTORS_HPP_
#define GETCONDUCTANCEFACTORS_HPP_

class DrugDataReader
{
	private:

	public:
		DrugDataReader(int drug_number)
		{
			std::ifstream drug_file;
			drug_file.open("projects/BethM/test/Drugs/drug_data_file.dat");
			if(!drug_file)
			{ // file couldn't be opened
				EXCEPTION("Couldn't open data file");
			}
		}
		double GetIC50Value(unsigned drugIndex, unsigned channelIndex)
		{
			assert(drugIndex < GetNumDrugs());
			assert(channelIndex<4u);
			double ic50 = mIc50values[drugIndex](channelIndex);

			/**
			 * If an IC50 value takes the value
			 *  -1 in the data file that means "affect unknown".
			 *  -2 in the data file means "known to have no affect".
			 */
			if (fabs(ic50 + 2) < 1e-9)
			{   // This drug is known to have no effect on this channel
				// We set the IC50 value to DBL_MAX for the block calculations.
				ic50 = DBL_MAX;
			}

			return ic50;
		}

		double GetClinicalDoseRange(unsigned drugIndex, unsigned lowOrHigh)
		{
			assert(lowOrHigh==0 || lowOrHigh==1);
			assert(drugIndex < GetNumDrugs());
			if (mClinicalDoseRange[drugIndex](0) < 0)
			{
				EXCEPTION("No data available on clinical dose for this drug.");
			}
			return mClinicalDoseRange[drugIndex](lowOrHigh);
		}

		void LoadDrugDataFromFile(std::string fileName)
		{

			std::ifstream indata; // indata is like cin
			indata.open(fileName.c_str()); // opens the file
			if(!indata)
			{ // file couldn't be opened
				EXCEPTION("Couldn't open data file: " + fileName);
			}

			bool line_is_not_blank = true;
			unsigned line_counter = 0;
			while (line_is_not_blank)
			{
				std::string this_line;
				getline(indata, this_line);
				std::stringstream line(this_line);

				unsigned category;
				std::string name;
				std::string in_redfern_figs;
				c_vector<double, 4> ic50s;
				c_vector<double, 2> doses;
				double grandi_measure;
				ic50s.clear();
				doses.clear();

				line >> name;
				line >> category;
				line >> ic50s(0);
				line >> ic50s(1);
				line >> ic50s(2);
				line >> ic50s(3);
				line >> doses(0);
				line >> doses(1);
				//line >> in_redfern_figs;
				//line >> grandi_measure;

				mDrugNames.push_back(name);
				mRedfernCategory.push_back(category);
				mIc50values.push_back(ic50s);
				mClinicalDoseRange.push_back(doses);
				mGrandiMeasure.push_back(grandi_measure);

				line_counter++;
				//std::cout << "line_counter = " <<line_counter << "\n" << std::flush;
				if (indata.eof())
				{
					break;
				}
			}

		}

		unsigned GetNumDrugs(void)
		{
			assert(mDrugNames.size()==mRedfernCategory.size());
			assert(mDrugNames.size()==mIc50values.size());
			assert(mDrugNames.size()==mClinicalDoseRange.size());
			return mDrugNames.size();
		}

		/**
		 * Calculate the probability of a channel being open given this drug, IC50 and hill coefficient.
		 *
		 * Note: A negative IC50 value is interpreted as "drug has no effect on this channel".
		 *
		 * @param rConc  concentration of the drug.
		 * @param rIC50  IC50 value for this drug and channel
		 * @param hill  Hill coefficient for this drug dependent inactivation curve (defaults to 1).
		 *
		 * @return proportion of channels which are still active at this drug concentration
		 */
		static double CalculateConductanceFactor(const double& rConc, const double& rIC50, double hill = 1.0)
		{
			if (rIC50 < 0)
			{
				return 1.0;
			}
			else
			{
				return 1.0/(1.0 + pow((rConc/rIC50), hill));
			}
		}
};

#endif // DRUGDATAREADER_HPP_
