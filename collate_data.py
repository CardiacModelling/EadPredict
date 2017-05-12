#! /usr/bin/python

# get risk categories and sort alphabetically
rf = open("redferns.tsv",'r')
redferns = rf.readlines()
redferns.sort()

a1f = open("curated_dataset_apd90.dat")
apd90 = a1f.readlines()
apd90.sort()

glf = open("curated_dataset_grandi_lancaster_sobie.dat")
grandils = glf.readlines()
grandils.sort()

olf = open("curated_dataset_ohara_lancaster_sobie.dat")
oharals = olf.readlines()
oharals.sort()

inaf = open("DruggedSteadyStateThresholds/CollatedINa.dat")
ina = inaf.readlines()
ina.sort()

icalf = open("DruggedSteadyStateThresholds/ohara_rudy_2011membrane_L_type_calcium_current_conductance_0")
ical = icalf.readlines()
ical.sort()

ikrf = open("DruggedSteadyStateThresholds/ohara_rudy_2011membrane_rapid_delayed_rectifier_potassium_current_conductance_0")
ikr = ikrf.readlines()
ikr.sort()

hergf = open("herg_block.tsv")
herg = hergf.readlines()
herg.sort()

output_file = open("collated_data.tsv",'w')

for i in range(0,len(redferns)):
	output_file.write(redferns[i].split()[0]+ "\t"+ redferns[i].split()[1]+ "\t")
	output_file.write(apd90[i].split()[1]+ "\t"+ grandils[i].split()[1]+ "\t"+ grandils[i].split()[2]+ "\t"+  oharals[i].split()[1]+ "\t")
	output_file.write(oharals[i].split()[2]+ "\t"+ ina[i].split()[1]+ "\t"+ ical[i].split()[1]+ "\t"+ ikr[i].split()[1]+ "\t"+ herg[i].split()[1])

