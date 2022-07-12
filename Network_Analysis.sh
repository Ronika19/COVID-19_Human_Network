cd HealthyVsICU
mkdir -p Output_Files
python3 DEG_NonDEG_Module_Pathway.py
cd ..
echo HealthyVsICU Done
cd HealthyVsModerate
mkdir -p Output_Files
python3 DEG_NonDEG_Module_Pathway.py
cd ..
echo HealthyVsModerate Done
cd HealthyVsSevere
mkdir -p Output_Files
python3 DEG_NonDEG_Module_Pathway.py
cd ..
echo HealthyVsSevere Done
cd ModerateVsSevere
mkdir -p Output_Files
python3 DEG_NonDEG_Module_Pathway.py
cd ..
echo ModerateVsSevere Done
cd SevereVsICU
mkdir -p Output_Files
python3 DEG_NonDEG_Module_Pathway.py
cd ..
echo SevereVsICU Done
