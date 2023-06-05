#! C:/Windows/System32/bash.exe

cd C:/Users/emmid/Desktop/Research/Merlin/merlin-1.1.2

#make sure path to merlin is in the path environment variable

echo $PATH > C:/Users/emmid/Desktop/Research/Merlin/path.txt

merlin -m C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/mc1r.map -d C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/mc1r_dorsum.dat -p C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/Ped_Files/ped_dorsum_mc1r.txt -xx --assoc --bits:30 > C:/Users/emmid/Desktop/Research/Merlin/mc1r_dorsum.txt
merlin -m C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/mc1r.map -d C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/mc1r_limb.dat -p C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/Ped_Files/ped_hindlimb_mc1r.txt -xx --assoc --bits:30 > C:/Users/emmid/Desktop/Research/Merlin/mc1r_hindlimb.txt

merlin -m C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/asip.map -d C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/asip_dorsum.dat -p C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/Ped_Files/ped_dorsum_asip.txt -xx --assoc --bits:30 > C:/Users/emmid/Desktop/Research/Merlin/asip_dorsum.txt
merlin -m C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/asip.map -d C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/asip_limb.dat -p C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/Ped_Files/ped_hindlimb_asip.txt -xx --assoc --bits:30 > C:/Users/emmid/Desktop/Research/Merlin/asip_hindlimb.txt

merlin -m C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/retsat.map -d C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/retsat_dorsum.dat -p C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/Ped_Files/ped_dorsum_retsat.txt -xx --assoc --bits:30 > C:/Users/emmid/Desktop/Research/Merlin/retsat_dorsum.txt
merlin -m C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/retsat.map -d C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/retsat_limb.dat -p C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/Ped_Files/ped_hindlimb_retsat.txt -xx --assoc --bits:30 > C:/Users/emmid/Desktop/Research/Merlin/retsat_hindlimb.txt

merlin -m C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/bsn2.map -d C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/bsn2_dorsum.dat -p C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/Ped_Files/ped_dorsum_bsn2.txt -xx --assoc --bits:30 > C:/Users/emmid/Desktop/Research/Merlin/bsn2_dorsum.txt
merlin -m C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/bsn2.map -d C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/bsn2_limb.dat -p C:/Users/emmid/Desktop/Spec_Analysis/Spec_Results/Ped_Files/ped_hindlimb_bsn2.txt -xx --assoc --bits:30 > C:/Users/emmid/Desktop/Research/Merlin/bsn2_hindlimb.txt
