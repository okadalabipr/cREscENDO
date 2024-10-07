import sys
import subprocess

args = sys.argv
samplename=str(args[1])
chanknum=10

for i in range(chanknum):
  command = ["python","script/7_2DNAembed.py",samplename,str(i)]
  subprocess.call(command)
  command = ["python","script/7_3RNAprediction_t2_raw_2.py",samplename,str(i)]
  subprocess.call(command)
  command = ["python","script/7_4RNAprediction_load_t2_raw_2.py",samplename,str(i)]
  subprocess.call(command)
  command = ["python","script/7_5RNAprediction_load_predrna_t2_raw_2.py",samplename,str(i)]
  subprocess.call(command)
  command = ["python","script/7_6DeepLift_full_split.py",samplename,str(1000),str(i)]
  subprocess.call(command)