import sys
import subprocess

args = sys.argv
samplename=str(args[1])
foldidx=int(str(args[1]))
chanknum=10

i=foldidx
command = ["python","script/7_6DeepLift_full_split_single.py",samplename,str(1000),str(i)]
subprocess.call(command)