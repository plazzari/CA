#! /bin/bash

myPCL1=/Users/plazzari/Documents/workspace/PCL/GET_POLYTOPES/build

## declare an array  of CA test
declare -a CATEST=("gold3-4-19" "goldlarg4-24" "fort.12-12" "fort.12-18"  "fort.12-24" "fort.12-25" )  
#declare -a CATEST=("goldlarg4-24")  

######################################################################## 
# Estimate the polytope feature TESTS
cd $myPCL1

rm *_b.out*
for ca in "${CATEST[@]}"
do
#  python ./create_border_cloud_funcXYZ_CA4.py -i ${ca}.xyz -o ${ca}_b.pcd
   ./detect_polytopes ${ca}_b.pcd ${ca}_b.pcd.inner  ${ca}_b.out
   python ./write_sheet.py -i ${ca}
done      

echo EOB
