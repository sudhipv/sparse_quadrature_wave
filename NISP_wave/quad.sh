dirname=d3l5

fstype=sparse

dim=3

level=5


echo "dir name is $dirname"

[ -d ./quadData/$dirname ] || mkdir ./quadData/$dirname

cd /Users/sudhipv/documents/UQTk_v3.0.4/UQTk-install/bin


./generate_quad -d $dim -p $level -g HG -x $fstype | tee /Users/sudhipv/documents/python_ni/NISP_wave/quadData/$dirname/output.txt



cp ./qdpts.dat /Users/sudhipv/documents/python_ni/NISP_wave/quadData/$dirname/
cp ./wghts.dat /Users/sudhipv/documents/python_ni/NISP_wave/quadData/$dirname/
cd -



















