#!/bin/bach

phn=$1
data=$2
iter=$3

cp models/phase$phn phases
cp $data unweighted

expect -c "
set timeout 20000
spawn docker run -it --rm -v $PWD:/app -w /app empht
expect \"Select 1 or 2:\"
send \"1\n\"
expect \"Select 1 or 2:\"
send \"2\n\"
expect \"Select (1-4):\"
send \"3\n\"
expect \"Number of phases of the PH-distribution to be fitted, (p):\"
send \"${phn}\n\"
expect \"Select(1-8):\"
send \"7\n\"
expect \"Number of EM-iterations:\"
send \"${iter}\n\"
expect \"Select 1 or 2\"
send \"1\n\"
expect \"$\"
exit 0
"
