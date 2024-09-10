# This script generates the GITHASH file which can be used to put the githash in the log
# by defining -DSETGITHASH using the cpp preprocesor
# This is not done by defualt in the makefile but is using FPM so you should make sure you run this

git describe --always --dirty > GITHASH.txt
sed '1s/./#define GITHASH "&/' GITHASH.txt > GITHASH2.txt
mv GITHASH2.txt GITHASH.txt
sed '1!b;s/$/\"/g' GITHASH.txt > GITHASH2.txt
mv GITHASH2.txt GITHASH.txt
