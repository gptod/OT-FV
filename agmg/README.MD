# Copy in this directory the files from the AGMG code.

agmg.m
dmtlagmg.mexa64
zmtlagmg.mexa64
dmtlagmg.mexmaci64
zmtlagmg.mexmaci64


# We include here a bash script "make_copy.sh" for generating copy of
# agmg solver when diferent matrices are solved at same time (see AGMG
# user guide for more details).
# 
# Example of use the following to create multiple copy

for i in `seq 1 129`; do ./make_copy.sh ${i}; done 

