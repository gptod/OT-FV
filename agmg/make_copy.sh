# program for creating multiple instance of agmg
# to preprocess and run different solver associate to different matrices.
# 
# See user guide.
#
# Use the following to create multiple copy
#
# for i in `seq 1 128`; do ./make_copy.sh ${i}; done 
#



mark=$1
echo ${mark}

cp agmg.m agmg${mark}.m
sed -i "s|agmg|agmg${mark}|g" agmg${mark}.m
cp dmtlagmg.mexa64 dmtlagmg${mark}.mexa64 
cp zmtlagmg.mexa64 zmtlagmg${mark}.mexa64
cp dmtlagmg.mexmaci64 dmtlagmg${mark}.mexmaci64
cp zmtlagmg.mexmaci64 zmtlagmg${mark}.mexmaci64
