Some functions used in this repository are based on the [AGMG](http://agmg.eu/) solver, that has to be downloaded separately (free under academic license).
Copy in this directory the files from the AGMG code.

```
agmg.m
dmtlagmg.mexa64
zmtlagmg.mexa64
dmtlagmg.mexmaci64
zmtlagmg.mexmaci64
```

Following the AGMG users' guide, in order to avoid useless
preprocess time, we need to create several copies of the
AGMG source code, one for each linear system that has to be
solved separately. We include here a bash script
"make_copy.sh" for doing it: for example the following command 

```
for i in `seq 1 129`; do ./make_copy.sh ${i}; done
```

will create 129 copies. Note that the number of copies has to be at least equal to the (number of time-steps) + 1. 



