cc -O3 -fomit-frame-pointer   -c -o phylip.o phylip.c
cc -O3 -fomit-frame-pointer   -c -o consense.o consense.c
cc -O3 -fomit-frame-pointer   -c -o cons.o cons.c
cc -O3 -fomit-frame-pointer consense.o phylip.o cons.o -lm -o consense
