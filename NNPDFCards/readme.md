To run the example on this folder, first install the code and then do:

```bash
../build/Vrap DYE605.dat DYE605_kinematics.dat
```

`vrap` will run and will produce a file `test.pineappl.lz4` containing the predictions per bin for the given kinematics:

```bash
pineappl convolute test.pineappl.lz4 NNPDF40_nnlo_as_01180
```

```
b       x1           x2        diff     scale uncertainty
        []           []         []             [%]       
-+-------+-------+----+----+-----------+--------+--------
0 8.50496 8.50496 -0.2 -0.2 1.6039350e2    -7.80     7.97
1  10.507  10.507 -0.1 -0.1 3.9390718e1    -8.95     8.77
2 8.50496 8.50496    0    0 1.7727698e2    -7.91     8.11
3  10.507  10.507  0.1  0.1 4.3949162e1    -9.09     8.89
4 8.50496 8.50496  0.2  0.2 1.8652510e2    -7.77     7.90
5  10.507  10.507  0.3  0.3 4.7241126e1    -9.13     8.78
6 8.90848 8.90848  0.4  0.4 1.3878889e2    -7.94     7.84 
```
