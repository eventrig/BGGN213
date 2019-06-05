Class 6: R Functions
================

Section 1: Writing Functions
============================

Today we will focus on **R functions** but we will start with a bit of **file reading**.

Example Plot
------------

``` r
plot(1:10, type="l", col="blue")
```

![](class06_rework_files/figure-markdown_github/unnamed-chunk-1-1.png)

Making a table for each of the three test files
-----------------------------------------------

We're keeping in mind what characters separate each variables in the columns and whether they have headers or not.

``` r
read.table("test1.txt", header=TRUE, sep=",") #commas separate the variables
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
read.table("test2.txt", header=TRUE, sep="$") #Dollar signs separate the variables
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
read.table("test3.txt") #the default sep= for read.table is tab, which means no input needed!
```

    ##   V1 V2 V3
    ## 1  1  6  a
    ## 2  2  7  b
    ## 3  3  8  c
    ## 4  4  9  d
    ## 5  5 10  e

You can even command R to draw data from a link to a file.

``` r
read.csv("https://bioboot.github.io/bggn213_S19/class-material/test1.txt")
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

Writing functions!
------------------

Here is the general outline/syntax of writing a function:

``` r
name.of.function <- function(argument1, argument2) {
 statements
 return(something)
}

square.it <- function(x) {  #making a function called square.it
 square <- x * x
 return(square)
}


# square a number
square.it(5)
```

    ## [1] 25

``` r
# square a vector
square.it(c(1, 4, 2))
```

    ## [1]  1 16  4

``` r
#even a matrix
matrix1 <- cbind(c(3, 10), c(4, 5))
square.it(matrix1)
```

    ##      [,1] [,2]
    ## [1,]    9   16
    ## [2,]  100   25

Here we make a custom function called "add" that requires an input x and a default input for y is provided (which is one). The function itself will add the x input given to the default y of 1.

``` r
add <- function(x,y=1) {
  # here is the body of the code
  x + y 
}
```

Testing our new function out, expecting to get the value of 4 + 1:

``` r
add (4)
```

    ## [1] 5

Just because the y have a default input, however, doesn't mean we can't change it later by providing the variable identity and the equal sign, followed by the new value we are temporarily inserting into the function. Here is an example where we expect the code to add 4 and 5:

``` r
add(4, y=5)
```

    ## [1] 9

You can vectorize this function, too, making your function have greater computational use!

``` r
add( c(1,2,3), 1)
```

    ## [1] 2 3 4

However, you can't vectorize with just any syntax. Using commas indicates to R that there are different arguments being input into a function, when in actuality you intend to input a vector. For example:

``` r
#add(1,2,3)
```

See, this doesn't work because it's not a vector but three separate arguments. Our function only even has an x and y component written for it.

Utilizing functions: Rescale!
=============================

``` r
rescale <- function(x) {
  rng <-range(x)
 (x - rng[1]) / (rng[2] - rng[1])
}
```

``` r
x <- c(1, 3, NA, 5, 10) 
 rng <-range(x, na.rm = TRUE)
 rng
```

    ## [1]  1 10

``` r
(x - rng[1]) / (rng[2] - rng[1])
```

    ## [1] 0.0000000 0.2222222        NA 0.4444444 1.0000000

Let's make another function: Rescale!
=====================================

``` r
rescale <- function(x, na.rm=TRUE) { #here we also have it neglect any na data!
 rng <-range(x, na.rm = TRUE)  #function within a function for cleanliness and clarity
 (x - rng[1]) / (rng[2] - rng[1])
}
```

Let's see it in action

``` r
rescale(c(1, 3, NA, 5, 10))
```

    ## [1] 0.0000000 0.2222222        NA 0.4444444 1.0000000

Rescale version 2
-----------------

Here we will implement if/else conditional statements! This allows R to compute whether one condition is met or not and will run unique functions for each.

``` r
rescale2 <- function(x, na.rm=TRUE, plot=FALSE) {
 if(na.rm) {
 rng <-range(x, na.rm=na.rm)
 } else {
 rng <-range(x)
 }
 print("Hello")
 answer <- (x - rng[1]) / (rng[2] - rng[1])
 print("is it me you are looking for?")
 if(plot) {
 plot(answer, typ="b", lwd=4)
 }
 print("I can see it in ...")
}
```

``` r
rescale2(1:10)
```

    ## [1] "Hello"
    ## [1] "is it me you are looking for?"
    ## [1] "I can see it in ..."

rescale version 3

``` r
rescale3 <- function(x, na.rm=TRUE, plot=TRUE) {
 if(na.rm) {
 rng <-range(x, na.rm=TRUE)
 } else {
 rng <-range(x)
 }
 print("Hello")
 answer <- (x - rng[1]) / (rng[2] - rng[1])
 return(answer)
 print("is it me you are looking for?")
 if(plot) {
 plot(answer, typ="b", lwd=4)
 }
 print("I can see it in ...")
}
```

``` r
rescale3(1:10)
```

    ## [1] "Hello"

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

Now on to *Simplyfing Functions*
================================

First we'll load up the package we've installed and will be using henceforth.

``` r
library(bio3d)
```

### Improving code

Here we will be working with this snippet of code to modify any mistakes (copy paste errors!) and etcetera.

``` r
s1 <- read.pdb("4AKE") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

``` r
s2 <- read.pdb("1AKE") # kinase no drug
```

    ##   Note: Accessing on-line PDB file
    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
s3 <- read.pdb("1E4Y") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

``` r
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s1, chain="A", elety="CA")
s1.b <- s1.chainA$atom$b
s2.b <- s2.chainA$atom$b
s3.b <- s3.chainA$atom$b
plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
```

![](class06_rework_files/figure-markdown_github/unnamed-chunk-22-1.png)

``` r
plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor")
```

![](class06_rework_files/figure-markdown_github/unnamed-chunk-22-2.png)

``` r
plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor")
```

![](class06_rework_files/figure-markdown_github/unnamed-chunk-22-3.png)

Here's the corrected version of code:

``` r
s1 <- read.pdb("4AKE") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 2q/h4gv5v9551j4vppw0svdstk40000gn/T//Rtmp6ksqJU/4AKE.pdb exists. Skipping
    ## download

``` r
s2 <- read.pdb("1AKE") # kinase no drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 2q/h4gv5v9551j4vppw0svdstk40000gn/T//Rtmp6ksqJU/1AKE.pdb exists. Skipping
    ## download

    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
s3 <- read.pdb("1E4Y") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 2q/h4gv5v9551j4vppw0svdstk40000gn/T//Rtmp6ksqJU/1E4Y.pdb exists. Skipping
    ## download

``` r
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s3, chain="A", elety="CA")
s1.b <- s1.chainA$atom$b
s2.b <- s2.chainA$atom$b
s3.b <- s3.chainA$atom$b
plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
```

![](class06_rework_files/figure-markdown_github/unnamed-chunk-23-1.png)

``` r
plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor")
```

![](class06_rework_files/figure-markdown_github/unnamed-chunk-23-2.png)

``` r
plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor")
```

![](class06_rework_files/figure-markdown_github/unnamed-chunk-23-3.png)

Answering questions about this code
-----------------------------------

### **Q1.** What type of object is returned from the read.pdb() function?

reads PDB files

### **Q2**. What does the trim.pdb() function do?

Trims a PDB object to a subset of atoms. It produces a smaller PDB object from a larger one. This was discovered by copy&pasting the relevant lines of code by themselves, to isolate their output.

### **Q3**. What input parameter would turn off the marginal black and grey rectangles in the plots and what do they represent in this case?

in the plotb3 function, override the default argument for top and for bot to equal FALSE, rather than TRUE. They represent the beta strands and alpha helices (2º structure annotations).

what the plots look like without these:

``` r
# without the secondary features plotted
library(bio3d)
s1 <- read.pdb("4AKE") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 2q/h4gv5v9551j4vppw0svdstk40000gn/T//Rtmp6ksqJU/4AKE.pdb exists. Skipping
    ## download

``` r
s2 <- read.pdb("1AKE") # kinase no drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 2q/h4gv5v9551j4vppw0svdstk40000gn/T//Rtmp6ksqJU/1AKE.pdb exists. Skipping
    ## download

    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
s3 <- read.pdb("1E4Y") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 2q/h4gv5v9551j4vppw0svdstk40000gn/T//Rtmp6ksqJU/1E4Y.pdb exists. Skipping
    ## download

``` r
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s3, chain="A", elety="CA")
s1.b <- s1.chainA$atom$b
s2.b <- s2.chainA$atom$b
s3.b <- s3.chainA$atom$b
plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor", top = FALSE, bot = FALSE)
```

![](class06_rework_files/figure-markdown_github/unnamed-chunk-24-1.png)

``` r
plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor", top = FALSE, bot = FALSE)
```

![](class06_rework_files/figure-markdown_github/unnamed-chunk-24-2.png)

``` r
plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor", top = FALSE, bot = FALSE)
```

![](class06_rework_files/figure-markdown_github/unnamed-chunk-24-3.png)

### **Q4**. What would be a better plot to compare across the different proteins?

a dendogram

### **Q5**. Which proteins are more similar to each other in their B-factor trends. How could you quantify this?

*hint: here's the code for it*

``` r
hc <- hclust( dist( rbind(s1.b, s2.b, s3.b) ) )
plot(hc)
```

![](class06_rework_files/figure-markdown_github/unnamed-chunk-25-1.png) \*\* s2.b and s3.b are most similar to one another.\*\*

rbind will "pair up" the matricies of separate data sets by their rows, whereas cbind would do columns.

``` r
rbind(s1.b, s2.b, s3.b)
```

    ##       [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11]
    ## s1.b 29.02 18.44 16.20 19.67 20.26 20.55 17.05 22.13 26.71 33.05 30.66
    ## s2.b 37.14 25.76 23.90 17.83 19.86 21.75 20.21 16.92 17.47 18.35 18.31
    ## s3.b 25.46 17.86 10.28  4.73  4.36  5.10  9.59 12.19 11.41  9.39  8.08
    ##      [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22]
    ## s1.b 32.73 25.61 33.19 41.03 24.09 16.18 19.14 29.19 14.79 19.63 28.54
    ## s2.b 20.57 14.56 17.87 11.87 24.63 21.29 35.13 29.68 23.96 32.34 35.34
    ## s3.b  9.01 11.77 12.15 12.72  9.62 12.18 19.95 19.59 15.73 22.51 25.87
    ##      [,23] [,24] [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33]
    ## s1.b 27.49 32.56 17.13 15.50  6.98 24.07 24.00 23.94 30.70 24.70 32.84
    ## s2.b 35.64 38.91 29.00 36.55 28.83 27.15 30.28 28.13 19.90 21.95 25.07
    ## s3.b 23.08 20.97 17.28 12.69 12.24 14.14 14.05  9.38  5.03  7.78 10.13
    ##      [,34] [,35] [,36] [,37] [,38] [,39] [,40] [,41] [,42] [,43] [,44]
    ## s1.b 34.60 33.01 44.60 50.74 57.32 47.04 67.13 81.04 75.20 59.68 55.63
    ## s2.b 16.15 18.35 21.19 27.13 28.55 21.10 38.88 33.63 29.51 29.21 33.01
    ## s3.b  8.96  7.50  5.48  2.97  2.73  3.23  7.81 10.40 10.67 12.79 17.90
    ##      [,45] [,46] [,47] [,48] [,49] [,50] [,51] [,52] [,53] [,54] [,55]
    ## s1.b 45.12 39.04 44.31 38.21 43.70 44.19 47.00 48.67 41.54 50.22 45.07
    ## s2.b 20.92 17.17 25.84 29.80 16.89 24.66 35.62 23.52 23.37 34.41 25.96
    ## s3.b 13.56 12.94 14.78 11.31  8.79 14.13 15.10  7.92  8.15 14.28 14.04
    ##      [,56] [,57] [,58] [,59] [,60] [,61] [,62] [,63] [,64] [,65] [,66]
    ## s1.b 49.77 52.04 44.82 39.75 35.79 38.92 37.93 27.18 26.86 27.53 31.16
    ## s2.b 16.79 20.20 23.72 23.29 25.23 19.81 19.00 20.21 22.62 21.40 23.47
    ## s3.b 12.42 11.84  6.57  9.59 11.84  9.61 12.18  7.89  5.74  5.31  7.67
    ##      [,67] [,68] [,69] [,70] [,71] [,72] [,73] [,74] [,75] [,76] [,77]
    ## s1.b 27.08 23.03 28.12 24.78 24.22 18.69 40.67 38.08 55.26 46.29 26.25
    ## s2.b 23.20 20.21 25.90 30.58 28.25 37.60 44.66 54.46 91.10 92.02 86.85
    ## s3.b  7.99  8.24 12.34 20.98 17.93 16.30 16.94 22.19 22.36 18.96 17.18
    ##      [,78] [,79] [,80] [,81] [,82] [,83] [,84] [,85] [,86] [,87] [,88]
    ## s1.b 37.14 27.50 16.86 27.76 19.27 22.22 26.70 25.52 21.22  15.9 15.84
    ## s2.b 80.21 68.72 42.01 27.69 23.06 21.98 18.60 20.17 15.06  14.2 23.07
    ## s3.b 18.99 16.65 13.39 11.61 10.10 11.03 13.31 12.66  9.44   6.6  5.20
    ##      [,89] [,90] [,91] [,92] [,93] [,94] [,95] [,96] [,97] [,98] [,99]
    ## s1.b 22.44 19.61 21.23 21.79 17.64 22.19 22.73 16.80 23.25 35.95 24.42
    ## s2.b 20.36 25.76 17.02 13.71 23.88 26.72 22.58 24.51 45.23 38.07 36.97
    ## s3.b  5.06  6.16  6.20  6.24  6.34  7.39  7.86 11.66 17.87 17.67 14.63
    ##      [,100] [,101] [,102] [,103] [,104] [,105] [,106] [,107] [,108] [,109]
    ## s1.b  20.96  20.00  25.99  24.39  17.19  12.16  17.35  24.97  14.08  22.01
    ## s2.b  35.17  37.83  43.69  29.14  24.56  25.20  19.27  20.88  18.27  16.96
    ## s3.b  14.30  16.98  19.84  13.36  10.93  11.52   7.56   8.85   7.07  10.08
    ##      [,110] [,111] [,112] [,113] [,114] [,115] [,116] [,117] [,118] [,119]
    ## s1.b  22.26  22.78  27.47  30.49  32.02  20.90  27.03  23.84  44.37  42.47
    ## s2.b  21.38  18.33  23.18  21.15  21.97  22.63   9.74  16.71  26.18  30.39
    ## s3.b  12.34  12.05  13.10  18.63  21.34  15.73  13.16  14.04  18.13  13.59
    ##      [,120] [,121] [,122] [,123] [,124] [,125] [,126] [,127] [,128] [,129]
    ## s1.b  33.48  44.56  56.67  60.18  66.62  59.95  70.81  88.63 100.11  86.60
    ## s2.b  22.95  25.51  20.28  16.86  21.94  20.59  21.64  27.42  35.72  23.47
    ## s3.b  12.12  13.37  10.57   6.60   7.73   7.91  11.31  14.38  14.60  12.25
    ##      [,130] [,131] [,132] [,133] [,134] [,135] [,136] [,137] [,138] [,139]
    ## s1.b  85.80  77.48  68.13  52.66  45.34  52.43  60.90  62.64  72.19  66.75
    ## s2.b  31.57  23.71  19.01  21.52  19.40  24.32  34.28  23.96  23.14  26.60
    ## s3.b  12.33  11.10  11.53  10.44   9.18  11.36  17.28  16.45  15.21  12.11
    ##      [,140] [,141] [,142] [,143] [,144] [,145] [,146] [,147] [,148] [,149]
    ## s1.b  58.73  74.57  79.29  79.53  76.58  66.40  64.76  70.48  74.84  70.11
    ## s2.b  24.94  28.49  28.18  41.64  23.85  28.67  28.76  35.16  35.46  28.74
    ## s3.b  12.12  14.10  14.94  21.72  16.82  12.61  13.40  12.64  12.24   9.13
    ##      [,150] [,151] [,152] [,153] [,154] [,155] [,156] [,157] [,158] [,159]
    ## s1.b  74.82  78.61  78.24  66.70  66.10  67.01  72.28  80.64  68.54  43.23
    ## s2.b  26.99  31.74  40.41  33.73  25.57  29.13  29.74  36.32  22.58  22.82
    ## s3.b  12.31  19.68  19.83  15.34  15.61  14.07  13.64  16.87  11.89  12.92
    ##      [,160] [,161] [,162] [,163] [,164] [,165] [,166] [,167] [,168] [,169]
    ## s1.b  51.24  45.72  61.60  45.61  42.57  41.03  41.02  33.34  19.48  34.38
    ## s2.b  46.67  29.44  25.40  17.27  20.38  21.55  19.19  15.89  18.37  30.51
    ## s3.b  19.93  23.72  23.13  13.35  11.51  18.51  17.24  11.92  12.36  13.42
    ##      [,170] [,171] [,172] [,173] [,174] [,175] [,176] [,177] [,178] [,179]
    ## s1.b  33.11  25.48  29.68  40.71  32.91  24.41  19.20  15.43  19.93  20.66
    ## s2.b  18.47  11.70  18.45  24.75  16.63  20.80  19.62  22.56  19.87  20.22
    ## s3.b  11.45  11.09  14.19  14.22  12.15  10.49  11.29  11.74   9.53   7.65
    ##      [,180] [,181] [,182] [,183] [,184] [,185] [,186] [,187] [,188] [,189]
    ## s1.b  12.72  21.40  18.21  26.68  34.50  25.77  26.52  36.85  31.05  39.84
    ## s2.b  21.16  22.13  20.66  22.82  32.86  26.04  20.60  44.44  35.28  38.03
    ## s3.b   7.21   7.56   8.14  11.07  16.93  11.12   8.79  16.03  18.87  17.72
    ##      [,190] [,191] [,192] [,193] [,194] [,195] [,196] [,197] [,198] [,199]
    ## s1.b  48.03  23.04  29.57  23.00  23.80  26.59  25.49  23.25  19.89  32.37
    ## s2.b  28.46  29.10  30.19  26.17  22.71  23.39  23.44  16.27  21.26  24.67
    ## s3.b  14.72  14.08  14.21   9.99   6.63  10.11  12.64  15.06  14.21  14.20
    ##      [,200] [,201] [,202] [,203] [,204] [,205] [,206] [,207] [,208] [,209]
    ## s1.b  30.97  42.16  29.64  29.69  33.15  26.38  23.17  29.35  32.80  25.92
    ## s2.b  19.12  23.26  21.75  24.59  27.26  22.63  26.40  31.60  29.57  30.90
    ## s3.b  16.39  16.31  16.07  17.83  20.24  14.28  17.10  17.00  18.88  17.13
    ##      [,210] [,211] [,212] [,213] [,214]
    ## s1.b  38.01  45.95  44.26  44.35  70.26
    ## s2.b  32.29  46.86  41.73  49.31  66.76
    ## s3.b  23.68  24.72  19.74  24.12  33.57

Below is the workout of my homework portion of this lab
=======================================================

### Q6

``` r
# without the secondary features plotted
library(bio3d)
s1 <- read.pdb("4AKE") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 2q/h4gv5v9551j4vppw0svdstk40000gn/T//Rtmp6ksqJU/4AKE.pdb exists. Skipping
    ## download

``` r
s2 <- read.pdb("1AKE") # kinase no drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 2q/h4gv5v9551j4vppw0svdstk40000gn/T//Rtmp6ksqJU/1AKE.pdb exists. Skipping
    ## download

    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
s3 <- read.pdb("1E4Y") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 2q/h4gv5v9551j4vppw0svdstk40000gn/T//Rtmp6ksqJU/1E4Y.pdb exists. Skipping
    ## download

``` r
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s3, chain="A", elety="CA")
s1.b <- s1.chainA$atom$b
s2.b <- s2.chainA$atom$b
s3.b <- s3.chainA$atom$b
plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor", top = FALSE, bot = FALSE)
```

![](class06_rework_files/figure-markdown_github/unnamed-chunk-27-1.png)

``` r
plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor", top = FALSE, bot = FALSE)
```

![](class06_rework_files/figure-markdown_github/unnamed-chunk-27-2.png)

``` r
plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor", top = FALSE, bot = FALSE)
```

![](class06_rework_files/figure-markdown_github/unnamed-chunk-27-3.png)

``` r
# . . . . . . . . . . . . . .

hc <- hclust( dist( rbind(s1.b, s2.b, s3.b) ) )
plot(hc)
```

![](class06_rework_files/figure-markdown_github/unnamed-chunk-27-4.png)

``` r
library(bio3d)
s1 <- read.pdb("4AKE") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 2q/h4gv5v9551j4vppw0svdstk40000gn/T//Rtmp6ksqJU/4AKE.pdb exists. Skipping
    ## download

``` r
s2 <- read.pdb("1AKE") # kinase no drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 2q/h4gv5v9551j4vppw0svdstk40000gn/T//Rtmp6ksqJU/1AKE.pdb exists. Skipping
    ## download

    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
s3 <- read.pdb("1E4Y") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 2q/h4gv5v9551j4vppw0svdstk40000gn/T//Rtmp6ksqJU/1E4Y.pdb exists. Skipping
    ## download

``` r
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s3, chain="A", elety="CA")
s1.b <- s1.chainA$atom$b
s2.b <- s2.chainA$atom$b
s3.b <- s3.chainA$atom$b
plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor", top = FALSE, bot = FALSE)
```

![](class06_rework_files/figure-markdown_github/unnamed-chunk-28-1.png)

``` r
plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor", top = FALSE, bot = FALSE)
```

![](class06_rework_files/figure-markdown_github/unnamed-chunk-28-2.png)

``` r
plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor", top = FALSE, bot = FALSE)
```

![](class06_rework_files/figure-markdown_github/unnamed-chunk-28-3.png)

``` r
# creating a function to assemble the hclust with the three PDB objects
hc_plot <- function (x,y,z) {
  hclust( dist( rbind(x, y, z)) )
}

# plotting out the created dendogram
  plot(
    hc_plot(s3.b, s2.b, s1.b)
  )
```

![](class06_rework_files/figure-markdown_github/unnamed-chunk-28-4.png)

``` r
s3 <- read.pdb("1E4Y") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 2q/h4gv5v9551j4vppw0svdstk40000gn/T//Rtmp6ksqJU/1E4Y.pdb exists. Skipping
    ## download

``` r
s3.chainA <- trim.pdb(s3, chain="A", elety="CA")

s3.b <- s3.chainA$atom$b

plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor", top = FALSE, bot = FALSE)
```

![](class06_rework_files/figure-markdown_github/unnamed-chunk-29-1.png) EDIT

``` r
s3 <- read.pdb("1E4Y") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 2q/h4gv5v9551j4vppw0svdstk40000gn/T//Rtmp6ksqJU/1E4Y.pdb exists. Skipping
    ## download

``` r
s3.chainA <- trim.pdb(s3, chain="A", elety="CA")


plotb3(s3.chainA$atom$b,
       sse=s3.chainA, typ="l", ylab="Bfactor", top = FALSE, bot = FALSE)
```

![](class06_rework_files/figure-markdown_github/unnamed-chunk-30-1.png)

EVEN MORE

``` r
s3 <- read.pdb("1E4Y") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 2q/h4gv5v9551j4vppw0svdstk40000gn/T//Rtmp6ksqJU/1E4Y.pdb exists. Skipping
    ## download

``` r
plotb3(trim.pdb(s3, chain="A", elety="CA")
       $atom$b,
       sse=s3.chainA, typ="l", ylab="Bfactor", top = FALSE, bot = FALSE)
```

![](class06_rework_files/figure-markdown_github/unnamed-chunk-31-1.png) EVEN MORE X2 \#\#\#this is the one!!

``` r
# entering in the PDB files you want to work with

b1 <- trim.pdb(read.pdb("1E4Y"), chain="A", elety="CA")$atom$b    # kinase with drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 2q/h4gv5v9551j4vppw0svdstk40000gn/T//Rtmp6ksqJU/1E4Y.pdb exists. Skipping
    ## download

``` r
b2 <- trim.pdb(read.pdb("4AKE"), chain="A", elety="CA")$atom$b    # kinase with drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 2q/h4gv5v9551j4vppw0svdstk40000gn/T//Rtmp6ksqJU/4AKE.pdb exists. Skipping
    ## download

``` r
b3 <-trim.pdb(read.pdb("1AKE"), chain="A", elety="CA")$atom$b     # kinase no drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 2q/h4gv5v9551j4vppw0svdstk40000gn/T//Rtmp6ksqJU/1AKE.pdb exists. Skipping
    ## download

    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
# self defined function to assemble the dendograph of each three PDB objects
hc_plot <- function (x,y,z) {
  hclust( dist( rbind(x, y, z)) )
}

# plotting out the created dendogram, provides a single graph providing a B-factor trends comparison
  plot( hc_plot(b1, b2, b3) )
```

![](class06_rework_files/figure-markdown_github/unnamed-chunk-32-1.png)
