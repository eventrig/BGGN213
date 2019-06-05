class07: R Functions and packages
================
Emilya Ventriglia
4/24/2019

More on fucntion writting
-------------------------

First we'll revisit our function from last class (class 6)

``` r
source("http://tinyurl.com/rescale-R")
```

Test the **rescale()** function

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
x <- (c(1:10, "string"))
!is.numeric(x)
```

    ## [1] TRUE

``` r
rescale2 <- function(x, na.rm=TRUE, plot=FALSE, ...) {
 if( !is.numeric(x) ) {
 stop("Input x should be numeric", call.=FALSE)
 }
 rng <-range(x, na.rm=na.rm)

 answer <- (x - rng[1]) / (rng[2] - rng[1])
 if(plot) {
 plot(answer, ...)
 }
 return(answer)
}

rescale2(c(1:10))
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

provides a much more useful error message! (call. is when the console prints out the line of input that created the error)

Function Practice
-----------------

Write a function to identify the NA elemetns in two vectors:

Start with a simpile example that is more conceptually manageable:

``` r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
```

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

but now design a function for when both vectors have simultaneous "TRUE"s:

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

calculate how many na matches there are using summation:

``` r
sum( is.na(x) & is.na(y) )
```

    ## [1] 1

### Compiling this snippet into a working function

This is my working snippet of code that I can use as the body of my first function:

``` r
both_na <- function(x, y) {
 sum( is.na(x) & is.na(y) )
}
```

test:

``` r
both_na(x,y)
```

    ## [1] 1

test2:

``` r
both_na(c(NA,NA,NA),c(NA,NA,1))
```

    ## [1] 2

test3:

``` r
both_na(c(NA,NA,NA),c(1,NA,NA))
```

    ## [1] 2

intentional break using two vectors of different lengths:

``` r
both_na(c(1,NA,NA),c(1,NA,NA,NA))
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 2

``` r
both_na(c(1,NA,NA),c(1,NA,NA,NA,NA,NA))
```

    ## [1] 4

Uh oh, R is recycling the inputs of the first vector in order to accomodate the mismatching length. In other words, it restarts from the beginning of the shorter vector to match it with the subsequent inputs in the longer vector. *we need to catch this before hand to prevent it!*

### Checking that the legnth of our inputs are equal

"!="" means not equal to

``` r
x <- c(1,NA,NA)
y <- c(1,NA,NA,NA,NA,NA)
length(x) != length(y)
```

    ## [1] TRUE

**failing early and loud** by stopping this process, therby preventing this erroneous recycling.

both\_na2 &lt;- function(x, y) { if(length(x) != length(y)) { stop("Input x and y should be the same length") } sum( is.na(x) & is.na(y) ) } \# trying it out both\_na2(c(1,NA,NA), c(1,NA,NA,NA,NA,NA))

now with more information feedback features:

``` r
both_na3 <- function(x, y) {
 if(length(x) != length(y)) {
 stop("Input x and y should be vectors of the same length")
 }

 na.in.both <- ( is.na(x) & is.na(y) )
 na.number <- sum(na.in.both)
 na.which <- which(na.in.both)
 message("Found ", na.number, " NA's at position(s):",
 paste(na.which, collapse=", ") )

 return( list(number=na.number, which=na.which) )
}

#trying it out
x <- c(1,NA,NA)
y <- c(1,NA,NA)
both_na3(x,y)
```

    ## Found 2 NA's at position(s):2, 3

    ## $number
    ## [1] 2
    ## 
    ## $which
    ## [1] 2 3

Writing a Function for grade()
------------------------------

grade() should determine an overall grade from a vector fo student homework scores and drop the lowest single alignment score

student 1
=========

``` r
student_1 <- c(100, 100, 100, 100, 100, 100, 100, 90)

student_1_sum <- sum(student_1) - min(student_1)

student_1_sum / 7
```

    ## [1] 100

student 2
=========

``` r
student_2 <- c(100, NA, 90, 90, 90, 90, 97, 80)

student_2_sum <- sum(student_2, na.rm=TRUE) - min(student_2,na.rm=TRUE)

student_2_sum / 7
```

    ## [1] 79.57143

you can't use the trim argument in mean() because it trims both ends

``` r
student_1 <- c(100, 100, 100, 100, 100, 100, 100, 90)

mean(student_1, trim=min(student_1))
```

    ## [1] 100

the concise one for class:

``` r
 grade <- function(x) {
   (sum(x, na.rm=TRUE) - min(x, na.rm=TRUE))/(length(x)-1)
 }
grade(student_1)
```

    ## [1] 100

``` r
grade(student_2)
```

    ## [1] 79.57143

Now to grade all students in an example class

``` r
url <- "http://tinyurl.com/gradeinput"
students <- read.csv(url, row.names = 1) 
head(students)
```

    ##           hw1 hw2 hw3 hw4 hw5
    ## student-1 100  73 100  88  79
    ## student-2  85  64  78  89  78
    ## student-3  83  69  77 100  77
    ## student-4  88  NA  73 100  76
    ## student-5  88 100  75  86  79
    ## student-6  89  78 100  89  77

``` r
grade <- function(x) {
   (sum(x, na.rm=TRUE) - min(x, na.rm=TRUE))/(length(x)-1)
}
```

``` r
grade(students[5,])
```

    ## [1] 88.25

``` r
ans <- apply(students,1,grade)
```

figuring out who the top students are by sorting the averages by high to low:

``` r
sort(ans, decreasing =TRUE)
```

    ##  student-7  student-8 student-13  student-1 student-12 student-16 
    ##      94.00      93.75      92.25      91.75      91.75      89.50 
    ##  student-6  student-5 student-17  student-9 student-14 student-11 
    ##      89.00      88.25      88.00      87.75      87.75      86.00 
    ##  student-3 student-19 student-20  student-2 student-18  student-4 
    ##      84.25      82.75      82.75      82.50      72.75      66.00 
    ## student-15 student-10 
    ##      62.50      61.00

Intersect data example
======================

### one last function example...

find the intersection of two vectors

``` r
# Start with a simple version of the problem
df1 <- data.frame(IDs=c("gene1", "gene2", "gene3"),
 exp=c(2,1,1),
 stringsAsFactors=FALSE)
df2 <- data.frame(IDs=c("gene2", "gene4", "gene3", "gene5"),
 exp=c(-2, NA, 1, 2),
 stringsAsFactors=FALSE)
# Simplify further to single vectors
x <- df1$IDs
y <- df2$IDs
# Now what do we do?
```

google: 'R intersection of two vectors' to find out how to get us started

``` r
x <- df1$IDs
y <- df2$IDs


intersect(x,y)
```

    ## [1] "gene2" "gene3"

``` r
x %in% y
```

    ## [1] FALSE  TRUE  TRUE

``` r
x
```

    ## [1] "gene1" "gene2" "gene3"

``` r
y
```

    ## [1] "gene2" "gene4" "gene3" "gene5"

this has given me the indices (the position) of where the things that overlap are!

and can show you which they are

``` r
x[x %in% y]
```

    ## [1] "gene2" "gene3"

``` r
gene_intersect <- function(x, y) {
  cbind( x[x %in% y], y[ y %in% x])
}
```

``` r
df1 <- data.frame(IDs=c("gene1", "gene2", "gene3"),
 exp=c(2,1,1),
 stringsAsFactors=FALSE)
df3 <- data.frame(IDs=c("gene2", "gene2", "gene5", "gene5"),
 exp=c(-2, NA, 1, 2),
 stringsAsFactors=FALSE)

merge(df1, df2, by="IDs")
```

    ##     IDs exp.x exp.y
    ## 1 gene2     1    -2
    ## 2 gene3     1     1
