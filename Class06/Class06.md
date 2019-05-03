Class 06: R Functions
================
Madeline R. Luth
4/19/2019

Overview
--------

Today we will focus on **R** functions but we will start with some **file reading**.

``` r
plot(1:10, type = "l", col = "blue")
```

![](Class06_files/figure-markdown_github/unnamed-chunk-1-1.png)

``` r
read.table("test1.txt", header = TRUE, sep = ",")
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

Read in first file which is delimited by ","

``` r
read.table("test2.txt", header = TRUE, sep = "$")
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

Read in second file which is delimited by "$"

``` r
read.table("test3.txt")
```

    ##   V1 V2 V3
    ## 1  1  6  a
    ## 2  2  7  b
    ## 3  3  8  c
    ## 4  4  9  d
    ## 5  5 10  e

Read in final file which has spaces and tabs

Note: Option + Command + i is shortcut to insert a new R chunk

``` r
read.csv("https://bioboot.github.io/bggn213_W19/class-material/test2.txt")
```

    ##   Col1.Col2.Col3
    ## 1          1$2$3
    ## 2          4$5$6
    ## 3          7$8$9
    ## 4          a$b$c

Can also read.csv from a web domain instead of downloading file to the working directory

Our first function
------------------

Add some numbers

``` r
add <- function(x, y=1) {
  #The body of the function
  x + y
}
```

``` r
add(4)
```

    ## [1] 5

``` r
add(4,5)
```

    ## [1] 9

Adding the 5 for y will override the default argument of y=1

``` r
add(c(1,3,5), 1)
```

    ## [1] 2 4 6

``` r
#add(1,3,5)
#returns error of unused argument
```

``` r
#add(x=1, y="barry")
#error because second value is non-numeric
```

``` r
rescale <- function(x) {
   rng <-range(x)
   (x - rng[1]) / (rng[2] - rng[1])
}
```

``` r
rescale( 1:10 )
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
rescale( c(1,3,NA, 5, 10) )
```

    ## [1] NA NA NA NA NA

``` r
rescale2 <- function(x, na.rm=TRUE) {
  rng <- range(x, na.rm=na.rm)
  (x - rng[1]) / (rng[2] - rng[1])
}
```

``` r
rescale2( c(1,3,NA, 10))
```

    ## [1] 0.0000000 0.2222222        NA 1.0000000

``` r
rescale3 <- function(x, na.rm=TRUE, plot=FALSE) {
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

This is another example

``` r
rescale3( 1:10 )
```

    ## [1] "Hello"
    ## [1] "is it me you are looking for?"
    ## [1] "I can see it in ..."

``` r
rescale3( 1:10, plot = TRUE)
```

    ## [1] "Hello"
    ## [1] "is it me you are looking for?"

![](Class06_files/figure-markdown_github/unnamed-chunk-19-1.png)

    ## [1] "I can see it in ..."
