Class07: R functions and packages
================
Madeline R. Luth
4/24/2019

More on function writing
------------------------

First we will revisit our functions from last week

``` r
source("http://tinyurl.com/rescale-R")
```

``` r
rescale <- function(x, na.rm=TRUE, plot=FALSE, ...) {
   rng <-range(x, na.rm=na.rm)
   answer <- (x - rng[1]) / (rng[2] - rng[1])
   if(plot) {
     plot(answer, ...)
}
   return(answer)
}
```

Test the rescale function This will fail

``` r
#rescale( c(1:10, "string") )
```

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
```

This new "if" statement will check if the input of x is numeric. If x is numeric, it will return TRUE. If x is not numeric, then a new warning message will indicate to the user that their input should be numeric.

Note: Using an "!" in front of a statement flips the logic of that statement.

``` r
#rescale2( c(1:10, "string") )
```

Note that there is now an error message that is more meaningful to the function user.

Function practice
-----------------

Write a function to identify NA elements in two vectors

Start with a simple example input where I know what the answer should be.

``` r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA,3,NA,3, 4)
```

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

I want to find the position where both entries are NA, i.e. the return value of is.na is TRUE for both. (will be the 3rd position)

Use an AND statement

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

Can calculate the sum to figure out how many times this occurs. We know that it will equal 1.

``` r
sum( is.na(x) & is.na(y) )
```

    ## [1] 1

This will be the working snippet of code. Now we can write the function.

``` r
both.na <- function(x,y) {
  sum( is.na(x) & is.na(y) )
}
```

``` r
both_na(x,y)
```

    ## [1] 1

Function works

``` r
both_na( c(NA,NA,NA), c(NA,NA,1) )
```

    ## [1] 2

``` r
both_na( c(NA,NA,NA), c(1,NA,NA,NA) )
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 3

Notice that an answer is returned, but there is also a warning message. Even though the vectors are different length, R will recycle through the x vector (the shorter vector) and put the first position (NA) into the fourth position.

We should check the lengths of each vector and ensure they are equal before we run the rest of the function.

Note that we used the "!" to return TRUE if x is NOT equal to y. We can then turn this into an if statement for our function to give an error message to the user.

``` r
x <- c(NA, NA, NA)
y <- c(1, NA, NA, NA, NA, NA)
length(x) != length(y)
```

    ## [1] TRUE

Let's write that function with the if statement:

``` r
both_na2 <- function(x, y) {
  if(length(x) != length(y)) {
   stop("Input x and y should be the same length")
}
  sum( is.na(x) & is.na(y) )
}
```

Check what the output is with above vectors for x and y.

``` r
#both_na2(x,y)
```

The error message tells us that the vectors are not equal in length.

We can make our function even more useful by adding extra output messages.

``` r
both_na3 <- function(x, y) {
if(length(x) != length(y)) {
stop("Input x and y should be vectors of the same length")
}
  na.in.both <- ( is.na(x) & is.na(y) )
  na.number  <- sum(na.in.both)
  na.which   <- which(na.in.both)
  message("Found ", na.number, " NA's at position(s):",
          paste(na.which, collapse=", ") )
  return( list(number=na.number, which=na.which) )
}
```

``` r
x <- c( 1, 2, NA, 3, NA)
y<-c(NA,3,NA,3, 4)
both_na3(x,y)
```

    ## Found 1 NA's at position(s):3

    ## $number
    ## [1] 1
    ## 
    ## $which
    ## [1] 3

Student grading example function
--------------------------------

Write a function grade() to determine an overall grade from a vector of student homework assignment scores,dropping the lowest single assignment score.

student 1 c(100, 100, 100, 100, 100, 100, 100, 90)

student 2 c(100, NA, 90, 90, 90, 90, 97, 80)

``` r
# I am going to try to get this to work with student 1 only first
x <- c(100, 100, 100, 100, 100, 100, 100, 90)
```

Need to find the minimum value from the vector and then subtract it from the sum of the overall value. Then divide by the length of the vector - 1.

``` r
min(x)
```

    ## [1] 90

``` r
# Remove the minimum value from the vector
x_drop <- x[x!= min(x)]
x_drop
```

    ## [1] 100 100 100 100 100 100 100

Next, sum the vector x\_drop

``` r
points_total <- sum(x_drop)
```

Now divide the total points (a single number, not a vector) by the length of the dropped score vector (again, a single number)

``` r
final_grade <- points_total/length(x_drop)
final_grade
```

    ## [1] 100

Alternatively: (sum(x) - min(x)) / (length(x) - 1)

Now let's write the full function: can just copy the snippet, go to Code, then click "Extract function", and name your function

``` r
grade <- function(x) {
  (sum(x, na.rm = TRUE) - min(x, na.rm = TRUE)) / (length(x) - 1)
}
```

Let's test the function

``` r
grade(x)
```

    ## [1] 100

We get the expected result, so check student 2.

``` r
y <- c(100, NA, 90, 90, 90, 90, 97, 80)
grade(y)
```

    ## [1] 79.57143

There is an NA, so the result is NA. Recall the argument na.rm Go back and add na.rm = TRUE to the function.

Now we will grade the students in an example class and determine the top 5 students.

url <https://tinyurl.com/gradeinput>

``` r
#read in the csv and assign student names to column 1
students <- read.csv("https://tinyurl.com/gradeinput", row.names = 1)
```

``` r
grade(students[1,])
```

    ## [1] 91.75

This graded the first student, but we need to run it across all 20 students. Use the apply function.

``` r
ans <- apply(students, 1, grade)
```

Now we want to know the top 5 students.

``` r
sort(ans, decreasing = TRUE)
```

    ##  student-7  student-8 student-13  student-1 student-12 student-16 
    ##      94.00      93.75      92.25      91.75      91.75      89.50 
    ##  student-6  student-5 student-17  student-9 student-14 student-11 
    ##      89.00      88.25      88.00      87.75      87.75      86.00 
    ##  student-3 student-19 student-20  student-2 student-18  student-4 
    ##      84.25      82.75      82.75      82.50      72.75      66.00 
    ## student-15 student-10 
    ##      62.50      61.00

Another function example, but more biologically relevant
--------------------------------------------------------

Find common genes in two data sets and return their associated data (from each data set)

``` r
df1
```

    ##     IDs exp
    ## 1 gene1   2
    ## 2 gene2   1
    ## 3 gene3   1

``` r
df2
```

    ##     IDs exp
    ## 1 gene2  -2
    ## 2 gene4  NA
    ## 3 gene3   1
    ## 4 gene5   2

Simplify the dataframes to single vectors by ID (gene 1, 2, etc.)

``` r
x <- df1$IDs
y <- df2$IDs
```

Try function intersect

``` r
intersect(x,y)
```

    ## [1] "gene2" "gene3"

We get the intersection points, but we want the rest of the data to be included

From the intersect help page "See Also": *%in%* is a more intuitive interface as a binary operator, which returns a logical vector indicating if there is a match or not for its left operand.

``` r
x %in% y
```

    ## [1] FALSE  TRUE  TRUE

Tells us gene 1 is not in y, but gene 2 and gene 3 are. Returns indices.

``` r
x[x %in% y]
```

    ## [1] "gene2" "gene3"

``` r
y %in% x
```

    ## [1]  TRUE FALSE  TRUE FALSE

We can now cbind the results

``` r
cbind( x[x %in% y], y[y %in% x] )
```

    ##      [,1]    [,2]   
    ## [1,] "gene2" "gene2"
    ## [2,] "gene3" "gene3"

Make this snippet into a function

``` r
gene_intersect <- function(x, y) {
  cbind( x[x %in% y], y[y %in% x] )
}
```

``` r
merge(df1, df2, by="IDs")
```

    ##     IDs exp.x exp.y
    ## 1 gene2     1    -2
    ## 2 gene3     1     1
