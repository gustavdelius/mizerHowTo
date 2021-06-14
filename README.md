# mizerHowTo
This package contains several tutorials aiming to learn how to handle mizer models and to be able to calibrate mizer models with empirical data.

How to install:

The `devtools` package is needed to install packages hosted on GitHub. 

Use `install.packages("devtools")` if necessary

Then `devtools::install_github("sizespectrum/mizerHowTo")`

Finally, load the newly installed package with

`library(mizerHowTo)`

To display any tutorial, use the `tutorial()` function.

Available tutorials are

- HTM1: How to parametrise a Mizer model

- HTM2: How to calibrate a Mizer model using averaged data over time

Examples: 

`tutorial("HTM1","html")` will display the online version of the first tutorial.

`tutorial("HTM2","Rmd")` will display the code of the second tutorial.

Several shiny apps (R code with user interface) are embedded in the tutorials. One can start them from the Rmarkdown files (.Rmd) by running the appropriate chunk of code.
All shiny apps functions follows this expression `shiny_XXX()`.

