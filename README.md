# mizerHowTo
This package contains several tutorials aiming to learn how to handle mizer models and to be able to calibrate mizer models with empirical data.

How to install:

The `remotes` package is needed to install packages hosted on GitHub. 

Use `install.packages("remotes")` if necessary

Then `remotes::install_github("sizespectrum/mizerHowTo")`

Finally, load the newly installed package with

`library(mizerHowTo)`

To display any tutorial, use the `tutorial()` function.

Available tutorials are
  
- HTM0: Examples on the use of Mizer

Other tutorials are in preparation but are currently incomplete and partially
outdated.

Examples: 

`tutorial("HTM0","html")` will display the online version of the first tutorial.

`tutorial("HTM0","Rmd")` will display the code of the first tutorial.

Several shiny apps (R code with user interface) are embedded in the tutorials. One can start them from the Rmarkdown files (.Rmd) by running the appropriate chunk of code.
All shiny apps functions follows this expression `shiny_XXX()`.

