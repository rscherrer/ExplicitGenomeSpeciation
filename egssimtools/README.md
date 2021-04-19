# egsimtools: R package for the analysis of the EGS simulation

Check out the vignette of the package in the `doc` folder. You can visualize the `vignette.pdf` version online from GitHub. Alternatively, download the repository and open the `vignette.html` version in your browser. You can also visit open the website of this package by opening `index.html` from the `docs` folder in your browser. This website contains the vignette as well as an overview of all the functions.

Installation instructions can be found in the vignette.

**Note:** the `tmp` folder contains old functions that have been removed from the package. These include functions to plot the data directly from simulation folders. I removed them from the package because they are too constraining with respect to the breadth of analyses we can perform when compared to regular tidyverse pipelines, which are way more flexible. For example, it is simple to expand a pipeline to plot multiple simulations instead of a single one, while this is cumbersome if we use the functions in `tmp`. See the vignette for examples of how to plot the data in multiple ways.
