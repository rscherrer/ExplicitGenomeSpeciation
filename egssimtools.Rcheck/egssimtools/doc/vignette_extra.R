#
# Functions are provided for (1) reading the data into R data frames, suitable for flexible manipulation, (2) plotting the data from individual simulations, and (3) combining data from multiple simulations into a single data frame.
#
# ## Plot single simulations
#
# We provide here a number of functions to extract relevant information from the simulations. First we set up the root folder of our simulation of interest.
#
# ```{r}
# root <- system.file("extdata", "example_1", package = "egssimtools")
# ```
#
# The EGS simulations are complex and can save many different kinds of data, with different units of observation: simulation-wise data (e.g. ecological divergence through time), individual-wise data (e.g. individual trait values), locus-wise data or edge-wise data. Again, please refer to the main page of the program for a detailed description.
#
# ### Simulation-wide variable
#
# To plot a simulation-wide variable through time, use
#
# ```{r, fig.width = 4, fig.height = 3}
# # Ecological isolation through time
# plot_sim_line(root, "EI")
# ```
#
# This function can also be used to plot variables against each other, for example:
#
#   ```{r, fig.width = 4, fig.height = 3}
# # Reproductive versus ecological isolation
# plot_sim_line(root, "RI", x = "EI")
# ```
#
# Even better, you can plot against each other variables that belong to the same `.dat` data file, for example:
#
#   ```{r, fig.width = 4, fig.height = 3}
# # Fst of the ecological trait versus that of the mating trait
# plot_sim_line(root, y = "Fst", x = "Fst", by = 3, j = 1, k = 2)
# ```
#
# ### Individual-specific variable
#
# To plot the distribution of an individual variable at a certain time point, use:
#
#   ```{r, fig.width = 4, fig.height = 3}
# # Plot the distribution of ecological trait values at the end of the simulation
# plot_pop_density(
#   root, "individual_trait", by = 3, j = 1, plot_type = "density"
# )
# ```
#
# And to get a distribution through time, use:
#
#   ```{r, fig.width = 4, fig.height = 3}
# plot_pop_bin2d(root, "individual_trait", by = 3, j = 1)
# ```
#
# It is also possible to plot the co-distribution of two individual-specific variables at a specific time point,
#
# ```{r, fig.width = 4, fig.height = 3}
# plot_pop_bin2d(
#   root, y = "individual_trait", x = "individual_trait", by = 3, j = 1, k = 2
# )
# ```
#
# ### Locus-specific variables
#
# For locus-specific variables, use the `plot_genome_*` family of functions, such as, for a genome scan at a specific time point,
#
# ```{r, fig.width = 4, fig.height = 2}
# plot_genome_scan(root, "genome_Fst")
# ```
#
# It is possible to map different variables of the genetic architecture to some aesthetics of the plot, like color:
#
#   ```{r, fig.width = 4, fig.height = 2}
# plot_genome_scan(root, "genome_Fst", mapping = list(color = "trait"))
# ```
#
# Genome scans through time can be viewed using:
#
#   ```{r, fig.width = 4, fig.height = 4}
# plot_genome_heatmap(root, "genome_Fst")
# ```
#
# To plot distributions over the genome, use:
#
#   ```{r, fig.width = 4, fig.height = 3}
# plot_genome_density(root, "genome_Fst", plot_type = "histogram", bins = 50)
# ```
#
# Which can itself take different aesthetic mappings and take different forms, such as violin plots:
#
#   ```{r, fig.width = 4, fig.height = 3}
# plot_genome_density(
#   root, "genome_Fst", x = "trait", mapping = list(fill = "trait"),
#   plot_type = "violin"
# )
# ```
#
# It is also possible to visualize distributions through time,
#
# ```{r, fig.width = 4, fig.height = 3, message = FALSE}
# plot_genome_ridges(
#   root, "genome_Fst", times = c(100, 200, 300), facet_by = "trait"
# )
# ```
#
# ### Plot the gene network
#
# (FIX THIS RAPH)
#
# To plot a variable mapped onto the gene regulatory network at a certain time point, use:
#
#   ```{r, fig.width = 4, fig.height = 3}
# # Plot Fst mapped onto the network (takes a while)
# plot_network(root, "genome_Fst")
# ```
#
# One can also plot a gene network without any simulated data through time, and instead map as the size aesthetics a variable from the genetic architecture (such as the degree or effect size of a node), or no mapping at all:
#
#   ```{r, fig.width = 4, fig.height = 3}
# plot_network_blank(root)
# ```
#
# ### Note
#
# ### Read a single simulation
#
# The most important function is `read_data`. This function allows the user to read the data from one simulation. The user provides the names of the variables to be read, optional parameters to be read too, and the output is returned in a data frame, suitable for further analyses or merging with other data frames. We can for example read speciation metrics through time for one simulation, and add a specific set of parameters to them:
#
#   ```{r}
# root <- "../inst/extdata/example_1"
# data <- read_data(
#   root, variables = c("time", "EI", "RI", "SI"),
#   parnames = c("ecosel", "hsymmetry")
# )
# head(data, 4)
# ```
#
# For these variables, there is one value per time point throughout the simulation, while the values of the two parameters are unique to the simulation. Parameter values are therefore duplicated as many times as there are rows in the resulting table.
#
# Here, all variables have the same length (one value per time point). Simulation output sometimes has different dimensions (see the main page for details). Genome-wide Fst, for example, is computed for each of the three traits at every time point, resulting in a linear, binary data file that is three times as long as `time.dat`. To append genome-wide Fst to time points, we need to split the `Fst.dat` file in three columns. We do this using the `by` argument:
#
#   ```{r}
# data <- read_data(root, variables = c("time", "Fst"), by = c(1, 3))
# head(data, 4)
# ```
#
# The argument `by` specifies how many columns to split each data vector into. We can also use it, for example, to split locus-specific Fst (calculated on a per locus per time point basis) into as many columns as there are loci:
#
#   ```{r}
# data <- read_data(root, variables = c("time", "genome_Fst"), by = c(1, 300))
# head(data[, 1:6], 4)
# ```
#
# But we may also want to add to this dataset other locus-specific data such as allele frequencies or locus-specific additive variance. Then, keeping the variables as single columns would allow to easily identify multiple values measured at a specific locus at a specific time by finding its row. To bind data recorded on a per locus, per time point basis to time point data, for example, we need to duplicate the time column by the number of loci. We do this using the `dupl` argument:
#
#   ```{r}
# data <- read_data(
#   root, variables = c("time", "genome_Fst", "genome_freq", "genome_varA"),
#   dupl = c(300, 1, 1, 1)
# )
# head(data, 4)
# ```
#
# The combination of `by` and `dupl` allows to perform all kinds of combinations of data computed at different frequencies and with different dimensions.
#
# One important note is that `dupl` can also take a character string as one of its elements. In this case, the function is expecting this character string to be the name of a binary file to read, in which are stored the numbers of times that each element of the corresponding variable must be repeated. This is useful when different values must be duplicated different numbers of times. A typical example is the case of individual trait values through time, because the number of individual varies from one generation to the next. The syntax may look like this:
#
#   ```{r}
# data <- read_data(
#   root, variables = c("time", "individual_trait"), by = c(1, 3),
#   dupl = list("population_size", 1)
# )
# head(data, 4)
# ```
#
# Here, `individual_trait` is split into three columns (one per trait), then each element in the `time` column (i.e. each time point) is repeated as many times as there are individuals in this time point, which is read in `population_size.dat`. Note that `dupl` is now a list because we cannot combine a character string and an integer into the same vector.
#
# ### Read parameters
#
# Simulation parameters can be added to the resulting dataset by specifying the names of the parameters to read in `parnames`. Parameter values will be read from a parameter file, assuming it exists in the simulation folder. By default the function looks for a file called `paramlog.txt`, but this can be changed. Not all parameter values are numbers, and for this reason the parameters are read as factors. You can specify in `as_numeric` the names of parameters to convert into numeric. Some parameters are composite and contain multiple values, such as `nvertices` which contains the number of loci for each trait. Use `combine = FALSE` to split these parameters into their constituent values and give one column to each.
#
# Note that it is possible to only read the parameters of a simulation by using the `read_parameters` function, which returns a list, e.g.:
#
#   ```{r}
# read_parameters(root, c("ecosel", "hsymmetry"))
# ```
#
# ## Combining simulations
#
# Use the `collect_sims` function to combine the data from multiple simulations. For example, we can combine speciation metrics across all simulations in our toy data folder:
#
#   ```{r, message = FALSE}
# data <- collect_sims(
#   root = "../data",
#   variables = c("time", "EI"),
#   parnames = c("ecosel",  "hsymmetry"),
#   id_column = "sim",
#   level = 1,
#   pattern = "example",
#   verbose = FALSE)
# head(data, 4)
# ```
#
# This will run the `read_data` function on the simulation folders found in the `root` directory, and bind the resulting data frames together by rows. Arguments of `read_data` can therefore be used in `collect_sims`. The argument `id_column` gives the name of the column that should identify the different simulations.
#
# The argument `root` can be a vector of paths to each of the simulation folders, but can also (as in the example) be one (or multiple) folder(s) where to find the simulations. If folders where simulations are stored are provided, the function will recursively look into these folders for a `pattern` characterizing simulation folders. The level of recursion can be specified by `level`, where zero (the default) means that `root` is assumed to be a vector of simulation folders. Use, for example, `level = 1` to access simulations from one or several parent directories. The recursive search for simulation folders is done with the `fetch_dirs` function, which may be useful in a number of cases where folder search is needed.
#
# The `collect_sims` function has a `check_extant` argument. If TRUE, the function will specifically look for non-missing, non-extinct, completed simulations. The two latter conditions rely on SLURM output log files being present in the folder, and so are applicable only for simulations run on the cluster. For simulations run locally, set this argument to FALSE.
#
# Note that in many of the functions dealing with multiple simulations, there is a `verbose` argument specifying whether to display messages, and a `pb` argument to display progress bars.
#
# ## Genetic architecture
#
# The genetic architecture can be retrieved from an architecture text file. Use `read_architecture` to do so and store the genetic architecture into a list with various fields, like this:
#
#   ```{r}
# arch <- read_architecture(root)
# str(arch)
# ```
#
# Note that by default this function looks for a file called "architecture.txt", but another name can be specified using the `filename` argument.
#
# Alternatively, the function `read_genome_architecture` returns a locus-based data frame version of the architecture, therefore omitting elements specific to network edges, but with additional information such as chromosome location or degree.
#
# ```{r}
# loci <- read_genome_architecture(root)
# head(loci, 4)
# ```
#
# Another one for networks specifically? I would have to look into network plotting functions to see what sort of stuff they can handle...
#
# ## Plotting functions
#
# Diagnostic plotting functions to eyeball single simulation results are provided.
#
# ### Locus-wise variables
#
# To visualize a genome scan of a given locus-specific variable at a specified time point (which defaults to the last time point if `t` is unspecified), for example, use:
#
#   ```{r, fig.width = 6, fig.height = 2, align = "center"}
# dplot_genome_scan(root, y = "genome_Fst")
# ```
#
# The output of these quick-plotting functions are all `ggplot` objects, so they can be further customized. For more details on how to make and customize plots using `ggplot2`, please refer to relevant online [documentation](https://ggplot2.tidyverse.org/).
#
# It is also possible to plot a genome scan through the entire duration of the simulation using a heatmap:
#
#   ```{r, fig.width = 4, fig.height = 5}
# dplot_genome_heatmap(root, y = "genome_Fst")
# ```
#
# Yet another way to explore genome scans is to look at densities of a given variable across groups of loci, e.g., loci underlying different traits:
#
#   ```{r, fig.width = 3, fig.height = 2}
# dplot_genome_violin(root, y = "genome_Fst", x = "trait")
# ```
#
# Densities across the genome can be viewed through time, using:
#
#   ```{r, fig.width = 4}
# dplot_genome_ridges(root, y = "genome_Fst", times = c(0, 1000, 5000, 10000, 15000))
# ```
#
# The same goes for locus-wise lines through time:
#
#   ```{r, fig.width = 4}
# dplot_genome_lines(root, y = "genome_Fst")
# ```
#
# which takes a `span` argument determining the span of the smoothing of the curves used.
#
# It is also possible to map locus-wise variables onto the gene network itself:
#
#   ```{r, fig.width = 6, fig.height = 5}
# dplot_network(root, y = "genome_Fst")
# ```
#
# ### Individual-wise variables
#
# We can also plot individual attributes such as trait values across the population at a certain time to make sure that speciation has happened:
#
#   ```{r, fig.width = 3, fig.height = 2}
# dplot_population_density(
#   root, y = "individual_trait", by = 3, j = 1, fill = "lightgreen"
# ) +
#   labs(x = "Ecological trait", y = "Density") +
#   xlim(c(-2, 2))
# ```
#
# Here, by setting `j = 1` we specify that we want to plot the density of the first trait. Because `individual_trait.dat` contains information for three traits for each individual, this long vector has to be split into three columns in `read_data`, hence the `by = 3`. Again, a `t` argument can be specified, and if not the last time point is shown.
#
# ```{r, fig.width = 5, fig.height = 3}
# dplot_population_bin2d(
#   root, y = "individual_trait", by = 3, j = 1, bins = 100
# ) +
#   labs(x = "Time (generations)", y = "Ecological trait", fill = "Count") +
#   scale_fill_continuous(type = "viridis")
# ```
#
# ### Simulation-wise variables
#
# It is possible, of course to plot simulation-wise observation through time or against each other. For example, we can plot the genome-wide Fst for all loci underlying the ecological trait through time (A), or we may want to show how reproductive isolation builds up with ecological divergence (B):
#
#   ```{r, fig.width = 6, fig.height = 2}
# p1 <- dplot_simulation_line(root, y = "Fst", by = 3, j = 1) +
#   labs(x = "Time (generations)", y = parse(text = "F[ST]"))
# p2 <- dplot_simulation_line(root, y = "RI", x = "EI") +
#   labs(x = "Ecological divergence", y = "Reproductive isolation")
# plot_grid(p1, p2, labels = c("A", "B"))
# ```
#
#
