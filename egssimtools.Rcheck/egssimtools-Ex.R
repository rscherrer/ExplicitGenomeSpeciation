pkgname <- "egssimtools"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "egssimtools-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('egssimtools')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("fetch_dirs")
### * fetch_dirs

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fetch_dirs
### Title: Fetch directories within directories
### Aliases: fetch_dirs

### ** Examples


root <- system.file("extdata", package = "egssimtools")
fetch_dirs(root, pattern = "example", level = 1)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fetch_dirs", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("find_completed")
### * find_completed

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: find_completed
### Title: Find completed simulations
### Aliases: find_completed

### ** Examples


root <- system.file("extdata", package = "egssimtools")
find_completed(root, pattern = "example_", level = 1)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("find_completed", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("find_extant")
### * find_extant

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: find_extant
### Title: Find extant simulations
### Aliases: find_extant

### ** Examples


root <- system.file("extdata", package = "egssimtools")
find_extant(root, pattern = "example_", level = 1)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("find_extant", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("find_extinct")
### * find_extinct

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: find_extinct
### Title: Find extinct simulations
### Aliases: find_extinct

### ** Examples


root <- system.file("extdata", package = "egssimtools")
find_extinct(root, pattern = "example_", level = 1)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("find_extinct", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("find_missing")
### * find_missing

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: find_missing
### Title: Find missing simulations
### Aliases: find_missing

### ** Examples


root <- system.file("extdata", package = "egssimtools")
find_missing(root, pattern = "example_", level = 1)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("find_missing", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_chromosomes")
### * get_chromosomes

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_chromosomes
### Title: Get locus chromosomal locations
### Aliases: get_chromosomes

### ** Examples


root <- system.file("extdata", "example_1", package = "egssimtools")
arch <- read_arch(root)
get_chromosomes(arch)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_chromosomes", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_degrees")
### * get_degrees

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_degrees
### Title: Get locus degrees
### Aliases: get_degrees

### ** Examples


root <- system.file("extdata", "example_1", package = "egssimtools")
arch <- read_arch(root)
get_degrees(arch)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_degrees", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_status")
### * get_status

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_status
### Title: SLURM status
### Aliases: get_status

### ** Examples


root <- system.file("extdata", "example_1", package = "egssimtools")
get_status(root)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_status", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("guess_nloci")
### * guess_nloci

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: guess_nloci
### Title: Guess the number of loci
### Aliases: guess_nloci

### ** Examples


root <- system.file("extdata", "example_1", package = "egssimtools")
guess_nloci(root)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("guess_nloci", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("is_extinct")
### * is_extinct

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: is_extinct
### Title: Did a simulation go extinct?
### Aliases: is_extinct

### ** Examples


root <- system.file("extdata", "example_1", package = "egssimtools")
is_extinct(root)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("is_extinct", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("is_missing")
### * is_missing

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: is_missing
### Title: Is a simulation missing?
### Aliases: is_missing

### ** Examples


root <- system.file("extdata", "example_1", package = "egssimtools")
is_missing(root)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("is_missing", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("read_arch")
### * read_arch

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: read_arch
### Title: Read architecture
### Aliases: read_arch

### ** Examples


root <- system.file("extdata", "example_1", package = "egssimtools")
read_arch(root)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("read_arch", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("read_arch_file")
### * read_arch_file

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: read_arch_file
### Title: Read an architecture file
### Aliases: read_arch_file

### ** Examples


f <- system.file("extdata", "example_1", package = "egssimtools")
f <- file.path(f, "architecture.txt")
read_arch_file(f)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("read_arch_file", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("read_arch_genome")
### * read_arch_genome

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: read_arch_genome
### Title: Read locus-specific genome architecture
### Aliases: read_arch_genome

### ** Examples


root <- system.file("extdata", "example_1", package = "egssimtools")
read_arch_genome(root)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("read_arch_genome", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("read_arch_network")
### * read_arch_network

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: read_arch_network
### Title: Read gene network architecture
### Aliases: read_arch_network

### ** Examples


root <- system.file("extdata", "example_1", package = "egssimtools")
read_arch_network(root)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("read_arch_network", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("read_binary")
### * read_binary

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: read_binary
### Title: Read binary file
### Aliases: read_binary

### ** Examples


f <- system.file("extdata", package = "egssimtools")
f <- file.path(f, "example_1/time.dat")
read_binary(f)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("read_binary", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("read_data")
### * read_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: read_data
### Title: Read simulation data
### Aliases: read_data

### ** Examples


root <- system.file("extdata", "example_1", package = "egssimtools")
read_data(root, c("time", "genome_Fst"), by = c(1, 300))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("read_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("read_genome")
### * read_genome

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: read_genome
### Title: Read locus-specific data through time
### Aliases: read_genome

### ** Examples


root <- system.file("extdata", "example_1", package = "egssimtools")

# Read Fst throughout the genome
read_genome(root, "genome_Fst")

# Read multiple metrics and attach architecture
read_genome(root, c("genome_Fst", "genome_Cst"), architecture = TRUE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("read_genome", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("read_param")
### * read_param

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: read_param
### Title: Read parameters
### Aliases: read_param

### ** Examples


root <- system.file("extdata", "example_1", package = "egssimtools")
read_param(root)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("read_param", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("read_param_file")
### * read_param_file

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: read_param_file
### Title: Read a parameter file
### Aliases: read_param_file

### ** Examples


f <- system.file("extdata", package = "egssimtools")
f <- file.path(f, "example_1/parameters.txt")
read_param_file(f)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("read_param_file", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("read_pop")
### * read_pop

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: read_pop
### Title: Read individual data through time
### Aliases: read_pop

### ** Examples


root <- system.file("extdata", "example_1", package = "egssimtools")

# Read individual trait values through time
read_pop(root, "individual_trait", by = 3)

# Read trait values, ecotypes and habitats
variables <- c("individual_trait", "individual_ecotype", "individual_habitat")
read_pop(root, variables, by = c(3, 1, 1))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("read_pop", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("read_sim")
### * read_sim

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: read_sim
### Title: Read simulation summary through time
### Aliases: read_sim

### ** Examples


root <- system.file("extdata", "example_1", package = "egssimtools")

# Read the ecological divergence through time
read_sim(root, "EI")

# Read multiple speciation metrics
read_sim(root, c("EI", "RI"))

# Read genome-wide metrics
read_sim(root, "Fst", by = 3)

# Or a combination of multiple metrics with different dimensions
read_sim(root, c("EI", "Fst"), by = c(1, 3))

# Even locus-specific data in the wide-format
read_sim(root, "genome_Fst", by = 300)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("read_sim", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rename_str")
### * rename_str

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rename_str
### Title: Function to rename the columns of a data frame given a character
###   vector
### Aliases: rename_str

### ** Examples


df <- data.frame(a = 1, b = 1)
rename_str(df, c("aa", "bb"))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rename_str", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("set_param_file")
### * set_param_file

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: set_param_file
### Title: Set up parameter file
### Aliases: set_param_file

### ** Examples


set_param_file(list(ecosel = 1.4, nvertices = "100 100 100"))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("set_param_file", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("smoothen_data")
### * smoothen_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: smoothen_data
### Title: Smoothen lines using LOESS
### Aliases: smoothen_data

### ** Examples


root <- system.file("extdata", "example_1", package = "egssimtools")
data <- read_genome(root, "genome_Fst", architecture = TRUE)
smoothen_data(data, "time", "genome_Fst", line = "locus")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("smoothen_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
