devtools::load_all()
setup_dir(all_options(main = "devs"))
all <- read_data(dirtree(main = "devs"), "all")

plot(all$newmoonnumber, all$total, type = "l")

y <- all$total
x <- all$newmoonnumber

abundances <- all
level <- "All"
metadata <- read_data(dirtree(main="devs"), "metadata")

mo <- LTAvg(all, metadata, level)

