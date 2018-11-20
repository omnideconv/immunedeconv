# Configuration file. All configuration options are stored in an environment.
config_env = new.env()

# xCell will fail if ran with more cores than available.
config_env$xcell_cores = min(4, detectCores())
