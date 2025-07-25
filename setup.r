# Copyright (C) 2025 Blair Shevlin <blairshevlin@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Record of Revisions
#
# Date            Programmers                         Descriptions of Change
# ====         ================                       ======================
# 2025/07/24    Blair Shevlin                           wrote original code

# List all required packages
required_packages <- c(
  # Core tidyverse and data manipulation
  "tidyverse",
  "readr",
  "purrr",
  "broom", 
  # File system management
  "fs",
  "here",
  # Plotting and visualization
  "ggpubr", 
  "patchwork",
  "cowplot", 
  "scales",
  "ggrepel", 
  "GGally",
  "ggeffects", 
  "ggnewscale", 
  # Colors and themes
  "RColorBrewer",
  # Statistical analysis
  "lmerTest",
  "rstatix",  
  "effectsize",
  "emmeans",
  "caTools",
  # Data import/export
  "readxl"
)

# Check which packages are missing
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

# Install missing packages
if(length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, dependencies = TRUE)
  cat("Installation complete!\n")
} else {
  cat("All required packages are already installed.\n")
}

# Load all packages and check for any issues
cat("Loading packages...\n")
failed_loads <- c()

for(pkg in required_packages) {
  result <- tryCatch({
    library(pkg, character.only = TRUE)
    TRUE
  }, error = function(e) {
    failed_loads <<- c(failed_loads, pkg)
    cat("Warning: Failed to load", pkg, "\n")
    FALSE
  })
}

if(length(failed_loads) > 0) {
  cat("Failed to load the following packages:\n")
  cat(paste(failed_loads, collapse = ", "), "\n")
  cat("You may need to install them manually or check for compatibility issues.\n")
} else {
  cat("All packages loaded successfully!\n")
}

# Display session info for reproducibility
cat("\nSession Information\n")
sessionInfo()

cat("\nPackage Versions\n")
package_versions <- installed.packages()[required_packages, "Version"]
names(package_versions) <- required_packages
print(package_versions)

cat("\nSetup complete\n")