# @file PackageMaintenance
#
# Copyright 2025 Observational Health Data Sciences and Informatics
#
# This file is part of EmpiricalCalibration
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Manually delete package from library. Avoids "Already in use" message when rebuilding
unloadNamespace("EmpiricalCalibration")
.rs.restartR()
folder <- system.file(package = "EmpiricalCalibration")
folder
unlink(folder, recursive = TRUE, force = TRUE)
file.exists(folder)

# Format and check code:
styler::style_pkg()
OhdsiRTools::checkUsagePackage("EmpiricalCalibration")
OhdsiRTools::updateCopyrightYearFolder()
devtools::spell_check()

# Create manual and vignettes:
unlink("extras/EmpiricalCalibration.pdf")
system("R CMD Rd2pdf ./ --output=extras/EmpiricalCalibration.pdf")

rmarkdown::render("vignettes/EmpiricalPCalibrationVignette.Rmd",
                  output_file = "../inst/doc/EmpiricalPCalibrationVignette.pdf",
                  rmarkdown::pdf_document(latex_engine = "pdflatex",
                                          toc = TRUE,
                                          number_sections = TRUE))
unlink("inst/doc/EmpiricalPCalibrationVignette.tex")

rmarkdown::render("vignettes/EmpiricalCICalibrationVignette.Rmd",
                  output_file = "../inst/doc/EmpiricalCiCalibrationVignette.pdf",
                  rmarkdown::pdf_document(latex_engine = "pdflatex",
                                          toc = TRUE,
                                          number_sections = TRUE))
unlink("inst/doc/EmpiricalCiCalibrationVignette.tex")

rmarkdown::render("vignettes/EmpiricalMaxSprtCalibrationVignette.Rmd",
                  output_file = "../inst/doc/EmpiricalMaxSprtCalibrationVignette.pdf",
                  rmarkdown::pdf_document(latex_engine = "pdflatex",
                                          toc = TRUE,
                                          number_sections = TRUE))
unlink("inst/doc/EmpiricalMaxSprtCalibrationVignette.tex")

pkgdown::build_site()
OhdsiRTools::fixHadesLogo()

# Release package:
devtools::check_win_devel()

rhub::rc_submit(platforms = "atlas")

devtools::release()
