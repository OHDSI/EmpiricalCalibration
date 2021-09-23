# Returns the list of reverse dependencies, and installs all dependencies of those packages
prepareForReverseDependencyCheck <- function(pkgdir = ".") {
  description <- scan(file.path(pkgdir, "DESCRIPTION"), what = character(), sep = "|", quiet = TRUE) 
  rootPackage <- gsub("^Package:[ \t]*", "", description[grepl("^Package:", description)])
  
  packageListUrl <- "https://raw.githubusercontent.com/OHDSI/Hades/master/extras/packages.csv"
  gitHubOrganization <- "ohdsi"
  hadesPackageList <- read.table(packageListUrl, sep = ",", header = TRUE) 
  
  dependencies <- lapply(hadesPackageList$name, getPackageDependenciesFromGitHub)
  dependencies <- do.call(rbind, dependencies)
  hadesDependencies <- dependencies[dependencies$dependency %in% hadesPackageList$name & 
                                      dependencies$type != "Suggests", ]
  
  recurseReverseDependencies <- function(package, maxDepth = 10) {
    # print(package)
    if (maxDepth == 0) {
      warning("Maximum recursion depth reached. Are there circular dependencies?")
      return(c())
    }
    reverseDependencies <- hadesDependencies$package[hadesDependencies$dependency == package]
    recursiveReverseDependencies <- lapply(reverseDependencies, recurseReverseDependencies, maxDepth = maxDepth - 1)
    recursiveReverseDependencies <- do.call(c, recursiveReverseDependencies)
    return(unique(c(reverseDependencies, recursiveReverseDependencies)))
  }
  reverseDependencies <- recurseReverseDependencies(rootPackage)
  
  if (length(reverseDependencies) == 0) {
    writeLines(sprintf("No reverse dependencies found for package '%s'", rootPackage))
    return(data.frame())
  }
  
  # Includes suggests of packages to check, but not suggests of deeper dependencies:
  packagesToInstall <- unique(c(reverseDependencies,
                                dependencies$dependency[dependencies$package %in% reverseDependencies]))
  packagesToInstall <- c(packagesToInstall, "formatR") # Required for some vignettes
  # Don't install packages that are already installed:
  packagesToInstall <- packagesToInstall[!packagesToInstall %in% rownames(installed.packages())]
  
  packagesToInstallFromGitHub <- packagesToInstall[packagesToInstall %in% hadesPackageList$name[!hadesPackageList$inCran]]
  packagesToInstallFromCran <- packagesToInstall[!packagesToInstall %in% packagesToInstallFromGitHub]
  
  if (length(packagesToInstallFromCran) > 0) {
    remotes::install_cran(packagesToInstallFromCran)
  }
  for (package in packagesToInstallFromGitHub) {
    remotes::install_github(sprintf("%s/%s", gitHubOrganization, package), upgrade = FALSE)
  }
  
  return(hadesPackageList[hadesPackageList$name %in% reverseDependencies, ])
}

#   for (package in reverseDependencies) {
#     if (hadesPackageList$inCran[hadesPackageList$name == package]) {
#       source <- "CRAN"
#     } else {
#       source <- "GitHub"
#     }
#     checkPackage(package, source)
#   }
# }

checkPackage <- function(package, inCran) {
  writeLines(sprintf("*** Checking package '%s' ***", package))
  gitHubOrganization <- "ohdsi"
  if (inCran) {
    sourcePackage <- remotes::download_version(package, type = "source")
    on.exit(unlink(sourcePackage))
  } else {
    sourcePackage <- remotes::remote_download(remotes::github_remote(sprintf("%s/%s", gitHubOrganization, package)))
    on.exit(unlink(sourcePackage))
  } 
  sourceFolder <- tempfile(pattern = package)
  dir.create(sourceFolder)
  on.exit(unlink(sourceFolder, recursive = TRUE), add = TRUE)
  untar(sourcePackage, exdir = sourceFolder)
  sourcePath <- list.dirs(sourceFolder, full.names = TRUE, recursive = FALSE)
  docDir <- file.path(sourcePath, "inst", "doc")
  if (dir.exists(docDir)) {
    unlink(docDir, recursive = TRUE)
  }
  rcmdcheck::rcmdcheck(path = sourcePath, args = c("--no-manual", "--no-multiarch"), error_on = "warning")
}

getPackageDependenciesFromGitHub <- function(package) {
  descriptionUrlTemplate <- "https://raw.githubusercontent.com/OHDSI/%s/master/DESCRIPTION"
  
  description <- scan(sprintf(descriptionUrlTemplate, package), what = character(), sep = "|", quiet = TRUE) 
  dependencies <- lapply(X = c("Depends", "Imports", "LinkingTo", "Suggests"), 
                         FUN = extractDependenciesFromDescriptionSection, 
                         description = description)
  dependencies <- do.call(rbind, dependencies)
  dependencies <- dependencies[dependencies$dependency != "R", ]
  coreRPackages <- rownames(installed.packages(priority = "base"))
  dependencies <- dependencies[!dependencies$dependency %in% coreRPackages, ]
  dependencies$package <- rep(package, nrow(dependencies))
  return(dependencies)
}

extractDependenciesFromDescriptionSection <- function(section, description) {
  tagsPos <- grep(":", description)
  sectionPos <- grep(sprintf("%s:", section), description)
  if (length(sectionPos) != 0) {
    endOfSection <- ifelse(sectionPos < max(tagsPos), min(tagsPos[tagsPos > sectionPos]), length(description) + 1)
    dependencies <- gsub("[\t ,]|(\\(.*\\))", "", gsub(sprintf("%s:", section), "", description[sectionPos:(endOfSection - 1)]))
    dependencies <- dependencies[dependencies != ""]
    if (length(dependencies) > 0) {
      return(data.frame(type = section,
                        dependency = dependencies))
    } 
  }
  return(NULL)
}
