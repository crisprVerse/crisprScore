crisprai_dependencies <- c("python==2.7.15",
                           "pip==20.1.1",
                           "wheel==0.36.2",
                           "setuptools==49.6.0",
                           "numpy==1.11.3",
                           "scipy==0.18.1",
                           "biopython==1.68",
                           "pysam=0.15.3",
                           "viennarna=2.4.18",
                           "scikit-learn==0.17.1",
                           "bx-python==0.7.3")

crisprai_dependencies_pip <- c("pandas==0.15.0")

#' @importFrom basilisk BasiliskEnvironment
env_crisprai <- BasiliskEnvironment(envname="crisprai_basilisk",
                                    pkgname="crisprScore",
                                    paths="python/crisprai",
                                    channels = c("bioconda", "conda-forge"),
                                    packages=crisprai_dependencies,
                                    pip=crisprai_dependencies_pip)

# env_azimuth <- BasiliskEnvironment(envname="azimuth_basilisk",
#                                    pkgname="crisprScore",
#                                    paths="python/azimuth",
#                                    packages=azimuth_dependencies,
#                                    channels = c("bioconda", "conda-forge"),
#                                    pip=azimuth_dependencies_pip)
# 
# env_lindel <- BasiliskEnvironment(envname="lindel_basilisk",
#                                   pkgname="crisprScore",
#                                   paths="python/lindel",
#                                   packages=lindel_dependencies,
#                                   channels=c("conda-forge", "bioconda"))



# if (.Platform$OS.type!="windows"){
#     env_deephf <- BasiliskEnvironment(envname="deephf_basilisk",
#                                       pkgname="crisprScore",
#                                       paths="python/deephf",
#                                       packages=deephf_dependencies,
#                                       channels=c("conda-forge", "bioconda"))
# 
# 
# 
#     env_enpamgb <- BasiliskEnvironment(envname="enpamgb_basilisk",
#                                        pkgname="crisprScore",
#                                        paths="python/enpamgb",
#                                        packages=enpamgb_dependencies,
#                                        channels=c("conda-forge", "bioconda"),
#                                        pip=enpamgb_dependencies_pip)
# 
#     env_deepcpf1 <- BasiliskEnvironment(envname="deepcpf1_basilisk",
#                                         pkgname="crisprScore",
#                                         paths="python/deepcpf1",
#                                         packages=deepcpf1_dependencies,
#                                         channels=c("conda-forge", "bioconda"))
# 
# } else {
#     env_deephf   <- NULL
#     env_enpamgb  <- NULL
#     env_deepcpf1 <- NULL
# }




