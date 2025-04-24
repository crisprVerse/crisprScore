azimuth_dependencies <- c("python==2.7",
                          "biopython==1.74",
                          "matplotlib==2.2.3",
                          "numpy==1.11.3",
                          "pandas==0.23.4",
                          "pip==20.1.1",
                          "scipy==1.2.1",
                          "six==1.16.0")
azimuth_dependencies_pip <- c("scikit-learn==0.17.1")

lindel_dependencies <- c("python==3.6",
                         "numpy==1.19.5",
                         "pip==21.0.1",
                         "scipy==1.5.3")

enpamgb_dependencies <- c("python==3.6",
                          "biopython==1.78",
                          "h5py==2.10.0",
                          "Keras-Preprocessing==1.1.2",
                          "numpy==1.19.5",
                          "pandas==1.1.5",
                          "pip==21.0.1",
                          "scikit-learn==0.21.2",
                          "scipy==1.5.3",
                          "setuptools==49.6.0",
                          "six==1.15.0")
enpamgb_dependencies_pip <- c("tensorflow==2.4.1")


deepspcas9_dependencies <- c("python==2.7",
                             "scipy==1.1.0",
                             "numpy==1.14.5")
deepspcas9_dependencies_pip <- c("tensorflow==1.4.1")


#3.7.7
rs3_dependencies <- c("python==3.7",
                      "pandas==1.3.5",
                      "scikit-learn==1.0.2",
                      "biopython==1.78",
                      "lightgbm==3.2.1")
rs3_dependencies_pip <- c("sglearn==1.2.5")



crisprai_dependencies <- c("python==2.7",
                           "biopython==1.68",
                           "bx-python==0.7.3",
                           "pysam==0.15.3",
                           "scikit-learn==0.17.1",
                           "scipy==0.18.1",
                           "viennarna==2.4.18")
crisprai_dependencies_pip <- c("pandas==0.15.0")



deephf_dependencies <- c("python==3.6",
                         "viennarna==2.4.5")
deephf_dependencies_pip <- c("biopython==1.71",
                             "h5py==2.9.0",
                             "dotmap==1.2.20",
                             "numpy==1.14.0",
                             "scipy==1.1.0",
                             "pandas==0.25.3",
                             "Keras==2.1.6",
                             "gpy==1.9.8",
                             "scikit-learn==0.19.1",
                             "matplotlib==3.1.1",
                             "tensorboard==1.8.0",
                             "tensorflow==1.8.0",
                             "Theano==1.0.5")


deepcpf1_dependencies <- c("python==2.7")
deepcpf1_dependencies_pip <- c("h5py==2.10.0",
                           "Keras==2.1.5",
                           "numpy==1.16.5",
                           "scipy==1.2.1",
                           "tensorflow==1.0.0",
                           "Theano==1.0.1")

#' @importFrom basilisk BasiliskEnvironment
env_azimuth <- BasiliskEnvironment(envname="azimuth_basilisk",
                                   pkgname="crisprScore",
                                   paths="python/azimuth",
                                   packages=azimuth_dependencies,
                                   channels = c("conda-forge", "bioconda"),
                                   pip=azimuth_dependencies_pip)

env_lindel <- BasiliskEnvironment(envname="lindel_basilisk",
                                  pkgname="crisprScore",
                                  paths="python/lindel",
                                  packages=lindel_dependencies,
                                  channels=c("conda-forge", "bioconda"))

env_crisprai <- BasiliskEnvironment(envname="crisprai_basilisk",
                                    pkgname="crisprScore",
                                    packages=crisprai_dependencies,
                                    channels = c("bioconda",
                                                 "conda-forge"),
                                    pip=crisprai_dependencies_pip)

env_deepspcas9 <- BasiliskEnvironment(envname="deepspcas9__basilisk",
                                      pkgname="crisprScore",
                                      packages=deepspcas9_dependencies,
                                      channels = c("bioconda", "conda-forge"),
                                      pip=deepspcas9_dependencies_pip)

env_rs3 <- BasiliskEnvironment(envname="rs3__basilisk",
                               pkgname="crisprScore",
                               packages=rs3_dependencies,
                               channels = c("bioconda", "conda-forge"),
                               pip=rs3_dependencies_pip)

if (.Platform$OS.type!="windows"){
    env_deephf <- BasiliskEnvironment(envname="deephf_basilisk",
                                      pkgname="crisprScore",
                                      paths="python/deephf",
                                      packages=deephf_dependencies,
                                      channels=c("conda-forge", "bioconda"),
                                      pip=deephf_dependencies_pip)



    env_enpamgb <- BasiliskEnvironment(envname="enpamgb_basilisk",
                                       pkgname="crisprScore",
                                       paths="python/enpamgb",
                                       packages=enpamgb_dependencies,
                                       channels=c("conda-forge", "bioconda"),
                                       pip=enpamgb_dependencies_pip)

    env_deepcpf1 <- BasiliskEnvironment(envname="deepcpf1_basilisk",
                                        pkgname="crisprScore",
                                        paths="python/deepcpf1",
                                        packages=deepcpf1_dependencies,
                                        channels=c("conda-forge", "bioconda"),
                                        pip=deepcpf1_dependencies_pip)

} else {
    env_deephf   <- NULL
    env_enpamgb  <- NULL
    env_deepcpf1 <- NULL
}








