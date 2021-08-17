crisprai_dependencies <- c("python==2.7.15",
                           "pip==20.1.1",
                           "wheel==0.36.2",
                           "setuptools==44.0.0",
                           "numpy==1.11.3",
                           "scipy==0.18.1",
                           "biopython==1.68",
                           "pysam=0.15.3",
                           "viennarna=2.4.18",
                           "scikit-learn==0.17.1",
                           "bx-python==0.7.3",
                           "libgomp==9.3.0",
                           "libgfortran5==9.3.0",
                           "libedit==3.1.20191231",
                           "sqlite==3.35.5")

crisprai_dependencies_pip <- c("pandas==0.15.0")


env_crisprai <- basilisk::BasiliskEnvironment(envname="crisprai_basilisk",
                                              pkgname="crisprScore",
                                              paths="python/crisprai",
                                              packages=crisprai_dependencies,
                                              channels = c("bioconda", "conda-forge"),
                                              pip=crisprai_dependencies_pip)


