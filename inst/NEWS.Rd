\name{crisprScoreNews}
\title{crisprScore News}
\encoding{UTF-8}


\section{Version 1.15.2}{\itemize{
\item Conda environment now have to be built outside of crisprScore. 
\item The tracer argument for RuleSet3 is now fixed; the previous version of crisprScore was now taking into account the tracer, and the Chen2013 tracer was always used. 
}}



\section{Version 1.15.1}{\itemize{
\item The scores developed in Python 2 are no longer available: DeepSpCas9, DeepCpf1, Azimuth, and the crisprai scores. 
}}


\section{Version 1.9.3}{\itemize{
\item Fixed documentation for enPamGB and DeepCpf1 algorithms. 
}}


\section{Version 1.9.2}{\itemize{
\item Fixed python environments based on the lastest changes in basilisk (using conda-forge installations).
}}

\section{Version 1.3.4}{\itemize{
\item Fixed unit testing for MIT scores.
}}

\section{Version 1.3.3}{\itemize{
\item One more fix to the MIT scoring. Forgot to save sysdata last time.
}}


\section{Version 1.3.1}{\itemize{
\item Fixed MIT formula. Previous calculations were erroneous. 
}}


\section{Version 1.0.0}{\itemize{
\item New package \pkg{crisprScore}, for on-target and off-target scoring of guide RNAs (gRNAs).
}}