# Lifetime
Statistics course exercise
The goal is to determine mass and lifetime of an unknown particle. Data consists of lifetime and mass measurements of particle candidates. masslife_signal.data only consists of signal events, while others have background as well. Background consists of a particle with no lifetime (immediate decay) and a particle with long (relatively) lifetime. Lifetime measurements are smeared with a Gaussian resolution (each particle type can have different resolution), this is reflected in lifetime PDF being an exponential decay convoluted with a gaussian.

To run the analysis, cern root is required with ROOFit installed.

//Loading the macros
root -l
.L analysis.C

//plots and saves a histogram of the data
plotDataSet()

//performs a lifetime fit on signal data only
signalLifetimeFit()

//performs a mass fit on signal data only
signalMassFit()

//performs a lifetime fit using lifetime data with signal + background
lifetimeFit()

//performs a mass fit using mass data with signal + background
//graphs a negative log likelihood dependance on fitting parameter (mass)
massFit()

//performs a lifetime fit using lifetime data with signal + background
//graphs a negative log likelihood dependance on fitting parameter (lifetime)
lifetimeFit()

//performs a combined fit of mass and lifetime using product of PDFs; this kind of fit is more "powerful" as it uses both lifetime and mass during fitting (not discarding one or the other as in previous fits)
//signalPDF=massGaussian*lifetimeConvolution; backgroundPDF=quadraticPolynome*(lifetimeConv1+fraction*lifetimeConv2)
//totalPDF = signalFraction*signalPDF + (1-signalFraction)*backgroundPDF
//also plots a confidence region as a contour plot
combinedFit()

//plots profile likelihood of fitted mass and lifetime () which follows a chi-squared(/2) distribution
profiles()
