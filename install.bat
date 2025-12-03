R -e "devtools::document(pkg = '/projects/tcrClustR')"
@REM devtools::reload(pkgload::inst("tcrClustR"))

R CMD build . && R CMD INSTALL --build tcrClustR_0.0.0.9000.tar.gz
@REM rm .\tcrClustR_0.0.0.9000.tar.gz