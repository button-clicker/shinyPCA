@echo off
echo Creating Conda environment and installing dependencies...

REM Ensure Conda is initialized
CALL C:\Users\%USERNAME%\miniconda3\Scripts\activate.bat

REM Create and activate the environment
CALL conda create -n shinyPCA -y
CALL conda activate shinyPCA

REM Install R and required packages
CALL conda install -c conda-forge r-base=4.4.3 r-essentials r-shiny r-plotly r-htmlwidgets r-biocmanager r-reactable -y

REM Install Bioconductor DESeq2 package in R
Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('DESeq2')"

set BIN_PATH=%CONDA_PREFIX%\bin
if not exist "%BIN_PATH%" mkdir "%BIN_PATH%"
copy pca.R "%BIN_PATH%\pca.R"
copy run_shinyPCA.bat "%BIN_PATH%\run_shinyPCA.bat"

echo Environment setup complete!
pause
