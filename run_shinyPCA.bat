@echo off
CALL C:\Users\%USERNAME%\miniconda3\Scripts\activate.bat

REM Run the Shiny app
CALL conda activate shinyPCA
cd %CONDA_PREFIX%\bin
R -e "shiny::runApp('C:\\Users\\%USERNAME%\\miniconda3\\envs\\shinyPCA\\bin\\pca.R', launch.browser = TRUE)"

REM Pause to keep the window open (optional)
pause