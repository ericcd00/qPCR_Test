# qPCR_Data_Analysis <img src="https://avatars.githubusercontent.com/u/157121172?s=400&u=b127d13c6c44a0f6e5368f79be98176105d1c30b&v=4" align="right" height="130"></a>

[![Lifecycle: experimental](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

## Overview

The LAB_qPCR_Data_Analysis repository aims to help you automatize the following processes:
* qPCR data transformation (from a pdf file to an xlsx file) --> `Read_pdf()` 
* Delta-delta ct analysis --> `ddct_analysis()` 
* Graphical visualization of the delta-delta ct analysis --> `ddct_plot()` 
* Data visualization and extreme values detection  --> `qPCR_report()` 
* Integration of different qPCR experiments --> `qPCR_experiments()` 

The main scripts contained in this repository are the following:

* `qPCR_Analysis_User_Script.R` - This file has to be downloaded. It contains all
the instructions necessary for the user to use the different functionalities of
the script.
* `Misc_qPCR_Analysis.R` - This file is stored in the "R/" subfolder and contains 
all the functions that can be used from this repository. The user will load the 
functions in the "Source" step, that will be explained in the __"Steps"__ section.

The first thing you will see when you open the script will be a brief summary of the 
distribution of the script. It is advised to read it thoroughly to avoid future 
errors during the execution of the functions.


## Steps

There are 3 main steps in this script. In the __Introduction__ step, the needed parameters 
will be loaded into the environment. In the __Source__ step, the functions will be loaded 
from the Miscelaneous script, and the libraries will also be installed and/or loaded. Lastly, 
in the __Execution__ step, the needed functions will be executed.

### 1. Introduction

The parameters that you will have to fill depend on the functions that you will execute on
the third step.

The first few parameters have to be always filled:

#### result_path
```{r}
result_path = "/home/user/Documents/Files/Projects/LAB_qPCR_Data_Analysis/Demo_Files/Results"
```

This parameter has to be the folder where you want your results to be stored. Make sure the 
path is written correctly. If you want to check if the path exists, run this code: `file.exists(result_path)`.
If "TRUE", it exists.

#### exp_name
```{r}
exp_name = "example_test"
```
This will be the first part of the namefiles generated by the different functions. For example, in the example above, 
the execution of the functions will generate files named _example_test_..._pdf_ (for example).

#### housekeeping_genes
```{r}
housekeeping_genes = c("GENE1")
```
The housekeeping genes are stably expressed in all cells and conditions, and therefore are used as control to compare 
the relative expression of other genes in different conditions. This parameter can be filled by one or more housekeeping
genes. If you want to use more than one housekeeping gene, run the line as shown below:

```{r}
housekeeping_genes = c("GENE1", "GENE2")
```

#### analized_groups
```{r}
analized_groups = c("KO", "WT")
```
This parameter contains the different groups (not replicates) that you will study and compare.


The CT data is obtained from a pdf file, and therefore the first step is to read this pdf. The
function you can use is the `Read_PDF()` function. If you intend to use this function, you will 
need to fill the following parameters:

#### pdf_path
```{r}
pdf_path = "/home/user/Documents/Files/Projects/LAB_qPCR_Data_Analysis/Demo_Files/test_files/test.pdf"
```
This parameter will tell the `Read_PDF()` function the path where the pdf file is stored. Notice 
the last thing in the parameter must be a ".pdf"

#### wells
```{r}
wells = 27
```
This parameter has to be filled with the number of rows that each table contains (it has to be the same
number for all the tables).

If you have already read this pdf and want to use the resulting excel as the input, you must fill this 
parameter instead:

#### Excel_path
```{r}
Excel_path = "/home/user/Documents/Files/Projects/LAB_qPCR_Data_Analysis/Demo_Files/Results/example_test_qPCR_data.xlsx"
```
That parameter will be used by another function, `read_excel()`. It can be used if you remove outliers or
unreliable labeled data.

After having the results of the `Read_PDF()` or `read_excel()` functions, the delta-delta ct analysis can be performed.
This analysis will need one more parameter to be filled:


#### Control_variable
```{r}
control_variable = "CT"
```
This variable will be used as a control for the ddct analysis.

The results of the ddct analysis can be used for the next function `ddct_plot()`. This function requires
between 3 and 5 variables:

#### Genes_of_interest
```{r}
Genes_of_interest = c("GENE1", "GENE2", "GENE3", "GENE4")
```
The genes in this parameter will be the ones to be plotted.

#### Groups_of_interest
```{r}
Groups_of_interest = c("WT", "KO")
```
The conditions in this parameter will be the ones to be plotted.

#### Title

```{r}
title = "DDCT test results"
```
This will be the title on top of the graph. Leave it empty "" if you do not want any title.

#### PairsTest

```{r}
PairsTest = "t.test"             # "t.test" or "wilcox.test"
```
This will be the test that will be used to compare each pair of conditions. It will only be performed
if there are 3 or more replicates of each condition.

#### GroupTest

```{r}
GroupTest = "kruskal.test"       # "anova" or "kruskal.test"
```
This will be the test that will be used to compare all conditions. It will only be performed if 
there are 3 or more conditions to compare.


If you want to compare the results of different ddct analysis of the same experiment (with the same conditions
and replicates), you can use the `qPCR_experiments()` function. This will require the 2 last parameters:

#### exp_paths
```{r}
exp_paths = list("/home/user/Documents/Files/Projects/LAB_qPCR_Data_Analysis/Demo_Files/Results/xlsx/example_test_qPCR_ddct.xlsx",
                 "/home/user/Documents/Files/Projects/LAB_qPCR_Data_Analysis/Demo_Files/test_files/example_test2_qPCR_ddct.xlsx", 
                 "/home/user/Documents/Files/Projects/LAB_qPCR_Data_Analysis/Demo_Files/test_files/example_test3_qPCR_ddct.xlsx")
```
This parameter has to be a list of paths where the ddct results will be stored.

#### exp_names
```{r}
exp_names = c("Test1", "Test2", "Test3")
```
This parameter will be used to diferentiate between experiments.

After filling all the needed parameters, you will go to the next step.

### 2. Source

This line of code will load all the functions from the `Misc_qPCR_Analysis.R` script and will
install and/or load the libraries that the functions will use. It is recommended to run this step 
after closing and opening RStudio even if the session was not closed, to ensure you are always running
the latest version of the script. 

```{r}
source("/home/user/Documents/Files/Projects/LAB_qPCR_Data_Analysis/R/Misc_qPCR_Analysis.R")
```

### 3. Execution

There are different functions that the script can execute:

#### Read_PDF

The Read_PDF function is used to extract the qPCR data from a pdf file and export it to a xlsx file.
To use this function it is required to have the "Java Developer Kit" (JDK) previously installed.  

The function follows this structure:

```{r}
Results <- Read_PDF(pdf_path = pdf_path, 
                    result_path = result_path, 
                    exp_name = exp_name, 
                    wells = wells, 
                    analized_groups = analized_groups, 
                    housekeeping_genes = housekeeping_genes)
```

#### read_excel

The read_excel function is used to recover the qPCR data from a previous pdf file that has already been
transformed into an xlsx file using the Read_PDF function. It can be very useful if outliers have been removed 
or to avoid refilling the pdf parameters.

The function follows this structure:

```{r}
Results <- read_excel(Excel_path = Excel_path)
```

Both Read_PDF and read_excel functions store their results in a variable called "Results". Do not change its name, 
because it is used by the following functions.


#### ddct_analysis

This function performs the delta-delta ct analysis using the control variable and housekeeping genes that have been
previously assigned. The resulting file is an xlsx file containing the final and inbetween values of the analyisis,
including the standard deviation of the different replicates.
```{r}
ddct_results <- ddct_analysis(data = Results,
                              result_path = result_path, 
                              exp_name = exp_name, 
                              housekeeping_genes = housekeeping_genes,
                              control_variable = control_variable)
```


#### qPCR_report

After obtaining the results from the delta-delta ct analysis, you can use this data to automatically create a report
of the data. This report creates different Principal Component Analysis, Heatmaps, Boxplots... and uses the Rosner 
test to find posible outliers within the data.

```{r}
qPCR_report(ct_data = Results,
            dd_data = ddct_results,
            pal = 1, # Number from 1 to 5
            output_dir = paste0(result_path),
            output_file = paste0(exp_name, "_qPCR_Analysis_Report"))
```

#### plot_result

This function will create a plot comparing the ddct values between the different conditions for 
each gene. It will consider if there are enough replicates or groups to perform the statistics.

```{r}
plot_result <- ddct_plot(data = Results,
                      ddct_values = ddct_results,
                      Genes_of_interest = Genes_of_interest,
                      title = title,
                      result_path = result_path,
                      Groups_of_interest = Groups_of_interest,
                      analized_groups = analized_groups,
                      GroupTest = GroupTest,
                      PairsTest = PairsTest,
                      exp_name = exp_name)
```

#### complete_ddct_analysis

If you want to do the Read_PDF, ddct_analysis, qPCR_report and ddct_plot functions, you can 
use this combined function:

```{r}
complete_ddct_analysis <- complete_ddct_analysis()
```

#### qPCR_experiments

The last function of this script is used to combine different experiments. It compares the means of the ddct values
of each group for each gene between the different experiments to observe if there are differences between experiments.

```{r}
combine_qPCRs <- qPCR_experiments(exp_paths = exp_paths,
                              exp_names = exp_names,
                              Genes_of_interest = Genes_of_interest,
                              title = title,
                              result_path = result_path,
                              Groups_of_interest = Groups_of_interest,
                              exp_name = exp_name)
```


## Example

The next table is an example of how the pdf tables are distributed:

#### Abs Quant/2nd Derivative Max for GENE1 (Abs Quant/2nd Derivative Max)
Results

| Inc | POS | Name |   Type  |   CP  |  Concentration |   Standard |      Status      |
|:---:|:---:|:----:|:-------:|:-----:|:--------------:|:----------:|:----------------:|
|     |  A1 |  CT1 | Unknown | 22.31 |                |            |                  |
|     |  A2 |  CT2 | Unknown | 22.35 |                |            |                  |
|     |  A3 |  CT3 | Unknown | 20.27 |                |            |         <        |
|     |  B1 |  CT4 | Unknown | 22.32 |                |            |                  |
|     |  B2 |  KO1 | Unknown | 21.29 |                |            |                  |
|     |  B3 |  KO2 | Unknown | 21.02 |                |            |                  |
|     |  C1 |  KO3 | Unknown | 20.98 |                |            |                  |
|     |  C1 |  KO4 | Unknown | 21.09 |                |            |                  |
|     |  C1 |  KO5 | Unknown | 21.21 |                |            |                  |
|     |  D1 |  H2O | Unknown |       |                |            |                  |
|     |  D2 |  H2O | Unknown |       |                |            |                  |
|     |  D3 |  H2O | Unknown | 25.27 |                |            |                  |

The title "Abs Quant..." includes the name of the gene that the table is referring to. If this gene
is housekeeping, it must be placed in the "housekeeping_genes" parameter.

The number of rows of this table is 12. That number must be placed in the "wells" parameter: `wells = 12`

The "Name" column contains the name of each __sample__ of the experiment. The conditions
of this study could be CT and KO. Therefore, the "analized_groups" parameter would be: `analized_groups = c("CT", "KO")`
