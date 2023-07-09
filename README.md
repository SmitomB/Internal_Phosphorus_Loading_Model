<h1 align="center">Characterizing and Projecting Internal Phosphorus Loading through Bayesian Mass-balance Modeling</h1> 

<h2 align="left">Description of the Repository</h2>
This repository contains the data and code used in the paper "Characterizing and Projecting Internal Phosphorus Loading through Bayesian Mass-balance Modeling". The details of the subrepositories are as follows:

| Sub-repositories | Description|
|:--------------------|:------------------------|
|DataFiles| Data files such as .csv, .xlsx, .txt are stored in this folder.|
|Pmodels| All model codes are stored in this folder. |
|RWorkSpaces| R workspaces required to run the model.|
|OutputFiles| Contains the output data from the P models.|

<h2 align="left">Prerequisites to run the model</h2>

<h3 align="left">GitHub Account</h3>

A _GitHub_ account is required to run the model. The details of creating a new _GitHub_ account can be easily found over the internet.

<h3 align="left">Installations and setup</h3>

Install _Git_, _R_, and _RStudio_ from respective websites. </br>

To link _Git_ to _RStudio_, open _RStudio_ and go to _Tools_ > _Global Options…_ click on _Git/SVN_. Check _Enable version control interface for RStudio projects_. Set the path to the _Git executable_ that was installed. Click _OK_ and restart _Rstudio_. If done correctly, the _Git_ icon will appear on the _RStudio_ toolbar. Next, configure _Git_ and set your user name and email (the email address you used to register on GitHub). You can directly open the _Git_ prompt from within _RStudio_. User name and email needs to be set only once. Go to _Tools_ > _Shell_ to open the Git Shell to tell Git your username and _GitHub_ email.

<h2 align="left">How to run the model</h2>

Step 1 : In the top-right corner of this page, click _Fork_. Under _Owner_, select the dropdown menu and click an owner for the forked repository. By default, forks are named the same as their upstream repositories. Optionally, to further distinguish your fork, in the _Repository name_ field, type a name. Click _Create fork_. A forked repository will be created in your GitHub account.</br> </br>
Step 2: Click on _Code_ icon in the forked repository to copy the HTTPS url. The HTTPS url will be something like https://github.com/SmitomB/Internal_Phosphorus_Loading_Model.git. </br></br>
Step 3 : Open RStudio and click _New Project..._ from _File_. Select _Version Control_ and _Git_. In the _Repository URL_, enter the HTTPS URL of the forked repository from your GitHUb account. Next, enter a _Project directory name_ and select the directory in which the you want to store the R project. Select _Create Project_</br></br>
Step 4 : Run the codes stored in _JLPMv016p001.Rmd_ to _JLPMv016p010.Rmd_ in _Pmodels_ repository.   



