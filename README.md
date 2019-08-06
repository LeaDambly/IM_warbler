README
================
Bob O’Hara & Lea Dambly
6/14/2019

## Introduction

This is the repository for the code for the integrated modelling example
for Isaac et al. (submitted), fitting a model for the black-throated
blue warbler (*Setophaga caerulescens*) in Pennsylvania, USA. This is an
extension of the analysis in Miller et al. (2019).

We use the following observation data:

  - Pennsylvania Breeding Bird Atlas point count data for part of
    Pennsylvania. This is taken straight from [the
    SI](https://besjournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2F2041-210X.13110&file=mee313110-sup-0001-supplementA.zip)
    of Miller et al. (2019)
  - eBird data from 2006-2009, down loaded from GBIF, using the spocc
    package (Chamberlain 2018)
  - North American BBS data (Pardieck and Hudson. 2018) for 2005-2009

These data are in the Data directory.

For environmental data we use elevation and canopy cover, both imported
directly into R

LEA, ADD DETAILS?

## Workflow

The workflow is as follows:

  - extract the observation data, with [ExtractData.R](ExtractData.R).
    This has already been done, and hte results are the .csv files in
    Data/
  - make the INLA stacks, using [MakeStacks.R](MakeStacks.R). This
    includes the covariate data, and creates some large objects that we
    don’t save in this repository.
  - Fit the model, with [FitWarblerModel.R](FitWarblerModel.R). This
    takes some time (overnight?), and creates another big .RData file.
  - Look at and plot the results, with
    [WarblerResults.R](WarblerResults.R). This produces some maps

## Files

This Directory

  - [ExtractData.R](ExtractData.R): Code to extract data from different
    sources, and save as .csvs.

  - [FitWarblerModel.R](FitWarblerModel.R): Code to, um, fit the warbler
    model

  - [IM\_warbler.Rproj](IM_warbler.Rproj): R project file

  - [MakeStacks.R](MakeStacks.R): Code to make INLA stacks from data.

  - [README.html](README.html): This file, unless you’re looking at…

  - [README.Rmd](README.Rmd): … this file (the Markdown file to make the
    html)

  - [References.bib](References.bib): BibTex file with references

  - [warblerfunctions.R](warblerfunctions.R): Lots of miscelleneous
    functions

  - [WarblerResults.R](WarblerResults.R): Code to plot & poke the
    results

  - [FitWarblerModelTwoRFs.R](FitWarblerModelTwoRFs.R): Code to fit
    model with a ransom field on eBird observation effort. This is in
    development, so doesn’t do anything useful

Data Folder

  - [BBA.csv](Data/BBA.csv): Pennsylvania breeding bird atlas data
  - [BBS.csv](Data/BBS.csv): North American breeding bird survey data
  - [eBird.csv](Data/eBird.csv): eBird data, downloaded from GBIF

Functions Folder

Some of these functions are from the current version of the PointedSDms
package. As that package is in early development, we cannot guarantee
the functions will be the same when you read this.

  - [Add2010Census.R](Functions/Add2010Census.R): Functions to add 2010
    census data to a data frame
  - [FitModel.R](Functions/FitModel.R): Function to fit a model with
    INLA (from Pointed SDMs)
  - [GetBBSData.R](Functions/GetBBSData.R): A couple of functions to
    update rBBS package
  - [GetNearestCovariate.R](Functions/GetNearestCovariate.R): Function
    to get covariate values nearest the data/integration point (from
    Pointed SDMs)
  - [MakeBinomStack.R](Functions/MakeBinomStack.R): Function to create
    stack for presence only points (from Pointed SDMs)
  - [MakeIntegrationStack.R](Functions/MakeIntegrationStack.R): Function
    to create stack for integration points (from Pointed SDMs)
  - [MakePointsStack.R](Functions/MakePointsStack.R): Function to create
    stack for presence only points (from Pointed SDMs)
  - [MakeProjectionGrid.R](Functions/MakeProjectionGrid.R): Function to
    create stack for predictions (from Pointed SDMs)
  - [MakeSpatialRegion.R](Functions/MakeSpatialRegion.R): Function to
    set up spatial structure for region (from Pointed SDMs)

## Acknowledgements

Most of the code was written by Lea Dambly, and then changed by Bob
O’Hara. Additional comments and help from …

## References

<div id="refs" class="references">

<div id="ref-Spocc">

Chamberlain, Scott. 2018. *Spocc: Interface to Species Occurrence Data
Sources*. <https://CRAN.R-project.org/package=spocc>.

</div>

<div id="ref-Miller2019">

Miller, David A. W., Krishna Pacifici, Jamie S. Sanderlin, and Brian J.
Reich. 2019. “The Recent Past and Promising Future for Data Integration
Methods to Estimate Species’ Distributions.” *Methods in Ecology and
Evolution* 10 (1): 22–37. <https://doi.org/10.1111/2041-210X.13110>.

</div>

<div id="ref-NABBS">

Pardieck, D. J. Ziolkowski Jr., K. L., and M.-A.R. Hudson. 2018. *North
American Breeding Bird Survey Dataset 1966 - 2017, Version 2017.0.* U.S.
Geological Survey, Patuxent Wildlife Research Center.
<https://doi.org/10.5066/F76972V8>.

</div>

</div>
