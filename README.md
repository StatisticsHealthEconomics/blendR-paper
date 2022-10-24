# Blended Survival Curves: A new Approach to Extrapolation for Time-to-Event Outcomes from Clinical Trials in Health Technology Assessment

[![arXiv](https://img.shields.io/badge/arXiv-2206.00154-f9f107.svg)](https://arxiv.org/abs/2206.00154)

This is the `GitHub` repository to host the `R` code used to run the example in the paper

> Che, Zhaojing and Green, Nathan and Baio, Gianluca (2022) Blended Survival Curves: A new Approach to Extrapolation for Time-to-Event Outcomes from Clinical Trials in Health Technology Assessment, https://doi.org/10.48550/arxiv.2206.00154

The [example dataset](Data/TA174.csv) is based on the [CLL-8 trial](https://doi.org/10.1016/S0140-6736(10)61381-5) data, which were also used in [NICE technology appraisal TA174](https://www.nice.org.uk/guidance/ta174). A detailed report [`case_study_blend.html`](Scripts/case_study_blend.html) in the `Scripts` folder explains how to perform the blending analysis, step-by-step. 

In addition, the corresponding [`Rmd`](Scripts/case_study_blend.Rmd) file is used to manipulate and restructure the information to produce the modelling and all `R` code from the `markdown` document are extracted to the script [`CODE_case_study.R`](Scripts/CODE_case_study.R).   

## Obtaining the files

The folder contents can be downloaded as a `.zip` file of the repository.
This can be done either by using the "Clone or download - Download ZIP" button in GitHub or from inside of the R environment using the following code.

```r
download.file(url = "https://github.com/StatisticsHealthEconomics/blendR-paper/archive/main.zip",
              destfile = "blendR-paper.zip")

# unzip the .zip file
unzip(zipfile = "blendR-paper.zip")
```

## 👂 Feedback

Please open an issue in our [issue tracker](https://github.com/statisticsHealthEconomics/blendR/issues).
